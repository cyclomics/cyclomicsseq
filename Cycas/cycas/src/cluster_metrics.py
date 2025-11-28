from abc import ABC, abstractmethod
from typing import Dict, List


class Metric(ABC):
    """Abstract base for cluster evaluation strategies."""

    @abstractmethod
    def calculate_metric(self, cluster_id: int, all_blocks: List[Dict]) -> float:
        """Calculate the strategy-specific metric."""
        pass


class SwitchFractionStrategy(Metric):
    """Calculates fraction of switches between different clusters along the concatemer."""

    def calculate_metric(self, cluster_id: int, all_blocks: List[Dict]) -> float:
        if len(all_blocks) < 2:
            return 0.0
        sorted_blocks = sorted(all_blocks, key=lambda b: int(b["read_pos"]))
        switches = sum(
            1
            for i in range(1, len(sorted_blocks))
            if sorted_blocks[i]["rel_id"] != sorted_blocks[i - 1]["rel_id"]
        )
        return switches / (len(sorted_blocks) - 1)


class FragmentationIndexStrategy(Metric):
    """Calculates the fraction of segments that a cluster is split into along the concatemer."""

    def calculate_metric(self, cluster_id: int, all_blocks: List[Dict]) -> float:
        sorted_all = sorted(all_blocks, key=lambda b: int(b["read_pos"]))

        num_segments = 0
        in_segment = False
        cluster_size = 0

        for b in sorted_all:
            if b["rel_id"] == cluster_id:
                cluster_size += 1
                if not in_segment:
                    num_segments += 1
                    in_segment = True
            else:
                in_segment = False

        if cluster_size <= 1:
            return 0.0

        return (num_segments - 1) / (cluster_size - 1)


class InterlacingRatioStrategy(Metric):
    """Calculates the fraction of interlaced positions of a cluster with other clusters"""

    def calculate_metric(self, cluster_id: int, all_blocks: List[Dict]) -> float:
        sorted_all = sorted(all_blocks, key=lambda b: int(b["read_pos"]))

        # Record positions of blocks belonging to target cluster by ID
        positions = [i for i, b in enumerate(sorted_all) if b["rel_id"] == cluster_id]
        if len(positions) <= 1:
            return 0.0

        # For each position, check if next block belongs to different cluster
        # Count how many times that happens
        interlaced_transitions = sum(
            1 for idx in positions[:-1] if sorted_all[idx + 1]["rel_id"] != cluster_id
        )
        return interlaced_transitions / (len(positions) - 1)


class OrientationSwitchStrategy(Metric):
    """Measures switches in orientation of a cluster_id based on longest run of consistent orientation."""

    def calculate_metric(self, cluster_id: int, all_blocks: List[Dict]) -> float:
        sorted_all = sorted(
            [b for b in all_blocks if b["rel_id"] == cluster_id],
            key=lambda b: int(b["read_pos"]),
        )

        if not sorted_all:
            return 0.0
        if len(sorted_all) == 1:
            return 1.0

        orientations = [b["ori"] for b in sorted_all]

        # Find longest run of consecutive same-orientation values
        longest_run = 1
        current_run = 1
        for i in range(1, len(orientations)):
            if orientations[i] == orientations[i - 1]:
                current_run += 1
                longest_run = max(longest_run, current_run)
            else:
                current_run = 1

        return longest_run / len(orientations)


class NumOrientationSwitches(Metric):
    """Calculate the number of orientation switches between forward and reverse

    Ignores unknown alignment segments

    Examples (F=forward, R=reverse, U=unknown):
        FFF -> 0
        RRR -> 0
        FFFRRR -> 1
        RFR -> 2
        FRFR -> 3
        FFFURRR -> 1
        FFUFF -> 0
    """

    def calculate_metric(self, cluster_id: int, all_blocks: List[Dict]) -> float:
        # get the alignments relevant to the cluster_id and sort them
        sorted_all = sorted(
            [b for b in all_blocks if b["rel_id"] == cluster_id],
            key=lambda b: int(b["read_pos"]),
        )
        oris = [b["ori"] for b in sorted_all] if sorted_all else []

        return self._calculate_metric(oris)

    @staticmethod
    def _calculate_metric(x: List[str]) -> float:
        if len(x) < 2:
            return 0.0

        num_changes = 0.0
        for i in range(len(x) - 1):
            if x[i] != x[i + 1]:
                num_changes += 1.0
        return num_changes


class NumSegments(Metric):
    """Calculate the number of aligned repeats in an alignment group regardless of
    orientation.
    """

    def calculate_metric(self, cluster_id: int, all_blocks: List[Dict]) -> float:
        return len([b for b in all_blocks if b["rel_id"] == cluster_id])


class MetricFactory:
    """ """

    def get_rules(self) -> List[Metric]:
        """ """
        return [
            SwitchFractionStrategy(),
            FragmentationIndexStrategy(),
            InterlacingRatioStrategy(),
            OrientationSwitchStrategy(),
            NumOrientationSwitches(),
            NumSegments(),
        ]


class ClusterMetrics:
    """Compute cluster metrics for a given cluster using multiple strategies."""

    def __init__(
        self,
        cluster_id: int,
        all_blocks: List[Dict],
        rule_factory_class=MetricFactory,
    ):
        """
        Args:
            cluster_id: The rel_id of the cluster to evaluate.
            all_blocks: List of all blocks for the read.
            rule_factory_class: Factory class providing fragmentation strategies.
        """
        self.cluster_id = cluster_id
        self.all_blocks = all_blocks
        self.rules = rule_factory_class().get_rules()

    def calculate_metrics(self) -> Dict[str, float]:
        """Apply all fragmentation strategies and return a dict of {strategy_name: metric}."""
        metrics = {}
        for rule in self.rules:
            metric = rule.calculate_metric(self.cluster_id, self.all_blocks)
            metrics[rule.__class__.__name__] = metric
        return metrics
