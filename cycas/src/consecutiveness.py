from abc import ABC, abstractmethod
from typing import Dict, List


class ConsecutivenessStrategy(ABC):
    """Abstract base for consecutiveness evaluation strategies."""

    @abstractmethod
    def calculate_metric(self, cluster_id: int, all_blocks: List[Dict]) -> float:
        """Calculate the strategy-specific metric."""
        pass

    @abstractmethod
    def is_consecutive(self, metric: float, threshold: float) -> bool:
        """Assess consecutiveness given a metric and a threshold."""
        pass


class NewOccurrenceStrategy(ConsecutivenessStrategy):
    """
    Checks number of new ID occurrences.
    If an ID stops appearing in the read, and then appears again,
    then it will summed as a new occurrence.
    """

    def calculate_metric(self, cluster_id: int, all_blocks: List[Dict]) -> float:
        sorted_all = sorted(all_blocks, key=lambda b: int(b["read_pos"]))
        num_segments = 0
        in_segment = False
        for block in sorted_all:
            if block["rel_id"] == cluster_id:
                if not in_segment:
                    num_segments += 1
                    in_segment = True
            else:
                in_segment = False
        return float(num_segments)

    def is_consecutive(self, metric: float, threshold: float) -> bool:
        return metric == threshold


class ConsecutivenessStrategyFactory:
    """Factory for retrieving consecutiveness strategies by name."""

    _strategies = {
        "new_occurrence": NewOccurrenceStrategy(),
    }

    @classmethod
    def get_strategy(cls, name: str) -> ConsecutivenessStrategy:
        try:
            return cls._strategies[name]
        except KeyError:
            raise ValueError(f"Unknown strategy: {name}")
