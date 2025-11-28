from typing import Dict, List

from .config import ReadClassificationConfig
from .structure_rules import StructureRuleFactory


class StructureClassifier:
    """
    Classifies a read structure string based on insert blocks using a set of classification rules.

    A structure string is parsed into insert blocks, grouped by their relative ID, and
    evaluated against multiple classification rules.
    """

    def __init__(
        self,
        config: ReadClassificationConfig,
        rule_factory_class=StructureRuleFactory,
    ):
        """
        Args:
            rule_factory_class: A factory class that returns rules via `.get_rules()`.
            config: Dictionary with necessary classification thresholds.
        """
        self.config = config or ReadClassificationConfig()
        self.rules = rule_factory_class(self.config.to_dict()).get_rules()

    def classify(
        self, cluster: tuple[int, List[Dict], bool], metrics: Dict[str, float]
    ) -> str:
        """
        Classify a cluster based on defined rules and calculated metrics.

        Args:
            cluster: Tuple of (relative id, List of insert block dicts, is_consecutive).
            metrics: Dict of calculated metrics {metric: value}.

        Returns:
            A classification label (e.g., "1D", "2D", "3D", etc.), or fallback labels
            such as "too_few_inserts" or "no_cluster_classification". Can be compound by
            joining with "+" (e.g. "3D_+_non-consecutive").
        """
        labels = set()

        if len(cluster[1]) < self.config.min_mapped_inserts:
            return "too_few_inserts"

        for rule in self.rules:
            label = rule.apply(cluster, metrics)
            if label:
                labels.add(label)

        if not any(labels):
            labels.add("no_cluster_classification")

        return self.resolve_mixed_labels(labels)

    @staticmethod
    def resolve_mixed_labels(labels: set) -> str:
        """
        Resolve multiple conflicting classification labels into a single unified label.

        This handles combinations of 1D/2D structures, as well as fallback classifications
        for chimeric or unknown outcomes.

        Args:
            labels: A set of classification labels resulting from rule application.

        Returns:
            A unified string label representing the final classification.
        """

        return "_+_".join(sorted(labels))
