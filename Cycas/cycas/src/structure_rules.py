from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Tuple

from .config import ReadClassificationConfig


class ClassificationRule(ABC):
    """
    Abstract base class for structure classification rules.
    """

    applicable = True  # Can be toggled per subclass

    def __init__(self, config: Optional[dict] = None) -> None:
        """
        Args:
            config: Optional dictionary of rule-specific thresholds or options.
        """
        self.config = config or {}

    @abstractmethod
    def apply(
        self, cluster: tuple[int, List[Dict], bool], metrics: Dict[str, float]
    ) -> Optional[str]:
        """
        Apply classification logic to a cluster of insert blocks.

        Args:
            cluster: Tuple of (relative id, List of insert block dicts, is_consecutive)
            metrics: Dict of calculated metrics {metric: value}.

        Returns:
            A classification label such as "1D", "2D", or "template_switch", or None.
        """
        pass

    @property
    def name(self) -> str:
        """
        Name of the rule, defaults to class name.
        """
        return self.__class__.__name__


class RuleConsecutive(ClassificationRule):
    """
    Classifies clusters as 'consecutive' if blocks appear consecutively in original read.
    """

    def apply(
        self, cluster: Tuple[int, List[Dict], bool], metrics: Dict[str, float]
    ) -> Optional[str]:
        _, _, is_consecutive = cluster
        if not is_consecutive:
            return "non-consecutive"
        return None


class Rule1D(ClassificationRule):
    """
    Classifies clusters as '1D' if all blocks share the same orientation,
    within a configurable orientation error threshold.
    """

    def apply(
        self, cluster: Tuple[int, List[Dict], bool], metrics: Dict[str, float]
    ) -> Optional[str]:
        if metrics["NumOrientationSwitches"] == 0:
            return "1D"

        return None


class Rule2D(ClassificationRule):
    """
    Classifies clusters as '2D' if both forward and reverse reads are present,
    and the orientation alternation rate is above a threshold.
    """

    def apply(
        self, cluster: Tuple[int, List[Dict], bool], metrics: Dict[str, float]
    ) -> Optional[str]:
        if metrics["NumOrientationSwitches"] == 1:
            return "2D"

        if metrics["NumOrientationSwitches"] > 1:
            if metrics["NumOrientationSwitches"] / metrics[
                "NumSegments"
            ] <= self.config.get("max_orientation_error_2D", 0.2):
                return "2D"
            else:
                return None

        return None


class Rule3D(ClassificationRule):
    """
    Classifies clusters as '3D' if both forward and reverse reads are present,
    and the orientation alternation rate is below a threshold.
    """

    def apply(
        self, cluster: Tuple[int, List[Dict], bool], metrics: Dict[str, float]
    ) -> Optional[str]:
        if metrics["NumOrientationSwitches"] > 1:
            if metrics["NumOrientationSwitches"] / metrics[
                "NumSegments"
            ] > self.config.get("max_orientation_error_2D", 0.2):
                return "3D"
            else:
                return None

        return None


class StructureRuleFactory:
    """
    Factory to construct a set of classification rules using a shared config.
    """

    def __init__(self, config: ReadClassificationConfig):
        """
        Args:
            config: Dictionary of configuration thresholds and options
                used by each rule, e.g.:
                {
                    "max_orientation_error_1D": 0.1,
                    "max_orientation_error_2D": 0.2,
                }
        """
        self.config = config

    def get_rules(self) -> List[ClassificationRule]:
        """
        Instantiates and returns all applicable classification rule objects.

        Returns:
            List of ClassificationRule instances.
        """
        return [
            RuleConsecutive(self.config),
            Rule1D(self.config),
            Rule2D(self.config),
            Rule3D(self.config),
        ]
