from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict


@dataclass(frozen=True)
class ReadClassificationConfig:
    max_orientation_error_1D: float = 0.1
    max_orientation_error_2D: float = 0.2
    min_mapped_inserts: int = 3
    backbone_prefix: str = "BB"
    consecutiveness_method: str = "new_occurrence"
    consecutiveness_thresholds: Dict[str, float] = field(
        default_factory=lambda: {
            "new_occurrence": 1.0,
            "switch_fraction": 0.5,
            "fragmentation_index": 0.5,
            "interlacing_ratio": 0.5,
        }
    )
    calibration_model: Path = Path("models/test/fraction")

    @property
    def consecutiveness_threshold(self) -> float:
        """Return the threshold for the currently selected method"""
        return self.consecutiveness_thresholds[self.consecutiveness_method]

    def to_dict(self) -> Dict[str, Any]:
        """Return the dataclass fields as a dictionary"""
        return asdict(self)
