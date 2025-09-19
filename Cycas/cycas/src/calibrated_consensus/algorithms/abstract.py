from abc import abstractmethod
from pathlib import Path
from typing import Any, Dict, Tuple


class AbstractCalibratedConsensus:
    """
    Abstract class for generation and calibration of consensus sequences
    """

    def __init__(self):
        pass

    @abstractmethod
    def prep_file(self, *args, **kwargs):
        raise NotImplementedError

    @abstractmethod
    def fit_model(self, *args, **kwargs):
        raise NotImplementedError

    @abstractmethod
    def make_consensus(self, *args, **kwargs) -> Dict[str, Any]:
        raise NotImplementedError

    @abstractmethod
    def calibrate(self, *args, **kwargs) -> Tuple[str, str]:
        raise NotImplementedError

    @abstractmethod
    def save(self, output_path: Path):
        raise NotImplementedError

    @abstractmethod
    def load(self, input_path: Path):
        raise NotImplementedError
