from abc import ABC, abstractmethod

from pysam import AlignedSegment


class Filter(ABC):
    """
    abstract filter class, the filter function should return True if it passes the filtering criteria
    """

    @abstractmethod
    def filter(self, alignment: AlignedSegment) -> bool:
        pass

    @property
    def name(self) -> str:
        return self.__class__.__name__


class FilterMinimunMappingQuality(Filter):
    """
    Check that the mapping quality is above a set value
    """

    def __init__(self, mapping_quality_min: int = 0):
        self.mapping_quality_min = mapping_quality_min

    def filter(self, alignment: AlignedSegment) -> bool:
        status = True
        if alignment.mapping_quality < self.mapping_quality_min:
            status = False
        return status


class FilterSecondaryAlignments(Filter):
    """
    We are not interested in secondary alignments, only in supplementary alignments.
    """

    def filter(self, alignment: AlignedSegment) -> bool:
        # flip the behaviour of is secondary to make it work in the framework.
        status = True
        if alignment.is_secondary:
            status = False
        return status


class FilterMinimunRawReadLength(Filter):
    """
    Check that the mapping quality is above a set value
    """

    def __init__(self, minimum_length: int = 100):
        self.minimum_length = minimum_length

    def filter(self, alignment: AlignedSegment) -> bool:
        status = True
        if alignment.infer_read_length() < self.minimum_length:
            status = False
        return status
