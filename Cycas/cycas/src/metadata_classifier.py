from abc import ABC, abstractmethod
from collections import Counter
from typing import Dict, Union
import random

from loguru import logger

from src.alignment_processor import AlignmentProcessor


METADATA = Dict[str, int]
UNKNOWN_CLASS_LABEL = "Unknown"


class BaseMetadataClass(ABC):
    def __init__(self) -> None:
        # Give all childern access to the default processor.
        self.alignment_processor = AlignmentProcessor()

    @property
    @abstractmethod
    def priority(self) -> int:
        """
        If two classifiers are both true, the priority defines the classification. Lower is better
        """
        pass

    @abstractmethod
    def inspect(self) -> Union[bool, str]:
        """
        return false if not classified as this type, return name if classified correctly
        """
        pass

    @abstractmethod
    def chop_blocks(self) -> bool:
        """
        Create a function that creates multiple blocks of alignment objects you want to call consensus on together.
        """
        pass

    def __lt__(self, other):
        """Used to sort classes by their priority"""
        return self.priority < other.priority

    def between(self, score, n1, n2) -> bool:
        """
        Check if score is larger than n1 and smaller than n2.
        """
        if score >= n1 and score <= n2:
            return True
        else:
            return False

    @property
    def color(self) -> str:
        random_number = random.randint(0, 16777215)
        hex_number = str(hex(random_number))
        hex_number = "#" + hex_number[2:]
        return hex_number

    @property
    def name(self):
        return self.__class__.__name__


class MetadataClassifier:
    def __init__(self, classifier_base_class=BaseMetadataClass) -> None:
        # TODO get and order all classifiers
        self.classifier_base_class = classifier_base_class
        self.classification_options = self.get_classes_available()
        self.result = {}

    def all_subclasses(self, cls):
        """
        Recursive function to get nested inheritance classes, usefull for when rules inherrit each other
        """
        return set(cls.__subclasses__()).union(
            [s for c in cls.__subclasses__() for s in self.all_subclasses(c)]
        )

    def get_classes_available(self):
        """Get all rules that are applicable."""
        rules = self.all_subclasses(self.classifier_base_class)
        rules = [x() for x in rules if x.applicable]
        # priority is used in sorting
        rules.sort()
        return rules

    def classify(self, metadata):
        for i in self.classification_options:
            outcome = i.inspect(metadata["rules"])
            # Report first passing class
            if outcome:
                self.result[metadata["readname"]] = outcome
                break

        return outcome

    def separate_alignments(self, classification, alignments):
        classifier = [
            x for x in self.classification_options if x.name == classification
        ]

        if not len(classifier) == 1:
            if len(classifier) > 1:
                raise RuntimeError(
                    f"Multiple classifiers ({len(classifier)}) found with the same name. Aborting..."
                )
            else:
                # empty list:
                raise RuntimeError(
                    f"No classifiers found for name ({classification}). Aborting..."
                )

        classifier = classifier[0]

        return classifier.chop_blocks(alignments)

    def show_counts(self):
        return Counter([v for v in self.result.values()])


class BackboneInsert(BaseMetadataClass):
    applicable = True
    priority = 1
    color = "#355E3B"

    def inspect(self, metadata: METADATA) -> str:
        outcome = False
        conditions = [
            metadata["Check1BackBoneInAlignments"] == 1,
            metadata["Check2ChromosomesPresent"] == 1
            and metadata["CheckTwoStartLocationsByOwnLength"] == 0,
        ]
        if all(conditions):
            return self.name
        else:
            return outcome

        return outcome

    def chop_blocks(self, alignments):
        return self.alignment_processor.create_paired_alignment_groups_per_chromosome_and_distance(
            alignments
        )


class BackboneDoubleInsert(BaseMetadataClass):
    applicable = True
    priority = 2
    color = "#228B22"

    def inspect(self, metadata: METADATA) -> str:
        outcome = False
        # No or 2 backbones is a fail
        if not metadata["Check1BackBoneInAlignments"]:
            return outcome
        # one backbone and two starting position
        if (
            metadata["CheckTwoStartLocationsByOwnLength"]
            and metadata["Check2ChromosomesPresent"]
        ):
            return self.name
        if metadata["Check3ChromosomesPresent"]:
            return self.name

        return outcome

    def chop_blocks(self, alignments):
        return self.alignment_processor.create_paired_alignment_groups_per_chromosome_and_distance(
            alignments
        )


class SingleBackboneUnalignedGaps(BaseMetadataClass):
    applicable = True
    priority = 9

    def inspect(self, metadata: METADATA) -> str:
        outcome = False
        conditions = [
            metadata["Check1BackBoneInAlignments"] == 1,
            metadata["CheckMinimumUnalignedGap"] == 1,
            metadata["CalculateAlternatingChromosomeRatio"] < 0.5,
            metadata["Check1ChromosomesPresent"] == 0,
        ]

        if all(conditions):
            return self.name
        else:
            return outcome

    def chop_blocks(self, alignments):
        return self.alignment_processor.create_paired_alignment_groups_per_chromosome_and_distance(
            alignments
        )


class SingleBackbone(BaseMetadataClass):
    applicable = True
    priority = 10
    color = "#dc7633"

    def inspect(self, metadata: METADATA) -> str:
        outcome = False
        if not metadata["Check1BackBoneInAlignments"]:
            return outcome
        if metadata["Check1ChromosomesPresent"]:
            return self.name
        return outcome

    def chop_blocks(self, alignments):
        return self.alignment_processor.create_paired_alignment_groups_per_chromosome_and_distance(
            alignments
        )


class LowAlignmentCount(BaseMetadataClass):
    applicable = True
    priority = 11

    def inspect(self, metadata: METADATA) -> str:
        outcome = False
        conditions = [metadata["CalculateSegmentCountPercentageTenSegments"] < 0.21]
        if all(conditions):
            return self.name
        else:
            return outcome

    def chop_blocks(self, alignments):
        return self.alignment_processor.create_paired_alignment_groups_per_chromosome_and_distance(
            alignments
        )


class SingleInsertUnalignedGaps(BaseMetadataClass):
    applicable = True
    priority = 19

    def inspect(self, metadata: METADATA) -> str:
        outcome = False
        conditions = [
            metadata["Check1ChromosomesPresent"] == 1,
            metadata["CheckTwoStartLocationsByOwnLength"] == 0,
            metadata["CalculateAlternatingChromosomeRatio"] < 0.5,
            metadata["CheckEqualStartLocations"] > 0.8,
            metadata["CheckMinimumUnalignedGap"] == 0,
        ]

        if all(conditions):
            return self.name
        else:
            return outcome

    def chop_blocks(self, alignments):
        return self.alignment_processor.create_paired_alignment_groups_per_chromosome_and_distance(
            alignments
        )


class SingleInsert(BaseMetadataClass):
    applicable = True
    priority = 20
    color = "#21618c"

    def inspect(self, metadata: METADATA) -> str:
        outcome = False
        if not metadata["Check1ChromosomesPresent"]:
            return outcome
        if metadata["CheckTwoStartLocationsByOwnLength"]:
            return outcome
        if metadata["CheckEqualStartLocations"] > 0.8:
            return self.name

        return outcome

    def chop_blocks(self, alignments):
        return self.alignment_processor.create_paired_alignment_groups_per_chromosome_and_distance(
            alignments
        )


class SingleInsertUncertain(BaseMetadataClass):
    applicable = True
    priority = 24

    def inspect(self, metadata: METADATA) -> str:
        outcome = False
        if metadata["CheckTwoStartLocationsByOwnLength"]:
            return outcome
        if metadata["Check1ChromosomesPresent"]:
            return self.name

        return outcome

    def chop_blocks(self, alignments):
        return self.alignment_processor.create_paired_alignment_groups_per_chromosome_and_distance(
            alignments
        )


class DoubleInsertUncertain(BaseMetadataClass):
    applicable = True
    priority = 25

    def inspect(self, metadata: METADATA) -> str:
        outcome = False
        if (
            not metadata["Check1BackBoneInAlignments"]
            and metadata["Check2ChromosomesPresent"]
        ):
            return self.name

        return outcome

    def chop_blocks(self, alignments):
        return self.alignment_processor.create_paired_alignment_groups_per_chromosome_and_distance(
            alignments
        )


class MessyAlignment(BaseMetadataClass):
    applicable = True
    priority = 1000

    def inspect(self, metadata: METADATA) -> bool:
        if metadata["Check4ChromosomesPresent"]:
            return self.name
        if metadata["Check3ChromosomesPresent"]:
            return self.name

        return False

    def chop_blocks(self, alignments):
        return []


class Unknown(BaseMetadataClass):
    name = UNKNOWN_CLASS_LABEL
    applicable = True
    priority = 9999
    color = "#D3D3D3"

    def inspect(self, metadata: METADATA) -> bool:
        return UNKNOWN_CLASS_LABEL

    def chop_blocks(self, alignments):
        return []
