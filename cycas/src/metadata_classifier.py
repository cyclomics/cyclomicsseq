from abc import ABC, abstractmethod
from collections import Counter
from email.mime import application
from importlib.metadata import metadata
from typing import Dict

from loguru import logger

from src.alignment_processor import AlignmentProcessor


METADATA = Dict[str, int]
UNKOWN_CLASS = "Unkown"


class BaseMetadataClass(ABC):
    def __init__(self) -> None:
        self.alignment_processor = AlignmentProcessor()

    @property
    @abstractmethod
    def priority(self) -> bool:
        pass

    @abstractmethod
    def inspect(self) -> bool:
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
    def name(self):
        return self.__class__.__name__


class BackboneInsert(BaseMetadataClass):
    applicable = True
    priority = 1

    def inspect(self, metadata: METADATA) -> str:
        outcome = False
        if not metadata["Check1BackBoneInAlignments"]:
            return outcome
        if (
            metadata["Check2ChromosomesPresent"]
            and not metadata["CheckTwoStartLocationsByOwnLength"]
        ):
            return self.__class__.__name__

        return outcome

    def chop_blocks(self, alignments):
        return self.alignment_processor.create_paired_alignment_groups_per_chromosome_and_distance(
            alignments
        )


class BackboneDoubleInsert(BaseMetadataClass):
    applicable = True
    priority = 1

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
            return self.__class__.__name__
        if metadata["Check3ChromosomesPresent"]:
            return self.__class__.__name__

        return outcome

    def chop_blocks(self, alignments):
        return self.alignment_processor.create_paired_alignment_groups_per_chromosome_and_distance(
            alignments
        )


class SingleBackbone(BaseMetadataClass):
    applicable = True
    priority = 10

    def inspect(self, metadata: METADATA) -> str:
        outcome = False
        if not metadata["Check1BackBoneInAlignments"]:
            return outcome
        if metadata["Check1ChromosomesPresent"]:
            return self.__class__.__name__
        return outcome

    def chop_blocks(self, alignments):
        return self.alignment_processor.create_paired_alignment_groups_per_chromosome_and_distance(
            alignments
        )


class SingleInsert(BaseMetadataClass):
    applicable = True
    priority = 20

    def inspect(self, metadata: METADATA) -> str:
        outcome = False
        if not metadata["Check1ChromosomesPresent"]:
            return outcome
        if ["CheckTwoStartLocationsByOwnLength"]:
            return outcome
        if metadata["CheckEqualStartLocations"] > 0.8:
            return self.__class__.__name__

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
            return self.__class__.__name__
        if metadata["Check3ChromosomesPresent"]:
            return self.__class__.__name__

        return False

    def chop_blocks(self, alignments):
        return []


class Unkown(BaseMetadataClass):
    applicable = True
    priority = 9999

    def inspect(self, metadata: METADATA) -> bool:
        return UNKOWN_CLASS

    def chop_blocks(self, alignments):
        return []


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
            raise RuntimeError(
                f"Multiple classifiers ({len(classifier)}) found with the same name. Aborting..."
            )

        classifier = classifier[0]

        return classifier.chop_blocks(alignments)

    def show_counts(self):
        return Counter([v for v in self.result.values()])
