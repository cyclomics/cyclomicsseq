from abc import ABC, abstractmethod
from copy import deepcopy
from typing import Dict, List
import random

import pysam
from loguru import logger

from .classification_rules import (
    CalculateAlternatingChromosomeRatio,
    CalculateRelativeChromosomalAbundance,
    CalculateContiguousStartPositions,
    CalculateTripletChromosomes,
    Check1BackBoneInAlignments,
    CheckEqualStartLocations,
    CheckMaximumCountAlignments,
    CheckMinimumCountAlignments,
    CheckNChromosomesPresent,
    CheckNDirectionalChromosomesPresent,
    CheckTwoStartLocations,
    ClassificationRule,
    CheckAllBackBoneInAlignments,
    Check1DirectionalChromosomesPresent,
    Check2DirectionalChromosomesPresent,
    Check3DirectionalChromosomesPresent,
)

# used for type hints
MinimumScores = Dict[ClassificationRule, int]


class BaseClassifier(ABC):
    @abstractmethod
    def __init__(self, priority: int, expected_minimum_results: MinimumScores):
        """
        A classifier is a collection of rules that determine the likelihood that an AlignmentGroup is of a certain read_type.
        To implement a Classifier all that is required is the property of self.rules set to several initiated ClassificationRules,
        as well as a set of minimum expected results set in self.expected_minimum_results, with the name as the key and the minimum value as the value.
        This is usually implemented by setting a dict in the child class valled expected_minimum_results, and calling the parent init class afterwards.

        The parent (this class) function of classify can than be used as a standard interface.

        priority: int
        rules: List[ClassificationRule]
        """
        self.priority = priority
        self.rules = expected_minimum_results.keys()
        self.expected_minimum_results = {}

        for k, v in expected_minimum_results.items():
            if k.name in self.expected_minimum_results.keys():
                error_text = (
                    f"{self.name} tries to apply a rule to itself twice! rule: {k.name}"
                )
                logger.critical(error_text)
                raise NotImplementedError(error_text)
            # add to the dict
            self.expected_minimum_results[k.name] = v

    @property
    def name(self) -> str:
        return self.__class__.__name__

    @property
    def color(self) -> str:
        random_number = random.randint(0, 16777215)
        hex_number = str(hex(random_number))
        hex_number = "#" + hex_number[2:]
        return hex_number

    def apply_rules(self, group):
        outcome = {}
        for rule in self.rules:
            outcome[rule.name] = rule.score(group)
        return outcome

    def manage_expectations(self, outcome):
        expectation_deficit = 0

        for rule, score in outcome.items():
            expected = self.expected_minimum_results[rule]
            # if we get less than expected, calculate the deficit
            if score < expected:
                logger.debug(
                    f"Rule {rule} failed to meet expectation of {expected} at {score}"
                )
                expectation_deficit += expected - score
        logger.debug(f"deficit is {expectation_deficit}")

        return expectation_deficit

    def classify(self, group) -> str:
        logger.debug(f"Classifying using {self.__class__.__name__}")

        outcome = self.apply_rules(group)
        class_score = self.manage_expectations(outcome)
        if class_score == 0:
            logger.debug("applying classification!")
            groups = self.create_paired_alignment_groups(group)

            return (self.name, groups, True, outcome)
        else:
            return (self.name, None, False, outcome)

    def create_paired_alignment_groups_per_chromosome(self, group):
        """
        Create new groups that only contain alignments against a single chromosome.
        TODO: make work without deepcopy as this takes 70% of the runtime
        """
        output = []
        for directional_chromosome in group.alignment_chromosomes_present():
            group_chr = deepcopy(group)
            group_chr.extract_alignments_by_directional_chromosome(
                directional_chromosome
            )
            output.append(group_chr)
        return output

    def _get_spaced_integers(self, integers, spacer):
        """
        given a list of integer return all integers that are at least spacer away from the last one
        so 1,2,3,6,7,8 with a spacer of 2 wil return 1,6

        """
        integers.sort()
        last = integers[0]
        result = [integers[0]]

        for i in integers:
            if i > last + spacer:
                result.append(i)
                last = i
        return result

    def create_paired_alignment_groups_per_chromosome_and_distance(
        self, group, distance=5000
    ):
        output = []

        for directional_chromosome in group.alignment_chromosomes_present():
            group_chr = deepcopy(group)
            group_chr.extract_alignments_by_directional_chromosome(
                directional_chromosome
            )
            start_positions = [x.alignment_chromosome_start for x in group.alignments]

            # multi occurence within chromosome
            if max(start_positions) - distance > min(start_positions):
                logger.debug("multi insert found in single chromosome ")
                spaced_start_positions = self._get_spaced_integers(
                    start_positions, distance
                )
                for i in spaced_start_positions:
                    current_group = deepcopy(group_chr)

                    current_group.extract_alignments_by_start_position(
                        i, margin=distance
                    )
                    output.append(current_group)
            else:
                # single group
                output.append(group_chr)
        return output

    @staticmethod
    def calculate_baseunit_copies(alns, baseunit_length):
        lengths = [x.alignment_length for x in alns]
        # calculate the relactive length and sum the result, ignore 0's
        return sum(
            [x / baseunit_length for x in lengths if x > 0 and baseunit_length > 0]
        )

    def create_metadata(self, group, consensus):
        """
        Basic function to create the metadata, specifics should, and could be altered by classifiers.
        """
        first_aln = group.alignments[0]

        metadata = {}
        # required info
        metadata[
            "id"
        ] = f"{first_aln.readname}_{first_aln.alignment_chromosome}_{first_aln.alignment_chromosome_start}"
        metadata["raw_length"] = first_aln.raw_read_length
        metadata["baseunit_copies"] = self.calculate_baseunit_copies(
            group.alignments, len(consensus)
        )
        metadata["baseunit_start_idx"] = first_aln.first_cigar_value
        metadata["baseunit_end_idx"] = group.alignments[-1].first_cigar_value
        metadata["baseunit_length"] = len(consensus)
        metadata["baseunit_certainty"] = 0
        metadata["baseunit_orientation"] = first_aln.alignment_direction
        metadata["baseunit_idx"] = ",".join(
            [str(x.first_cigar_value) for x in group.alignments]
        )
        # aditional info
        metadata["chromosome"] = first_aln.alignment_chromosome
        metadata["classification"] = self.name
        return metadata


class SingleInsert(BaseClassifier):
    """Find reads that consist of only an insert that start contiguously."""

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")

        expected_minimum_results = {
            CheckEqualStartLocations(): 1,
            CalculateRelativeChromosomalAbundance(): 0.9,
            CalculateContiguousStartPositions(): 0.7,
            CheckMinimumCountAlignments(): 1,
        }
        super().__init__(priority, expected_minimum_results)

    @property
    def color(self):
        return "#21618c"

    def create_paired_alignment_groups(self, group):
        return self.create_paired_alignment_groups_per_chromosome_and_distance(group)

    def create_metadata(self, group, consensus):
        metadata = super().create_metadata(group, consensus)
        return metadata


class SingleInsertGaps(BaseClassifier):
    """Find reads that only contain an insert, but there are gaps in between the alignments."""

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")

        expected_minimum_results = {
            CheckEqualStartLocations(): 1,
            CalculateRelativeChromosomalAbundance(): 0.9,
            # CalculateContiguousStartPositions(): 0.2,
            CheckMinimumCountAlignments(): 1,
        }
        super().__init__(priority, expected_minimum_results)

    @property
    def color(self):
        return "#2980b9"

    def create_paired_alignment_groups(self, group):
        return self.create_paired_alignment_groups_per_chromosome_and_distance(group)

    def create_metadata(self, group, consensus):
        metadata = super().create_metadata(group, consensus)
        return metadata


class SingleBackbone(BaseClassifier):
    """Find reads that only contain a backbone and that start contiguously."""

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")

        expected_minimum_results = {
            CheckEqualStartLocations(): 1,
            CalculateRelativeChromosomalAbundance(): 0.9,
            CalculateContiguousStartPositions(): 0.7,
            CheckMinimumCountAlignments(): 1,
            Check1BackBoneInAlignments(): 1,
        }
        super().__init__(priority, expected_minimum_results)

    @property
    def color(self):
        return "#dc7633"

    def create_paired_alignment_groups(self, group):
        return self.create_paired_alignment_groups_per_chromosome_and_distance(group)

    def create_metadata(self, group, consensus):
        metadata = super().create_metadata(group, consensus)
        return metadata


class DoubleBackbone(BaseClassifier):
    """Find reads that only contain a backbone and that start contiguously."""

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")

        expected_minimum_results = {
            CheckMinimumCountAlignments(): 1,
            CheckAllBackBoneInAlignments(): 1,
            CheckNDirectionalChromosomesPresent(directional_chromosomes=2): 1,
        }
        super().__init__(priority, expected_minimum_results)

    @property
    def color(self):
        return "#d35400"

    def create_paired_alignment_groups(self, group):
        return self.create_paired_alignment_groups_per_chromosome_and_distance(group)

    def create_metadata(self, group, consensus):
        metadata = super().create_metadata(group, consensus)
        return metadata


class TripleBackbone(BaseClassifier):
    """Find reads that only contain a backbone and that start contiguously."""

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")

        expected_minimum_results = {
            CheckMinimumCountAlignments(): 1,
            CheckAllBackBoneInAlignments(): 1,
            CheckNDirectionalChromosomesPresent(directional_chromosomes=3): 1,
        }
        super().__init__(priority, expected_minimum_results)

    @property
    def color(self):
        return " #ba4a00"

    def create_paired_alignment_groups(self, group):
        return self.create_paired_alignment_groups_per_chromosome_and_distance(group)

    def create_metadata(self, group, consensus):
        metadata = super().create_metadata(group, consensus)
        return metadata


class QuadtripleBackbone(BaseClassifier):
    """Find reads that only contain a backbone and that start contiguously."""

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")

        expected_minimum_results = {
            CheckEqualStartLocations(): 1,
            CalculateRelativeChromosomalAbundance(): 0.9,
            CalculateContiguousStartPositions(): 0.7,
            CheckMinimumCountAlignments(): 1,
            CheckAllBackBoneInAlignments(): 1,
            CheckNDirectionalChromosomesPresent(directional_chromosomes=3): 1,
        }
        super().__init__(priority, expected_minimum_results)

    @property
    def color(self):
        return " #a04000"

    def create_paired_alignment_groups(self, group):
        return self.create_paired_alignment_groups_per_chromosome_and_distance(group)

    def create_metadata(self, group, consensus):
        metadata = super().create_metadata(group, consensus)
        return metadata


class SingleBackboneGaps(BaseClassifier):
    """Find reads that only consist of backbone alignments, but there are gaps in between the alignments."""

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")

        expected_minimum_results = {
            CheckEqualStartLocations(): 1,
            CalculateRelativeChromosomalAbundance(): 0.9,
            CalculateContiguousStartPositions(): 0.2,
            CheckMinimumCountAlignments(): 1,
            Check1BackBoneInAlignments(): 1,
        }
        super().__init__(priority, expected_minimum_results)

    @property
    def color(self):
        return "#e59866"

    def create_paired_alignment_groups(self, group):
        return self.create_paired_alignment_groups_per_chromosome_and_distance(group)

    def create_metadata(self, group, consensus):
        metadata = super().create_metadata(group, consensus)
        return metadata


class TwoInsertsChromosomal(BaseClassifier):
    """
    There are two ways we can know that there are two inserts: they map to different chromosomes or they map on two separate locations.
    This class checks the first case where they are on two different chromosomes.
    """

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")

        expected_minimum_results = {
            CheckMinimumCountAlignments(): 1,
            CheckEqualStartLocations(): 1,
            CalculateContiguousStartPositions(): 0.7,
            CheckNDirectionalChromosomesPresent(directional_chromosomes=2): 1,
        }
        # raise NotImplementedError
        super().__init__(priority, expected_minimum_results)

    @property
    def name(self):
        return "DoubleInsert"

    @property
    def color(self):
        # Note that this is defined twice!
        return "#85c1e9"

    def create_paired_alignment_groups(self, group):
        return self.create_paired_alignment_groups_per_chromosome_and_distance(group)

    def create_metadata(self, group, consensus):
        metadata = super().create_metadata(group, consensus)
        return metadata


class TwoInsertsChromosomalGaps(BaseClassifier):
    """
    There are two ways we can know that there are two inserts: they map to different chromosomes or they map on two separate locations.
    This class checks the first case where they are on two different chromosomes. Allows for gaps in between the alignments.
    """

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")

        expected_minimum_results = {
            CheckMinimumCountAlignments(): 1,
            CheckEqualStartLocations(): 1,
            CalculateAlternatingChromosomeRatio(): 0.2,
            CheckNDirectionalChromosomesPresent(directional_chromosomes=2): 1,
        }
        super().__init__(priority, expected_minimum_results)

    @property
    def name(self):
        return "DoubleInsertGaps"

    @property
    def color(self):
        return "#2980b9"

    def create_paired_alignment_groups(self, group):
        return self.create_paired_alignment_groups_per_chromosome_and_distance(group)

    def create_metadata(self, group, consensus):
        metadata = super().create_metadata(group, consensus)
        return metadata


class TwoInsertsStartLocations(BaseClassifier):
    """
    There are two ways we can know that there are two inserts: they map to different chromosomes or they map on two separate locations.
    This class checks the start positions.
    """

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")
        # group start locations
        expected_minimum_results = {
            CheckMinimumCountAlignments(): 1,
            CalculateContiguousStartPositions(): 0.7,
            CheckTwoStartLocations(): 1,
        }
        super().__init__(priority, expected_minimum_results)

    @property
    def name(self):
        return "DoubleInsert"

    @property
    def color(self):
        # Note that this is defined twice!
        return "#85c1e9"

    def create_paired_alignment_groups(self, group):
        return self.create_paired_alignment_groups_per_chromosome_and_distance(group)

    def create_metadata(self, group, consensus):
        metadata = super().create_metadata(group, consensus)
        return metadata


class BackboneInsert(BaseClassifier):
    """Find reads that consist of a backbone and an insert that alternate, without major non mapping regions"""

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")

        expected_minimum_results = {
            CheckMinimumCountAlignments(): 1,
            CheckEqualStartLocations(): 1,
            CalculateContiguousStartPositions(): 0.7,
            CalculateAlternatingChromosomeRatio(): 0.7,
            Check1BackBoneInAlignments(): 1,
            CheckNDirectionalChromosomesPresent(directional_chromosomes=2): 1,
        }
        super().__init__(priority, expected_minimum_results)

    @property
    def color(self):
        return "#355E3B"

    def create_paired_alignment_groups(self, group):
        return self.create_paired_alignment_groups_per_chromosome_and_distance(group)

    def create_metadata(self, group, consensus):
        metadata = super().create_metadata(group, consensus)
        return metadata


class BackboneInsertGaps(BaseClassifier):
    """
    Find reads that consist of a backbone and an insert.
    Less strict on the alternation pattern since missing mappings can cause this metric to drop.
    """

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")

        expected_minimum_results = {
            CheckMinimumCountAlignments(): 1,
            CheckEqualStartLocations(): 1,
            Check1BackBoneInAlignments(): 1,
            CalculateAlternatingChromosomeRatio(): 0.2,
            CheckNDirectionalChromosomesPresent(directional_chromosomes=2): 1,
        }
        super().__init__(priority, expected_minimum_results)

    @property
    def color(self):
        return "#4F7942"

    def create_paired_alignment_groups(self, group):
        return self.create_paired_alignment_groups_per_chromosome_and_distance(group)

    def create_metadata(self, group, consensus):
        metadata = super().create_metadata(group, consensus)
        return metadata


class ShortContiguous(BaseClassifier):
    """Find reads that are short, but have relatively consistant alignments throughout the read"""

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")

        expected_minimum_results = {
            CheckMaximumCountAlignments(): 1,
            CheckEqualStartLocations(): 1,
            CalculateContiguousStartPositions(): 0.7,
        }
        super().__init__(priority, expected_minimum_results)

    @property
    def color(self):
        return "#bb8fce"

    def create_paired_alignment_groups(self, group):
        return self.create_paired_alignment_groups_per_chromosome_and_distance(group)

    def create_metadata(self, group, consensus):
        metadata = super().create_metadata(group, consensus)
        return metadata

    @property
    def name(self):
        return "LowAlignmentCount"


class ShortGaps(BaseClassifier):
    """Find reads that are short."""

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")

        expected_minimum_results = {
            CheckMaximumCountAlignments(): 1,
            # CheckEqualStartLocations(): 1,
            # CalculateContiguousStartPositions(): 0.2,
        }
        super().__init__(priority, expected_minimum_results)

    @property
    def color(self):
        return "#9b59b6"

    def create_paired_alignment_groups(self, group):
        return self.create_paired_alignment_groups_per_chromosome_and_distance(group)

    def create_metadata(self, group, consensus):
        metadata = super().create_metadata(group, consensus)
        return metadata

    @property
    def name(self):
        return "LowAlignmentCountGaps"


class BackboneDoubleInsert(BaseClassifier):
    """Find reads that have a double insert concatermerized with a single backbone"""

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")

        expected_minimum_results = {
            CheckMinimumCountAlignments(): 1,
            CheckEqualStartLocations(): 1,
            CalculateContiguousStartPositions(): 0.7,
            Check1BackBoneInAlignments(): 1,
            CalculateTripletChromosomes(): 0.7,
        }
        super().__init__(priority, expected_minimum_results)

    @property
    def color(self):
        return "#228B22"

    def create_paired_alignment_groups(self, group):
        # TODO: split these out by chromosome if they are both the same orientation
        # raise NotImplementedError
        return self.create_paired_alignment_groups_per_chromosome_and_distance(group)

    def create_metadata(self, group, consensus):
        metadata = super().create_metadata(group, consensus)
        return metadata


class TripleInsert(BaseClassifier):
    """Find reads that have 3 inserts, statistically very rare"""

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")

        expected_minimum_results = {
            CheckMinimumCountAlignments(): 1,
            CheckEqualStartLocations(): 1,
            CalculateContiguousStartPositions(): 0.7,
            CheckNDirectionalChromosomesPresent(directional_chromosomes=3): 1,
        }
        super().__init__(priority, expected_minimum_results)

    @property
    def color(self):
        return "#3498db"

    def create_paired_alignment_groups(self, group):
        # TODO: Pick the longest strech that looks good.
        return self.create_paired_alignment_groups_per_chromosome_and_distance(group)

    def create_metadata(self, group, consensus):
        metadata = super().create_metadata(group, consensus)
        return metadata


class DirectionalFlip(BaseClassifier):
    """Find reads that flip their directionality patern trougout the read"""

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")

        #  check if we have 2 chromosomes but 4 directional chromosomes
        expected_minimum_results = {
            CheckMinimumCountAlignments(): 1,
            CheckNDirectionalChromosomesPresent(directional_chromosomes=4): 1,
            CheckNChromosomesPresent(chromosome_count=2): 1,
            CalculateAlternatingChromosomeRatio(): 0.6,
        }
        super().__init__(priority, expected_minimum_results)

    @property
    def color(self):
        return "#98FB98"

    def create_paired_alignment_groups(self, group):
        # TODO: Pick the longest strech that looks good.
        return self.create_paired_alignment_groups_per_chromosome_and_distance(group)

    def create_metadata(self, group, consensus):
        metadata = super().create_metadata(group, consensus)
        return metadata


class SingleDirectionalFlip(BaseClassifier):
    """Find reads that flip their directionality patern trougout the read"""

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")

        #  check if we have 2 chromosomes but 4 directional chromosomes
        expected_minimum_results = {
            CheckMinimumCountAlignments(): 1,
            CheckNDirectionalChromosomesPresent(directional_chromosomes=2): 1,
            CheckNChromosomesPresent(chromosome_count=1): 1,
            CalculateContiguousStartPositions(): 0.9,
        }
        super().__init__(priority, expected_minimum_results)

    @property
    def color(self):
        return " #e59866"

    def create_paired_alignment_groups(self, group):
        # TODO: Pick the longest strech that looks good.
        return self.create_paired_alignment_groups_per_chromosome_and_distance(group)

    def create_metadata(self, group, consensus):
        metadata = super().create_metadata(group, consensus)
        return metadata


class SingletonBackbone(BaseClassifier):
    """Find reads that only have a single aligned object and that are fairly short."""

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")

        expected_minimum_results = {
            CheckMinimumCountAlignments(1): 1,
            CheckMaximumCountAlignments(1): 1,
            CheckNDirectionalChromosomesPresent(directional_chromosomes=1): 1,
            Check1BackBoneInAlignments(): 1,
        }
        super().__init__(priority, expected_minimum_results)

    @property
    def color(self):
        return "#f6ddcc"

    def create_paired_alignment_groups(self, group):
        return self.create_paired_alignment_groups_per_chromosome_and_distance(group)


class Singleton(BaseClassifier):
    """Find reads that only have a single aligned object and that are fairly short."""

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")

        expected_minimum_results = {
            CheckMinimumCountAlignments(1): 1,
            CheckMaximumCountAlignments(1): 1,
            CheckNDirectionalChromosomesPresent(directional_chromosomes=1): 1,
        }
        super().__init__(priority, expected_minimum_results)

    @property
    def color(self):
        return "#d6eaf8"

    def create_paired_alignment_groups(self, group):
        return self.create_paired_alignment_groups_per_chromosome_and_distance(group)


class ComplexConcatemer(BaseClassifier):
    """
    All reads containing more than 5 mapping directional chromosomes that are not short.
    """

    def __init__(self, priority: int) -> None:
        logger.debug(f"Classifier initiated: {self.__class__.__name__}")

        expected_minimum_results = {
            CheckMinimumCountAlignments(): 1,
            Check1DirectionalChromosomesPresent(): 1,
            Check2DirectionalChromosomesPresent(): 1,
            Check3DirectionalChromosomesPresent(): 1,
        }
        # raise NotImplementedError
        super().__init__(priority, expected_minimum_results)

    @property
    def color(self):
        return "#6c3483"

    def create_paired_alignment_groups(self, group):
        return []

    def create_metadata(self, group, consensus):
        first_aln = group.alignments[0]

        metadata = {}
        # required info
        metadata["id"] = f"{first_aln.readname}_complex"
        metadata["raw_length"] = first_aln.raw_read_length
        metadata["baseunit_copies"] = 0
        metadata["baseunit_start_idx"] = 0
        metadata["baseunit_end_idx"] = 0
        metadata["baseunit_length"] = 0
        metadata["baseunit_certainty"] = 0
        metadata["baseunit_orientation"] = "="
        metadata["baseunit_idx"] = "0"
        # aditional info
        metadata["chromosome"] = first_aln.alignment_chromosome
        metadata["classification"] = self.name
        return metadata
