import statistics
from abc import ABC, abstractmethod
from collections import Counter
from functools import lru_cache
from unittest import result

from loguru import logger


class RuleFactory:
    """
    Class to create all `ClassificationRule`'s that can be applied to Alignment groups
    """

    def all_subclasses(self, cls):
        """
        Recursive function to get nested inheritance classes, usefull for when rules inherrit each other
        """
        return set(cls.__subclasses__()).union(
            [s for c in cls.__subclasses__() for s in self.all_subclasses(c)]
        )

    def get_rules(self):
        """Get all rules that are applicable."""
        rules = [x() for x in self.all_subclasses(ClassificationRule) if x.applicable]
        return rules


class ClassificationRule(ABC):
    """
    rules used for classification,
    Checks are 1 or 0
    Calculates are a range between 1 an 0

    lru cache is used for speed since every classification calls rules.
    Rules that are applicable will the used executed by the factory.
    """

    @abstractmethod
    def score(self, group) -> int:
        pass

    @property
    def name(self) -> str:
        return self.__class__.__name__

    def __repr__(self) -> str:
        return f"instance of Rule {self.name}"

    @property
    @abstractmethod
    def applicable(self) -> bool:
        pass


class CheckEqualStartLocations(ClassificationRule):
    """
    Check if all alignment (without the first and last) start locations are within margin of the median, looks per chromosome.
    """

    applicable = True

    def __init__(self, margin: int = 150):
        logger.debug(f"Rule {self.__class__.__name__} created, margin: {margin}")
        self.margin = margin

    @lru_cache(maxsize=1)
    def score(self, group):
        result = 1
        # if there is only one, we already know the awnser
        if len(group) < 2:
            return result

        # remove the first and last alignment, since these are more effected by cuts than actual differences
        # TODO: handle this smarter
        alignments = group.alignments.copy()
        alignments.pop(0)
        alignments.pop(-1)

        if len(alignments) <= 1:
            return 0
        #  set compresension, since we only need every chromosome once

        chromosomes_observed = {
            (x.alignment_chromosome, x.alignment_direction) for x in alignments
        }

        for chromosome_direction in chromosomes_observed:
            chromosome = chromosome_direction[0]
            direction = chromosome_direction[1]
            relevant_reads = [
                x for x in alignments if x.alignment_chromosome == chromosome
            ]
            relevant_reads = [
                x for x in relevant_reads if x.alignment_direction == direction
            ]

            positions = [x.alignment_chromosome_start for x in relevant_reads]
            median = statistics.median(positions)
            within_margin = [
                median - self.margin <= x <= median + self.margin for x in positions
            ]
            if not all(within_margin):
                result = 0

        logger.debug(f"Rule {self.__class__.__name__} is scored at {result}")
        return result


class CalculateRelativeChromosomalAbundance(ClassificationRule):
    """
    Calculate the fraction that aligns to the most occuring (directional) chromosome
    """

    applicable = True

    def __init__(self):
        logger.debug(f"Rule {self.__class__.__name__}")

    @lru_cache(maxsize=1)
    def score(self, group):
        directional_chromosomes = [
            f"{x.alignment_chromosome}_{x.alignment_direction}"
            for x in group.alignments
        ]
        if len(set(directional_chromosomes)) == 1:
            result = 1
        else:
            # calculate the relative occurance of the most common chromosome
            directional_chromosome_count = Counter(directional_chromosomes)
            total = sum(directional_chromosome_count.values())
            most_common_count = directional_chromosome_count.most_common(1)[0][1]

            result = most_common_count / total

        logger.debug(f"Rule {self.__class__.__name__} is scored at {result}")
        if result > 1:
            logger.critical("{self.__class__.__name__} Over 1!")
        return result


class CalculateAlternatingChromosomeRatio(ClassificationRule):
    """
    Calculate the amount of locations where the previous and next chromosome are the same, whilst not the same as the current chromosome
    """

    applicable = True

    def __init__(self, margin: int = 10):
        logger.debug(f"Rule {self.__class__.__name__} created, margin: {margin}")

    @lru_cache(maxsize=1)
    def score(self, group):
        alignments = group.alignments
        if len(group) <= 2:
            logger.debug("skipping alternation, group to short")
            return 0

        alternation_events = 0
        for i, aln in enumerate(alignments):
            # skip first and last
            if i == 0 or i == len(group) - 1:
                continue
            previous = f"{alignments[i-1].alignment_chromosome}_{alignments[i-1].alignment_direction}"
            current = f"{alignments[i].alignment_chromosome}_{alignments[i].alignment_direction}"
            following = f"{alignments[i+1].alignment_chromosome}_{alignments[i+1].alignment_direction}"
            if previous == following and previous != current:
                alternation_events += 1
        #  calculate the fraction, substraction since we cannot check position 1 and last.
        result = alternation_events / (len(group) - 2)
        if result > 1:
            logger.critical("{self.__class__.__name__} Over 1!")
        return result


class CheckNDirectionalChromosomesPresent(ClassificationRule):
    """
    Check if there are `directional_chromosomes` each occuring at least `minimum_occurrence` times within the alignment group.
    """

    applicable = False

    def __init__(self, directional_chromosomes=1, minimum_occurrence: int = 2):
        self.minimum_occurrence = minimum_occurrence
        self.directional_chromosomes = directional_chromosomes
        logger.debug(
            f"Rule {self.__class__.__name__} created, minimum_occurrence: {minimum_occurrence}"
        )

    @lru_cache(maxsize=1)
    def score(self, group):
        alignments = group.alignments
        result = 1
        chroms = [
            f"{x.alignment_chromosome}_{x.alignment_direction}" for x in alignments
        ]
        counts = Counter(chroms)
        # Filter on min occurrence
        counts = {k: v for k, v in counts.items() if v >= self.minimum_occurrence}

        if len(counts.keys()) == self.directional_chromosomes:
            result = 1
        else:
            result = 0
        return result


class Check1DirectionalChromosomesPresent(CheckNDirectionalChromosomesPresent):
    applicable = True

    def __init__(self) -> None:
        super().__init__(directional_chromosomes=1)

    def score(self, group):
        super_score = super().score(group)
        return super_score


class Check2DirectionalChromosomesPresent(CheckNDirectionalChromosomesPresent):
    applicable = True

    def __init__(self) -> None:
        super().__init__(directional_chromosomes=2)

    def score(self, group):
        super_score = super().score(group)
        return super_score


class Check3DirectionalChromosomesPresent(CheckNDirectionalChromosomesPresent):
    applicable = True

    def __init__(self) -> None:
        super().__init__(directional_chromosomes=3)

    def score(self, group):
        super_score = super().score(group)
        return super_score


class Check4DirectionalChromosomesPresent(CheckNDirectionalChromosomesPresent):
    applicable = True

    def __init__(self) -> None:
        super().__init__(directional_chromosomes=4)

    def score(self, group):
        super_score = super().score(group)
        return super_score


class CheckNChromosomesPresent(ClassificationRule):
    """
    Check if there are two chromosomes within the alignment group.
    """

    applicable = False

    def __init__(self, chromosome_count=2, minimum_occurrence: int = 1):
        self.minimum_occurrence = minimum_occurrence
        self.chromosome_count = chromosome_count
        logger.debug(
            f"Rule {self.__class__.__name__} created, minimum_occurrence: {minimum_occurrence}"
        )

    @lru_cache(maxsize=1)
    def score(self, group):
        alignments = group.alignments
        result = 1
        chroms = [x.alignment_chromosome for x in alignments]
        counts = Counter(chroms)
        # Filter on min occurrence
        counts = {k: v for k, v in counts.items() if v >= self.minimum_occurrence}
        if len(counts.keys()) == self.chromosome_count:
            result = 1
        else:
            result = 0
        return result


class Check1ChromosomesPresent(CheckNChromosomesPresent):
    """
    Check if there is one chromosome in the alignments
    """

    applicable = True

    def __init__(self) -> None:
        super().__init__(chromosome_count=1)


class Check2ChromosomesPresent(CheckNChromosomesPresent):
    """
    Check if there is one chromosome in the alignments
    """

    applicable = True

    def __init__(self) -> None:
        super().__init__(chromosome_count=2)


class Check3ChromosomesPresent(CheckNChromosomesPresent):
    """
    Check if there is one chromosome in the alignments
    """

    applicable = True

    def __init__(self) -> None:
        super().__init__(chromosome_count=3)


class Check4ChromosomesPresent(CheckNChromosomesPresent):
    """
    Check if there is one chromosome in the alignments
    """

    applicable = True

    def __init__(self) -> None:
        super().__init__(chromosome_count=4)


class CheckMinimumCountAlignments(ClassificationRule):
    """
    Check if there at least `minimum` alignments in the group. return 1 if this is true
    """

    applicable = True

    def __init__(self, minimum: int = 4):
        self.minimum = minimum
        logger.debug(f"Rule {self.__class__.__name__} created, margin: {minimum}")

    @property
    def name(self):
        return f"{self.__class__.__name__}_{self.minimum}"

    @lru_cache(maxsize=1)
    def score(self, group):
        result = 1
        if len(group) < self.minimum:
            result = 0
        logger.debug(f"Rule {self.__class__.__name__} is scored at {result}")
        return result


class CheckMaximumCountAlignments(ClassificationRule):
    """
    Check if there at most `maximum` alignments in the group. return 1 if this is true
    """

    applicable = True

    def __init__(self, maximum: int = 10):
        self.maximum = maximum
        logger.debug(f"Rule {self.__class__.__name__} created, margin: {maximum}")

    @property
    def name(self):
        return f"{self.__class__.__name__}_{self.maximum}"

    @lru_cache(maxsize=1)
    def score(self, group):
        result = 1
        if len(group) >= self.maximum:
            result = 0
        logger.debug(f"Rule {self.__class__.__name__} is scored at {result}")
        return result


class CalculateContiguousStartPositions(ClassificationRule):
    """
    Check that all start positions follow each other up with at maximum `margin` space between them.
    This uses the length to determine the next expected start site.
    """

    applicable = True

    def __init__(self, margin: int = 50):
        self.margin = margin
        logger.debug(f"Rule {self.__class__.__name__} created, margin: {margin}")

    @lru_cache(maxsize=1)
    def score(self, group):
        if len(group) < 2:
            return 0

        contiguous_events = 0
        expected_start = group.alignments[0].first_cigar_value
        for x in group.alignments:
            start = x.first_cigar_value
            # check if within margin
            if start < expected_start - self.margin:
                pass
            elif start > expected_start + self.margin:
                pass
            else:
                contiguous_events += 1
            # set params for next loop
            end = start + x.alignment_length
            expected_start = end

        result = contiguous_events / (len(group))
        logger.debug(f"Rule {self.__class__.__name__} is scored at {result}")
        if result > 1:
            logger.critical("{self.__class__.__name__} Over 1!")
        return result


class Check1BackBoneInAlignments(ClassificationRule):
    """
    Check if there is any mapping against a referent where the name starts with `BB`.

    """

    applicable = True

    def __init__(self):
        logger.debug(f"Rule {self.__class__.__name__} created")

    @lru_cache(maxsize=1)
    def score(self, group):
        backbone_in_chromosome = [
            x.startswith("BB")
            for x in group.alignment_chromosomes_present(directional=False)
        ]
        if sum(backbone_in_chromosome) == 1:
            result = 1
        else:
            result = 0
        logger.debug(f"Rule {self.__class__.__name__} is scored at {result}")
        return result


class CheckAllBackBoneInAlignments(ClassificationRule):
    """
    Check if there is any mapping against a referent where the name starts with `BB`.

    """

    applicable = True

    def __init__(self):
        logger.debug(f"Rule {self.__class__.__name__} created")

    @lru_cache(maxsize=1)
    def score(self, group):
        backbone_in_chromosome = [
            x.alignment_chromosome.startswith("BB") for x in group.alignments
        ]
        if all(backbone_in_chromosome):
            result = 1
        else:
            result = 0
        logger.debug(f"Rule {self.__class__.__name__} is scored at {result}")
        return result


class CheckTwoStartLocations(ClassificationRule):
    """Check if we have starting positions that are at least `minimum distance` apart in a single chromosome."""

    applicable = True

    def __init__(self, minimum_distance=1000):
        logger.debug(f"Rule {self.__class__.__name__} created")
        self.minimum_distance = minimum_distance

    def score(self, group):
        alignments = group.alignments
        #  set compresension, since we only need every chromosome once
        chromosomes_observed = {
            (x.alignment_chromosome, x.alignment_direction) for x in alignments
        }
        result = 0

        for chromosome_direction in chromosomes_observed:
            chromosome = chromosome_direction[0]
            direction = chromosome_direction[1]
            relevant_reads = [
                x for x in alignments if x.alignment_chromosome == chromosome
            ]
            relevant_reads = [
                x for x in relevant_reads if x.alignment_direction == direction
            ]

            positions = [x.alignment_chromosome_start for x in relevant_reads]
            pos_min, pos_max = min(positions), max(positions)

            if pos_max - pos_min > self.minimum_distance:
                result = 1

        return result


class CheckTwoStartLocationsByOwnLength(ClassificationRule):
    """Check if we have starting positions that are at least `minimum distance` apart in a single chromosome."""

    applicable = True

    def __init__(self, minimum_distance=1000):
        logger.debug(f"Rule {self.__class__.__name__} created")
        self.minimum_distance = minimum_distance

    def score(self, group):
        alignments = group.alignments
        #  set compresension, since we only need every chromosome once
        chromosomes_observed = {
            (x.alignment_chromosome, x.alignment_direction) for x in alignments
        }
        result = 0

        for chromosome_direction in chromosomes_observed:
            chromosome = chromosome_direction[0]
            direction = chromosome_direction[1]
            relevant_reads = [
                x for x in alignments if x.alignment_chromosome == chromosome
            ]
            relevant_reads = [
                x for x in relevant_reads if x.alignment_direction == direction
            ]
            mean_length = statistics.mean([x.alignment_length for x in relevant_reads])

            positions = [x.alignment_chromosome_start for x in relevant_reads]
            pos_min, pos_max = min(positions), max(positions)

            if pos_max - pos_min > mean_length:
                result = 1

        return result


class CalculateTripletChromosomes(ClassificationRule):
    """
    Check if the pattern corresponds with a Triplicate (eg Babo-Babo-Ins or Babo-Ins-Ins)

    So let say we have:
    BBIBBIBBI
    123456789

    1: skip
    2: case 1
    3: case 3
    4: case 2
    5: case 1
    6: case 3
    """

    applicable = True

    def score(self, group):
        alignments = group.alignments

        if len(group) <= 2:
            logger.debug(f"skipping {self.__class__.__name__}, group to short")
            return 0

        alternation_events = 0
        for i, aln in enumerate(alignments):
            # skip first and last
            if i == 0 or i == len(group) - 1:
                continue
            previous = f"{alignments[i-1].alignment_chromosome}_{alignments[i-1].alignment_direction}"
            current = f"{alignments[i].alignment_chromosome}_{alignments[i].alignment_direction}"
            following = f"{alignments[i+1].alignment_chromosome}_{alignments[i+1].alignment_direction}"
            try:
                following_two = f"{alignments[i+1].alignment_chromosome}_{alignments[i+2].alignment_direction}"
            except IndexError:
                following_two = None

            if previous == current and previous != following:
                # case 1
                alternation_events += 1
            elif following == current and previous != following:
                # case 2
                alternation_events += 1

            elif following_two:
                if current != following and following == following_two:
                    # case 3
                    alternation_events += 1

        result = alternation_events / (len(group) - 2)
        #  calculate the fraction, substraction since we cannot check position 1 and last.
        if result > 1:
            logger.critical("{self.__class__.__name__} Over 1!")
        return result


# class CalculateRelativeStartPositions(ClassificationRule):
#     """
#     Calculate the uniformaty of the starting positions, to make more fine grained decisions
#     """

#   applicable = True

#     def __init__(self, margin: int = 10):
#         self.margin = margin
#         logger.debug(f"Rule {self.__class__.__name__} created, margin: {margin}")
#         raise NotImplementedError

#     @lru_cache(maxsize=1)
#     def score(self, group):
#         pass


# class CheckNoAlignmentOverlaps(ClassificationRule):
#     """
#     Check that the Alignments do not overlap each other. overlaps less than `margin` are allowed.
#     """

#   applicable = True

#     def __init__(self, margin: int = 10):
#         logger.debug(f"Rule {self.__class__.__name__} created, margin: {margin}")

#     @lru_cache(maxsize=1)
#     def score(self, group):
#         logger.debug(f"Rule {self.__class__.__name__} is scored at {result}")
#         raise NotImplementedError
#         return 1


class CheckMinimumUnalignedGap(ClassificationRule):
    """
    Check if the minimum unaligned gap between segments is smaller than the target value
    """

    applicable = True

    def __init__(self, target_value: int = 30):
        logger.debug(
            f"Rule {self.__class__.__name__} created, target_value: {target_value}"
        )
        self.target_value = target_value

    @lru_cache(maxsize=1)
    def score(self, group):
        # result = 1 if group.find_median_unaligned_region() > self.target_value else 0
        # get all gaps between segments (start and end might be hard to align)
        inter_segment_gaps = group.find_unaligned_regions()[1:-1]
        if len(inter_segment_gaps) > 0:
            minimum_segment_gap = min(inter_segment_gaps)
        else:
            minimum_segment_gap = 0
        if minimum_segment_gap < self.target_value:
            result = 1
        else:
            result = 0
        logger.debug(f"Rule {self.__class__.__name__} is scored at {result}")
        return result


class CalculateSegmentCountPercentageTenSegments(ClassificationRule):
    """
    Check that the Alignments do not overlap each other. overlaps less than `margin` are allowed.
    """

    applicable = True

    def __init__(self, segments: int = 10):
        self.segments = segments
        logger.debug(f"Rule {self.__class__.__name__} created, segments: {segments}")

    @lru_cache(maxsize=1)
    def score(self, group):
        if len(group) >= self.segments:
            result = 1
        else:
            result = len(group) / self.segments
        return result
