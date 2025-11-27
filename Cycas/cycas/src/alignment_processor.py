from copy import deepcopy
from statistics import median
from loguru import logger


class AlignmentProcessor:
    def _get_spaced_integers(self, integers, spacer):
        """
        given a list of integers, this will return all integers that are at least `spacer` away from the last one
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
        self, group, distance=20
    ):
        """
        Find all groups of alignment per chromosome that are at least `distance` seperated on said chromosome.
        """
        output = []

        for chromosome in group.alignment_chromosomes_present(directional=False):
            group_chr = deepcopy(group)
            group_chr.extract_alignments_by_chromosome(chromosome)
            start_positions = [
                x.alignment_chromosome_start for x in group_chr.alignments
            ]
            lengths = [x.alignment_length for x in group_chr.alignments]

            # if we can set the distance automatically we use 2x the median
            # 2x caused issues with the multi amplicon panels and is now deemed to convervative
            if distance == "auto":
                distance = round(max(lengths))
            distance = int(distance)

            # if max(start_positions) - min(start_positions) > median(lengths):
            #     print('potential overlapping amplicon')
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
                    current_group.consensus_chromosome = chromosome
                    output.append(current_group)
            else:
                # single group
                group_chr.consensus_chromosome = chromosome
                output.append(group_chr)
        return output

    def create_paired_alignment_groups_per_directional_chromosome_and_distance(
        self, group, distance="auto"
    ):
        output = []

        for directional_chromosome in group.alignment_chromosomes_present():
            group_chr = deepcopy(group)
            group_chr.extract_alignments_by_directional_chromosome(
                directional_chromosome
            )
            start_positions = [x.alignment_chromosome_start for x in group.alignments]
            # if we can set the distnace automatically we use 2x the median
            if distance == "auto":
                lengths = [x.alignment_length for x in group_chr.alignments]
                distance = round(2 * median(lengths))
            distance = int(distance)

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
                    current_group.consensus_chromosome = directional_chromosome
                    output.append(current_group)
            else:
                # single group
                group_chr.consensus_chromosome = directional_chromosome
                output.append(group_chr)

        return output
