import re
from collections import defaultdict
from functools import lru_cache
from math import ceil, floor
from pprint import pprint
from statistics import mean, median
from typing import Dict, List, Tuple

import pysam
from loguru import logger

from .config import ReadClassificationConfig
from .consecutiveness import ConsecutivenessStrategyFactory
from .filters import Filter
from .plots.spagettogram import create_spagettogram


class Alignment:
    def __init__(self, pysam_obj, keep_obj=True):
        """
        Alternative to the pysam AlignedSegment, stores all info about a read segment
        """
        # logger.debug(f"{self.__class__.__name__} created")
        if keep_obj:
            self.pysam_obj = pysam_obj

        self.readname = pysam_obj.query_name
        self.alignment_length = pysam_obj.query_alignment_length
        self.alignment_chromosome = pysam_obj.reference_name
        self.alignment_chromosome_start = pysam_obj.reference_start
        self.alignment_direction = "-" if pysam_obj.is_reverse else "+"
        self.cigarstring = pysam_obj.cigarstring
        self.cigars = self.get_cigars(allow_flip=True)  # flips if -
        self.raw_read_length = pysam_obj.infer_read_length()
        self.mapping_quality = pysam_obj.mapping_quality
        self.sequence = pysam_obj.query_alignment_sequence
        self.nm_tag = pysam_obj.get_cigar_stats()[0][-1]
        self.mapq = pysam_obj.mapq

        self.extended_cigar = self.extend_cigar(self.cigars)
        self.extended_seq = self.extend_sequence(
            seq=self.sequence, extended_cigar=self.extended_cigar
        )
        self.extended_qual = self.extend_qual(
            pysam_obj.query_alignment_qualities, self.extended_cigar
        )

    def __str__(self) -> str:
        """String representation."""
        print(f"{self.readname} - {self.cigarstring} ({self.alignment_direction})")

    def to_pysam(self) -> pysam.AlignedSegment:
        seg = self.pysam_obj
        seg.query_name = self.readname  # e.g. with _0 suffix if renamed
        seg.reference_name = self.alignment_chromosome  # updated reference id
        seg.reference_start = self.alignment_chromosome_start  # updated start
        seg.mapping_quality = self.mapping_quality

        seg.flag = seg.flag  # updated flag
        seg.cigar = seg.cigar  # updated cigar, like with hard clips etc.

        # TODO: modify query sequence and qualities
        # seg.query_sequence = seg.query_sequence  # processed sequence if changed
        # seg.query_qualities = seg.query_qualities  # processed qualities if changed

        # TODO: modify supplementary alignments?
        # NOTE: if read is split, only 1 primary will remain. Recalculate flag based on MAPQ?

        # Update any optional or auxiliary tags (NM, AS, etc.)
        # for tag, value in seg.tags:
        #     if tag == "SA":
        #         sup_alignments = value.split(";")
        #     seg.set_tag(tag, value)
        return seg

    def get_cigars(self, allow_flip=True) -> List[Tuple[str, int]]:
        """
        Converts a cigarstring into a workable format:
        '721H149M686H' -> [('H', 721), ('M', 149), ('H', 686)]
        """
        # Find all alternating letters and numbers
        cigars = list(
            zip(
                re.findall(r"\D+", self.cigarstring),
                re.findall(r"\d+", self.cigarstring),
            )
        )
        # convert numerical part to an int
        cigars = [(x[0], int(x[1])) for x in cigars]

        # reverse the array if the read is actually backwards
        if self.alignment_direction == "-" and allow_flip:
            cigars.reverse()

        return cigars

    @property
    def directional_chromosome(self) -> str:
        """return the chromosome and direction as a string eg chr17- when it maps to chr17 and is mapping in reverse"""
        return self.alignment_chromosome + self.alignment_direction

    @property
    def first_cigar_value(self) -> int:
        """Extract the first integer from the cigars tuple created by `get_cigars`"""
        try:
            value = int(self.cigars[0][1])
        except IndexError:
            logger.critical("Could not retrieve cigar value for alignment")
            value = 0

        return value

    def extend_sequence(self, seq, extended_cigar, reverse=False):
        """Given a sequence and a extended cigar string it will expand the sequence to match the cigarstring"""
        pos = 0
        ext_seq = ""
        previous_deletion = False
        skip_next_M = False
        for i, cig in enumerate(extended_cigar):
            if cig in ["M", "I"]:
                if cig == "M" and skip_next_M:
                    skip_next_M = False
                    continue
                ext_seq += seq[pos]
                pos += 1
                previous_deletion = False

            elif cig == "D":
                if previous_deletion:
                    ext_seq += "-"
                    continue

                if reverse:
                    ext_seq += seq[pos]
                    pos += 1
                    skip_next_M = True

                ext_seq += "-"
                previous_deletion = True

        return ext_seq

    def extend_qual(self, qual, ext_cigar: List[str]) -> List[int]:
        """Create an extended representation for the quality string. It will add 0 quality for each location that is missing in the original read."""
        pos = 0
        ext_qual = []

        for cig in ext_cigar:
            if cig in ["M", "I"]:
                ext_qual.append(qual[pos])
                pos += 1

            if cig == "D":
                ext_qual.append(0)

        return ext_qual

    # Reverse the Cigars!!
    def extend_cigar(self, cigars, reverse_cigar=False):
        """
        Extend the cigar string into a full representation.
        so 3H2M1I5H becomes:
        HHHMMIHHHHH
        """
        cigar = ""

        for cig in cigars:
            match = cig[0]
            amount = int(cig[1])
            if match in ["H", "S"]:
                continue

            cigar += match * amount

        if reverse_cigar or self.alignment_direction == "-":
            cigar = cigar[::-1]
        return cigar


class AlignmentGroup:
    def __init__(
        self,
        read_name: str,
        bam: pysam.IndexedReads,
        filters: List[Filter],
        keep_per_read_objects: bool = True,
        config: ReadClassificationConfig = None,
    ):
        """
        Class to store all alignments of a given read in a bam file
        """
        logger.debug(f"{self.__class__.__name__} created")
        self.config = config or ReadClassificationConfig()
        self.read_name = read_name
        self.read_length_raw = 0  # gets updated by the get_alignments
        self.filters = filters
        self.keep_per_read_objects = keep_per_read_objects
        self.high_reject_ratio = False

        (self.alignments, self.alignment_count_prefilter) = self.get_alignments(bam)

        self.sort_alignments()
        self.assign_alignment_ids()
        self.structure = self.create_intermediate_read_structure()
        self.clusters = self.parse_structure()
        self.consensus_chromosome = None

    def __len__(self):
        """
        The length of the group is the number of alignments in it.
        """
        return len(self.alignments)

    def __repr__(self) -> str:
        return f"AlignmentGroup for {self.read_name} containing {len(self)} Alignments on {self.alignment_chromosomes_present()}"

    @classmethod
    def from_existing(
        cls, existing_group: "AlignmentGroup", id: str, alignments: List
    ) -> "AlignmentGroup":
        """
        Create a new AlignmentGroup instance based on an existing one,
        but replacing the alignments list.
        """
        new_group = cls.__new__(cls)  # create an uninitialized instance

        # Copy relevant attributes from existing_group
        new_group.read_name = existing_group.read_name + f"_{id}"
        new_group.read_length_raw = existing_group.read_length_raw
        new_group.filters = existing_group.filters
        new_group.keep_per_read_objects = existing_group.keep_per_read_objects
        new_group.consensus_chromosome = existing_group.consensus_chromosome

        # Override the alignments and related attributes
        new_group.alignments = alignments
        new_group.alignment_count_prefilter = len(alignments)

        return new_group

    def get_start_position(self):
        "return the median start position of all alignments in the group, note that it does not take chromosomes into account."
        start = min([x.alignment_chromosome_start for x in self.alignments])
        return floor(start)

    def get_end_position(self):
        "return the median end position of all alignments in the group, note that it does not take chromosomes into account."
        end = max(
            [x.alignment_chromosome_start + x.alignment_length for x in self.alignments]
        )
        return ceil(end)

    def count_aligned_bases(self) -> int:
        """
        Count the amount of bases that align against the reference genome for all alignments in the group.
        """
        return sum([aln.alignment_length for aln in self.alignments])

    def find_longest_unaligned_region(self) -> int:
        """
        Find the longest region that is unaligned
        """
        regions = self.find_unaligned_regions()
        if regions:
            return max(regions)
        else:
            return 0

    def find_median_unaligned_region(self) -> int:
        regions = self.find_unaligned_regions()
        return median(regions)

    def find_mean_unaligned_region(self) -> int:
        regions = self.find_unaligned_regions()
        return mean(regions)

    @lru_cache(1)
    def find_unaligned_regions(self) -> List[int]:
        result = []
        # add initial unaligned segment
        result.append(self.alignments[0].first_cigar_value)

        # calculate the gaps between the alignments
        for x, y in zip(self.alignments, self.alignments[1:]):
            gap = y.first_cigar_value - (x.first_cigar_value + x.alignment_length)
            result.append(gap)

        # add last unaligned segment
        result.append(self.alignments[-1].cigars[-1][1])

        return result

    def create_intermediate_read_structure(self):
        """
        Create an alternative representation of a read with the structure:
        bases_on_read:type:orient:assembly:alignment_position:read_position:mapping_length:ID
        """
        structure = []
        gaps = self.find_unaligned_regions()

        for aln, gap in zip(self.alignments, gaps):
            if gap > 0:
                structure.append(f"{gap}:U")
            elif gap < 0:
                structure.append(f"{gap}:O")

            len_bases = aln.alignment_length
            len_mapping = len(aln.sequence)

            type = "BB" if aln.alignment_chromosome.startswith("BB") else "I"
            orientation = "F" if aln.alignment_direction == "+" else "R"
            assembly = aln.alignment_chromosome
            aln_position = aln.alignment_chromosome_start
            read_position = aln.first_cigar_value

            structure.append(
                f"{len_bases}:{type}:{orientation}:{assembly}:{aln_position}:{read_position}:{len_mapping}:{aln.id}"
            )

        # If concatemer ends in a gap, add it
        last_aln_cigars = self.alignments[-1].cigars
        if last_aln_cigars[-1][0] in ["H", "S"]:
            structure.append(f"{last_aln_cigars[-1][1]}:U")

        return ",".join(structure)

    def parse_structure(self) -> List[Tuple[int, List[Dict], bool]]:
        """
        Parse a structure string into clusters with consecutiveness info.

        Args:
            structure: Comma-separated segment string.

        Returns:
            List of tuples: (cluster_id, cluster_blocks, is_consecutive)
        """
        blocks = self.structure.strip().split(",")
        insert_blocks = self.parse_insert_blocks(blocks)
        # NOTE: TEMPORARY
        self.all_blocks = insert_blocks
        clusters = self.group_by_rel_id(insert_blocks)
        return clusters

    def parse_insert_blocks(self, blocks: List[str]) -> List[Dict]:
        """
        Extract insert blocks and their metadata from parsed structure strings.

        Args:
            blocks: List of colon-delimited segment strings.

        Returns:
            A list of dictionaries containing parsed insert block fields:
            - "ori": orientation
            - "ref": reference name
            - "pos": genomic position
            - "read_pos": position within read
            - "rel_id": relative segment ID
        """
        insert_blocks = []

        for block in blocks:
            parts = block.split(":")
            if len(parts) < 8 or parts[1] != "I":
                continue
            _, btype, ori, ref, aln_pos, read_pos, _, rel_id = parts
            if ref.startswith(self.config.backbone_prefix):
                continue  # Skip backbone
            try:
                pos = int(aln_pos)
            except ValueError:
                continue
            insert_blocks.append(
                {
                    "ori": ori,
                    "ref": ref,
                    "pos": pos,
                    "read_pos": int(read_pos),
                    "rel_id": int(rel_id),
                }
            )

        return insert_blocks

    def group_by_rel_id(self, blocks: List[Dict]) -> List[Tuple[int, List[Dict], bool]]:
        """
        Cluster insert blocks by their relative ID (rel_id) and determine consecutiveness.

        Args:
            blocks: A list of parsed insert block dictionaries.

        Returns:
            List of tuples: (cluster_id, cluster_blocks, is_consecutive)
        """
        grouped_blocks = defaultdict(list)
        for block in blocks:
            grouped_blocks[int(block["rel_id"])].append(block)

        clusters = []
        for cluster_id, cluster in grouped_blocks.items():
            is_consec = self.is_cluster_consecutive(cluster_id, blocks)
            clusters.append((cluster_id, cluster, is_consec))

        return clusters

    def is_cluster_consecutive(self, cluster_id: int, all_blocks: List[Dict]) -> bool:
        """
        Determine if a cluster is consecutive in read order using the configured strategy.

        Args:
            cluster_id: The relative ID of the cluster.
            all_blocks: All blocks for this read.

        Returns:
            True if consecutive according to the selected strategy.
        """
        method = self.config.consecutiveness_method
        threshold = self.config.consecutiveness_threshold

        strategy = ConsecutivenessStrategyFactory().get_strategy(method)
        metric = strategy.calculate_metric(cluster_id, all_blocks)
        return strategy.is_consecutive(metric, threshold)

    def alignment_chromosomes_present(self, directional=True) -> Tuple[str, str]:
        """Create a tuple from the chromosomes where the is alignment."""
        if directional:
            return {
                (x.alignment_chromosome, x.alignment_direction) for x in self.alignments
            }
        else:
            return {(x.alignment_chromosome) for x in self.alignments}

    def extract_alignments_by_chromosome(self, chromosome):
        """Get alignments that align against a certain chromosome."""
        self.alignments = [
            x for x in self.alignments if x.alignment_chromosome == chromosome
        ]

    def extract_alignments_by_directional_chromosome(self, chromosome):
        self.alignments = [
            x for x in self.alignments if x.alignment_chromosome == chromosome[0]
        ]
        self.alignments = [
            x for x in self.alignments if x.alignment_direction == chromosome[1]
        ]

    def extract_alignments_by_start_position(self, position, margin):
        """Find all alignments that align within `margin` of `position` and only keep those in the object."""
        result = []
        for x in self.alignments:
            if abs(x.alignment_chromosome_start - position) < margin:
                # read close by
                result.append(x)
        self.alignments = result
        logger.debug("done filtering by position")

    def _alignments_rejected_warn(
        self, alignments_rejected, alignments, warn_threshold=1.8
    ):
        if (
            len(alignments_rejected) / (len(alignments_rejected) + len(alignments))
            > warn_threshold
        ):
            self.high_reject_ratio = True

    @staticmethod
    def _update_alignments_rejected(cigar_string, alignments_rejected, filter_status):
        """
        Helper function to update the failed alignments and why
        """
        failed_filters = [key for key, value in filter_status.items() if not value]
        for filt in failed_filters:
            alignments_rejected[cigar_string].add(filt)
        return alignments_rejected

    def get_alignments(self, bam: pysam.IndexedReads) -> List[Alignment]:
        """
        Collect all pysam alignmentsegment objects in a wrapper class in the propper way.
        This is done since the cigarstring of reversely mapping reads needs to be reversed.
        """

        logger.debug(f"generating alignment objects for {self.read_name}")
        pysam_alignments = bam.find(self.read_name)

        alignments = []
        alignments_rejected = defaultdict(set)

        # need to do a manual count due to the fact that this is a generator
        alignment_count = 0

        for alignment in pysam_alignments:
            alignment_count += 1
            # set the raw read length if its missing, prevents double lookup.
            if self.read_length_raw == 0:
                self.read_length_raw = alignment.infer_read_length()

            # we use the cigar string as the ID for the alignments
            cigar = alignment.cigarstring
            logger.debug(f"\t{cigar}")

            # apply all filters applied to this group
            # TODO: set filter status on object
            filter_status = {x.name: x.filter(alignment) for x in self.filters}

            #  If a filter returns False we ignore the alignment
            if not all(filter_status.values()):
                self._update_alignments_rejected(
                    cigar, alignments_rejected, filter_status
                )
            else:
                # All reads passing filtering
                alignments.append(Alignment(alignment, self.keep_per_read_objects))

        # Trow some warnings if we reject stuff
        self._alignments_rejected_warn(alignments_rejected, alignments)
        return (
            alignments,
            alignment_count,
        )

    def assign_alignment_ids(self, position_margin: int = 250) -> None:
        """
        Assigns an integer ID to each alignment in self.alignments.
        IDs are reused for alignments mapped to similar regions within the same chromosome.
        """
        visited: defaultdict[str, dict[int, int]] = defaultdict(dict)
        next_id = 0

        for aln in self.alignments:
            assembly = aln.alignment_chromosome
            aln_position = aln.alignment_chromosome_start

            ID = None
            for pos, existing_id in visited[assembly].items():
                if abs(pos - aln_position) <= position_margin:
                    ID = existing_id
                    break

            if ID is None:
                ID = next_id
                next_id += 1

            visited[assembly][aln_position] = ID
            aln.id = ID  # Assign ID to the alignment object

    def sort_alignments(self):
        logger.debug("sorting AlignmentGroup")
        self.alignments.sort(key=lambda x: x.first_cigar_value)

    def create_spagettogram(self, title=None):
        create_spagettogram(self.alignments, title=title)

    def _print_group_properties(self, cigar_length=6) -> None:
        """
        QOL debug function to print info about the ReadGroup
        """
        pprint(
            [
                (
                    x.alignment_chromosome,
                    x.alignment_direction,
                    x.cigars[:cigar_length],
                    x.mapping_quality,
                )
                for x in self.alignments
            ]
        )
