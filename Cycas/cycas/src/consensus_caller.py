import math
from abc import ABC, abstractmethod
from collections import Counter
from typing import Any, Dict

from Bio.Seq import Seq
from loguru import logger

from .barcode_extractor import BarcodeExtrator


class ConsensusAlignment:
    def __init__(self, extended_seq, extended_cigar, extended_qual, direction=None):
        self.consensus_cigar = ""
        self.consensus_seq = ""
        self.consensus_quality = ""
        self.current_pos = 0
        self.seq = extended_seq
        self.cigar = extended_cigar
        self.qual = extended_qual
        self.direction = direction
        # self.qual = np.array(extended_qual)

    def __repr__(self):
        return f"{self.seq}\n{self.cigar}\n{self.convert_to_ascii()}"

    def __len__(self):
        return len(self.seq)

    def get_cigar_pos(self, pos):
        try:
            cigar = self.cigar[pos]
        except IndexError:
            cigar = "-"
        return cigar

    def get_seq_pos(self, pos):
        try:
            nucleotide = self.seq[pos]
        except IndexError:
            nucleotide = "-"
        return nucleotide

    def get_qual_pos(self, pos):
        try:
            phred_quality = self.qual[pos]
        except IndexError:
            phred_quality = 0
        return phred_quality

    def add_empty_location(self, pos, end_start_char="_"):
        """
        Add information to the consensus object by checking its status for that position and extending the required data.
        """
        # if the read has started there is more than just _ in the sequence
        started = set(self.seq[:pos]) != set(end_start_char)
        # if the read has ended there is only _ left
        ended = set(self.seq[pos:]) == set(end_start_char)
        if started and not ended:
            extender = "-"
        else:
            extender = end_start_char
        self.seq = self.seq[:pos] + extender + self.seq[pos:]
        self.cigar = self.cigar[:pos] + extender + self.cigar[pos:]
        self.qual.insert(pos, 0)

    def convert_to_ascii(self, phred_offset=33):
        """
        convert all qualities to its representative base (+33)
        """
        return "".join([chr(x + phred_offset) for x in self.qual])

    def get_aligned_position(self, position, offset):
        """
        Used to make the barcode, offset for when alignment doesnt start correctly.
        Also corrects for inserts etc since we are looking for a specific location in the referent.

        eg position 32 in below structure will return A at index 36:
        >BBCR
        GGGCGGTATGTCA-TGCACAC--G-AATCCCGAAGANTGTTGTCCATTCATTGAATATGAGATCTCNATGGTATGATCAATATNCGGATGCGATATTGATANCTGATAAATCATATATGCATAATCTCACATTATATTTATTATAATAAATCATCGTAGATATACACAATGTGAATTGTATACAATGGATAGTATAACTATCCAATTTCTTTGAGCATTGGCCTTGGTGTAGATTGCATGACATACCGCCC
        >structure:
        GGGCGG--TGTCA-TG---AC--G-AATGCCGAAGAATGTTG-TC-CATTCATTGAATAT--A-ATCT-GATGGTATGATCAATATGCGGATGCGATATTGATA-TCTAAT-AATCACTATATGCATAATCTCA-CATTATATTTATTATAAT-AAAT-CATCGTAGAT-ATACACAATGTGAATTGTATACAATGGAT-AGTA-TAACTATCCAATTTCTTT-GAGCATTGGCCTTGGTGTAGATTGCA--TGACAT-CCGCCC_
        MMMMMMDDMMMMM-MMDDDMM--M-MMMMMMMMMMMMMMMMM-MM-MMMMMMMMMMMMMM-DMDMMMMDMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM-MMMMMMDMMMMMIMMMMMMMMMMMMMMMM-MMMMMMMMMMMMMMMMMM-MMMM-MMMMMMMMMM-MMMMMMMMMMMMMMMMMMMMMMMMMMMMM-MMMM-MMMMMMMMMMMMMMMMMM-MMMMMMMMMMMMMMMMMMMMMMMMMM--MMMMMMDMMMMMM_
        GGGCGGTATG--G-TGCGCAC----AATGCCGAAGAATGTTG-TC-CATTCATTGAATATAGAGATCTCGATGGTATGATCAATATGCGGATGCGATATTGATA-TCTGATAAATCA-TATATGCATAATCTCA-CATTATATTTATTATAAT-AAAT-CATCGTAGAT-ATACACAATGTGAATTGTATACAATGGAT-AGTA-TAACTATGCAA-ATCTTT-GAGCATTGGCCTTGGTGTAGATTGCA--TGACGTACCGCCC_
        MMMMMMMMMMDDM-MMMMMMM--D-MMMMMMMMMMMMMMMMM-MM-MMMMMMMMMMMMMMIMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM-MMMMMMMMMMMM-MMMMMMMMMMMMMMMM-MMMMMMMMMMMMMMMMMM-MMMM-MMMMMMMMMM-MMMMMMMMMMMMMMMMMMMMMMMMMMMMM-MMMM-MMMMMMMMMMMDMMMMMM-MMMMMMMMMMMMMMMMMMMMMMMMMM--MMMMMMMMMMMMM_
        GGGCGGTATGTAG-TGCACAC----AATCCCGAAGAATGTTG-TC-CATTCATTGAATAT-GAGATCTCGATGGTATCATCAATATGCGGATGCGATATTGATA-TCTGATAAATCA-TATATGCATAATCTCA-CATTATATTTATTATAAT-AAAT-CATCGTAGAT-ATACACAATGTGAATTGTATACAATGGAT-AGTA-TAACTATCCAATTTCTTT-GAGCATTGGCCTTGGTGTAGATTGCA--TGACAT-CCGCCC_
        MMMMMMMMMMMMM-MMMMMMM--D-MMMMMMMMMMMMMMMMM-MM-MMMMMMMMMMMMMM-MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM-MMMMMMMMMMMM-MMMMMMMMMMMMMMMM-MMMMMMMMMMMMMMMMMM-MMMM-MMMMMMMMMM-MMMMMMMMMMMMMMMMMMMMMMMMMMMMM-MMMM-MMMMMMMMMMMMMMMMMM-MMMMMMMMMMMMMMMMMMMMMMMMMM--MMMMMMDMMMMMM_
        GGGCGGTATGTCA-TGCACAC--A-AA--CCGAAGAATGTTG-TC-CATTCATTGAATAT-GAGATCTCGATGGTATGATCAATATGCGGATGCGATATTGATA-TCTGATAAATCA-TATATGCATAATCTCA-CATTATATTTATTATAAT-AAAT-CATCGTAGAT-ATACACAATGTGAATTGTATACAATGGAT-AGTA-TAACTATCCAATTTCTTT-GAGCATTGGCCTTGGTGTAGATTGCA--TGACATACCGCCC_
        MMMMMMMMMMMMM-MMMMMMM--M-MMDDMMMMMMMMMMMMM-MM-MMMMMMMMMMMMMM-MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM-MMMMMMMMMMMM-MMMMMMMMMMMMMMMM-MMMMMMMMMMMMMMMMMM-MMMM-MMMMMMMMMM-MMMMMMMMMMMMMMMMMMMMMMMMMMMMM-MMMM-MMMMMMMMMMMMMMMMMM-MMMMMMMMMMMMMMMMMMMMMMMMMM--MMMMMMMMMMMMM_
        GGGCGGTATGTAG-TGCACACGAG-AATCCCGAAGAATGTTG-TC-CATTCATTGAATAT-GAGATCTCGATGGTATGATCAATATGCGGATGCGATATTGATA-TCTGATAAATCA-TATATGCATAATCTCA-CATTATATTTATTATAAT-AAAT-AGTCGTAG-T-ATACACAATGTGAATTGTATACAATGGAT-AGTATTAACTATCCAATTTCTTT-GAGCATTGGCCTTGGTGTAGATTGCA--TGACATACCGCCC_
        MMMMMMMMMMMMM-MMMMMMMIIM-MMMMMMMMMMMMMMMMM-MM-MMMMMMMMMMMMMM-MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM-MMMMMMMMMMMM-MMMMMMMMMMMMMMMM-MMMMMMMMMMMMMMMMMM-MMMM-MMMMMMMMDM-MMMMMMMMMMMMMMMMMMMMMMMMMMMMM-MMMMIMMMMMMMMMMMMMMMMMM-MMMMMMMMMMMMMMMMMMMMMMMMMM--MMMMMMMMMMMMM_

        """

        # use offset to indicate start, ignore all - and insertions, but keep deletions
        target_pos = position - offset
        target_pos += sum([x == "-" or x == "I" for x in self.cigar[:target_pos]])

        return self.get_seq_pos(target_pos)


class BaseConsensusCaller(ABC):
    @abstractmethod
    def call(self, block):
        pass

    def determine_consensus_direction_backwards(self, block) -> bool:
        """
        Check if the most common direction is backwards
        """
        directions = Counter([x.alignment_direction for x in block.alignments])
        if directions.most_common()[0][0] == "-":
            return True
        else:
            return False

    def _calculate_best_nucleotide(self, nucleotides, probabilities):
        """
        Calculate the probability that a nucleotide is incorrect and return the lowest probability nucleotide with its probability
        """
        nucleotide_likelihood = {}
        for nuc in set(nucleotides):
            nucleotide_likelihood[nuc] = self.calculate_nucleotide_likelihood(
                nuc, nucleotides, probabilities
            )

        best_nuc = min(nucleotide_likelihood, key=nucleotide_likelihood.get)
        best_nuc_probability = nucleotide_likelihood[best_nuc]

        return best_nuc, best_nuc_probability

    def _limit_quality_values(self, quality):
        if quality > self.max_qual:
            quality = self.max_qual
        if quality < self.min_qual:
            quality = self.min_qual
        return quality

    def create_consensus_alignment_objects(self, relevant_alns, start_extender="_"):
        # TODO square this off the matrix to increase performance is last part
        # create a mapping from alignment to ConsensusAlignment
        consensus = {}
        first_start = min([x.alignment_chromosome_start for x in relevant_alns])

        # Create the objects
        for i in relevant_alns:
            start_extension = i.alignment_chromosome_start - first_start

            # add zeros to the quality array
            extended_qual = [0] * start_extension
            extended_qual.extend(i.extended_qual)

            consensus[i] = ConsensusAlignment(
                extended_seq=start_extender * start_extension + i.extended_seq,
                extended_cigar=start_extender * start_extension + i.extended_cigar,
                extended_qual=extended_qual,
                direction=i.alignment_direction,
            )
        return consensus

    def _add_insert_gap(self, alignment_objects, i):
        """
        Add an empty location for all non insert alignments
        """
        return [
            x.add_empty_location(i)
            for x in alignment_objects.values()
            if x.get_cigar_pos(i) != "I"
        ]

    def square_off_alignment(self, alignment_object, desired_length, extender="_"):
        """
        Add _ at the start and end of COnsensus Alignment to make them all the same length
        """
        required_additions = desired_length - len(alignment_object)

        alignment_object.seq += extender * required_additions
        alignment_object.cigar += extender * required_additions
        alignment_object.qual += [0] * required_additions

    def update_insert_locations(self, alignment_objects):
        max_length = max([x.alignment_length for x in alignment_objects.keys()])
        max_insert_n = max(
            [x.extended_cigar.count("I") for x in alignment_objects.keys()]
        )
        maximal_consensus_length = max_length + max_insert_n

        for i in range(maximal_consensus_length):
            cigar = [x.get_cigar_pos(i) for x in alignment_objects.values()]
            if "I" in cigar:
                self._add_insert_gap(alignment_objects, i)

        max_length = max([len(x) for x in alignment_objects.values()])
        # add _ to the ends of the shorter reads
        [
            self.square_off_alignment(aln, max_length)
            for aln in alignment_objects.values()
        ]

    @staticmethod
    def calculate_nucleotide_likelihood(desired_nuc, nucs, probs, start_likelihood=1):
        likelihood = start_likelihood
        for nuc, prob in zip(nucs, probs):
            if nuc == desired_nuc:
                likelihood *= prob
            else:
                if prob == 0:
                    continue
                likelihood /= prob
        return likelihood

    @staticmethod
    def phred_to_prob(phred):
        return 10 ** (-phred / 10)

    def prob_to_phred(self, prob):
        """
        convert a probability to a phred score, takes care of zero and infinite values
        """
        # if we are sure (and python made the float 0), return a high phred
        if prob == 0:
            result = 100
        else:
            try:
                result = int(-10 * math.log10(prob))
            # caused by infinity by large disagreements on a lot of nucleotides
            except OverflowError:
                result = self.min_qual
        return result

    def convert_to_ascii(self, quality, phred_offset=33):
        """
        convert all qualities to its representative base (+33)
        """
        return "".join([chr(x + phred_offset) for x in quality])


class ConsensusCallerMetadata(BaseConsensusCaller):
    """
    Return a consensus sequence and quality in the following format
    backbone0: {
        seq: ACGT,
        qual: ~~~~
        len: 123
        bases before: 10000
        },
    insert0: {
        seq: ACGT,
        qual: ~~~~
        len: 123
        bases before: 10000
    },
    insert1: {
        seq: ACGT,
        qual: ~~~~
        len: 123
        bases before: 10000
    }
    """

    def __init__(self, max_qual=93, minimal_deletion_ratio=0.55):
        logger.debug(f"Consensus caller: {self.__class__.__name__}")
        # lowest ascii char that is human readable, converter takes care of the offset
        self.min_qual = 0
        self.max_qual = max_qual
        self.minimal_deletion_ratio = minimal_deletion_ratio

        if self.max_qual > 93:
            logger.warn(
                """max qual is set to a valuehigher than 93, 
                this might cause invalid fastq outputs"""
            )

    def call(
        self, blocks, minimum_per_block_alignment_count: int = 3
    ) -> Dict[str, Dict[str, Any]]:
        """
        returns a dict with all information about the consensus calling.
        {
            "backbone0": {
            "seq": "ACGT,",
            "qual": "~~~~",
            "len": "123",
            "bases before": "10000",
            },
            "insert0": {
            "seq": "ACGT,",
            "qual": "~~~~",
            "len": "123",
            "bases before": "10000",
            },
        }
        """
        result = {}
        backbone_id = 0
        insert_id = 0

        for block in blocks:
            # Assign tag and increase appropriate counter
            if block.consensus_chromosome.startswith("BB"):
                tag = f"backbone{backbone_id}"
                backbone_id += 1
            else:
                tag = f"insert{insert_id}"
                insert_id += 1
            # filter short blocks
            filter_on_count = False
            if len(block) < minimum_per_block_alignment_count:
                logger.debug(
                    f"rejected consensus creation for {block.read_name} - {block.consensus_chromosome} ({tag}) alns: {len(block)} < required ({minimum_per_block_alignment_count})"
                )
                filter_on_count = True

            # create consensus objects and get values,
            (
                cons,
                qual,
                barcode,
                barcode_arrays,
                alignment_start,
                best_nuc_support,
                flipped,
                consensus_structure,
            ) = self.create_block_consensus(block)

            full_consensus_structure = "|".join(
                [str(len(consensus_structure)), str(len(cons))] + consensus_structure
            )
            # add to the reporting for the output
            result[tag] = {
                "filtered": filter_on_count,
                "seq": cons,
                "qual": qual,
                "len": len(cons),
                "barcode": barcode,
                "barcode_array": barcode_arrays,
                "alignment_count": len(block),
                "aligned_bases_before_consensus": block.count_aligned_bases(),
                "alignment_position": f"{block.consensus_chromosome}:{alignment_start}:{alignment_start+len(cons)}",
                "alignment_orientation": "R" if flipped else "F",
                "original_read_positions": self.generate_block_positions(block),
                "nt_repeat_support": best_nuc_support,
                "consensus_structure": full_consensus_structure,
            }
        return result

    def generate_block_positions(self, block):
        [
            (x.first_cigar_value, x.raw_read_length - x.cigars[-1][1])
            for x in block.alignments
        ]
        output = []
        for x in block.alignments:
            start = x.first_cigar_value
            end = (
                x.raw_read_length - x.cigars[-1][1]
                if x.cigars[-1][0] == "H"
                else x.raw_read_length
            )
            output.append((start, end))
        return output

    def create_block_consensus(self, block):
        """
        Convert a block(grouped alignments) to consensus.
        """
        logger.debug(f"creating consensus using: {self.__class__.__name__}")

        consensus_alignments = self.create_consensus_alignment_objects(block.alignments)
        self.update_insert_locations(consensus_alignments)

        # needs to be 148 for standard cycloseq
        max_length = max([len(x.seq) for x in consensus_alignments.values()])
        # TODO add inserts!
        if max_length > 10000:
            logger.critical(
                f"consensus region is over 10k! Trying to construct consensus over {max_length} nt"
            )

        consensus = ""
        quality = []
        alignment_start = 0
        aligned_segment_count = len(consensus_alignments)
        best_nuc_support = []
        consensus_structure = []
        consensus_structure_seperator = ","
        # loop over all posible positions
        for i in range(max_length):
            nucs = [x.get_seq_pos(i) for x in consensus_alignments.values()]
            quals = [x.get_qual_pos(i) for x in consensus_alignments.values()]
            started_or_ended_count = len([x for x in nucs if x == "_"])
            directions = [x.direction for x in consensus_alignments.values()]

            # remove non-participating reads, check if we have enough participation to make consensus
            filtered_nucs = []
            filtered_quals = []
            filtered_directions = []

            for n, q, d in zip(nucs, quals, directions):
                if n != "_":
                    filtered_nucs.append(n)
                    filtered_quals.append(q)
                    filtered_directions.append(d)

            reprs = []

            for n, q, d in zip(filtered_nucs, filtered_quals, filtered_directions):
                if n == "-":
                    base = "D"
                else:
                    base = n.capitalize() if d == "+" else n.lower()
                string_repr = f"{base}{q}"
                reprs.append(string_repr)

            nucs = filtered_nucs
            quals = filtered_quals

            # not more than 2 nucs? or at least 10% of the reads, rounded up
            if (
                len(nucs) < 2
                or math.ceil(aligned_segment_count * 0.1) < started_or_ended_count
            ):
                consensus_structure_entry = consensus_structure_seperator.join(
                    ["D" + str(len(filtered_nucs))] + reprs
                )
                consensus_structure.append(consensus_structure_entry)
                logger.debug("skipping consensus base creation due to low coverage")
                continue
            logger.debug("starting consensus creation")

            # if we detect a real delition we dont have to do anything
            deletion_ratio = nucs.count("-") / len(nucs)
            if deletion_ratio >= self.minimal_deletion_ratio:

                consensus_structure_entry = consensus_structure_seperator.join(
                    ["D" + str(len(filtered_nucs))] + reprs
                )
                consensus_structure.append(consensus_structure_entry)
                logger.debug(
                    f"skipping consensus base creation due to deletion ratio {deletion_ratio}"
                )
                continue

            filtered_nucs = []
            filtered_quals = []
            # Remove deltions from the descision making
            for n, q in zip(nucs, quals):
                if n != "-":
                    filtered_nucs.append(n)
                    filtered_quals.append(q)

            nucs = filtered_nucs
            quals = filtered_quals
            # now we can calculate probs
            probs = [self.phred_to_prob(x) for x in quals]

            best_nuc, best_nuc_prob = self._calculate_best_nucleotide(nucs, probs)
            best_nuc_support.append(tuple([nucs.count(best_nuc), len(nucs)]))

            # If we start the consensus string we store the position
            if alignment_start == 0 and best_nuc != "-":
                alignment_start = i + block.get_start_position()

            # TODO: check if two nucleotides have the exact same likelihood.
            determined_qual = self.prob_to_phred(best_nuc_prob)

            # limit the quality
            determined_qual = self._limit_quality_values(determined_qual)

            consensus += best_nuc
            quality.append(determined_qual)
            consensus_structure_entry = consensus_structure_seperator.join(
                [best_nuc + str(len(filtered_nucs))] + reprs
            )
            consensus_structure.append(consensus_structure_entry)
            pass
        # get barcode, we need to do this here since we need positional awareness within the barcode alignment
        barcode_extractor = BarcodeExtrator(
            alignment_blocks=consensus_alignments, offset=block.get_start_position()
        )

        barcode, barcode_arrays = barcode_extractor.extract_barcode()
        logger.debug("barcode extractor made")

        # flip the consensus and quality to reflect real directionality of the read if needed
        flipped = False
        if self.determine_consensus_direction_backwards(block):
            consensus = str(Seq(consensus).reverse_complement())
            quality = quality[::-1]
            best_nuc_support = best_nuc_support[::-1]
            consensus_structure = consensus_structure[::-1]
            flipped = True

        return (
            consensus,
            self.convert_to_ascii(quality),
            barcode,
            barcode_arrays,
            alignment_start,
            best_nuc_support,
            flipped,
            consensus_structure,
        )
