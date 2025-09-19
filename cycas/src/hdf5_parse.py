import json
from typing import Dict, List

import h5py

from .alignment import Alignment


class Hdf5Writer:
    def __init__(self, filepath) -> None:
        self.file_path = filepath
        self.hdf5 = h5py.File(self.file_path, "w")

        # only include single insert reads
        # possible cases
        #   BackboneDoubleInsert
        #   BackboneInsert
        #   DoubleInsertUncertain
        #   LowAlignmentCount
        #   SingleBackbone
        #   SingleBackboneUnalignedGaps
        #   SingleInsert
        #   SingleInsertUnalignedGaps
        #   SingleInsertUncertain'
        self.allowed_types = {"SingleInsert", "BackboneInsert"}

    def close(self):
        self.hdf5.close()

    def save_raw_blocks(self, metadata: Dict, block: List[Alignment]) -> None:
        """Saves raw alignment blocks of the same locus in the style of MSA into an HDF5 file.

        Args:
            metadata (Dict): Dictionary with obtained metadata, including classification and metrics.
            block (List[Alignment]): List of Alignment objects representing the MSA.

        Note:
            Writes data to an HDF5 file.
        """

        if "structure_classification" not in metadata:
            return None

        read_uuid = metadata["readname"]
        grp = self.hdf5.create_group(read_uuid)
        grp.attrs["readname"] = read_uuid

        for k, v in metadata.items():
            try:
                grp.attrs[k] = v
            except TypeError:
                grp.attrs[k] = json.dumps(v)

        block_ref = block[0].alignment_chromosome
        block_first_pos = min([x.alignment_chromosome_start for x in block])

        grp.attrs["ref"] = block_ref
        grp.attrs["first_pos"] = block_first_pos

        block_data = []
        for aln in block:
            sam = aln.pysam_obj
            fwd_direction = 1 if sam.is_forward else 0
            # Group by key name
            nuc_array = list(aln.extended_seq)  # ['A', 'A', 'C', 'A', ....
            cigar_array = list(
                aln.extended_cigar
            )  # ['_', '_', 'M', 'M', 'D', '-', ....
            qual_array = list(aln.extended_qual)  # [ 40, 23, 12, 5, ....

            # tuple of mapq, sam_flag, fwd, read_start_position, query_start, query_end, query_start, query_end
            aln_data = (
                sam.mapping_quality,
                sam.flag,
                fwd_direction,
                aln.first_cigar_value,
                sam.query_alignment_start,
                sam.query_alignment_end,
                block_ref,
                sam.reference_start,
                sam.reference_end,
            )

            block_data.append([nuc_array, cigar_array, qual_array, aln_data])

        for i, (nuc_array, cigar_array, qual_array, aln_data) in enumerate(block_data):
            data_group = grp.create_group(f"data_{i}")
            data_group.create_dataset("nuc_array", data=nuc_array, dtype="S1")
            data_group.create_dataset("cigar_array", data=cigar_array, dtype="S1")
            data_group.create_dataset("qual_array", data=qual_array, dtype="int8")
            data_group.create_dataset("aln_data", data=aln_data, dtype="S128")

    def save_consensus_blocks(self, metadata, blocks):
        if "classification" not in metadata:
            return None

        if metadata["classification"] not in self.allowed_types:
            return None

        read_uuid = list(blocks[0].keys())[0].readname

        grp = self.hdf5.create_group(read_uuid)

        grp_attrs = metadata.copy()
        [grp_attrs.pop(x) for x in ("rules", "consensus_reads")]
        for k, v in grp_attrs.items():
            grp.attrs[k] = v

        for (i, block), (consensus_key, consensus_data) in zip(
            enumerate(blocks), metadata["consensus_reads"].items()
        ):
            if consensus_data["filtered"]:
                continue
            if "readname_unique" not in consensus_data:
                continue
            if len(block) == 0:
                continue

            block_ref = list(block.keys())[0].alignment_chromosome
            block_first_pos = min([x.alignment_chromosome_start for x in block.keys()])

            block_data = []
            for aln, spaced_aln in block.items():
                sam = aln.pysam_obj
                fwd_direction = 1 if sam.is_forward else 0
                # Group by key name
                nuc_array = list(spaced_aln.seq)  # ['A', 'A', 'C', 'A', ....
                cigar_array = list(
                    spaced_aln.cigar
                )  # ['_', '_', 'M', 'M', 'D', '-', ....
                qual_array = spaced_aln.qual  # [ 40, 23, 12, 5, ....

                # tuple of mapq, sam_flag, fwd, read_start_position, query_start, query_end, query_start, query_end
                aln_data = (
                    sam.mapping_quality,
                    sam.flag,
                    fwd_direction,
                    aln.first_cigar_value,
                    sam.query_alignment_start,
                    sam.query_alignment_end,
                    block_ref,
                    sam.reference_start,
                    sam.reference_end,
                )

                block_data.append([nuc_array, cigar_array, qual_array, aln_data])

            group_name = f"block_{i}"
            # make each block a subgroup of the read
            subgrp = grp.create_group(group_name)

            # store block consensus data:
            subgrp.attrs["consensus_seq"] = consensus_data["seq"]
            subgrp.attrs["consensus_qual"] = consensus_data["qual"]
            subgrp.attrs["consensus_alignment_position"] = consensus_data[
                "alignment_position"
            ]
            subgrp.attrs["consensus_readname_unique"] = consensus_data[
                "readname_unique"
            ]
            subgrp.attrs["consensus_orientation"] = consensus_data[
                "alignment_orientation"
            ]

            # Store block attributes
            subgrp.attrs["ref"] = block_ref
            subgrp.attrs["first_pos"] = block_first_pos

            # Store block data as datasets
            for j, (nuc_array, cigar_array, qual_array, aln_data) in enumerate(
                block_data
            ):
                data_group = subgrp.create_group(f"data_{j}")
                data_group.create_dataset("nuc_array", data=nuc_array, dtype="S1")
                data_group.create_dataset("cigar_array", data=cigar_array, dtype="S1")
                data_group.create_dataset("qual_array", data=qual_array, dtype="int8")
                data_group.create_dataset("aln_data", data=aln_data, dtype="S128")
