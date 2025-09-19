"""Module on read simulation, for documentation on all the different parameters read
the main function `simulate_reads`.
"""

import os
import uuid
from copy import deepcopy
from typing import Any, Dict, List, Tuple
from warnings import warn

import numpy as np
from tqdm import tqdm

UUID_NAMESPACE = uuid.UUID("12345678-1234-5678-1234-567812345678")
DNA = {"A", "C", "G", "T"}
INT_TO_DNA = {0: "A", 1: "C", 2: "G", 3: "T"}
DNA_TO_INT = {v: k for k, v in INT_TO_DNA.items()}
DNA_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
DEBUG = False


def printdbg(*args):
    if DEBUG:
        print("DEBUG: ", *args)


def simulate_reads(
    reference_file: str,
    output_dir: str,
    num_reads: int,
    read_type: List[str],
    insert_length: Tuple[int, int],
    num_repeats: Tuple[int, int],
    error_rate_mismatch: float,
    error_rate_deletion: float,
    error_rate_insertion: float,
    tdt_length: Tuple[int, int],
    num_template_switches: Tuple[int, int],
    num_pore_chimeras: Tuple[int, int],
    num_library_chimeras: Tuple[int, int],
    jagged_end_left: Tuple[int, int],
    jagged_end_right: Tuple[int, int],
    seed: int,
):
    """Main function to simulate Cyclomics Reads from a reference genome

    Parameters that have a range (min and max), values will be sampled uniformly from
    the given range with both values included.

    Error rates are independent, so consider the overall error profile when setting them.

    To not add TDT, template switching, or chimeras set both values to 0.

    Args:
        reference_file: path to a reference genome file to simulate data from. A
            position is selected at random, if the selected DNA contains lowercase
            letters, indicating a repetitive region another position is selected.
        output_dir: path to an output directory to save the simulated reads to.
        read_type: type of read to simulate (1D, 3D).
        num_reads: number of reads to simulate.
        insert_length: min and max insert length
        num_repeats: min and max number of repeats
        error_rate_mismatch: rate of mismatch errors.
        error_rate_deletion: rate of deletion errors.
        error_rate_insertion: rate of insertion errors.
        tdt_length: min and max random TDT length.
        num_template_switches: min and max number of template switches.
        num_pore_chimeras: min and max number of chimeras due to poor pore splitting.
        num_library_chimeras: min and max number of chimeras due to library preparation.
        jagged_end_left: min and max jagged end left length.
        jagged_end_right: min and max jagged end right length.
        seed: random seed for reproducibility.
    """

    # Assertion checks for input files and output directory
    assert os.path.isfile(
        reference_file
    ), f"Reference file '{reference_file}' does not exist or is not a file."
    assert os.access(
        reference_file, os.R_OK
    ), f"Reference file '{reference_file}' is not readable."
    assert os.path.isdir(
        output_dir
    ), f"Output directory '{output_dir}' does not exist or is not a directory."
    assert os.access(
        output_dir, os.W_OK
    ), f"Output directory '{output_dir}' is not writable."

    reference = read_reference_file(reference_file)

    output_file_fastq = os.path.join(output_dir, "simulated_reads.fastq")
    # output_file_read_metadata = os.path.join(output_dir, 'simulated_reads_metadata.csv')
    # output_file_insert_metadata = os.path.join(output_dir, 'simulated_inserts_metadata.csv')
    output_file_fastq_con = open(output_file_fastq, "w")

    total_valid_reads = 0
    for _ in tqdm(range(num_reads), disable=DEBUG):
        valid_read, fastq_string, read_metadata = simulate_read(
            reference=reference,
            read_type=read_type,
            insert_length=insert_length,
            num_repeats=num_repeats,
            error_rate_mismatch=error_rate_mismatch,
            error_rate_deletion=error_rate_deletion,
            error_rate_insertion=error_rate_insertion,
            tdt_length=tdt_length,
            num_template_switches=num_template_switches,
            num_pore_chimeras=num_pore_chimeras,
            num_library_chimeras=num_library_chimeras,
            jagged_end_left=jagged_end_left,
            jagged_end_right=jagged_end_right,
            seed=seed,
        )
        seed += 1

        if not valid_read:
            continue

        total_valid_reads += 1

        output_file_fastq_con.write(fastq_string)

    output_file_fastq_con.close()

    print(f"Total valid reads simulated: {total_valid_reads}/{num_reads}")


def simulate_read(
    reference: Dict[str, str],
    read_type: List[str],
    insert_length: Tuple[int, int],
    num_repeats: Tuple[int, int],
    error_rate_mismatch: float,
    error_rate_deletion: float,
    error_rate_insertion: float,
    tdt_length: Tuple[int, int],
    num_template_switches: Tuple[int, int],
    num_pore_chimeras: Tuple[int, int],
    num_library_chimeras: Tuple[int, int],
    jagged_end_left: Tuple[int, int],
    jagged_end_right: Tuple[int, int],
    seed: int,
) -> Tuple[bool, str, Dict[str, Any]]:
    """Simulate a single read.

    For info on the args check `simulate_reads`

    Returns:
        A tuple with:
        - A bool, if the read is valid or not
        - A string with line breaks, ready to be written to a fastq file, so it contains
        the read_id, sequence, direction and quality scores.
        - A list of dictionaries with the insert metadata, ready to be converted to a
        pandas data.frame. If there are no chimeras, then the list is of length 1,
        otherwise of length the number of chimeras sampled.
    """

    rng = np.random.default_rng(seed=seed)
    read_id = uuid.uuid5(UUID_NAMESPACE, str(seed))

    printdbg(f"Read id: {read_id}")

    # inserts can be 1D or 3D, since 2D refers to the read
    possible_insert_types = []
    for it in ["1D", "3D"]:
        if it in read_type:
            possible_insert_types.append(it)

    # Sample the number of independent reads that will be joined together
    # into a single read, simulating reads that are not correctly split by the
    # sequencer.
    num_pore_chimeras = rng.integers(
        low=num_pore_chimeras[0],
        high=num_pore_chimeras[1],
        endpoint=True,
    )

    printdbg(f"Num pore chimeras: {num_pore_chimeras}")

    read_metadata = {
        "read_id": str(read_id),
        "length": None,
        "num_pore_chimeras": num_pore_chimeras,
        "num_inserts": 0,
        "seed": seed,
    }

    insert_index = 0
    final_read = ""
    # pc_idx = pore_chimera_index
    read_insert_metadata_rows = list()
    for pc_idx in range(num_pore_chimeras + 1):
        num_library_chimeras = rng.integers(
            num_library_chimeras_min, num_library_chimeras_max, endpoint=True
        )
        read_metadata["num_inserts"] += num_library_chimeras + 1

        printdbg(f"Num library chimeras: {num_library_chimeras}")

        insert_metadata_rows = list()
        # lc_idx = library_chimera_index
        for lc_idx in range(num_library_chimeras + 1):
            read_type = rng.choice(possible_insert_types)

            printdbg(f"Read type: {read_type}")

            # establish metadata of each insert
            insert_metadata = {
                "read_id": str(read_id),
                "insert_id": insert_index,
                "read_type": read_type,
                "mapping_fwd": "",
                "seq_fwd": "",
                "tdt_fwd": "",
                "mapping_rev": "",
                "seq_rev": "",
                "tdt_rev": "",
                "template_switches": 0,
            }

            seq, contig, start, end, reverse = random_genomic_sequence(
                reference=reference,
                length_min=insert_length_min,
                length_max=insert_length_max,
                rng=rng,
            )
            direction = "rev" if reverse else "fwd"

            printdbg(f"Mapping: {contig}:{start}-{end}")

            if read_type == "1D":
                mapping = f"{contig}:{start}-{end}"
                insert_metadata[f"mapping_{direction}"] = mapping
                insert_metadata[f"seq_{direction}"] = seq
                insert_metadata[f"tdt_{direction}"] = tdt(
                    length_min=tdt_length_min, length_max=tdt_length_max, rng=rng
                )

            elif read_type == "3D":
                if jagged_end_max_left > 0:
                    jagged_end_left = rng.integers(
                        jagged_end_min_left, jagged_end_max_left
                    )
                else:
                    jagged_end_left = 0

                if jagged_end_max_right > 0:
                    jagged_end_right = rng.integers(
                        jagged_end_min_right, jagged_end_max_right
                    )
                else:
                    jagged_end_right = 0

                st_fwd = start
                nd_fwd = end
                st_rev = start
                nd_rev = end

                if rng.choice([0, 1]) == 1:
                    st_rev += jagged_end_left
                else:
                    st_fwd += jagged_end_left

                if rng.choice([0, 1]) == 1:
                    nd_rev -= jagged_end_right
                else:
                    nd_fwd -= jagged_end_right

                if direction == "fwd":
                    insert_metadata["mapping_fwd"] = f"{contig}:{st_fwd}-{nd_fwd}"
                    insert_metadata["mapping_rev"] = f"{contig}:{st_rev}-{nd_rev}"
                    insert_metadata["seq_fwd"] = reference[contig][st_fwd:nd_fwd]
                    insert_metadata["seq_rev"] = reverse_complement(
                        reference[contig][st_rev:nd_rev]
                    )
                else:
                    insert_metadata["mapping_rev"] = f"{contig}:{st_fwd}-{nd_fwd}"
                    insert_metadata["mapping_fwd"] = f"{contig}:{st_rev}-{nd_rev}"
                    insert_metadata["seq_rev"] = reverse_complement(
                        reference[contig][st_fwd:nd_fwd]
                    )
                    insert_metadata["seq_fwd"] = reference[contig][st_rev:nd_rev]

                insert_metadata["tdt_fwd"] = tdt(
                    length_min=tdt_length_min, length_max=tdt_length_max, rng=rng
                )
                insert_metadata["tdt_rev"] = tdt(
                    length_min=tdt_length_min, length_max=tdt_length_max, rng=rng
                )

            else:
                raise ValueError("Insert type can only be 1D or 3D")

            insert_metadata_rows.append(insert_metadata)
            insert_index += 1

        # concatenate inserts into circle
        final_insert = list()
        inserts_3d = list()
        for row in insert_metadata_rows:
            if row["read_type"] == "1D":
                if row["seq_fwd"] != "":
                    final_insert.append(row["seq_fwd"] + row["tdt_fwd"])
                else:
                    final_insert.append(row["seq_rev"] + row["tdt_rev"])
            else:
                final_insert.append(row["seq_fwd"] + row["tdt_fwd"])
                inserts_3d.append((row["seq_rev"] + row["tdt_rev"]))

        inserts_3d = inserts_3d[::-1]
        final_insert = "".join(final_insert + inserts_3d)
        insert_len = len(final_insert)
        printdbg("Final insert length: ", insert_len)
        # RCA
        ## repeat
        num_repeats = rng.integers(num_repeats_min, num_repeats_max, endpoint=True)
        printdbg(f"Num repeats: {num_repeats}")

        final_insert = final_insert * num_repeats
        printdbg("Final insert length after repeats: ", len(final_insert))
        ## chunk
        right_chop = rng.integers(0, insert_len // 2)
        printdbg("Right chopped bases: ", right_chop)
        left_chop = rng.integers(0, insert_len // 2)
        printdbg("Left chopped bases: ", left_chop)
        final_insert = final_insert[left_chop:-right_chop]

        printdbg("Final insert length after chunking: ", len(final_insert))

        num_template_switches = rng.integers(
            num_template_switches_min, num_template_switches_max, endpoint=True
        )
        current_template = deepcopy(final_insert)
        for nts in range(num_template_switches):
            if len(current_template) <= insert_len * 2:
                break

            template_copy_len = rng.integers(insert_len * 2, len(current_template))
            current_template = reverse_complement(current_template)[:template_copy_len]
            final_insert += current_template

        # concatenate with other pore chimeras
        final_read += final_insert

        for row in insert_metadata_rows:
            row["template_switches"] = num_template_switches
            read_insert_metadata_rows.append(row)

    # final_read = add_errors(final_read)

    valid_read = True
    if len(final_read) == 0:
        valid_read = False
        warn(
            "At least one simulated read has a length of zero, probably insert length is too short, or too few repeats"
        )

    final_qual = "&" * len(final_read)
    # make string
    fastq_string = f"@{read_id}\n{final_read}\n+\n{final_qual}\n"
    read_metadata["length"] = len(final_read)

    return valid_read, fastq_string, read_metadata, read_insert_metadata_rows


def random_genomic_sequence(
    reference: Dict[str, str],
    length_min: int,
    length_max: int,
    rng: np.random.Generator,
) -> Tuple[str, str, int, int, bool]:
    """Get a random genomic sequence from a reference genome. Retry if there are
    non ACGT bases or lowercase.
    """

    contigs = [f"chr{n}" for n in range(1, 23)]
    contig = rng.choice(contigs)
    # contig = rng.choice(list(reference.keys()))
    contig_length = len(reference[contig])

    start = rng.integers(0, contig_length - length_min)
    length = rng.integers(length_min, length_max)
    reverse = bool(rng.integers(0, 2))
    end = start + length

    seq = reference[contig][start:end]
    for base in seq:
        if base not in DNA:
            return random_genomic_sequence(reference, length_min, length_max, rng)

    if reverse:
        seq = reverse_complement(seq)

    return seq, contig, start, end, reverse


def tdt(length_min: int, length_max: int, rng: np.random.Generator) -> str:
    if length_min == 0 and length_max == 0:
        return ""

    size = rng.integers(low=length_min, high=length_max, size=1, endpoint=True)
    if size == 0:
        return ""

    dna_ints = rng.integers(low=0, high=3, size=size, endpoint=True)

    return "".join([INT_TO_DNA[i] for i in dna_ints])


def reverse_complement(seq: str) -> str:
    """Reverse complement a sequence"""

    return "".join(DNA_COMPLEMENT[base] for base in seq[::-1])


def read_reference_file(file: str) -> Dict[str, str]:
    reference = dict()
    with open(file, "r") as f:
        for line in f:
            if line.startswith(">"):
                contig = line[1:].split(" ")[0].strip()
                reference[contig] = list()
            else:
                reference[contig].append(line.strip())

    for contig in reference:
        reference[contig] = "".join(reference[contig])
    return reference
