#!/usr/bin/env python

import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import List
from collections import Counter 

import mappy
import parasail

NAME = str
SEQ = str
QUAL = str or None
FX_READ = (NAME,SEQ,QUAL,str) 

@dataclass
class Backbone:
    name:str
    seq:str

    def __repr__(self) -> str:
        return f"Backbone {self.name}"

    def __name__(self) -> str:
        return self.name

    def map(self, sequence, show_aln=True):
        mapping = parasail.sg_trace(sequence,self.seq, open=2, extend =1, matrix= parasail.dnafull)
        if show_aln:
            print(f"score: {mapping.score}")
            print(mapping.traceback.ref)
            print(mapping.traceback.comp)
            print(mapping.traceback.query)

        return mapping

    def check_for_barcode(self, sequence, min_score=500, show_aln=False):
        
        # Parasail does not automatically check the sense and anti-sense strand, thus we do it manually.
        mapping = self.map(sequence, show_aln=show_aln)
        sequence_rc = sequence[::-1].translate(str.maketrans("ATCG", "TAGC"))
        mapping_rc = self.map(sequence_rc, show_aln=show_aln)
        best_mapping = max(mapping, mapping_rc, key=lambda x: x.score)

        # If the best mapping is mapping poorly, no need to continue.
        if best_mapping.score < min_score:
            return (best_mapping.score, "Unknown")
        
        # If the best barcode contains N bases, we extract these from the mapping and report them as the name
        if self.has_ambiguous_bases:
            ambigous_positions = self.ambigous_positions_from_other(best_mapping.traceback.ref)
            # Look up the barcode for each position that we have an n.
            barcode = "".join([best_mapping.traceback.query[x] if best_mapping.traceback.comp[x] in [".","|"] else "N" for x in ambigous_positions])
        else:
            barcode = ""

        return (best_mapping.score, barcode)

    @property
    def has_ambiguous_bases(self) -> bool:
        return "N" in self.seq
    
    def ambigous_positions(self) -> List[int]:
        return [i for i, char in enumerate(self.seq) if char == "N"]

    def ambigous_positions_from_other(self, other) -> List[int]:
        """
        Find the indexes for the N's in other sequences
        """
        return [i for i, char in enumerate(other) if char == "N"]


@dataclass
class BackboneCollection:
    file_path: Path
    

    def __post_init__(self):
        # print("initializing Backbone class")
        
        if not self.file_path.exists():
            print("backbone file does not exist.")
            exit(1)
        read_gen = mappy.fastx_read(str(self.file_path), read_comment=True)
        # extract the backbone from the file, check that its alone in the file.
        self.backbones = []        
        for r in read_gen:
            read = r
            self.backbones.append(Backbone(read[0], read[1]))
        # lookup dict based on the name
        self.backbone_lookup = {x.name: x for x in self.backbones}


    def demux_read(self, read:FX_READ, length_err_preference_factor:float = 0.5):
        # print("demuxing read")
        results = {}
        for bb in self.backbones:
            mapped = list(bb.check_for_barcode(read[1]))
            results[bb.name] = mapped if mapped[0] > 0 else None

        best_result = max(results, key=lambda x: results[x][0])
        best_backbone = self.backbone_lookup[best_result]
        if results[best_result][1] == "Unknown":
            classification = "Unknown"
        else:
            if best_backbone.has_ambiguous_bases:
                classification = f"{best_result}_{results[best_result][1]}"
            else:
                classification = f"{best_result}"
        # print(classification)
        return classification

    def get_barcode_from_cs(self,read_mapped, backbone) -> str:
        ambigous = backbone.ambigous_positions()
        seq = self.seq_from_cs(backbone.seq, read_mapped)
        barcode = [seq[x] for x in ambigous]
        barcode = "".join(barcode)
        return barcode

    def seq_from_cs(self, ref_seq, read_mapped):
        cs_elements = re.findall(r'[:*\+\-~]?\w+', read_mapped.cs)
        pos = read_mapped.r_st
        seq = 'N' * pos
        # print(f"starting at {pos}")
        for elem in cs_elements:
            # print(elem)
            if elem[0] == ':':
                amount = int(elem[1:])
                seq += ref_seq[pos:pos+amount]
                pos += amount 
            elif elem[0] == '*':
                pos += 1
                seq += elem[-1].upper()
            elif elem[0] == '-':
                amount = len(elem[1:])
                seq += "N" * amount
                pos += amount
            elif elem[0] == '+':
                pass

            # print(seq)
        seq += "N" * (read_mapped.ctg_len - len(seq))
        return seq
   

def process_cyclomics_read_name(read_name:str) -> str:
    """
    Get the universal part of the read name from the sequencing read,
    this is neccecery since some consensus algorithms add info to the name.
    """

    return read_name.split("_")[0]

def main(backbone_fasta_path: Path,
         read_fastq_path: Path,
         write_fastq_path: Path,
         output_by_readname: bool,
         multi_read:bool):
    # check if input fastq for demuxing exists
    if not read_fastq_path.exists():
        print(f"{read_fastq_path} doest not exist. exiting...")
        exit(1)
    # create and check the output folder
    write_fastq_path.mkdir(exist_ok=True)
    if len([x for x in write_fastq_path.glob('*')]) != 0:
        print(f"{write_fastq_path.resolve()} exists and is not empty. exiting...")
        exit(1)

    # create the objects neccessary for demuxing
    backbones = BackboneCollection(backbone_fasta_path)
    reads = mappy.fastx_read(str(read_fastq_path), read_comment=True)
    demux_results = []
    # demux each read.
    for  i, read in enumerate(reads):
        # print(read[0])
        barcode = backbones.demux_read(read)
        # print(barcode)
        demux_results.append((barcode, read))


    demux_results_multi_read = {}
    if multi_read:
        # First check each outcome so that we have a mapper for the original readname
        for result in demux_results:
            read = result[1]
            conserved_read_name = process_cyclomics_read_name(read[0])
            conserved_read_name_in_results = conserved_read_name in demux_results_multi_read.keys()
            # Overwrite the unknown barcodes.
            if conserved_read_name_in_results and demux_results_multi_read[conserved_read_name] == "Unknown":
                demux_results_multi_read[conserved_read_name] = result[0]
            # Barcode collision
            elif conserved_read_name_in_results:
                demux_results_multi_read[conserved_read_name] = "Collision"
            # Simple addition 
            else:
                demux_results_multi_read[conserved_read_name] = result[0]
        new_results = []
        # Than update all readnames with the mapper
        for result in demux_results:
            conserved_read_name = process_cyclomics_read_name(result[1][0])
            new_bc = demux_results_multi_read[conserved_read_name]
            new_results.append((new_bc, result[1]))
        demux_results = new_results
    

    print(Counter([x[0] for x in demux_results]).most_common())
    
    open_files = {}
    for result in demux_results:
        if result[0] not in open_files:
            open_files[result[0]] = (write_fastq_path / (result[0] + ".fastq")).open(mode="w")
        read = f"@{result[1][0]} {result[1][3]}\n{result[1][1]}\n+\n{result[1][2]}\n"
        open_files[result[0]].write(read)

    [x.close() for x in open_files.values()]

    print('DOne')

if __name__ == "__main__":
    import argparse

    dev = True
    if not dev:
        parser = argparse.ArgumentParser(
            description=("Detect indels in BAM alignment file.")
        )

        parser.add_argument("backbone_fasta", type=Path)
        parser.add_argument("read_fastq", type=Path)
        parser.add_argument("by_name", type=bool, default=False)
        args = parser.parse_args()

        main(args.backbone_fasta, args.read_fastq, args.by_name)

    if dev:
        # EGFR
        backbone_fasta = Path(
            "/home/dami/Software/0220/cyclomicsseq/backbones/BB42_dummy_bc.fasta"
        )
        read_fastq = Path(
            "/home/dami/Software/0220/cyclomicsseq/test4_000128/consensus/FAW01764_pass_fc459b5d_61ee14c3_123_filtered.consensus.fastq"
        )
        barcode_folder = Path("./barcode_demux")

        main(backbone_fasta, read_fastq, barcode_folder,False,True)

