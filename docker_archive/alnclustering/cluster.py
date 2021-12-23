import sys, os, glob
import pysam
from Bio.Cluster import kcluster 
import numpy as np 
import argparse

def parse_fasta(fasta_file):
    """file_path => dict
    Return a dict of id:sequences.
    """
    d = {}
    _id = False
    seq = ""
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith("\n"):
                continue
            if line.startswith(">"):
                if not _id:
                    _id = line.strip()[1:]
                elif _id and seq:
                    d.update({_id: seq})
                    _id = line.strip()[1:]
                    seq = ""
            else:
                seq += line.strip()
        d.update({_id: seq})
    return d

def main(msa_file ,output1,output2, cluster):
    sequences = parse_fasta(msa_file).values()
    matrix = np.asarray([np.frombuffer(s.encode(), dtype=np.uint8) for s in sequences]) 
    clusterid, error, found = kcluster(matrix, nclusters=cluster)
    # print(clusterid) #[1 0 0 1]

    with open(output1, 'w') as hap0:
        with open(output2, 'w') as hap1:
            for seq, clust_id in zip(parse_fasta(msa_file).items(), clusterid):
                target = hap1 if clust_id ==1 else hap0
                # remove all - and
                clean_seq = seq[1].replace('-','').replace(">", "")
                target.write(f">{seq[0]}\n{clean_seq}\n")


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Cluster MSA reads into n groups, writes files 0 and 1')
    parser.add_argument('-i','--input_msa', type=str)
    parser.add_argument('-o','--output_1', type=str)
    parser.add_argument('-O','--output_2', type=str)
    parser.add_argument('-c','--clusters', type=int)
    
    args = parser.parse_args()
    main(args.input_msa, args.output_1, args.output_2, args.clusters)