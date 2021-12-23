#!/usr/bin/env python3
from dataclasses import dataclass, field
from pathlib import Path
import argparse

@dataclass
class Cluster:
    name: str
    read_names: list[str] = field(default_factory=list)
    reads: dict[str,str] = field(default_factory=dict)

    def __repr__(self):
        return self.name
    
    def __hash__(self) -> int:
        return hash(self.name)

    
def parse_fasta(fasta_file):
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

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Split clusters into seperate fastas')
    parser.add_argument('-i','--input_fasta', type=str)
    parser.add_argument('-s','--sample', type=str)
    parser.add_argument('-m','--min_reads', type=int, default = 10)
    parser.add_argument('-r','--read_ratio', type=float, default = 0.25)

    args = parser.parse_args()

    clusters = []
    sample = args.sample
    fa_file = args.input_fasta
    clstr_file = f"{sample}.clstr"
    min_reads_cluster = args.min_reads
    read_ratio =    args.read_ratio

    print(sample)
    print(fa_file)
    print(clstr_file)

    with open(clstr_file, "r") as f:
    # initialize the first cluster
        cur_clstr = Cluster(next(f).strip())
        for line in f.readlines():
            if line.startswith(">"):
                clusters.append(cur_clstr)
                cur_clstr = Cluster(line.strip())
            else:
                read_name = line.split(",")[1].split("...")[0].strip()
                cur_clstr.read_names.append(read_name)
    print(f"found clusters: {len(clusters)}")

    fa_reads = parse_fasta(fa_file)

    for read_name, read_seq in fa_reads.items():
        for cluster in clusters:
            search_name = f">{read_name[:19]}"
            if search_name in cluster.read_names:
                cluster.reads[read_name] = read_seq
                break

    total_reads = 0
    for clst in clusters:
        total_reads += len(clst.read_names)
    print(f"clustered reads: {total_reads}")

    filtered_clusters = [x for x in clusters if len(x.read_names) > min_reads_cluster]
    print(f"found clusters: {len(filtered_clusters)}")
    print([x.name for x in filtered_clusters])
    min_th_clusters = [x for x in clusters if len(x.read_names) > total_reads * read_ratio]
    print(f"found clusters th: {len(min_th_clusters)}")
    print([x.name for x in min_th_clusters])

    clusters_accepted = set(filtered_clusters + min_th_clusters)
    clusters_rejected = set(clusters).difference(clusters_accepted)
    print(f"accepted clusters: {len(clusters_accepted)}")
    print(f"rejected clusters: {len(clusters_rejected)}")


    for cluster in clusters_accepted:
        file = Path(fa_file).with_suffix(f".C_{str(cluster).split(' ')[-1]}.fa")
        with open(file, "w") as f:
            print(f"{cluster} found {len(cluster.reads)} read, writing into {file}")
            for read_name, read_seq in cluster.reads.items():
                f.write(f">{read_name}\n{read_seq}\n")

    for cluster in clusters_rejected:
        file = Path(fa_file).with_suffix(f".C_rejected.fa")
        with open(file, "w") as f:
            for read_name, read_seq in cluster.reads.items():
                f.write(f">{read_name}\n{read_seq}\n")