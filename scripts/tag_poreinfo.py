
import pysam
import sys

def load_seqsum(seqsum_path):
    print('Loading sequencing summary')
    with open(seqsum_path,'r') as seqsum:
        header = seqsum.readline()
        for i,name in enumerate(header.split()):
            print('\t',i,name)
        lines = [line.split() for line in seqsum.readlines()]
        readnames = [x[2] for x in lines]
    print('- Done')
    return header,dict(zip(readnames, lines))



def read_getter_bam_f5(seqsum_path, in_bam_path, out_bam_path):
    # Initialize file handles
    in_bam_file   = pysam.AlignmentFile(in_bam_path, "rb")
    hss, seqsum   = load_seqsum(seqsum_path)
    out_bam_file   = pysam.AlignmentFile(out_bam_path+'.bam', "wb", header=in_bam_file.header)


    for aln in in_bam_file.fetch():
        if aln.query_name in seqsum:
            readsum = seqsum[aln.query_name]
            # print(readsum)
            aln.tags = aln.tags+[
                ('CH',int(readsum[4])), # CHannel
                ('MU',int(readsum[5])), # MUx
                ('CM',readsum[4]+'.'+readsum[5]), # Channel Mux
                ('ST',float(readsum[6])), # Start Time
                ('DU',float(readsum[7])), # DUration
                ('AD',float(readsum[10])-float(readsum[6])), # Adapter Duration?
                ('SL',int(readsum[13])), # Sequence_Length_template
                ('MQ',float(readsum[14])), # Mean_Qscore_template
                ('MT',float(readsum[16])), # Median_Template
                ('MA',float(readsum[17])), # Mad_Template
                ('ER',readsum[21]) # End_Reason
            ]
            # print(aln.query_name,aln.tags)
            out_bam_file.write(aln)
        # exit()

    in_bam_file.close()
    out_bam_file.close()

file_seqsum = sys.argv[1]
file_bam = sys.argv[2]
file_out = sys.argv[3]
pore_time_aln = read_getter_bam_f5(
    file_seqsum,
    file_bam,
    file_out
)

pysam.sort("-o", file_out+'.sort.bam', file_out+'.bam')
pysam.index(file_out+'.sort.bam')
