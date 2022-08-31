#!/usr/bin/env python

import os
import pysam
from collections import namedtuple
from typing import List, Tuple
from array import array
import pandas as pd
import re
import numpy as np

def split_cigar(cigar) -> List[Tuple[int, str]]:
    """
    Converts a cigarstring into a workable format:
    '721H149M686H' -> [(721, 'H'), (149, 'M'), (686, 'H')]

    Adapted from Cycas' alignment.py
    """
    # Find all alternating letters and numbers
    cigars = list(
        zip(
            re.findall("\d+", cigar), re.findall("\D", cigar)
        )
    )
    # convert numerical part to an int
    cigars = [(int(x[0]), x[1]) for x in cigars]

    return cigars

def correct_cigar(pysam_obj) -> List[Tuple[int, str]]:
    """
    Corrects reverse CIGAR strings to be readable from a fowrard orientation:
    If '721H149M686H' is reverse -> '686H149M721H'
    """
    cigar = pysam_obj.cigarstring
    if not pysam_obj.is_reverse:
        return split_cigar(cigar)
    else:
        return split_cigar(cigar)[::-1]

def extract_section(full_seq, full_seq_quality, start, end, min_length=30) -> Tuple[int, int, str, str]:
    """
    Extracts the sequence and per-base quality of a segment.
    Only considers outputs segments that are larger than the set minimum length 'min_length'.
    """
    if end - start > min_length:
        return tuple([start, end, full_seq[start:end], full_seq_quality[start:end]])
    else:
        return

def get_realignment_df(realignment_index, reads, df) -> List[Tuple[str, int, float]]:
    """
    From the alignment index of each realignment BAM file,
    create a dataframe with recalculated mappability.
    """

    realignment_statistics = []
    for read in reads:
        try:
            realignment_index.find(read)
            
        except:
            continue

        segments = [x for x in realignment_index.find(read)]
        alignment_lengths = [x.query_alignment_length for x in segments]
        raw_length = df.loc[df['read'] == read]['raw_length']
        new_aligned_bases = df.loc[df['read'] == read]['aligned_bases'] + sum(alignment_lengths)
        new_mappability = new_aligned_bases/raw_length
        realignment_statistics.append(tuple([read, int(new_aligned_bases), float(new_mappability)]))
        
    return realignment_statistics

def get_unaligned_segments(bam_filename, unaligned_fastq, reads_df):
    """ Get unaligned read segments from initial read mapping.
    """
    # Load BAM file into pysam object
    bam_fh = pysam.AlignmentFile(bam_filename, 'rb')
    index = pysam.IndexedReads(bam_fh)
    index.build()

    # Initialize a Pandas DataFrame as a list. Will be turned to df later
    df = []

    # Fetch list of read names
    readnames = [read.query_name for read in bam_fh.fetch()]
    readnames = np.unique(readnames).tolist()

    Segment = namedtuple("Segment", ['rotated_cigar', 'pysam_obj'])
    fastq_reads = []
    rejected_reads = []
    for read in readnames:
        # Get list of aligned segments for each read
        segments = [x for x in index.find(read)]
        # Get primary alignment
        try:
            primary_aln = [x for x in segments if not x.is_secondary and not x.is_supplementary][0]
        except:
            # No primary alignment found for this read.
            # This could be because the read is SingleBackbone. Log and skip.
            rejected_reads.append(read)
            continue
        
        # We need read sequence and quality to write out new reads from unaligned segments
        full_seq = primary_aln.query_sequence
        full_seq_quality = ''.join(map(lambda x: chr( x+33 ), primary_aln.query_qualities))
        # Get length of primary alignment
        raw_length = primary_aln.query_length

        # Initiate sorted list of segments for each read
        segments = [Segment(correct_cigar(x),x) for x in segments]
        segments.sort(key = lambda x: x[0][0])

        aligned_bases = 0
        start_previous = 0
        unaligned_segments = []
        for segment in segments:
            # Calculate nr of aligned bases from alignment length of each segment
            aln_length = segment.pysam_obj.query_alignment_length
            aligned_bases += aln_length

            # Iterate over CIGAR of this segment and append unaligned subsegments
            start_segment = segment.rotated_cigar[0][0]
            if start_segment > start_previous:
                section = extract_section(full_seq, full_seq_quality, start_previous, start_segment, min_length=30)
                unaligned_segments.append(section) if section else None
            
            start_previous = aln_length + start_segment
        
        # Calculate original mappability and append to DataFrame
        mappability = aligned_bases/raw_length
        df.append({'read': read, 'raw_length': raw_length, 'aligned_bases': aligned_bases, 'mappability': mappability})

        # Get unaligned sequence from the ends of CIGAR strings of each segment
        last_unaligned_segment = segments[-1].rotated_cigar
        last_unaligned_segment = start_previous + last_unaligned_segment[-1][0] if last_unaligned_segment[-1][1] in ["H", "S"] else None

        if last_unaligned_segment:
            section = extract_section(full_seq, full_seq_quality, start_previous, start_segment, min_length=30)
            unaligned_segments.append(section) if section else None

        fastq_reads.append(tuple([read, unaligned_segments]))

    df = pd.DataFrame(df)
    df.to_csv(reads_df, index=False)

    # write out all reads   
    with open(unaligned_fastq, 'w+') as fastq:
        for read, sections in fastq_reads:
            for i, section in enumerate(sections):
                start, end, seq, qual = section[0], section[1], section[2], section[3]
                fastq.write(f"@{read}_{i+1}_{start}-{end}\n{seq}\n+\n{qual}\n")

    return unaligned_fastq, reads_df

def merge_segment_alignments(original_bam, realigned_bam, df_path, merged_bamfile):
    """ Retrieve segment realignments, recalculate their CIGAR strings and merge with original alignments.
    """

    text="""@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:248387328
@SQ	SN:chr2	LN:242696752
@SQ	SN:chr3	LN:201105948
@SQ	SN:chr4	LN:193574945
@SQ	SN:chr5	LN:182045439
@SQ	SN:chr6	LN:172126628
@SQ	SN:chr7	LN:160567428
@SQ	SN:chr8	LN:146259331
@SQ	SN:chr9	LN:150617247
@SQ	SN:chr10	LN:134758134
@SQ	SN:chr11	LN:135127769
@SQ	SN:chr12	LN:133324548
@SQ	SN:chr13	LN:113566686
@SQ	SN:chr14	LN:101161492
@SQ	SN:chr15	LN:99753195
@SQ	SN:chr16	LN:96330374
@SQ	SN:chr17	LN:84276897
@SQ	SN:chr18	LN:80542538
@SQ	SN:chr19	LN:61707364
@SQ	SN:chr20	LN:66210255
@SQ	SN:chr21	LN:45090682
@SQ	SN:chr22	LN:51324926
@SQ	SN:chrX	LN:154259566
@SQ	SN:chrY	LN:62460029
@SQ	SN:chrM	LN:16569
@SQ	SN:BB41G	LN:220
@SQ	SN:BB22	LN:247
@SQ	SN:BB24	LN:247
@SQ	SN:BB25	LN:247
@SQ	SN:BBCR	LN:247
@SQ	SN:PJET	LN:2974
"""
    # open output file as pysam object
    merged_bowtie_bam = pysam.AlignmentFile(merged_bamfile, 'wb', text=text)

    # load in the original bam file
    bam_fh = pysam.AlignmentFile(original_bam, 'rb')
    index = pysam.IndexedReads(bam_fh)
    index.build()

    # load in the realignment file
    bowtie_fh = pysam.AlignmentFile(realigned_bam, 'rb')
    bowtie_index = pysam.IndexedReads(bowtie_fh)
    bowtie_index.build()

    # get list of readnames from previously created dataframe
    #reads = list(df['read'])
    # fetch all realigned reads
    bowtie_reads = [x for x in bowtie_fh.fetch()]

    # recalculate mappability statistics and turn it into a pandas dataframe
    #bowtie_statistics = get_realignment_df(bowtie_index, reads, df)
    #bowtie_statistics_df = pd.DataFrame(bowtie_statistics, columns=['read', 'bowtie_aligned_bases', 'bowtie_mappability'])

    # Load the reads dataframe into a Pandas DataFrame
    df = pd.read_csv(df_path)

    for aln in bowtie_reads:
        # parse info from alignment query names
        read_name, segment_nr, segment_pos = aln.query_name.split('_')
        segment_start, segment_end = segment_pos.split('-')
        aln.query_name = read_name
        
        # get total read length from dataframe
        length = int(df.loc[df['read']==read_name]['raw_length'])

        # get corrected cigars in list format,
        # so that we can update the clipped segments
        # to correspond with the full sequence length
        cigar = correct_cigar(aln)
        cigar = [list(ash) for ash in cigar]

        # first clipped edge
        if cigar[0][1] in ['H', 'S']:
            cigar[0][1] = 'H' if cigar[0][1] == 'S' else 'H'
            cigar[0][0] += int(segment_start)
        else:
            cigar.insert(0, [segment_start, 'H'])

        # second clipped edge
        if cigar[-1][1] in ['H', 'S']:
            cigar[-1][1] = 'H' if cigar[-1][1] == 'S' else 'H'
            cigar[-1][0] += length - int(segment_end)
        else:
            cigar.append([length - int(segment_end), 'H'])

        ######
        # From original BAM
        segments = [x for x in index.find(read_name)]
        # Get primary alignment
        try:
            primary_aln = [x for x in segments if not x.is_secondary and not x.is_supplementary][0]
        except:
            continue
        
        full_seq = primary_aln.query_sequence
        full_seq_quality = primary_aln.query_qualities
        ######

        if aln.is_reverse:
            # 2048 + 16
            aln.flag = 2064
            cigar = cigar[::-1]
        else:
            aln.flag = 2048

        # get nucleotide sequence and quality for this realigned region
        aln.query_sequence = full_seq[int(cigar[0][0]):length-int(cigar[-1][0])]
        aln.query_qualities = full_seq_quality[int(cigar[0][0]):length-int(cigar[-1][0])]

        # update cigarstring
        cigarstring = [[str(x) for x in ash] for ash in cigar]
        cigarstring = ''.join(sum(cigarstring, []))
        aln.cigarstring = cigarstring

        # write to opened output file
        merged_bowtie_bam.write(aln)

    # write all original mappings to the opened output file
    for aln in bam_fh.fetch():
        merged_bowtie_bam.write(aln)
    
    return merged_bamfile

if __name__ == "__main__":
    from sys import argv

    if argv[1] == "get_unaligned_segments":
        original_bam, unaliged_fastq, reads_df = argv[2], argv[3], argv[3]
        get_unaligned_segments(bam=original_bam)

    elif argv[1] == "merge_segment_alignments":
        original_bam, realigned_bam, df_path, merged_bamfile = argv[2], argv[3], argv[4], argv[5]
        merge_segment_alignments(original_bam=original_bam, realigned_bam=realigned_bam, df_path=df_path, merged_bamfile=merged_bamfile)

    else:
        raise ValueError("Invalid function name. Please run one of two functions: \
        'recover_segmens.py get_unaligned segments...' or 'recover_segments.py merge_segment_alignments ...'")
