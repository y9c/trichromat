#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2025 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2025-05-14 15:00

import sys
import os

def revcomp(seq):
    # str.translate is much faster for complementation
    trans = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(trans)[::-1]

class FastaIndex:
    def __init__(self, fasta_file):
        self.fasta_file = open(fasta_file, 'r')
        self.index = {}
        fai_file = fasta_file + ".fai"
        if not os.path.exists(fai_file):
            raise FileNotFoundError(f"FASTA index {fai_file} not found.")
        
        with open(fai_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                name = fields[0]
                length = int(fields[1])
                offset = int(fields[2])
                linebases = int(fields[3])
                linewidth = int(fields[4])
                self.index[name] = (length, offset, linebases, linewidth)

    def fetch(self, chrom, start, end):
        if chrom not in self.index:
            return "N" * (end - start + 1)
        
        length, offset, linebases, linewidth = self.index[chrom]
        
        # Clip to chromosome boundaries
        fetch_start = max(1, start)
        fetch_end = min(length, end)
        
        if fetch_start > fetch_end:
            return "N" * (end - start + 1)
            
        # Calculate start and end offsets in the file
        start_line = (fetch_start - 1) // linebases
        start_col = (fetch_start - 1) % linebases
        start_offset = offset + start_line * linewidth + start_col
        
        end_line = (fetch_end - 1) // linebases
        end_col = (fetch_end - 1) % linebases
        end_offset = offset + end_line * linewidth + end_col
        
        # Read the range including newlines
        self.fasta_file.seek(start_offset)
        raw_seq = self.fasta_file.read(end_offset - start_offset + 1)
        
        # Remove newlines and join
        seq = raw_seq.replace("\n", "").replace("\r", "").upper()
        
        # Pad with Ns if we clipped
        if start < 1:
            seq = "N" * (1 - start) + seq
        if end > length:
            seq = seq + "N" * (end - length)
            
        return seq

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Annotate motifs for hisat3n-table output.")
    parser.add_argument("-r", "--ref", required=True, help="Reference FASTA file.")
    parser.add_argument("-k", "--motif-flanking", type=int, default=2, help="Number of flanking bases.")
    args = parser.parse_args()

    fa = FastaIndex(args.ref)
    k = args.motif_flanking

    for line in sys.stdin:
        if line.startswith("#") or not line.strip():
            continue
        
        fields = line.strip().split('\t')
        if len(fields) < 5:
            continue
        
        chrom = fields[0]
        pos = int(fields[1])
        strand = fields[2]
        converted = fields[3]
        unconverted = fields[4]
        
        # hisat3n-table is 1-indexed
        motif = fa.fetch(chrom, pos - k, pos + k)
        if strand == "-":
            motif = revcomp(motif)
        
        print(f"{chrom}\t{pos}\t{motif}\t{converted}\t{unconverted}")

if __name__ == "__main__":
    main()
