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
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return "".join(complement.get(base, base) for base in reversed(seq))

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
        
        seq = []
        for p in range(start, end + 1):
            if p < 1 or p > length:
                seq.append("N")
                continue
            
            byte_offset = offset + ((p-1) // linebases) * linewidth + ((p-1) % linebases)
            self.fasta_file.seek(byte_offset)
            seq.append(self.fasta_file.read(1))
        
        return "".join(seq).upper()

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
