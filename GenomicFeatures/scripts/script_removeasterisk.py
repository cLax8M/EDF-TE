#!/usr/bin/env python3
import sys

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]

with open(input_fasta) as inp, open(output_fasta, 'w') as out:
    for line in inp:
        if line.startswith('>'):
            # keep only first token (before space)
            header = line.strip().split()[0]
            out.write(header + '\n')
        else:
            # remove stop codon symbols
            out.write(line.replace('*', ''))