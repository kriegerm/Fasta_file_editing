#!/usr/bin/env python3

'''
findseqinfastafile_v1.py
UPDATED: 1/5/2021

This code takes in 1. input file of genome(s), fasta format; 2. Genome name (i.e. NC_); 3. Starting coordinate of sequence; 4. Ending coordinate of sequence; 5. Strand, either plus or minus
and outputs the descriptor and sequence.

'''

import sys
from Bio import SeqIO
from Bio.Seq import Seq

def return_fasta_seq (inputgenome, genome, start, end, strand):
    recordtracker = 0

    start = int(start)
    end = int(end)
    strands = ["plus", "minus"]

    if strand not in strands:
        print("ENTER CORRECT VALUE FOR STRAND: EITHER PLUS OR MINUS")

    for seq_record in SeqIO.parse(inputgenome, "fasta"):

        if genome == seq_record.id:
            recordtracker = recordtracker + 1
            myseq = seq_record.seq
            seqslice = myseq[start:end]

            description = seq_record.description
            coordinates = str(start) + ".." + str(end)


            if strand == "minus":
                RCseqslice = seqslice.reverse_complement()
                print(">" + description + ", c" + coordinates + "\n" + RCseqslice)


            if strand == "plus":
                print(">" + description + ", " + coordinates + "\n" + seqslice)


    if recordtracker == 0:
        print("ERROR: GENOME NAME NOT FOUND IN INPUT FASTA FILE.")

    if recordtracker > 1:
        print("ERROR: MORE THAN ONE GENOME WITH THAT NAME FOUND IN INPUT FASTA FILE.")
        exit()


if __name__ == '__main__':
    if len(sys.argv) == 6:
         return_fasta_seq(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    else:
         print("Usage: 1. input file of genome(s), fasta format; 2. Genome name (i.e. NC_); 3. Starting coordinate of sequence; 4. Ending coordinate of sequence; 5. Strand, either plus or minus ")
         sys.exit(0)
