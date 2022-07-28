#!/usr/bin/env python3

import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq

'''
getseqfromfasta_v1.py
UPDATED: v1 1/27/2021

INPUT: Just one genome sequence in a FASTA format.

sRNA list:
accession	sRNA	seq_from	seq_to	strand
NC_020211.1	6S	    4362482	    4362664	   +

Make sure that if it's on the - strand then the seq from is greater than the seq to (if you do seq_to - seq_from you get a negative number, where you'd get a postivie number if it was on the + strand)

OUTPUT: Whatever file name you want for your output fasta file.
'''


def getseqfromFASTA(inputsRNAlist, inputFASTA, outputFASTA):

    input_accession = []
    input_sRNA = []
    input_start = []
    input_end = []
    input_strand = []

    FASTAseq = []

    with open(inputsRNAlist) as fh:
        fhcsv = csv.reader(fh, delimiter=',')
        next(fhcsv)
        for entry in fhcsv:
            input_accession.append(entry[0])
            input_sRNA.append(entry[1])
            input_start.append(entry[2])
            input_end.append(entry[3])
            input_strand.append(entry[4])



    for value in range (0, len(input_sRNA)):
        recordtracker = 0

        for seq_record in SeqIO.parse(inputFASTA, "fasta"):

            recordtracker = recordtracker + 1

            start = int(input_start[value])
            end = int(input_end[value])

            myseq = seq_record.seq
            description = seq_record.description
            coordinates = str(start) + ".." + str(end)

            if start > end:
                # if strand == "minus":
                seqslice = myseq[end:start]
                RCseqslice = seqslice.reverse_complement()
                record = (">" + input_sRNA[value] + " " + description + ", c" + coordinates + "\n" + RCseqslice)
                FASTAseq.append(record)

            if end > start:
                # if strand == "plus":
                seqslice = myseq[start:end]
                record = (">" + input_sRNA[value] + " " + description + ", " + coordinates + "\n" + seqslice)
                FASTAseq.append(record)

        if recordtracker > 1:
            print("ERROR in code block 4: MORE THAN ONE GENOME WITH THAT NAME FOUND IN INPUT FASTA FILE.")
            exit()

        if recordtracker == 0:
            record = (">ERROR: NO SEQ FOUND")
            FASTAseq.append(record)

    if len(input_sRNA) != len(FASTAseq):
        print("ERROR: Different number of input and results.")


    with open(outputFASTA, 'a') as outputfh:
        for line in FASTAseq:
            outputfh.write(str(line) + "\n")


if __name__ == '__main__':
    if len(sys.argv) == 4:
        getseqfromFASTA(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print("")
        sys.exit(0)

