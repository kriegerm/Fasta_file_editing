#!/usr/bin/env python3

'''
remove_repeptitive_nts.py

This script scans a .fastq file and removes strings of 10 or more redudant nt's.
Example:
The sequence CCCCGGTTCGAAAAAAAAAAAAAAAAAAAAAAGACATGATAATGTGTCGAAACACACTGGGTTTCCCCATTCGGGTATCGCCGGTTATAACGGTTCATATC
would be removed.

It returns an edited fastq file with only the non-redudant sequences remainin.
'''

import csv
import sys

def removerepetitive(inputfastq, outputfastq):
    with open(inputfastq, 'r') as fh:

        totalcounter = 0
        Tcounter = 0
        Acounter = 0
        Gcounter = 0
        Ccounter = 0

        for line in fh:

            line = str(line)
            if line.startswith("@"):

                totalcounter += 1
                seq = next(fh)

                if "TTTTTTTTTT" in seq:
                    Tcounter += 1
                    with open(outputfastq, 'a') as output:
                        output.write(line)
                        output.write(seq)
                        output.write(next(fh))
                        output.write(next(fh))

                if "AAAAAAAAAA" in seq:
                    Acounter += 1
                    with open(outputfastq, 'a') as output:
                        output.write(line)
                        output.write(seq)
                        output.write(next(fh))
                        output.write(next(fh))

                if "GGGGGGGGGG" in seq:
                    Gcounter += 1
                    with open(outputfastq, 'a') as output:
                        output.write(line)
                        output.write(seq)
                        output.write(next(fh))
                        output.write(next(fh))

                if "CCCCCCCCCC" in seq:
                    Ccounter += 1
                    with open(outputfastq, 'a') as output:
                        output.write(line)
                        output.write(seq)
                        output.write(next(fh))
                        output.write(next(fh))

    print("%d total sequences in the file" %(totalcounter))
    print("%d strings of 10 or more T's found" %(Tcounter))
    print("%d strings of 10 or more A's found" % (Acounter))
    print("%d strings of 10 or more G's found" % (Gcounter))
    print("%d strings of 10 or more C's found" % (Ccounter))

    return

if __name__ == '__main__':
    if len(sys.argv) == 3:
        removerepetitive(sys.argv[1], sys.argv[2])
    else:
        print("Usage: input.fasq output.fastq ")
        sys.exit(0)
