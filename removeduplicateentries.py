#!/usr/bin/env python3

'''
removeduplicateentries.py
UPDATES: v1 1/10/21
This code removes duplicate entires from a CSV genome file.

How to prepare for running this code:
Compile a large .txt file of all the Feature Tables from the NCBI FTP website. This can be done by running a shell script with the FTP addresses to these feature tables, as below:
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/267/845/GCF_001267845.1_ASM126784v1/GCF_001267845.1_ASM126784v1_feature_table.txt.gz
One you download and unzip all these .txt feature tables, concatenate them into one big file.

This is your input file, as a .txt. If it is provided as a .csv you will need to change the delimitor to "," in the first open statement.
'''

import csv
import sys
import pandas as pd

def removeduplicates(input, outputCSV):
    with open(input, 'r') as fh:
        fhcsv = csv.reader(fh, delimiter='\t')

        db_genomename = []
        db_featuretype = []
        db_genestrand = []
        db_genename = []
        db_genestart =[]
        db_geneend = []
        db_geneshortname = []

        for entry in fhcsv:
            ##There are some random entries at the bottom of the .tab file but they start with #, so we are going to ignore those
            if entry[0] is not "#":
                allinfo = []
                featuretype = entry[0]  # gene, CDS, tRNA, etc
                genomeassembly = entry[6]  # name of organism, need to match this to the genomename from the .tab output
                genestart = entry[7]
                geneend = entry[8]
                genestrand = entry[9]
                genename = str(entry[13])
                geneshortname = str(entry[14])
                db_genename.append(genename)
                db_genomename.append(genomeassembly)
                db_featuretype.append(featuretype)
                db_genestrand.append(genestrand)
                db_genestart.append(genestart)
                db_geneend.append(geneend)
                db_geneshortname.append(geneshortname)

    data = [db_genomename, db_genename, db_geneshortname, db_featuretype, db_genestrand, db_genestart, db_geneend]
    df = pd.DataFrame(data)
    df = df.transpose()
    df.columns = (['Genome', "Gene", "Symbol", "FeatureType", "Strand", "Start", "End"])

    sorted = df.drop_duplicates(['Genome', 'Start', 'End'], keep='last')
    sorted.to_csv(outputCSV, index=False)

if __name__ == '__main__':
    if len(sys.argv) == 3:
         removeduplicates(sys.argv[1], sys.argv[2])
    else:
         print("Usage: 1) Input .txt file of Feature Tables; 2) Output file name CSV ")
         sys.exit(0)
