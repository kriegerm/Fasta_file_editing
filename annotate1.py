#!/usr/bin/env python3

'''
annotate1.py - USE WITH ONE KNOWN GENOME LOCATION
UPDATED 8/26/2020

This code will take in a csv gb file and create a dictionary with location (start, end) as the key and the entire description as the value.

It will then take an input of positions within the genome and, using the gb dictionary, it will annotate the gene that the
position falls in.

'''

import csv
import sys
def annotate (inputcsv, locationscsv, outputcsv):

    #input csv file structure by column (get this by exporting a .gff file, converting to csv, and then arranging the columns by hand)
    #1 - feature type
    #2 - start
    #3 - end
    #4 - strand
    #5 - gene
    #6 - locus tag
    #7 - product
    #8 - protein ID

##Creates a dictionary with location (start, end) as the key and the entire description as the value
    gbdict = {}
    with open(locationscsv, 'r') as fh:
        fhcsv = csv.reader(fh, delimiter=',')
        field_names_list = next(fhcsv)

        for entry in fhcsv:
            locus_tag = entry[5]
            locus_tag = locus_tag.rstrip()
            locus_tag = locus_tag.lstrip()

            if entry[0] == "CDS": #if I did include gene then all the CDS would be overwritten
                locustag = []
                locustag.append(locus_tag)
                gbdict[str(locustag)] = entry

            if entry[0] != "CDS":  # Adding back in all the entires that did not have CDS records
                locustag = []
                locustag.append(locus_tag)
                if str(locustag) not in gbdict.keys():
                    gbdict[str(locustag)] = entry


    ##Takes an input file of a list of locations
    with open(inputcsv, 'r', encoding='utf-8-sig') as fh:
        fhcsv = csv.reader(fh, delimiter=',')
        for entry in fhcsv:
            coordinate = int(entry[0])
            coordinatecheck = 0

            for location, description in gbdict.items():
                location = list(location.split(","))
                start = int(location[0].lstrip("[").rstrip())
                end = int(location[1].lstrip().rstrip("]"))
                description = str(description)
                description = description.rstrip("]")
                description = description.lstrip("[")

                results = []
                writeline = []

                if coordinate in range(start, end):
                    results.append(description)
                    writeline.append(coordinate)
                    writeline.append(results)
                    coordinatecheck = coordinatecheck + 1

                    writeline = str(writeline)
                    writeline = writeline.replace("'", "")
                    writeline = writeline.replace('"', "")
                    writeline = writeline.replace("]", "")
                    writeline = writeline.replace("[", "")
                    writeline = writeline.replace("/", "")
                    writeline = writeline.replace("\\", "")

                    with open(outputcsv, 'a') as output:
                        writer = csv.writer(output)
                        writer.writerow(writeline)

            if coordinatecheck == 0:
                with open(outputcsv, 'a') as output:
                    writer = csv.writer(output)
                    writer.writerow(writeline)



if __name__ == '__main__':
    if len(sys.argv) == 4:
         annotate(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
         print("Usage: input csv with gene locations, input csv of genome (see documentation for arragement), output file name .csv ")
         sys.exit(0)
