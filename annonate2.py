#!/usr/bin/env python3

'''
annonate2.py - USE WITH KNOWN LOCUS TAG INPUT
UPDATED: 8/27/20

This code will take in a csv gb file and create a dictionary with location (start, end) as the key and the entire description as the value.
It will then take an input gene names, using the gb dictionary, it will annotate the gene that the
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

            if entry[0] == "CDS": #if I did include gene then all the CDS would be overwritten (which has the gene description)
                locustag = []
                locustag.append(locus_tag)
                gbdict[str(locustag)] = entry

            if entry[0] != "CDS":  # Adding back in all the entires that did not have "gene" records
                locustag = []
                locustag.append(locus_tag)
                if str(locustag) not in gbdict.keys():
                    gbdict[str(locustag)] = entry


    ##Takes an input file of a list of locations
    with open(inputcsv, 'r', encoding='utf-8-sig') as fh:

        inputcheck = 0
        outputcheck = 0

        fhcsv = csv.reader(fh, delimiter=',')
        for entry in fhcsv:
            inputcheck += 1

            new_locus = entry[0]
            new_locus = new_locus.rstrip()
            new_locus = new_locus.lstrip()

            coordinatecheck = 0

            for locus, description in gbdict.items():
                locus = str(locus)
                locus = locus.rstrip(" ']")
                locus = locus.lstrip("['")
                locus = locus.replace(" ", "")

                description = str(description)
                description = description.rstrip("]")
                description = description.lstrip("[")

                results = []
                writeline = []

                if new_locus == locus:
                    outputcheck += 1

                    results.append(description)
                    writeline.append(locus)
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
                        writer.writerow([writeline])

            if coordinatecheck == 0:
                outputcheck += 1

                with open(outputcsv, 'a') as output:
                    writer = csv.writer(output)
                    writer.writerow([writeline])

    if inputcheck != outputcheck:
        print("DIFFERENT NUMBER OF INPUT AND OUTPUT RESULTS!!! CHECK RESULTS CAREFULLY.")


if __name__ == '__main__':
    if len(sys.argv) == 4:
         annotate(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
         print("Usage: csv of locus tags for input, csv with genome information (see code documentation for arrangement), output csv file name ")
