__author__ = 'dgaston'

import csv
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', help="Input manifest file to parse [Required]")
    parser.add_argument('-o', '--outfile', help="Output FASTA file name [Required]")

    args = parser.parse_args()

    with open(args.outfile, 'w') as outfile:
        with open(args.infile, 'rU') as infile:
            reader = csv.reader(infile, dialect='excel-tab')
            reader.next()
            for row in reader:
                outfile.write(">%s-%s\n" % (row[0], "ULSO"))
                outfile.write("%s\n" % row[7])

                outfile.write(">%s-%s\n" % (row[0], "DLSO"))
                outfile.write("%s\n" % row[8])
