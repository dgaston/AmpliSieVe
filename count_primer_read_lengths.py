__author__ = 'dgaston'

import csv
import sys
import argparse
from collections import defaultdict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', help="Input manifest file to parse [Required]")

    args = parser.parse_args()

    counts = defaultdict(int)
    with open(args.infile, 'rU') as infile:
        reader = csv.reader(infile, dialect='excel-tab')
        reader.next()
        for row in reader:
            counts[len(row[7])] += 1
            counts[len(row[8])] += 1

    sys.stdout.write("Length\tNumber\n")
    for num in counts:
        sys.stdout.write("%s\t %s\n" % (num, counts[num]))