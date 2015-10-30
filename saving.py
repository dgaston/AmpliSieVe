__author__ = 'dgaston'

import csv
import sys


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)

    return bases


def read_manifest(infile):
    primers = dict()
    with open(infile, 'rb') as infh:
        reader = csv.reader(infh, dialect='excel-tab')
        reader.next()
        for row in reader:
            primer_dict = dict()

            primer_dict['ulso'] = row[1]
            primer_dict['dlso'] = reverse_complement(row[2])

            primer_dict['ulso_length'] = len(primer_dict['ulso'])
            primer_dict['dlso_length'] = len(primer_dict['dlso'])

            primers[row[0]] = primer_dict

            sys.stdout.write("%s - ULSO: %s, DLSO: %s, RevComp: %s\n" % (row[0], primer_dict['ulso'], row[2],
                                                                         primer_dict['dlso']))
    return primers


if __name__ == "__main__":
    primer_sets = read_manifest("/data/Resources/TruSight-Myeloid-Primers.txt")