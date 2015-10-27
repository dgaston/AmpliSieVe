__author__ = 'dgaston'

import argparse


def read_manifest(infile):
    primers = dict()
    with open(infile, 'rb') as infh:
        reader = csv.reader(infh, dialect='excel-tab')
        reader.next()
        for row in reader:
            primer_dict = dict()
            primer_dict['ulso'] = row[1]
            primer_dict['dlso'] = row[2]
            primers[row[0]] = primer_dict
    return primers

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_root', help="Input root from fastq files per amplicon [Required]")
    parser.add_argument('-m', '--manifest', help="Manifest CSV file with ULSO and DLSO primer sequences per amplicon")
    parser.add_argument('-t', '--threads', help="Number of processor threads to use")

    args = parser.parse_args()

    amplicons = read_manifest(args.manifest)

    for amplicon in amplicons:
        fastq1 = "%s-%s_R1.fastq" % (args.input_root, amplicon)
        fastq2 = "%s-%s_R2.fastq" % (args.input_root, amplicon)
        output = "%s-%s-merged.fastq" % (args.input_root, amplicon)

        command = ("pear -f %s -r %s -o %s -j %s" % (fastq1, fastq2, output, args.threads))
