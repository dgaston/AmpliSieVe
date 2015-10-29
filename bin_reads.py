__author__ = 'dgaston'

import sys
import csv
import argparse
import itertools
import HTSeq

from multiprocessing import Pool
from fuzzywuzzy import fuzz
from collections import defaultdict
from nltk.metrics import edit_distance


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


def match_read_primers(read1, read2, primer_sets):
    amplicon_assignment = defaultdict(lambda: 0.0)

    for amplicon in primer_sets:
        read1_primer_region = read1.seq[0:len(primer_sets[amplicon]['dlso'])]
        read2_primer_region = read2.seq[0:len(primer_sets[amplicon]['ulso'])]

        dlso = HTSeq.Sequence(primer_sets[amplicon]['dlso'], "DLSO")

        ratio1 = fuzz.ratio(read1_primer_region, dlso.get_reverse_complement())
        ratio2 = fuzz.ratio(read2_primer_region, primer_sets[amplicon]['ulso'])

        if ratio1 > amplicon_assignment['ratio1']:
            amplicon_assignment['ratio1_amplicon'] = amplicon
            amplicon_assignment['ratio1'] = ratio1
            amplicon_assignment['read1_primer_seq'] = read1_primer_region
        if ratio2 > amplicon_assignment['ratio2']:
            amplicon_assignment['ratio2_amplicon'] = amplicon
            amplicon_assignment['ratio2'] = ratio2
            amplicon_assignment['read2_primer_seq'] = read2_primer_region

    return amplicon_assignment


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infiles', help="Input fastq files, comma-separated with read 1 first [Required]")
    parser.add_argument('-m', '--manifest', help="Manifest CSV file with ULSO and DLSO primer sequences per amplicon")
    parser.add_argument('-t', '--threads', help="Number of processor threads to use")
    parser.add_argument('-o', '--output', help='Output root for naming output fastq files per amplicon')

    args = parser.parse_args()
    read_files = args.infiles.split(',')

    sys.stdout.write("Reading primer manifest file\n")
    primer_sets = read_manifest(args.manifest)

    fastq_file1 = HTSeq.FastqReader(read_files[0])
    fastq_file2 = HTSeq.FastqReader(read_files[1])

    read_comparisons = list()
    amplicon_reads = dict()
    sys.stdout.write("Sorting reads into appropriate amplicons\n")
    num_processed = 0
    for read1, read2 in itertools.izip(fastq_file1, fastq_file2):
        num_processed += 1
        if (num_processed % 10000) == 0:
            sys.stdout.write("Processed %s reads\n" % num_processed)
        result = match_read_primers(read1, read2, primer_sets)
        if result['ratio1_amplicon'] != result['ratio2_amplicon']:
            read1_fastq_file = open("%s-mismatched_R1.fastq" % args.output, "a")
            read2_fastq_file = open("%s-mismatched_R2.fastq" % args.output, "a")

            read1.write_to_fastq_file(read1_fastq_file)
            read2.write_to_fastq_file(read1_fastq_file)

            read1_fastq_file.close()
            read2_fastq_file.close()
            continue

        read1_fastq_name = "%s-%s_R1.fastq" % (args.output, result['ratio1_amplicon'])
        read2_fastq_name = "%s-%s_R2.fastq" % (args.output, result['ratio2_amplicon'])

        read1_fastq_file = open(read1_fastq_name, "w")
        read2_fastq_file = open(read2_fastq_name, "w")

        read1.write_to_fastq_file(read1_fastq_file)
        read2.write_to_fastq_file(read2_fastq_file)

        read1_fastq_file.close()
        read2_fastq_file.close()

    sys.stdout.write("Writing results to file\n")
