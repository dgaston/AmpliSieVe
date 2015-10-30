__author__ = 'dgaston'

import sys
import csv
import argparse
import itertools
import timeit

import readfq

from multiprocessing import Pool
from fuzzywuzzy import fuzz
from collections import defaultdict
from nltk.metrics import edit_distance


def simple_reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    bases = reversed([complement.get(base, base) for base in bases])
    bases = ''.join(bases)

    return bases


def read_manifest(infile):
    primers = dict()
    dlso_dict = dict()
    ulso_dict = dict()

    with open(infile, 'rb') as infh:
        reader = csv.reader(infh, dialect='excel-tab')
        reader.next()
        for row in reader:
            primer_dict = dict()

            primer_dict['ulso'] = row[1]
            primer_dict['dlso'] = simple_reverse_complement(row[2])

            primer_dict['ulso_length'] = len(primer_dict['ulso'])
            primer_dict['dlso_length'] = len(primer_dict['dlso'])

            primers[row[0]] = primer_dict

            ulso_dict[primer_dict['ulso']] = row[0]
            dlso_dict[primer_dict['dlso']] = row[0]
    return primers, ulso_dict, dlso_dict


def match_read_primers(data):
    read1, read2, primers, ulso, dlso = data

    amplicon_assignment = defaultdict(int)

    for amplicon in primers:
        read1_primer_region = read1[1][0:primers[amplicon]['dlso_length']]
        read2_primer_region = read2[1][0:primers[amplicon]['ulso_length']]

        # If we get a perfect match to the same amplicon for both primer regions return immediately
        if read1_primer_region == primers[amplicon]['dlso'] and read2_primer_region == primers[amplicon]['ulso']:
            amplicon_assignment['ratio1_amplicon'] = amplicon
            amplicon_assignment['ratio1'] = 100
            amplicon_assignment['read1_primer_seq'] = read1_primer_region
            amplicon_assignment['read1'] = read1

            amplicon_assignment['ratio2_amplicon'] = amplicon
            amplicon_assignment['ratio2'] = 100
            amplicon_assignment['read2_primer_seq'] = read2_primer_region
            amplicon_assignment['read2'] = read2

            return amplicon_assignment

        ratio1 = fuzz.ratio(read1_primer_region, primers[amplicon]['dlso'])
        ratio2 = fuzz.ratio(read2_primer_region, primers[amplicon]['ulso'])

        if ratio1 > amplicon_assignment['ratio1']:
            amplicon_assignment['ratio1_amplicon'] = amplicon
            amplicon_assignment['ratio1'] = ratio1
            amplicon_assignment['read1_primer_seq'] = read1_primer_region
            amplicon_assignment['read1'] = read1
        if ratio2 > amplicon_assignment['ratio2']:
            amplicon_assignment['ratio2_amplicon'] = amplicon
            amplicon_assignment['ratio2'] = ratio2
            amplicon_assignment['read2_primer_seq'] = read2_primer_region
            amplicon_assignment['read2'] = read2

    return amplicon_assignment


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infiles', help="Input fastq files, comma-separated with read 1 first [Required]")
    parser.add_argument('-m', '--manifest', help="Manifest CSV file with ULSO and DLSO primer sequences per amplicon")
    parser.add_argument('-t', '--threads', help="Number of processor threads to use")
    parser.add_argument('-o', '--output', help='Output root for naming output fastq files per amplicon')

    args = parser.parse_args()
    read_files = args.infiles.split(',')

    start_time = timeit.default_timer()

    sys.stdout.write("Reading primer manifest file\n")
    primer_sets, ulso_dict, dlso_dict = read_manifest(args.manifest)

    fastq_file1 = open(read_files[0], 'r')
    fastq_file2 = open(read_files[1], 'r')

    pool = Pool(processes=int(args.threads))

    read_comparisons = list()
    amplicon_reads = dict()
    amplicon_fastq_file_pairs = dict()

    num_processed = 0
    num_perfect = 0
    num_mismatched = 0

    sys.stdout.write("Sorting reads into appropriate amplicons\n")
    for result in pool.imap(match_read_primers, itertools.izip(readfq.readfq(fastq_file1), readfq.readfq(fastq_file2),
                                                               itertools.repeat(primer_sets),
                                                               itertools.repeat(ulso_dict),
                                                               itertools.repeat(dlso_dict)), 100000):
        num_processed += 1
        if (num_processed % 10000) == 0:
            sys.stdout.write("Processed %s reads\n" % num_processed)
        if result['ratio1_amplicon'] != result['ratio2_amplicon']:
            num_mismatched += 1
            read1_fastq_file = open("%s-mismatched_R1.fastq" % args.output, "a")
            read2_fastq_file = open("%s-mismatched_R2.fastq" % args.output, "a")

            read1_fastq_file.write("%s\n%s\n+\n%s\n" % (result['read1'][0], result['read1'][1], result['read1'][2]))
            read2_fastq_file.write("%s\n%s\n+\n%s\n" % (result['read2'][0], result['read2'][1], result['read2'][2]))

            read1_fastq_file.close()
            read2_fastq_file.close()
            continue

        if result['ratio1'] == 100 and result['ratio2'] == 100:
            num_perfect += 1

        read1_fastq_name = "%s-%s_R1.fastq" % (args.output, result['ratio1_amplicon'])
        read2_fastq_name = "%s-%s_R2.fastq" % (args.output, result['ratio2_amplicon'])

        if result['ratio1_amplicon'] not in amplicon_fastq_file_pairs.keys():
            amplicon_files_dict = dict()
            amplicon_files_dict['read1'] = read1_fastq_name
            amplicon_files_dict['read2'] = read2_fastq_file
            amplicon_fastq_file_pairs[result['ratio1_amplicon']] = amplicon_files_dict

        read1_fastq_file = open(read1_fastq_name, "a")
        read2_fastq_file = open(read2_fastq_name, "a")

        read1_fastq_file.write("%s\n%s\n+\n%s\n" % (result['read1'][0], result['read1'][1], result['read1'][2]))
        read2_fastq_file.write("%s\n%s\n+\n%s\n" % (result['read2'][0], result['read2'][1], result['read2'][2]))

        read1_fastq_file.close()
        read2_fastq_file.close()

    fastq_file1.close()
    fastq_file2.close()

    elapsed = timeit.default_timer() - start_time
    sys.stdout.write("********************************************\n")
    sys.stdout.write("Total elapsed time for pipeline: %s sec\n" % elapsed)
    sys.stdout.write("Number of mismatched read pairs: %s\n" % num_mismatched)
    sys.stdout.write("Number of perfectly matching reads: %s" % num_perfect)
    sys.stdout.write("Completed pipeline\n")
