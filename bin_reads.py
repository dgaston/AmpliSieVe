__author__ = 'dgaston'

import sys
import csv
import argparse
import itertools
import HTSeq

from multiprocessing import Pool
from fuzzywuzzy import fuzz
from collections import defaultdict


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
    # (read1, read2, primer_sets) = comparison
    amplicon_assignment = defaultdict(int)

    for amplicon in primer_sets:
        read1_primer_region = read1.get_reverse_complement()[0:len(primer_sets[amplicon]['dlso'])]
        read2_primer_region = read2.seq[0:len(primer_sets[amplicon]['ulso'])]

        ratio1 = fuzz.ratio(read1_primer_region, primer_sets[amplicon]['dlso'])
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
    results = list()
    sys.stdout.write("Constructing comparisons for multiprocessing\n")
    for read1, read2 in itertools.izip(fastq_file1, fastq_file2):
        result = match_read_primers(read1, read2, primer_sets)
        sys.stdout.write("Read1 Amplicon: %s (%s), Read2 Amplicon: %s (%s)\n" % (result['ratio1_amplicon'],
                                                                                 result['ratio1'],
                                                                                 result['ratio2_amplicon'],
                                                                                 result['ratio2']))
        # read_comparisons.append((read1, read2, primer_sets))

    # sys.stdout.write("Performing comparisons\n")
    # pool = Pool(processes=int(args.threads))
    # result = pool.map_async(match_read_primers, itertools.izip(fastq_file1, fastq_file2, itertools.repeat(primer_sets)))
    # results = result.get()
    # pool.close()
    # pool.join()

    # for result in results:
    #     sys.stdout.write("Read1 Amplicon: %s, Read2 Amplicon: %s\n" % (result['ratio1_amplicon'],
    #                                                                    result['ratio2_amplicon']))

    sys.stdout.write("Writing results to file\n")
