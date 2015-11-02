__author__ = 'dgaston'

import os
import sys
import csv
import argparse
import itertools
import timeit
import subprocess as sub

import readfq

from multiprocessing import Pool
from fuzzywuzzy import fuzz
from collections import defaultdict


def simple_reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    bases = reversed([complement.get(base, base) for base in bases])
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
            primer_dict['dlso'] = simple_reverse_complement(row[2])

            primer_dict['ulso_length'] = len(primer_dict['ulso'])
            primer_dict['dlso_length'] = len(primer_dict['dlso'])

            primers[row[0]] = primer_dict
    return primers


def match_read_primers(data):
    read1, read2, primers  = data

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


def write_fastq_pairs(file1, file2, read_result):
    read1_fastq_file = open(file1, "a")
    read2_fastq_file = open(file2, "a")

    read1_fastq_file.write("@%s\n%s\n+\n%s\n" % (read_result['read1'][0], read_result['read1'][1],
                                                 read_result['read1'][2]))
    read2_fastq_file.write("@%s\n%s\n+\n%s\n" % (read_result['read2'][0], read_result['read2'][1],
                                                 read_result['read2'][2]))

    read1_fastq_file.close()
    read2_fastq_file.close()


def read_bed(infile):
    amplicons = dict()
    with open(infile, 'rb') as infh:
        reader = csv.reader(infh, dialect='excel-tab')
        reader.next()
        for row in reader:
            amplicon_dict = dict()
            amplicon_dict['chr'] = row[0]
            amplicon_dict['start'] = int(row[1])
            amplicon_dict['end'] = int(row[2])
            amplicon_dict['size'] = int(row[2]) - int(row[1])
            amplicons[row[3]] = amplicon_dict

    return amplicons

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infiles', help="Input fastq files, comma-separated with read 1 first [Required]")
    parser.add_argument('-m', '--manifest', help="Manifest CSV file with ULSO and DLSO primer sequences per amplicon")
    parser.add_argument('-b', '--bed', help="BED file for target amplicons")
    parser.add_argument('-l', '--length', help="Read length")
    parser.add_argument('-t', '--threads', help="Number of processor threads to use")
    parser.add_argument('-o', '--output', help='Output root for naming output fastq files per amplicon')

    args = parser.parse_args()
    read_files = args.infiles.split(',')

    start_time = timeit.default_timer()

    sys.stdout.write("Reading primer manifest file\n")
    primer_sets = read_manifest(args.manifest)
    amplicon_bed = read_bed(args.bed)

    fastq_file1 = open(read_files[0], 'r')
    fastq_file2 = open(read_files[1], 'r')

    pool = Pool(processes=int(args.threads))

    read_comparisons = list()
    amplicon_reads = dict()
    amplicon_fastq_file_pairs = dict()
    no_read_amplicons = list()

    num_processed = 0
    num_perfect = 0
    num_mismatched = 0

    sys.stdout.write("Sorting reads into appropriate amplicons\n")
    for result in pool.imap(match_read_primers, itertools.izip(readfq.readfq(fastq_file1), readfq.readfq(fastq_file2),
                                                               itertools.repeat(primer_sets)), 100000):
        num_processed += 1
        if (num_processed % 10000) == 0:
            sys.stdout.write("Processed %s reads\n" % num_processed)

        if result['ratio1_amplicon'] != result['ratio2_amplicon']:
            num_mismatched += 1
            read1_filename = "%s-mismatched_R1.fastq" % args.output
            read2_filename = "%s-mismatched_R2.fastq" % args.output
            write_fastq_pairs(read1_filename, read2_filename, result)
            continue

        read1_fastq_name = "%s-%s_R1.fastq" % (args.output, result['ratio1_amplicon'])
        read2_fastq_name = "%s-%s_R2.fastq" % (args.output, result['ratio2_amplicon'])
        write_fastq_pairs(read1_fastq_name, read2_fastq_name, result)

        if result['ratio1'] == 100 and result['ratio2'] == 100:
            num_perfect += 1

        if result['ratio1_amplicon'] not in amplicon_fastq_file_pairs.keys():
            amplicon_files_dict = dict()
            amplicon_files_dict['read1'] = read1_fastq_name
            amplicon_files_dict['read2'] = read2_fastq_name
            amplicon_fastq_file_pairs[result['ratio1_amplicon']] = amplicon_files_dict

    fastq_file1.close()
    fastq_file2.close()

    elapsed = timeit.default_timer() - start_time
    sys.stdout.write("********************************************\n")
    sys.stdout.write("Elapsed time for binning amplicons: %s sec\n" % elapsed)
    sys.stdout.write("Number of mismatched read pairs: %s\n" % num_mismatched)
    sys.stdout.write("Number of perfectly matching reads: %s\n" % num_perfect)

    sys.stdout.write("Merging read pairs for each amplicon using Pear\n")
    for amplicon in primer_sets:
        output = "%s-%s" % (args.output, amplicon)
        logfile = "%s-%s-pear-merge.log" % (args.output, amplicon)
        sys.stdout.write("Merging read pairs for amplicon: %s\n" % amplicon)

        read1 = "%s-%s_R1.fastq" % (args.output, amplicon)
        read2 = "%s-%s_R2.fastq" % (args.output, amplicon)

        read1gz = "%s.gz" % read1
        read2gz = "%s.gz" % read2

        bgzip_command1 = ("bgzip %s" % read1)
        bgzip_command2 = ("bgzip %s" % read2)

        command = ("pear -f %s -r %s -o %s -j %s" % (read1gz, read2gz, output, args.threads))

        if os.path.isfile(read1):
            with open(logfile, "wb") as err:
                if not os.path.isfile(read1gz):
                    p = sub.Popen(bgzip_command1, stdout=sub.PIPE, stderr=err, shell=True)
                    output = p.communicate()
                    code = p.returncode
                    if code:
                        sys.stdout.write("An error occurred. Please check %s for details\n" % logfile)
                        sys.stdout.write("%s\n" % output)
                        sys.stderr.write("An error occurred. Please check %s for details\n" % logfile)

                if not os.path.isfile(read2gz):
                    p = sub.Popen(bgzip_command2, stdout=sub.PIPE, stderr=err, shell=True)
                    output = p.communicate()
                    code = p.returncode
                    if code:
                        sys.stdout.write("An error occurred. Please check %s for details\n" % logfile)
                        sys.stdout.write("%s\n" % output)
                        sys.stderr.write("An error occurred. Please check %s for details\n" % logfile)

                if not os.path.isfile(output):
                    p = sub.Popen(command, stdout=sub.PIPE, stderr=err, shell=True)
                    output = p.communicate()
                    code = p.returncode
                    if code:
                        sys.stdout.write("An error occurred. Please check %s for details\n" % logfile)
                        sys.stdout.write("%s\n" % output)
                        sys.stderr.write("An error occurred. Please check %s for details\n" % logfile)
        else:
            no_read_amplicons.append(amplicon)

    sys.stdout.write("Calculating expected and observed overlaps for amplicons\n")
    for amplicon in primer_sets:
        merged_file = "%s-%s.assembled.fastq" % (args.output, amplicon)
        overlap_dist = list()
        overlap_diff_dist = list()
        expected_overlap = (amplicon_bed[amplicon]['size'] - (2 * int(args.length)) * -1)

        if os.path.isfile(merged_file):
            sys.stdout.write("Reading file: %s\n" % merged_file)
            for name, seq, qual in readfq.readfq(merged_file):
                obs_overlap = (len(seq) - (2 * int(args.length)) * -1)
                overlap_dist.append(obs_overlap)
                overlap_diff_dist.append(obs_overlap - expected_overlap)

            sys.stdout.write("Amplicon %s\tExpected Overlap: %s\tMin Obs Overlap: %s\tMax ObsOverlap:%s\n " %
                             (amplicon, expected_overlap, min(overlap_dist), max(overlap_dist)))
        elif amplicon in no_read_amplicons:
            sys.stdout.write("Amplicon %s\tNo Data\n")
        else:
            sys.stderr.write("WARNING: Merged file %s for Amplicon %s does not exist but was not logged in no coverage"
                             "amplicons\n" % (merged_file, amplicon))

    sys.stdout.write("Completed pipeline\n")
