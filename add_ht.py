#!/usr/bin/env python3

import pysam
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-h", "--haplotype_bam", required=True,
            help='....')
    parser.add_argument("-m", "--methylation_tsv")
    parser.add_argument("-o", "--output_file_name")
    args = parser.parse_args()

    haplo_dic = {}
    samfile = pysam.AlignmentFile(args.haplotype_bam)

    for read in samfile:
        if read.flag & 3844 == 0:
            hpt = str(read.get_tag('HP')) if read.has_tag('HP') else None
            rname = str(read.qname)
            haplo_dic[rname] = hpt

    with open(args.methylation_tsv, 'r') as mc, open(args.output_file_name, 'w') as mh:
        mh.write(mc.readline().strip() + '\tHaplotype\n')
        for line in mc:
            read = line.split('\t')[4]
            mh.write(line.strip() + '\t' + str(haplo_dic[read] if read in haplo_dic else None) + '\n')


if __name__ == '__main__':
    main()

