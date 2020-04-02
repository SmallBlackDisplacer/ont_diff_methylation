import pysam
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-hb", "--haplotype_bam", type=str)
parser.add_argument("-mt", "--methylation_tsv", type=str)
parser.add_argument("-o", "--output_file_name", type=str)
args = parser.parse_args()

haplo_dic = {}
samfile = pysam.AlignmentFile(args.haplotype_bam, "rb")

for read in samfile:
    if read.flag & 3844 == 0:
        hpt = str(read.get_tag('HP')) if read.has_tag('HP') else None
        rname = str(read.qname)
        haplo_dic[rname] = hpt

with open (args.methylation_tsv, 'r') as mc, open (args.output_file_name, 'w') as mh:
    mh.write(mc.readline().strip() + '\tHaplotype\n')
    for line in mc:
        read = line.split('\t')[4]
        mh.write(line.strip() + '\t' + str(haplo_dic[read] if read in haplo_dic else None) + '\n')
