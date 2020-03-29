import pysam
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-hb", "--haplotype_bam", type=str)
parser.add_argument("-mt", "--methylation_tsv", type=str)
parser.add_argument("-o", "--output_file_name", type=str)
args = parser.parse_args()

haplo_dic = {}
samfile = pysam.AlignmentFile(args.haplotype_bam, "rb")

for read in samfile:
    if read.is_secondary:
        continue
    hpt = str(read.get_tag('HP')) if read.has_tag('HP') else None
    rname = str(read.qname)
    if rname in haplo_dic:
        if not hpt:
            continue
        if not haplo_dic[rname]:
            haplo_dic[rname] = hpt
        elif haplo_dic[rname] != hpt:
            haplo_dic[rname] = None
    else:
        haplo_dic[rname] = hpt

met_calls = pd.read_csv(args.methylation_tsv, sep="\t")

Haplotype = []
for read in met_calls.read_name:
    Haplotype.append(haplo_dic[read] if read in haplo_dic else None)
met_calls['Haplotype'] = Haplotype

met_calls.to_csv(args.output_file_name, sep='\t')
