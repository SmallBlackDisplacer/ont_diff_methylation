import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file_name", type=str)
parser.add_argument("-o", "--output_file", type=str)
args = parser.parse_args()

haplo_df = pd.read_csv(args.file_name, sep='\t')

haplo_df['is_methylated'] = list(map(lambda x: 1 if x > 2.5 else 0, haplo_df.log_lik_ratio))

minimalistic = haplo_df[['chromosome', 'start', 'end', 'Haplotype', 'is_methylated']].query(
    'Haplotype != "None"').groupby(['chromosome', 'start', 'end', 'Haplotype', 'is_methylated']).size().reset_index(
    name='Count')

minimalistic['Haplotype'] = minimalistic['Haplotype'].replace('1', 'hp_1').replace('2', 'hp_2')
minimalistic['is_methylated'] = minimalistic['is_methylated'].replace(1, 'met_').replace(0, 'unm_')

minimalistic_pivot = minimalistic.pivot_table(values="Count", index=['chromosome', 'start', 'end'],
                                              columns=['is_methylated', 'Haplotype'], fill_value=0).reset_index()

ind = pd.Index([str(e[0]) + str(e[1]) for e in minimalistic_pivot.columns.tolist()])
minimalistic_pivot.columns = ind

minimalistic_pivot.to_csv(args.output_file, sep='\t', index=False)
