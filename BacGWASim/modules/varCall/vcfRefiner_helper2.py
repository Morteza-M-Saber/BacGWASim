import pandas as pd

input_vcf = snakemake.input.vcf
output_pickle = snakemake.output.pickle

info_line_c = 0
with open(input_vcf, 'r') as file:
    line = file.readline()
    while line[:2] == '##':
        info_line_c += 1
        line = file.readline()
        
df = pd.read_csv(input_vcf, sep='\t', skiprows=info_line_c)
df = df.drop(['#CHROM', 'POS', 'REF', 'ALT', 'QUAL',
              'FILTER', 'INFO', 'FORMAT'], axis=1)
df.set_index('ID', inplace=True)
df = df.transpose()
df.to_pickle(output_pickle)
