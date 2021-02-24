import pandas as pd

input_file = snakemake.input.gcta_phen
output_file = snakemake.output.pickle

df = pd.read_csv(input_file, sep=' ', header=None, index_col=0)
df.columns = ['Sample', 'phenotype', 'nan']
df.drop(['Sample', 'nan'], axis=1, inplace=True)
df.to_pickle(output_file)
