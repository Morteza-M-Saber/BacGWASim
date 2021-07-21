import pandas as pd

input_phen = snakemake.input.phen
output_pickle = snakemake.output.pickle

#Convert phenotypes to matrix
df = pd.read_csv(input_phen, sep=' ', header=None, index_col=0)
df.columns = ['Sample','phenotype','nan']
df.drop(['Sample','nan'], axis=1, inplace=True)
df.to_pickle(output_pickle)