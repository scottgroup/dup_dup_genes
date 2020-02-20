import numpy as np
import pandas as pd

score = int(snakemake.wildcards.score)/100

# Loading data
data = pd.read_csv(snakemake.input.data, sep='\t')

# Loading GTF
gtf = pd.read_csv(
    snakemake.input.gtf, sep='\t', skiprows=[0,1,2,3,4],
    usecols=[2, 8], names=['feature', 'attributes'],
)
gtf = gtf[gtf['feature'] == 'gene']
gtf['gene_id'] = gtf['attributes'].str.extract('gene_id \"(.*?)\"').astype("string")
gtf['biotype'] = gtf['attributes'].str.extract('biotype \"(.*?)\"').astype("string")
gtf['gene_name'] = gtf['attributes'].str.extract('gene_name \"(.*?)\"').astype("string")

# Creating mean score edge
data['edge_score'] = data[['qseq_plen', 'sseq_plen']].mean(axis=1)

# Adding nn liaison columns
data['link'] = 'nn'

# Creating node details
node_df = gtf[['gene_id', 'gene_name', 'biotype']]

# Creating edges details
edge_df = data[['qgene', 'link', 'sgene', 'edge_score']]
edge_df = edge_df[edge_df['edge_score'] >= score]

print(edge_df)

# Saving to files
node_df.to_csv(snakemake.output.node, sep='\t', index=False)
edge_df.to_csv(snakemake.output.edge, sep='\t', index=False)
