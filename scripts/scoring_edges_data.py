import matplotlib.pyplot as plt
import palettable as pal
import numpy as np
import pandas as pd

score = int(snakemake.wildcards.score)/100

# Loading data
data = pd.read_csv(snakemake.input.data, sep='\t')
data['edge_score'] = data[['qseq_plen', 'sseq_plen']].mean(axis=1)
data = data[data['edge_score'] >= score]

# Loading GTF
gtf = pd.read_csv(
    snakemake.input.gtf, sep='\t', skiprows=[0,1,2,3,4],
    usecols=[2, 8], names=['feature', 'attributes'],
)
gtf = gtf[gtf['feature'] == 'gene']
gtf['gene_id'] = gtf['attributes'].str.extract('gene_id \"(.*?)\"').astype("string")
gtf['biotype'] = gtf['attributes'].str.extract('biotype \"(.*?)\"').astype("string")
gtf['gene_name'] = gtf['attributes'].str.extract('gene_name \"(.*?)\"').astype("string")

# Merging similar biotypes
biotype = {bio:bio for bio in gtf['biotype'].values}
as_pseudogenes = [
    "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene",
    "transcribed_unprocessed_pseudogene",  "translated_unprocessed_pseudogene",
    "translated_processed_pseudogene", "unprocessed_pseudogene",
    "unitary_pseudogene", "polymorphic_pseudogene", "processed_pseudogene"
]
as_IG_gene = [
    "IG_C_gene", "IG_C_pseudogene", "IG_D_gene", "IG_J_gene",
    "IG_J_pseudogene", "IG_V_gene", "IG_V_pseudogene", "IG_pseudogene"
]
as_TR_gene = [
    "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_J_pseudogene", "TR_V_gene", "TR_V_pseudogene"
]
for pseudo in as_pseudogenes:
    biotype[pseudo] = "pseudogene"
for pseudo in as_IG_gene:
    biotype[pseudo] = "IG_gene"
for pseudo in as_TR_gene:
    biotype[pseudo] = "TR_gene"
gtf['biotype'] = gtf['biotype'].map(biotype)

# Checkings
qgene = data.groupby('qgene')['sgene'].count().to_dict()
sgene = data.groupby('sgene')['qgene'].count().to_dict()

gtf['qgene'] = gtf['gene_id'].map(qgene)
gtf['sgene'] = gtf['gene_id'].map(sgene)
gtf.fillna(0, inplace=True)
gtf['hits'] = gtf['qgene'] + gtf['sgene']



print('Genes per biotypes')
print(gtf.groupby('biotype')['hits'].count())

print('Mean hit per biotypes')
print(gtf.groupby('biotype')['hits'].mean())



# Biotype hit repartition
colors = pal.colorbrewer.sequential.GnBu_8.mpl_colors
bins = [-0.1, 0, 1, 2, 5, 10, 20, 50, 1000]
biotypes = ['protein_coding', 'lncRNA', 'pseudogene', 'snRNA', 'miRNA', 'snoRNA', 'rRNA_pseudogene', 'rRNA', 'misc_RNA']
biotypes_bins = pd.DataFrame(columns=biotypes)
for biot in biotypes:
    biot_gtf = gtf[gtf['biotype'] == biot]
    biot_len = len(biot_gtf)
    biot_hist = biot_gtf.groupby('hits')['gene_id'].count()#/biot_len
    biot_binned = biot_hist.groupby(pd.cut(biot_hist.index, bins)).sum()
    biotypes_bins[biot] = biot_binned

# Plotting
sum_bins = np.zeros(len(biotypes))
for i, (index, row) in enumerate(biotypes_bins.iterrows()):
    plt.bar(biotypes, row, bottom=sum_bins, color=colors[i], edgecolor='k')
    sum_bins += row
legend = biotypes_bins.index.tolist()
for i in range(3):
    legend[i] = i
plt.legend(legend)
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()
