import os

import matplotlib.pyplot as plt
import palettable as pal
import numpy as np
import pandas as pd
from snakemake import shell

from plotting_utils import mm2inch, rcParams, biotype_colors
for k, v in rcParams.items():
    plt.rcParams[k] = v


def getting_biotype_hits(x):
    hits = list()
    for biotype in sorted(list(set(x))):
        hits.append(
            str(list(x).count(biotype)) + '-' + biotype
        )
    return '|'.join(hits)


# Loading parameters
score = int(snakemake.wildcards.score)/100
graph_path = snakemake.params.graph_path.format(**snakemake.wildcards)
biotypes = [
    'protein_coding', 'lncRNA', 'pseudogene', 'snRNA', 'miRNA',
    'snoRNA', 'rRNA_pseudogene', 'rRNA', 'misc_RNA'
]

# Creating folder
shell("mkdir {dir}".format(dir=os.path.dirname(graph_path)))

# Loading data
data = pd.read_csv(snakemake.input.data, sep='\t')
data['edge_score'] = data[['qseq_plen', 'sseq_plen']].mean(axis=1)
data = data[data['edge_score'] >= score]

# Loading GTF
gtf = pd.read_csv(
    snakemake.input.gtf, sep='\t', skiprows=[0, 1, 2, 3, 4],
    usecols=[2, 8], names=['feature', 'attributes'],
)
gtf = gtf[gtf['feature'] == 'gene']
gtf['gene_id'] = gtf['attributes'].str.extract('gene_id \"(.*?)\"')
gtf['biotype'] = gtf['attributes'].str.extract('biotype \"(.*?)\"')
gtf['gene_name'] = gtf['attributes'].str.extract('gene_name \"(.*?)\"')

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

# Illustrating connections between things
biotype_dict = gtf.set_index('gene_id')['biotype'].to_dict()
name_dict = gtf.set_index('gene_id')['gene_name'].to_dict()
data['qgene_biotype'] = data['qgene'].map(biotype_dict)
data['sgene_biotype'] = data['sgene'].map(biotype_dict)
data['qgene_name'] = data['qgene'].map(name_dict)
data['sgene_name'] = data['sgene'].map(name_dict)

hits_biotypes = data[['qgene', 'sgene_biotype']].groupby('qgene')['sgene_biotype'].apply(getting_biotype_hits).reset_index()
hits_biotypes['qgene_biotype'] = hits_biotypes['qgene'].map(biotype_dict)

# Calculating global links between biotypes
per_bio = data[data['qgene_biotype'].isin(biotypes) & data['sgene_biotype'].isin(biotypes)]
per_bio = per_bio.groupby(['qgene_biotype', 'sgene_biotype']).count()['sgene'].reset_index()
per_bio_dict = per_bio.groupby('qgene_biotype').sum()['sgene'].to_dict()
per_bio['norm'] = per_bio['qgene_biotype'].map(per_bio_dict)
per_bio['perc_sgene'] = per_bio['sgene'] / per_bio['norm'] * 100

# Creating barchart datastructure
_barc = dict()
for biotype in biotypes:
    _barc[biotype] = list()

for biotype in biotypes:
    for sub_biotype in biotypes:
        _hit = per_bio[
            (per_bio['qgene_biotype'] == sub_biotype) &
            (per_bio['sgene_biotype'] == biotype)
        ]

        if len(_hit) == 0:
            _barc[biotype].append(0)
        else:
            _barc[biotype].append(_hit['perc_sgene'].tolist()[0])


def filtering_hits_biotypes(df_bio, biotype):
    df_bio = df_bio[df_bio['qgene_biotype'] == biotype]
    df_bio = df_bio.groupby(['qgene_biotype', 'sgene_biotype']).count().reset_index()
    df_bio.sort_values(by=['qgene_biotype', 'qgene'], inplace=True)
    df_bio['qgene'] = df_bio['qgene'] / np.sum(df_bio['qgene'])*100
    return df_bio


for biotype in set(hits_biotypes['qgene_biotype']):
    _df = filtering_hits_biotypes(hits_biotypes, biotype)
    _df.to_csv(graph_path.format(biotype=biotype))


# ---------------- GRAPHING -------------------#
# Biotype hit repartition
colors = pal.colorbrewer.sequential.GnBu_8.mpl_colors
bins = [-0.1, 0, 1, 2, 5, 10, 20, 50, 1000]
biotypes_bins = pd.DataFrame(columns=biotypes)
for biot in biotypes:
    biot_gtf = gtf[gtf['biotype'] == biot]
    biot_len = len(biot_gtf)
    biot_hist = biot_gtf.groupby('hits')['gene_id'].count()/biot_len
    biot_binned = biot_hist.groupby(pd.cut(biot_hist.index, bins)).sum()
    biotypes_bins[biot] = biot_binned

# Plotting axis 0
fig, axes = plt.subplots(ncols=2, nrows=1, figsize=mm2inch((178,100)))
sum_bins = np.zeros(len(biotypes))
for i, (index, row) in enumerate(biotypes_bins.iterrows()):
    axes[0].barh(biotypes, row, left=sum_bins, color=colors[i], edgecolor='k')
    sum_bins += row
legend = biotypes_bins.index.tolist()
for i in range(3):
    legend[i] = i
axes[0].legend(legend)

# Plotting axis 1
sum_bins = np.zeros(len(biotypes))
for biotype in biotypes:
    axes[1].barh(biotypes, _barc[biotype], left=sum_bins, color=biotype_colors[biotype])
    sum_bins += _barc[biotype]

axes[0].invert_yaxis()
axes[1].invert_yaxis()
plt.tight_layout()
plt.savefig(snakemake.output.bar_plot)
