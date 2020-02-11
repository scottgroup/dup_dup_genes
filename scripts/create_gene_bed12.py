import pandas as pd

"""
BED12 format
chrom | chromStart | chormEnd | name | score | strand |
thickStart | thickEnd | itemRgb | blockCounts | blockSizes | blockStarts
"""

col_names_gtf = [
    'seqname', 'source', 'feature', 'start',
    'end', 'score', 'strand', 'frame', 'attributes'
]

col_names_gene_bed12 = [
    'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart',
    'thickEnd', 'itemRgb', 'blockCounts', 'blockSizes', 'blockStarts'
]

col_dtypes = {
    'seqname': 'string',
    'source': 'string',
    'feature': 'string',
    'start': 'Int64',
    'end': 'Int64',
    'score': 'string',
    'strand': 'string',
    'frame': 'string',
    'attributes': 'string',
}

# Reading gtf_gene file
gtf = pd.read_csv(
    snakemake.input.gtf, sep='\t', skiprows=[0,1,2,3,4],
    names=col_names_gtf, dtype=col_dtypes
)

# # Correcting for offset
# gtf['start'] -= 1

# Extracting gene_id
gtf['gene_id'] = gtf['attributes'].str.extract('gene_id \"(.*?)\"').astype("string")
gtf['transcript_id'] = gtf['attributes'].str.extract('transcript_id \"(.*?)\"').astype("string")

# Extracting transcripts ID
exon_gtf = gtf[gtf['feature'] == 'exon']
gene_tr = exon_gtf.groupby(['gene_id', 'transcript_id']).count()['seqname'].reset_index()
gene_tr = gene_tr[gene_tr.groupby(['gene_id'])['seqname'].transform(max) == gene_tr['seqname']].drop_duplicates(subset='gene_id', keep='first')
gene_tr = gene_tr[gene_tr['seqname'] > 1]
transcripts = gene_tr['transcript_id'].tolist()

# GTF with selected transcripts and exon
ex_tr_gtf = gtf[gtf['transcript_id'].isin(transcripts) & gtf['feature'].isin(['exon', 'transcript'])]

# Split into two gtf
tr_gtf = ex_tr_gtf[ex_tr_gtf['feature'] == 'transcript']
ex_gtf = ex_tr_gtf[ex_tr_gtf['feature'] == 'exon']

# Correcting ex_gtf to remove the transcript start
tr_start = tr_gtf.set_index('transcript_id')['start'].to_dict()
ex_gtf['tr_start'] = ex_gtf['transcript_id'].map(tr_start)
ex_gtf['corrected_start'] = ex_gtf['start'] - ex_gtf['tr_start']

# Extract list of values for exons
ex_gtf.loc[:, 'ex_len'] = ex_gtf['end'] - (ex_gtf['start'] - 1)
exon_count_dict = ex_gtf.groupby('transcript_id')['start'].count()
# Splitting for strand
pos_ex_gtf = ex_gtf[ex_gtf['strand'] == '+']
neg_ex_gtf = ex_gtf[ex_gtf['strand'] == '-']
exon_starts_dict = {
    **pos_ex_gtf.groupby('transcript_id')['corrected_start'].apply(list).apply(lambda x: ','.join([str(num) for num in x])).to_dict(),
    **neg_ex_gtf.groupby('transcript_id')['corrected_start'].apply(list).apply(lambda x: ','.join([str(num) for num in x[::-1]])).to_dict()
}
exon_lens_dict = {
    **pos_ex_gtf.groupby('transcript_id')['ex_len'].apply(list).apply(lambda x: ','.join([str(num) for num in x])).to_dict(),
    **neg_ex_gtf.groupby('transcript_id')['ex_len'].apply(list).apply(lambda x: ','.join([str(num) for num in x[::-1]])).to_dict()
}

# tr_gtf.set_index('transcript_id', inplace=True)
tr_gtf.loc[:, 'ex_start'] = tr_gtf['transcript_id'].map(exon_starts_dict)
tr_gtf.loc[:, 'ex_len'] = tr_gtf['transcript_id'].map(exon_lens_dict)
tr_gtf.loc[:, 'ex_count'] = tr_gtf['transcript_id'].map(exon_count_dict)

# Creating BED12 dataframe for genes
gtf_gene = gtf[gtf['feature'] == 'gene']
gene_bed = pd.DataFrame(index=gtf_gene.index)
gene_bed['chrom'] = gene_bed.index.map(gtf_gene['seqname'].to_dict())
gene_bed['chromStart'] = gene_bed.index.map(gtf_gene['start'].to_dict()) - 1
gene_bed['chromEnd'] = gene_bed.index.map(gtf_gene['end'].to_dict())
gene_bed['name'] = gene_bed.index.map(gtf_gene['gene_id'].to_dict())
gene_bed['score'] = 0
gene_bed['strand'] = gene_bed.index.map(gtf_gene['strand'].to_dict())
gene_bed['thickStart'] = 0
gene_bed['thickEnd'] = 0
gene_bed['itemRgb'] = 0
gene_bed['blockCounts'] = 1
gene_bed['blockSizes'] = gene_bed['chromEnd'] - gene_bed['chromStart']
gene_bed['blockStarts'] = 0

# Creating BED12 DataFrame for transcripts
tr_bed = pd.DataFrame(index=tr_gtf.index)
tr_bed['chrom'] = tr_bed.index.map(tr_gtf['seqname'].to_dict())
tr_bed['chromStart'] = tr_bed.index.map(tr_gtf['start'].to_dict()) - 1
tr_bed['chromEnd'] = tr_bed.index.map(tr_gtf['end'].to_dict())
tr_bed['name'] = tr_bed.index.map(tr_gtf['transcript_id'].to_dict())
tr_bed['score'] = 0
tr_bed['strand'] = tr_bed.index.map(tr_gtf['strand'].to_dict())
tr_bed['thickStart'] = 0
tr_bed['thickEnd'] = 0
tr_bed['itemRgb'] = 0
tr_bed['blockCounts'] = tr_bed.index.map(tr_gtf['ex_count'].to_dict())
tr_bed['blockSizes'] = tr_bed.index.map(tr_gtf['ex_len'].to_dict())
tr_bed['blockStarts'] = tr_bed.index.map(tr_gtf['ex_start'].to_dict())

# Saving to file
bed = pd.concat([tr_bed, gene_bed])
bed.to_csv(snakemake.output.bed12, sep='\t', index=False, header=False)
