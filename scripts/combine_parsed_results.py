import pandas as pd

def remove_duplicates(df, columns):
    """ """
    df['duplicated'] = df[columns].apply(lambda x: '_'.join(sorted(list(set(x)))), axis=1)
    df.drop_duplicates(subset=['duplicated'], inplace=True)
    df.drop(columns='duplicated', inplace=True)
    return df

# Loading GTF and creating dictionary
gtf = pd.read_csv(
    snakemake.input.gtf, sep='\t', skiprows=[0,1,2,3,4],
    usecols=[2,8], names=['feature', 'attributes']
)
gtf = gtf[gtf['feature'] == 'transcript']
gtf['gene_id'] = gtf['attributes'].str.extract('gene_id \"(.*?)\"').astype("string")
gtf['transcript_id'] = gtf['attributes'].str.extract('transcript_id \"(.*?)\"').astype("string")
tr_dict = {
    **gtf.set_index('transcript_id')['gene_id'].to_dict(),
    **{gene:gene for gene in gtf['gene_id'].tolist()}
}

# Combine all DF
data = pd.DataFrame(columns=['qseqid', 'sseqid', 'qseq_plen', 'sseq_plen'])
for file in snakemake.input.files:
    _df = pd.read_csv(file, sep='\t')
    data = pd.concat([data, _df])


# Removing duplicates A-B B-A from qseqid and sseqid
data = remove_duplicates(data, ['qseqid', 'sseqid'])

# Remove if gene is with itself
data['qgene'] = data['qseqid'].map(tr_dict)
data['sgene'] = data['sseqid'].map(tr_dict)
data = data[data['qgene'] != data['sgene']]

# Removing duplicates from transcripts
data = remove_duplicates(data, ['qgene', 'sgene'])

# To file
data.to_csv(snakemake.output.results, sep='\t', index=False)
