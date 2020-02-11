import pandas as pd

col_names_blast = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
    'qend', 'sstart', 'send', 'evalue', 'bitscore'
]

# Importing sequence lengths
seq_len = pd.read_csv(
    snakemake.input.seq_len, sep='\t', names=['seq', 'len'], index_col=0
)['len'].to_dict()

# Importing BLAST results
df = pd.read_csv(snakemake.input.blast, sep='\t', names=col_names_blast)
df['qseqid'] = df['qseqid'].str.split(':').str[0]
df['sseqid'] = df['sseqid'].str.split(':').str[0]

# Removing hit on the same gene
df = df[df['qseqid'] != df['sseqid']]

# Adding adding percentage of overlap
df['qseq_plen'] = df['length'] / df['qseqid'].map(seq_len) * df['pident'] / 100
df['sseq_plen'] = df['length'] / df['sseqid'].map(seq_len) * df['pident'] / 100

# Keeping minimum threshold each side
threshold = 0.33
df = df[
    (df['qseq_plen'] > threshold) &
    (df['sseq_plen'] > threshold)
]

# Removing duplicates A-B B-A from qseqid and sseqid
df['duplicated'] = df[['qseqid','sseqid']].apply(lambda x: '_'.join(sorted(list(set(x)))), axis=1)
df.drop_duplicates(subset=['duplicated'], inplace=True)
df.drop(columns='duplicated', inplace=True)

# Keeping only pertinent columns
df = df[['qseqid', 'sseqid', 'qseq_plen', 'sseq_plen']]

# Saving to file
df.to_csv(snakemake.output.parsed, sep='\t', index=False)
