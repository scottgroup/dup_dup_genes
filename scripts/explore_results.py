import pandas as pd
import numpy as np
import json

def sets2intersecting_iter(data, n=1000):

    hits = list()
    i=0
    for i in range(len(data)//n):
        hits.extend(data[i*n:(i+1)*n])
        hits = sets2intersecting(hits)
    hits.extend(data[(i+1)*n:])
    hits = sets2intersecting(hits)
    return hits


def sets2intersecting(data):
    # Items
    items = sorted(list(set([g for hit in data for g in hit])))

    # Finding M and N
    M = len(data)
    N = len(items)
    print(M, N)

    # Items dictionary
    idict = {k:v for k, v in zip(items, range(N))}
    iidict = {v:k for k, v in idict.items()}

    # Creating the MN bool array
    array = np.zeros((M, N), dtype=np.bool)

    # Filling the MN array with original sets
    for i in range(M):
        for j in data[i]:
            array[i, idict[j]] = True

    # Combining rows
    for i in range(N):
        filt = array[:, i]
        if sum(filt) > 1:
            array[filt, :] = np.any(array[filt], axis=0)

    # Keep unique lines
    bool_patterns = np.unique(array, axis=0)

    # Convert to list of integers
    ints = np.arange(0, N)
    patterns = list()
    for i in range(bool_patterns.shape[0]):
        inds = list(set(ints[bool_patterns[i, :]]))
        patterns.append([iidict[item] for item in inds])
    return patterns


# Loading DATA
data = pd.read_csv(snakemake.input.data, sep='\t')

gene_pairs = [tuple(x) for x in data[['qgene', 'sgene']].to_numpy()]

patterns = sets2intersecting_iter(gene_pairs)
dict_patterns = {k:v for k, v in zip(range(len(patterns)), patterns)}

json = json.dumps(dict_patterns)
f = open("dict.json","w")
f.write(json)
