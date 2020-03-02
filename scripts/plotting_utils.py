def mm2inch(*tupl):
    """
    https://stackoverflow.com/questions/14708695/specify-figure-size-in-centimeter-in-matplotlib
    """
    inch = 25.4
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

# Figure params
rcParams = {
    'svg.fonttype': 'none',
    'font.size': 6,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial']
}

# Colors
biotype_colors = {
    'protein_coding': '#f64f40ff',
    'lncRNA': '#6d9fc9ff',
    'pseudogene': '#7ccc6aff',
    'snRNA': '#ab7bbbff',
    'miRNA': '#fdab39ff',
    'snoRNA': '#ffff37ff',
    'rRNA_pseudogene': '#c6aa70ff',
    'rRNA': '#f5419eff',
    'misc_RNA': '#9b9b9bff'
}
