import pandas

groups = [
    ('gene-coordinate-change-data.tsv', 'comments', 'reference'),
    ('new-gene-data.tsv', 'comment_addition', 'reference_addition'),
    ('removed-gene-data.tsv', 'comment_removal', 'reference_removal'),
]

for file_name, comment_column, reference_column in groups:

    data = pandas.read_csv(file_name, sep='\t', na_filter=False)

    with_solexa = data[comment_column].str.contains('Solexa|solexa',regex=True)


    def formatting_function(ref):
        if 'PMID:18488015' in ref:
            return ref
        if ref != '':
            return ref + ',PMID:18488015'
        return 'PMID:18488015'

    data.loc[with_solexa, reference_column] = data.loc[with_solexa, reference_column].apply(formatting_function)
    data.to_csv(file_name, sep='\t', index=False)