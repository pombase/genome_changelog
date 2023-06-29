"""
Create new coordinate changes files where changes to type of RNA are ignored
"""
import pandas


input_files = ['results/all_coordinate_changes_file_comments.tsv',
                'results/pre_svn_coordinate_changes_file_comments.tsv',
                'results/only_modified_coordinates_comments.tsv',
]

chromosome_dictionary = {
    'chromosome1': 'I',
    'chromosome2': 'II',
    'chromosome3': 'III',
    'mating_type_region': 'mating_type_region',
    'pMIT': 'mitochondrial',
    'telomeric': 'telomeric',
}

genome_changes = pandas.read_csv('results/genome_sequence_changes.tsv',sep='\t',na_filter=False)
genome_changes['chromosome'] = genome_changes['chromosome'].apply(lambda x : chromosome_dictionary[x])
genome_changes['combined_column'] = genome_changes.apply(lambda r: [r['revision'], r['chromosome']], axis=1)

for input_file in input_files:
    output_file = input_file.split('.')[0] + '_no_type_change.tsv'
    data = pandas.read_csv(input_file,sep='\t',na_filter=False)
    data['original_index'] = data.index
    rna_data = data[data['feature_type'].str.contains(r'(?i)RNA')].copy()
    rna_data_where_type_changed = rna_data[rna_data[['systematic_id', 'revision', 'value']].duplicated(keep=False)].copy()
    rna_data_where_type_changed['combined_column'] = rna_data_where_type_changed.apply(lambda r: [r['revision'], r['chromosome']], axis=1)
    # We exclude revisions where genome changed
    rna_data_where_type_changed = rna_data_where_type_changed[~rna_data_where_type_changed['combined_column'].isin(genome_changes['combined_column'])]

    data = data[~data['original_index'].isin(set(rna_data_where_type_changed['original_index']))]
    data.drop(columns='original_index', inplace=True)
    data.to_csv(output_file, sep='\t', index=False)