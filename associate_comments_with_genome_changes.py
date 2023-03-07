import pandas


## 1. GENE STRUCTURE CHANGES ############################################

# Load data ===================================================

# Data from pombase comments
changelog_pombase = pandas.read_csv('gene_changes_comments_and_pmids/gene-coordinate-change-data.tsv',sep='\t',na_filter=False)

# Remove known repetitions (see https://github.com/pombase/genome_changelog/issues/7)
changelog_pombase = changelog_pombase[~((changelog_pombase['systematic_id'] == 'SPBC16E9.16c') & (changelog_pombase['date'] == '2007-01-03'))]
changelog_pombase = changelog_pombase[~((changelog_pombase['systematic_id'] == 'SPAC23D3.08') & (changelog_pombase['date'] == '2007-02-05'))]
changelog_pombase = changelog_pombase[~((changelog_pombase['systematic_id'] == 'SPBC4.02c') & (changelog_pombase['date'] == '2008-05-02'))]

# Remove changes that have not been included in the genome (Chr_I:682993!TC->T)
changelog_pombase = changelog_pombase[~((changelog_pombase['systematic_id'] == 'SPAC22F3.11c') & (changelog_pombase['date'] == '2017-04-04'))]

# Fill in empty dates with the closest
changelog_pombase['date'][changelog_pombase['date'] == ''] = pandas.NaT
changelog_pombase['date'] = changelog_pombase['date'].fillna(method='bfill')

# Data generated from get_modifications_on_main_features_only.py
changelog_script = pandas.read_csv('results/only_modified_coordinates.tsv',sep='\t',na_filter=False)
changelog_script['original_index'] = changelog_script.index

# Data about changes in qualifiers
db_xref_script = pandas.concat([
    pandas.read_csv(f, delimiter='\t', na_filter=False, dtype=str) for f in ['results/all_qualifier_changes_file.tsv','gene_changes_comments_and_pmids/pre_svn_qualifier_changes_file.tsv']
])

# Only the relevant columns, only db_xref qualifiers
db_xref_script = db_xref_script[db_xref_script['qualifier_type'] == 'db_xref']
db_xref_script = db_xref_script[['systematic_id', 'revision', 'feature_type', 'value','added_or_removed']]
# Rename value column to db_xref, that will be the final column in the output file
db_xref_script.rename(columns={'value': 'db_xref'}, inplace=True)

def formatting_references(x):
    sorted_list = sorted(list(set([i for i in x if i != ''])))
    return ','.join(sorted_list)


# Combine multiple db_xref changes in one line (comma separated)
unique_identifier_cols = ['systematic_id', 'revision', 'feature_type','added_or_removed']
db_xref_script = db_xref_script.groupby(unique_identifier_cols).agg({'db_xref': formatting_references})

changelog_with_comments = changelog_script.merge(db_xref_script, on=unique_identifier_cols, how='left')

changelog_with_comments['date'] = pandas.to_datetime(changelog_with_comments['date'],utc=True)
changelog_pombase['date'] = pandas.to_datetime(changelog_pombase['date'],utc=True)

changelog_with_comments = changelog_with_comments.sort_values(['date','systematic_id', 'added_or_removed'],ascending=[True, False, True])
changelog_pombase = changelog_pombase.sort_values(['date'],ascending=[True])

# We associate the comments and references of pombase with both addition and removal
changelog_pombase['added_or_removed'] = 'added,removed'
changelog_pombase['added_or_removed'] = changelog_pombase['added_or_removed'].apply(str.split, args=[','])
changelog_pombase = changelog_pombase.explode('added_or_removed')

# A note on this: merge_asof(left, right) finds for every row in left the nearest row in right, so to have a single match between
# comments in changelog_pombase to output_data rows, we have to do it like this
temp_data = pandas.merge_asof(changelog_pombase[['systematic_id','reference','comments', 'date', 'added_or_removed']],changelog_with_comments, by=['systematic_id', 'added_or_removed'], on=['date'], direction='nearest')
# Remove the orphan lines
temp_data = temp_data[~temp_data['revision'].isna()]
# Merge back
changelog_with_comments = pandas.merge(changelog_with_comments,temp_data[['original_index', 'reference','comments']], on='original_index',how='outer')

changelog_with_comments['date'] = changelog_with_comments['date'].dt.date
changelog_with_comments.rename(columns={'reference': 'pombase_reference'}, inplace=True)
changelog_with_comments.rename(columns={'comments': 'pombase_comments'}, inplace=True)
changelog_with_comments = changelog_with_comments.drop_duplicates()
changelog_with_comments.sort_values(['original_index'], inplace=True)
changelog_with_comments.drop(columns=['original_index'], inplace=True)
changelog_with_comments.to_csv('results/only_modified_coordinates_comments.tsv',sep='\t', index=False)



## 2. GENE REMOVAL / ADDITION on genome_changes_summary ############################################

# Associate with comments from pombase for removal / addition

genome_changes_summary = pandas.read_csv('results/genome_changes_summary.tsv', sep='\t', na_filter=False)

comments_new_genes = pandas.read_csv(
    'gene_changes_comments_and_pmids/new-gene-data.tsv', sep='\t', na_filter=False
    ).drop(columns=['date']
    ).groupby(['systematic_id'], as_index=False).agg({'comment_addition': lambda x: '|'.join(i for i in x if i != ''), 'reference_addition': formatting_references})
comments_removed_genes = pandas.read_csv(
    'gene_changes_comments_and_pmids/removed-gene-data.tsv', sep='\t', na_filter=False
    ).drop(columns=['date']
    ).groupby(['systematic_id'], as_index=False).agg({'comment_removal': lambda x: '|'.join(i for i in x if i != ''), 'reference_removal': formatting_references})

# Remove known orphan comments
known_orphans = pandas.read_csv('gene_changes_comments_and_pmids/known_orphan_comments_gene_add_remove.tsv', sep='\t')
orphans_new_genes = set(known_orphans.loc[known_orphans['type']!='removed_only','systematic_id'])

comments_new_genes = comments_new_genes[~comments_new_genes['systematic_id'].isin(orphans_new_genes)].copy()
comments_removed_genes = comments_removed_genes[~comments_removed_genes['systematic_id'].isin(known_orphans['systematic_id'])].copy()


genome_changes_summary = genome_changes_summary.merge(comments_new_genes, on='systematic_id', how='outer')
genome_changes_summary = genome_changes_summary.merge(comments_removed_genes, on='systematic_id', how='outer')

genome_changes_summary.fillna('', inplace=True)

genome_changes_summary.to_csv('results/genome_changes_summary_comments.tsv', sep='\t', index=False)

## 3. Comments on coordinate changes ############################################

for input_file in ['results/all_coordinate_changes_file.tsv', 'results/pre_svn_coordinate_changes_file.tsv']:
    output_file_name = input_file.split('.')[0] + '_comments.tsv'
    all_coordinate_changes = pandas.read_csv(input_file,sep='\t',na_filter=False)
    all_coordinate_changes['revision'] = all_coordinate_changes['revision'].astype(str)
    unique_identifier_cols = ['systematic_id', 'revision', 'feature_type','added_or_removed']
    all_coordinate_changes_with_comments = all_coordinate_changes.merge(db_xref_script, on=unique_identifier_cols, how='left')

    # We include "value" not to duplicate rows in the rare cases where a repeated feature was added or removed
    all_coordinate_changes_with_comments = all_coordinate_changes_with_comments.merge(changelog_with_comments[unique_identifier_cols + ['value', 'pombase_reference', 'pombase_comments']], on=unique_identifier_cols + ['value'], how='left')
    all_coordinate_changes_with_comments.to_csv(output_file_name, sep='\t', index=False)
