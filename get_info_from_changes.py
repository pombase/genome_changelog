import pandas
import argparse
from genome_functions import make_synonym_dict

class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
parser.add_argument('--coordinate_changes_files', nargs='+', default=['results/all_coordinate_changes_file.tsv','results/pre_svn_coordinate_changes_file.tsv'], help='coordinate changes file (in order that you want them concatenated)')
parser.add_argument('--qualifier_changes_files', nargs='+', default=['results/all_qualifier_changes_file.tsv','gene_changes_comments_and_pmids/pre_svn_qualifier_changes_file.tsv'], help='Qualifier changes tsv files')
parser.add_argument('--genes_in_wrong_chromosomes',default='valid_ids_data/genes_in_wrong_chromosomes.tsv', help='file with wrong chromosomes')
parser.add_argument('--systematic_ids_with_two_CDS',default='valid_ids_data/systematic_ids_associated_with_two_CDS.tsv', help='file with wrong chromosomes')
parser.add_argument('--output_modified_coordinates',default='results/only_modified_coordinates.tsv', help='output tsv file')
parser.add_argument('--output_summary',default='results/genome_changes_summary.tsv', help='output tsv file')

args = parser.parse_args()

# Load coordinate changes
data = pandas.concat([
    pandas.read_csv(f, delimiter='\t', na_filter=False, dtype=str) for f in args.coordinate_changes_files
])

# Main features only (mRNA was only used a few times, gene is a special one used only in the mating type region)
main_features_data = data[~data['feature_type'].isin(["5'UTR","3'UTR",'intron','promoter','LTR', 'misc_feature', 'mRNA','CDS_before','CDS_BEFORE','gene']) & data.systematic_id.str.startswith('SP')].copy()

# Make sure they are sorted properly (below, earliest_modification relies on 'added' rows being on top of 'removed' ones)
main_features_data['date'] = pandas.to_datetime(main_features_data['date'],utc=True)
# We add a special column for all types of RNA, because often in a revision both coordinates and feature_type would be changed
main_features_data['feature_type_temp'] = main_features_data['feature_type'].apply(lambda x: 'RNA' if 'RNA' in x else x)
main_features_data.sort_values(['date','revision','chromosome','systematic_id','feature_type_temp','added_or_removed'],inplace=True, ascending=[False, False,True, True, True, True])
main_features_data['date'] = main_features_data['date'].dt.date

# See the columns that only differ in value and added_removed > coordinates were modified
d = main_features_data.drop(columns=['value', 'added_or_removed', 'primary_name', 'feature_type'])
modification_logical_index = d.duplicated(keep=False)

# See the columns that only differ in value and added_removed > coordinates were modified
modification_data = main_features_data[modification_logical_index].copy().drop(columns=['feature_type_temp'])
modification_data.to_csv(args.output_modified_coordinates, sep='\t', index=False)

## Classify genes by type of change =====================================================================================

# There were some cases in which features were removed accidentally and re-inserted some revisions
# later, so we cannot rely only on those which are added and removed in the same revision.
# Also, some are added and then edited, some are edited and then deleted.

data_subset = main_features_data[['systematic_id','chromosome', 'primary_name', 'added_or_removed', 'date', 'value']].copy()

# We remove known errors
error_data = pandas.read_csv(args.genes_in_wrong_chromosomes, delimiter='\t', na_filter=False, dtype=str)
error_data['known_error'] = True
data_subset = data_subset.merge(error_data,on=['systematic_id', 'chromosome'], how='outer')
data_subset = data_subset[data_subset['known_error'] != True]
data_subset.drop(columns=['known_error'])

data_subset['net_change'] = 0
data_subset['nb_changes'] = 1
data_subset.loc[data_subset['added_or_removed'] == 'added','net_change'] = 1
data_subset.loc[data_subset['added_or_removed'] == 'removed','net_change'] = -1
data_subset['earliest_modification'] = data_subset['net_change'].copy()

# Columns containing latest and earliest change
data_subset['latest_change'] = data_subset['date'].copy()
data_subset['latest_coords'] = data_subset['value'].copy()
data_subset['earliest_change'] = data_subset['date'].copy()

data_subset2 = data_subset.groupby('systematic_id', as_index=False).agg({'net_change': sum, 'nb_changes': sum, 'earliest_modification': lambda x: list(x)[-1], 'chromosome': lambda x: list(x)[0], 'primary_name': lambda x: list(x)[0], 'latest_change': lambda x: list(x)[0], 'earliest_change': lambda x: list(x)[-1], 'latest_coords': lambda x: list(x)[0]})
data_subset2['category'] = ''

data_subset2.loc[(data_subset2['net_change'] == 1) & (data_subset2['nb_changes'] == 1), 'category'] = 'added'
data_subset2.loc[(data_subset2['net_change'] == 1) & (data_subset2['nb_changes'] > 1), 'category'] = 'added_and_changed'
data_subset2.loc[(data_subset2['net_change'] == -1) & (data_subset2['nb_changes'] == 1), 'category'] = 'removed'
data_subset2.loc[(data_subset2['net_change'] == -1) & (data_subset2['nb_changes'] > 1), 'category'] = 'changed_and_removed'

data_subset2.loc[(data_subset2['net_change'] == 0) & (data_subset2['earliest_modification'] == -1), 'category'] = 'changed'
data_subset2.loc[(data_subset2['net_change'] == 0) & (data_subset2['nb_changes'] == 2) & (data_subset2['earliest_modification'] == 1), 'category'] = 'added_and_removed'
data_subset2.loc[(data_subset2['net_change'] == 0) & (data_subset2['nb_changes'] > 3) &  (data_subset2['earliest_modification'] == 1), 'category'] = 'added_changed_and_removed'

# Control for uneven changes
data_subset2['abs_change'] = abs(data_subset2['net_change'])
uneven_changes = data_subset2[data_subset2['abs_change'] > 1]
if not uneven_changes.empty:
    affected_ids = set(uneven_changes['systematic_id'])
    raise ValueError(f'The following ids have too many changes {affected_ids}')

data_subset2.drop(columns=['net_change', 'earliest_modification','nb_changes', 'abs_change'], inplace=True)

## Find merges as removal of feature + added as synonym of another  =====================================================================================

synonym_dict = make_synonym_dict('valid_ids_data/gene_IDs_names.tsv', 'valid_ids_data/obsoleted_ids.tsv', 'valid_ids_data/missing_synonyms.tsv')

merged_logi = data_subset2['category'].str.contains('removed') & data_subset2['systematic_id'].isin(synonym_dict)

# Rename categories where merged happened
data_subset2.loc[merged_logi,'category'] = data_subset2.loc[merged_logi,'category'].apply(lambda x: x.replace('remove','merge'))
data_subset2['merged_into'] = ''
data_subset2.loc[merged_logi,'merged_into'] = data_subset2.loc[merged_logi,'systematic_id'].apply(lambda x: ','.join(sorted(synonym_dict[x])))

# Some systematic ids have been associated with two CDSs in the past, and this messes up the count. This is only relevant if data from the pre-svn folder is used.
if 'svn_2' in main_features_data['revision'].values:
    with open(args.systematic_ids_with_two_CDS) as f:
        multi_cds_ids = [line.rstrip('\n') for line in f]
    data_subset2.loc[data_subset2.systematic_id.isin(multi_cds_ids),'category'] = 'multi_CDS'

# Add extra column indicating what types of feature the systematic_id has ever contained
extra_column_dataset = main_features_data[['systematic_id', 'feature_type']].groupby('systematic_id', as_index=False).agg({'feature_type': lambda x: ','.join(sorted(list(set(x))))})
data_subset2 = data_subset2.merge(extra_column_dataset, on='systematic_id')
data_subset2 = data_subset2[['systematic_id','chromosome','primary_name','feature_type','category','merged_into', 'earliest_change', 'latest_change', 'latest_coords']]
data_subset2.to_csv(args.output_summary, sep='\t', index=False)
