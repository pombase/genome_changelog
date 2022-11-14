import pandas

# Load coordinate changes
data = pandas.concat([
    pandas.read_csv('all_coordinate_changes_file.tsv', delimiter='\t', na_filter=False),
    pandas.read_csv('pre_svn_coordinate_changes_file.tsv', delimiter='\t', na_filter=False)
])
# We only consider CDS features in genes that have alleles with sequence errors
output_data = data[~data['feature_type'].isin(["5'UTR","3'UTR",'intron','promoter','LTR', 'misc_feature'])]

# See the columns that only differ in value and added_removed > coordinates were modified
d = output_data.drop(columns=['value', 'added_or_removed'])
logi = d.duplicated(keep=False)

output_data[logi].to_csv('only_modified_coordinates.tsv', sep='\t', index=False)
