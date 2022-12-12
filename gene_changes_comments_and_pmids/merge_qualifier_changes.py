"""
Merge the qualifier_changes data from the release and the repository.
"""

import pandas

data_release = pandas.read_csv('qualifier_changes.tsv',sep='\t',na_filter=False)
data_repository = pandas.read_csv('../all_qualifier_changes_file.tsv',sep='\t',na_filter=False)

out_data = pandas.concat([data_repository,data_release])
out_data.drop_duplicates(inplace=True)
out_data.to_csv('qualifier_changes.tsv',sep='\t', index=False)