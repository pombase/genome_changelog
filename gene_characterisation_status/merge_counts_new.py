import pandas

svn_counts = pandas.read_csv('svn_counts.tsv', delimiter='\t', header=None, names=['revision', 'counts'])
revisions_monthly = pandas.read_csv('revisions_monthly.tsv', delimiter='\t')

svn_counts = svn_counts.merge(revisions_monthly, on='revision', how='left')
existing_table = pandas.read_csv('results/merged_counts.tsv', sep='\t')

all_data = pandas.concat([existing_table, svn_counts]).sort_values('date')

all_data.to_csv('results/merged_counts.tsv', sep='\t', index=False)