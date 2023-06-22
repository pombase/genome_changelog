import pandas


pre_svn_counts = pandas.read_csv('pre_svn_counts.tsv', delimiter='\t', header=None, names=['revision', 'counts'])
pre_svn_counts['date'] = pre_svn_counts['revision'].apply(lambda x: f'{str(x)[:4]}-{str(x)[4:6]}-{str(x)[6:8]}')

svn_counts = pandas.read_csv('svn_counts.tsv', delimiter='\t', header=None, names=['revision', 'counts'])
revisions_monthly = pandas.read_csv('revisions_monthly.tsv', delimiter='\t')

svn_counts = svn_counts.merge(revisions_monthly, on='revision')

all_data = pandas.concat([pre_svn_counts, svn_counts]).sort_values('date')

all_data.to_csv('results/formatted_counts.tsv', sep='\t', index=False)