import pandas


all_data = pandas.read_csv('latest_revision/counts.txt', delimiter='\t', names=['counts'], header=None)

all_data['counts'] = all_data['counts'].apply(lambda x: [int(i) for i in x.split(' ')])

all_data['biological_role_published'] = all_data['counts'].apply(lambda x: x[2])
all_data['biological_role_inferred'] = all_data['counts'].apply(lambda x: x[7])
all_data['conserved_unknown'] = all_data['counts'].apply(lambda x: x[10])
all_data['fission_yeast_unknown'] = all_data['counts'].apply(lambda x: x[8] + x[12])

all_data.drop(columns=['counts'], inplace=True)
all_data.to_csv('latest_revision/formatted_counts.tsv', sep='\t', index=False)
