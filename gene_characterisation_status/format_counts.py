import pandas
import matplotlib.pyplot as plt
import datetime

all_data = pandas.read_csv('results/merged_counts.tsv', delimiter='\t')
all_data['date'] = pandas.to_datetime(all_data['date'], utc=True).dt.date

all_data['counts'] = all_data['counts'].apply(lambda x: [int(i) for i in x.split(' ')])

all_data['published'] = all_data['counts'].apply(lambda x: x[2])
all_data['inferred_role'] = all_data['counts'].apply(lambda x: x[7])
all_data['conserved_unknown'] = all_data['counts'].apply(lambda x: x[10])
all_data['fission_yeast_unknown'] = all_data['counts'].apply(lambda x: x[8] + x[12])

all_data.drop(columns=['counts'], inplace=True)
all_data.to_csv('results/formatted_counts.tsv', sep='\t', index=False)

# Plot a stackplot of the counts
plt.figure()
plt.stackplot(all_data['date'], all_data['published'], all_data['inferred_role'], all_data['conserved_unknown'], all_data['fission_yeast_unknown'], labels=['Published', 'Inferred role', 'Conserved unknown', 'Fission yeast unknown'])

plt.xlim([datetime.date(2003,6,1), datetime.date.today()])
# place the legend outside the plot
plt.legend(loc='lower center', bbox_to_anchor=(0.5, 1), ncol=2)
plt.tight_layout()
plt.savefig('results/figure.svg')