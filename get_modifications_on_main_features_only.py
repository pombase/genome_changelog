import pandas
import argparse

class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
parser.add_argument('--input_files', nargs='+', help='coordinate changes file (in order that you want them concatenated)')
parser.add_argument('--output_file', help='output tsv file')
args = parser.parse_args()

# Load coordinate changes
data = pandas.concat([
    pandas.read_csv(f, delimiter='\t', na_filter=False) for f in args.input_files
])

# We only consider CDS features in genes that have alleles with sequence errors
output_data = data[~data['feature_type'].isin(["5'UTR","3'UTR",'intron','promoter','LTR', 'misc_feature'])]

# See the columns that only differ in value and added_removed > coordinates were modified
d = output_data.drop(columns=['value', 'added_or_removed'])
logi = d.duplicated(keep=False)

output_data = output_data[logi].copy()
output_data.to_csv(args.output_file, sep='\t', index=False)
