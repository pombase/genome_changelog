import pandas
import datetime
import subprocess
import os

data = pandas.read_csv('revisions.tsv', sep='\t', na_filter=False, names=['revision', 'date'])[::-1]
# Exclude first revision
data = data.iloc[1:]

lower_limit = datetime.date(2011, 8, 1)

first_days_of_month = list()

for year in range(2011, datetime.date.today().year + 1):
    for month in range(1, 13):
        first_day_of_month = datetime.date(year, month, 1)
        if first_day_of_month < lower_limit:
            continue
        if first_day_of_month > datetime.date.today():
            break
        first_days_of_month.append(first_day_of_month)
    else:
        continue
    break

data.date = pandas.to_datetime(data.date)

# Format first_day_of_month to string
data2 = pandas.DataFrame(list(map(lambda x: x.strftime('%Y-%m-%d'), first_days_of_month)), columns=['date'])
data['date'] = pandas.to_datetime(data['date'],utc=True)
data2['date'] = pandas.to_datetime(data2['date'],utc=True)

data = pandas.merge_asof(data2, data, on=['date'], direction='nearest')
data['date'] = data['date'].dt.date

data.to_csv('revisions_monthly.tsv', sep='\t', index=False)

# Download them if they are new
existing_dates = set(pandas.read_csv('results/formatted_counts.tsv', sep='\t')['date'])

if not os.path.isdir('data'):
    os.mkdir('data')

for date, revision in zip(data['date'], data['revision']):
    # We check for the date and not the revision, because mergeasof can give different revisions for same date
    if str(date) in existing_dates:
        continue
    if not os.path.isdir(f'data/{revision}'):
        os.mkdir(f'data/{revision}')

    for contig in 'chromosome1 chromosome2 chromosome3 mating_type_region pMIT'.split(' '):
        if os.path.isfile(f'data/{revision}/{contig}.contig'):
            continue
        print(f'downloading data/{revision}/{contig}...')

        url = f'https://curation.pombase.org/pombe-embl-repo/trunk/{contig}.contig'
        print(url)
        proc = subprocess.Popen([f'svn export --force -r {revision} {url} data/{revision}/{contig}.contig'], shell=True)
        proc.wait()
        proc.terminate()
