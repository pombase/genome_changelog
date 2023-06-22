import requests
import os
import glob

main_url = 'https://www.pombase.org/data/genome_sequence_and_features/artemis_files/OLD/'

contigs = ['chromosome1','chromosome2','chromosome3','mating_type_region','pMIT', 'telomeric']

with open('results/pre_svn_folder_list.tsv') as ins:

    ftp_folders = [l.strip() for l in ins.readlines()]

revision_info = list()

for ftp_folder in ftp_folders:
    date = ftp_folder[:4] + '-' + ftp_folder[4:6] + '-' + ftp_folder[6:8]
    revision_info.append(f'{ftp_folder} - {date}')
    print(f'\033[0;32m ftp folder {ftp_folder}\033[0m')
    for contig in contigs:
        out_file = f'pre_svn_data/{contig}/{ftp_folder}.contig'
        if os.path.exists(out_file):
            print(f'{out_file} exists, skipping')
            continue
        resp = requests.get(f'https://www.pombase.org/data/genome_sequence_and_features/artemis_files/OLD/{ftp_folder}/{contig}.contig')
        if resp.status_code == 404:
            print(f'no {contig} in {ftp_folder}')
        else:
            print(f'writing {out_file}')
            with open(out_file, 'wb') as out:
                out.write(resp.content)

revision_info = ['svn_2 kmr44 2011-08-22'] + sorted(revision_info,reverse=True)
# Generate the revisions.txt files
for contig in contigs:

    # Only if the contig file exists
    filtered_revision_info = [r for r in revision_info if os.path.isfile(f'pre_svn_data/{contig}/{r.split(" ")[0]}.contig')]

    with open(f'pre_svn_data/{contig}/revisions.txt','w') as out:
        out.write('\n'.join(filtered_revision_info))

# Fix known errors

# File that has iD instead of ID on first line
with open('pre_svn_data/chromosome1/20070312.contig') as ins:
    lines = ins.readlines()
lines[0] = lines[0].replace('iD','ID')
with open('pre_svn_data/chromosome1/20070312.contig', 'w') as out:
    out.writelines(lines)

# Weird first line in mating_type_region files
for f in glob.glob('pre_svn_data/mating_type_region/*.contig'):
    with open(f) as ins:
        lines = ins.readlines()

    lines[0] = 'ID   mating_type_region standard; DNA; FUN;  20128BP.\n'

    with open(f, 'w') as out:
        out.writelines(lines)

# Wrong first line
f = 'pre_svn_data/chromosome3/20070620.contig'
with open(f) as ins:
        lines = ins.readlines()
lines[0] = 'ID   chromosome_3      standard; DNA; FUN; 2452883 BP.\n'
with open(f, 'w') as out:
    out.writelines(lines)