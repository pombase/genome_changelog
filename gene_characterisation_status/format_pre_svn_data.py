import pandas
import os
import shutil

# This script requires having run get_ftp_site_files.py first

contigs = ['chromosome1','chromosome2','chromosome3','mating_type_region','pMIT', 'telomeric']

with open('../results/pre_svn_folder_list.tsv') as ins:
    ftp_folders = [l.strip() for l in ins.readlines()]

if not os.path.isdir('pre_svn_data'):
    os.mkdir('pre_svn_data')

for revision in ftp_folders:
    if not os.path.isdir(f'pre_svn_data/{revision}'):
        os.mkdir(f'pre_svn_data/{revision}')

    for contig in contigs:
        target_file = f'../pre_svn_data/{contig}/{revision}.contig'
        if os.path.exists(target_file):
            # Copy the file
            shutil.copy(target_file, f'pre_svn_data/{revision}/{contig}.contig')
