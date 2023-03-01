import os
import glob

revision_files = glob.glob('./data/*/revisions.txt')

for f in revision_files:
    output_dir = os.path.dirname(f)
    contig = output_dir.split('/')[-1]
    with open(f'{output_dir}/revisions.txt') as ins:
        first_line = ins.readline()
        if not first_line:
            continue
        contig = output_dir.split('/')[-1]
        old_revision = first_line.strip().split()[0]
        for line in ins:
            new_revision = old_revision
            old_revision = line.strip().split()[0]

            outfile = f'{output_dir}/diff/diff_{new_revision}_{old_revision}.txt'
            # if os.path.isfile(outfile):
            #     continue
            print(f'downloading diff {new_revision}_{old_revision} at {output_dir}')
            os.system(f'svn diff -r {new_revision}:{old_revision} https://curation.pombase.org/pombe-embl-repo/trunk/{contig}.contig > {outfile}')

