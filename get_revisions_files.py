import os
import glob

revision_files = glob.glob('./data/*/revisions.txt')

for f in revision_files:
    output_dir = os.path.dirname(f)

    with open(f'{output_dir}/revisions.txt') as ins:
        for line in ins:
            revision = line.strip().split()[0]
            contig = output_dir.split('/')[-1]
            outfile = f'{output_dir}/{revision}.contig'
            if os.path.isfile(outfile):
                continue
            print(f'downloading {revision} at {output_dir}')
            os.system(f'svn cat -r {revision} svn+ssh://manu@curation.pombase.org/var/svn-repos/pombe-embl/trunk/{contig}.contig > {outfile}')
