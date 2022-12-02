import os
import glob
import argparse
class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
parser.add_argument('--revision_files', nargs='+', default=glob.glob('./data/*/revisions.txt'), help='revision text file')
args = parser.parse_args()

for f in args.revision_files:
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
