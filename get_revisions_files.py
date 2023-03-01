import os
import glob
import argparse
from multiprocessing import Pool
import subprocess

class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


parser = argparse.ArgumentParser(description=__doc__, formatter_class=Formatter)
parser.add_argument('--nb_processes', default=1, help='number of processes to download in parallel')
parser.add_argument('--revision_files', nargs='+', default=glob.glob('./data/*/revisions.txt'), help='revision text file')
args = parser.parse_args()


def download(arg_tuple):
    revision, output_dir, contig = arg_tuple
    outfile = f'{output_dir}/{revision}.contig'
    if os.path.isfile(outfile):
        return
    print(f'downloading {revision} at {output_dir}')
    proc = subprocess.Popen([f'svn export --force -r {revision} https://curation.pombase.org/pombe-embl-repo/trunk/{contig}.contig {outfile}'], shell=True)
    proc.wait()
    proc.terminate()

if __name__ == '__main__':

    args4download = list()

    for f in args.revision_files:
        output_dir = os.path.dirname(f)

        with open(f'{output_dir}/revisions.txt') as ins:
            for line in ins:
                revision = line.strip().split()[0]
                contig = output_dir.split('/')[-1]
                args4download.append((revision,output_dir, contig))

        if int(args.nb_processes) > 1:
            with Pool(processes=5) as pool:
                pool.map(download,args4download)
        else:
            for arg_tuple in args4download:
                download(arg_tuple)