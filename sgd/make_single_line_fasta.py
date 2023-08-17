import sys
from Bio.SeqIO import parse
import subprocess

def main(fasta_file):
    no_extension = fasta_file.split('.')[0]
    out_file = f'{no_extension}_single_line.tsv'
    with open(out_file, 'w') as out:
        for record in parse(fasta_file, 'fasta'):
            out.write(f'{record.id}\t{record.seq}\n')

    # sort -o filename
    subprocess.call(['sort', '-o', out_file, out_file])

if __name__ == '__main__':
    for arg in sys.argv[1:]:
        main(arg)