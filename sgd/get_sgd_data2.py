import requests
import re
import os
import subprocess
import glob
from Bio.SeqIO import read

webpage = requests.get('http://sgd-archive.yeastgenome.org.s3-us-west-2.amazonaws.com/?delimiter=/&prefix=sequence/S288C_reference/orf_protein/archive/').text
release_info = re.findall(r'<Key>(sequence/S288C_reference/orf_protein/archive/orf_trans_all.*?(\d{6}).fasta.gz)</Key>', webpage)

# Download files
for path, date in release_info:
    file_name = 'data/{}.fasta.gz'.format(date)
    # If the file does not exist, download it:
    url = 'http://sgd-archive.yeastgenome.org/{}'.format(path)
    if not os.path.exists(file_name):
        print('Downloading {}'.format(file_name))
        with open(file_name, 'wb') as out:
            out.write(requests.get(url).content)

# # Download current release
with open('current.fasta.gz', 'wb') as out:
    out.write(requests.get('http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans_all.fasta.gz').content)

subprocess.call(['gzip', '-fd', 'current.fasta.gz'])

for uncompressed_file in glob.glob('data/*fasta.gz'):
    print('Uncompressing {}'.format(uncompressed_file))
    subprocess.call(['gzip', '-fd', uncompressed_file])

