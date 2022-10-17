#%%

from calculate_difference import build_seqfeature_dict
from Bio import SeqIO



d = build_seqfeature_dict(SeqIO.read('revision_files/chromosome1_7936.txt','embl'))

print(d['SPAC27F1.09c']['intron'])
