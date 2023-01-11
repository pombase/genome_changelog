from genome_functions import genome_dict_diff, build_seqfeature_dict
from Bio import SeqIO

with open('pre_svn_data/chromosome1/20090904.contig', errors='replace') as ins:
    new_genome_dict = build_seqfeature_dict(SeqIO.read(ins,'embl'), False)


print(len(new_genome_dict['SPAC8E11.03c']['CDS']))
print(len(new_genome_dict['SPAC8E11.02c']['CDS']))