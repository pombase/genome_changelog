from Bio.SeqIO import read

pombase_file = read('pre_svn_data/chromosome1/svn_2.contig', 'embl')
genbank_file = read('pre_svn_data/chromosome1/20110603.contig', 'embl')

print(pombase_file.seq == genbank_file.seq)


