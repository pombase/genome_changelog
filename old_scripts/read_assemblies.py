import Bio.SeqIO as SeqIO

data1 = SeqIO.parse('GCA_000002945.1.txt', 'embl')
data2 = SeqIO.parse('GCA_000002945.2.txt', 'embl')

while True:
    try:
        record1 = next(data1)
        record2 = next(data2)
        print(record1.seq == record2.seq)
    except StopIteration:
        break

data2 = SeqIO.parse('GCA_000002945.2.txt', 'embl')

current_contigs = ['chromosome1.contig', 'chromosome2.contig', 'chromosome3.contig', 'pMIT.contig', 'mating_type_region.contig', 'telomeric.contig']
for record1 in data2:
    contig = current_contigs.pop(0)
    record2 = SeqIO.read('../latest_genome/' + contig, 'embl')
    print(contig)
    print(record1.id, record2.id)
    print(record1.seq == record2.seq)
