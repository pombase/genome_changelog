# Fix known errors

# File that has iD instead of ID on first line
with open('data/chromosome1/3219.contig', errors='replace') as ins:
    lines = ins.readlines()
lines[0] = lines[0].replace('disorders;','')
with open('data/chromosome1/3219.contig', 'w') as out:
    out.writelines(lines)

for f in ['data/chromosome2/5338.contig', 'data/chromosome2/5334.contig']:
    with open(f, errors='replace') as ins:
        lines = ins.readlines()
    if 'ID' not in lines[0]:
        lines = ['ID   CU329671    standard; DNA; FUN; 4539804 BP.\n'] + lines

        with open(f, 'w') as out:
            out.writelines(lines)

# A bunch of revisions of chromosome 3 that do not have header
for rev in [28, 30, 31, 33, 34, 35, 39, 40]:
    f = f'data/chromosome3/{rev}.contig'

    with open(f, errors='replace') as ins:
        lines = ins.readlines()
    if 'ID' not in lines[0]:
        missing_header = '''ID   chromosome_3    standard; DNA; FUN; 2452883 BP.\nXX\nAC   chromosome_3;\nXX\n'''
        lines = [missing_header] + lines
        with open(f, 'w') as out:
            out.writelines(lines)

