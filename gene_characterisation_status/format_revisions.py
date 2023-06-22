import re

with open('revisions.tsv', 'w') as out:
    with open('svn_revisions.txt', 'r') as ins:
        line = ins.readline()
        while line:
            line = line.strip()
            if not line.startswith('----'):
                line = ins.readline()
                continue
            info_line = ins.readline().strip()
            if not info_line.startswith('r'):
                break
            revision, date = re.search(r'r(\d+) \| [^\s]+ \| ([^\s]+)', info_line).groups()
            out.write(f'{revision}\t{date}\n')
            line = ins.readline()

