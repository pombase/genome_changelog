import os

with open('data/chromosome1/revisions.txt') as ins:
    newer_revision = ins.readline().strip().split()[0]
    older_revision = ins.readline().strip().split()[0]
    i = 0
    while older_revision:
        i += 1
        if i > 100:
            break
        outfile = f'diffs/diff_{newer_revision}_{older_revision}.txt'
        if os.path.isfile(outfile):
            newer_revision = older_revision
            older_revision = ins.readline().strip().split()[0]
            continue
        print(f'downloading diff between {newer_revision} & {older_revision}')
        os.system(f'svn diff -r {newer_revision}:{older_revision} svn+ssh://manu@curation.pombase.org/var/svn-repos/pombe-embl/trunk/chromosome1.contig > {outfile}')
        newer_revision = older_revision
        older_revision = ins.readline().strip().split()[0]



