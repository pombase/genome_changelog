import os

with open('revision_numbers.txt') as ins:
    i = 0
    for line in ins:
        revision = line.strip()
        i += 1
        if i > 100:
            break

        outfile = f'revision_files/chromosome1_{revision}.txt'
        if os.path.isfile(outfile):
            continue
        print(f'downloading {revision}')
        os.system(f'svn cat -r {revision} svn+ssh://manu@curation.pombase.org/var/svn-repos/pombe-embl/trunk/chromosome1.contig > {outfile}')



