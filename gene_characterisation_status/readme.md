# Gene characterisation status

See https://github.com/pombase/website/issues/2048 for details on the meaning of the colours of CDS features.

For pre-svn revisions (from the ftp site), all revisions were used. For svn revisions, we take the revision that is closest to each first of the month, to have one revision chosen per month (see `results/formatted_counts.tsv`). This is done in the file `get_monthly_revisions.py`, using `pandas` `merge_asof`.

To update the files and graph, run `update_characterisation.sh`, which will only update files if a new first of the month has been reached since the last date recorded in `results/formatted_counts.tsv`.

For a given revision in svn, the script creates a folder `data/revision-number` that contains the contig files
in that revision. Then it runs the script `../colours_from_genome.py` that generates the file `data/revision-number/counts.txt`
with a single row that looks like this:

```
0 0 2502 0 13 0 64 1912 109 0 375 0 154 31 0 0
```

Each value corresponds to the number of occurrences of colour == 0 to colour == 15 in CDS features
that have a systematic_id (some uncertain CDS in some revisions do not have systematic ids).

For this particular example, there was 2502 qualifiers of colour=2, 13 of colour=4, etc.

For all revisions, these are gathered in the file `results/merged_counts.tsv` by the script `merge_counts_new.py`.

Then the file `format_counts.py` interprets the colour values and writes the file `results/formatted_counts.tsv`. It also makes the plot `results/figure.svg`.

