# Errors in files

## Missing space in first line

Some of the files are missing a space between the number and BP in the first line:

```
# first line is this
ID   chromosome_1    standard; DNA; FUN;  5579133BP.
# should be this
ID   chromosome_1    standard; DNA; FUN;  5579133 BP.
```

To support this, we override the function `EmblScanner._feed_seq_length` from `BioPython` in `pombe_svn_full_diff.py`.

## Other errors

Some of the files need fixing. This was done manually:
* Remove `disorders;` from first line in chromosome1 in revision 3219.
* chomosome2: 5338, 5334, 26 missing header.
* chromosome3: revisions 40-28 are missing the header.

## Incomplete files

Some files have a premature EOF, likely because Artemis was closed before finishing saving them. They are excluded from the analysis, listed in the dictionary `skip_files` in `pombe_svn_full_diff.py`.
