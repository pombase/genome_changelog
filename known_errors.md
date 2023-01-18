# Errors in files

## Skipping versions

Sometimes some version has formatting errors that make it impossible to read the contig files. For example, some old svn files have a premature EOF, likely because Artemis was closed before finishing saving them. You can simply exclude those files from the analysis. For that, add them to the dictionary in `skip_files` in the script `pombe_svn_diff.py`.

## Some fixed errors

### Missing space in first line

Some of the files are missing a space between the number and BP in the first line:

```
# first line is this
ID   chromosome_1    standard; DNA; FUN;  5579133BP.
# should be this
ID   chromosome_1    standard; DNA; FUN;  5579133 BP.
```

To support this, we override the function `EmblScanner._feed_seq_length` from `BioPython` in `pombe_svn_full_diff.py`.

### Other errors

Some of the old files from svn need fixing. This can be done by running the script:

```
python fix_known_errors_svn_files.py
```

Some of the pre-svn files are directly fixed in `get_ftp_site_files.py`
