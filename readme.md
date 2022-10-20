# Genome Changelog

This repository contains scripts for:

* Downloading previous versions of `.contig` files in the svn repository `curation.pombase.org/var/svn-repos/pombe-embl`.
* Summarise the differences between subsequent versions of the contig files, namely:
  * Coordinates of removed/added features.
  * Changes in coordinates of features that are present in both versions.
  * Changes in qualifiers of features that are present in both versions.
* It also contains a script to perform the diff for any two embl files (`two_genomes_diff.py`)

## TL;DR; to update diff files

```bash
# install dependencies
poetry install

# activate venv
poetry shell

# run this script
bash update_file.sh

# Commit changes
```

## Installing dependencies

To install the dependencies, we used poetry (see [poetry installation instructions](https://python-poetry.org/docs/)).

In the source directory run:

```
poetry install
```

This should create a folder `.venv` with the python virtual environment. To activate the virtual environment, then run:

```
poetry shell
```

Now when you call `python`, it will be the one from the `.venv`.

## Calculating differences between two genomes

Most of the repository was made for a pombe genome analysis, in which many diffs of pombe genome were compared. If you are here only to quickly compare two genomes, you can run the script:

```
python two_genomes_diff.py --new_genome data/chromosome1/8485.contig --old_genome data/chromosome1/8338.contig --output_locations_file 'a.tsv' --output_qualifiers_file 'b.tsv'
```

Arguments:

* `--new_genome` and `--old_genome`: files to be compared, in embl format.
* `--revision_string`: a string with 3 space-separated revision-related values (revision number, user, date). If not provided, it is not printed.
* `--output_locations_file`: the file where the diff in locations will be stored.
* `--output_qualifiers_file`: the file where the diff in qualifiers will be stored.

From now on the instructions, are to perform the analysis for pombase svn revisions.

## Getting the data

> **WARNING:** Downloading all revisions and generating the full diffs will require ~100GB of space.

### Access to the svn repository

You will need to make many calls to the svn server to download the files, you should set up ssh access:

```bash
# In the local machine - connect to the server:
ssh 'username@curation.pombase.org'

# In the remote machine
# Generate a key for your user
ssh-keygen

# exit the server
exit

# In the local machine - create an ssh key so you are not prompted for password each time
ssh-copy-id username@curation.pombase.org
```

### Getting the revisions that you want to analyse

```bash
# If you haven't, activate the local python environment
poetry shell

# Create the basic folder structure, and download the information about revisions

# If you want to download last revisions (since last revision mentioned in either all_qualifier_changes_file or all_coordinate_changes_file.tsv)
bash get_revisions_where_contigs_changed.sh last

# If you want to download ALL revisions (~100 GB)
bash get_revisions_where_contigs_changed.sh all
```

This will create a directory structure:

```
data
├── chromosome1
│   ├── change_log
│   │   ├── locations
│   │   ├── qualifiers
│   │   └── diff
│   └── revisions.txt
├── chromosome2
│   ├── change_log
│   │   ├── locations
│   │   ├── qualifiers
│   │   └── diff
│   └── revisions.txt
├── chromosome3
│   ├── change_log
│   │   ├── locations
│   │   ├── qualifiers
│   │   └── diff
│   └── revisions.txt
├── mating_type_region
│   ├── change_log
│   │   ├── locations
│   │   ├── qualifiers
│   │   └── diff
│   └── revisions.txt
└── pMIT
    ├── change_log
    │   ├── locations 
    │   ├── qualifiers
    │   └── diff
    └── revisions.txt
```

In each `revisions.txt` file there is information about the revisions that affect a given chromosome space separated (revision user date), e.g.:
```
8485 vw253 2022-10-08
```

After this, you can download all versions each contig file where changes were made by running:

```
python get_revisions_files.py
```

The generated folder tree looks like this:

```
data
├── chromosome1
│   ├── 10.contig
│   ├── 1000.contig
│   ├── 1002.contig
│   ├── 1007.contig
...
```
Where `data/chromosome1/10.contig` is the file `chromosome1.contig` at revision 10, etc.

You can also download the svn diffs by running `python get_svn_diff.py`

### Known errors

There are known erros in the downloaded contig files, some of which will require manual fixing if you are to run the analysis pipeline. See [`known_errors.md`](known_errors.md).

## Running the analysis

For running the analysis:

```bash
# If you haven't activated the environment
poetry shell

python pombe_svn_diff.py
```

This will generate the following output:

```
data
├── chromosome1
│   ├── change_log
│   │   ├── locations
│   │   │   ├── 10.tsv
│   │   │   ├── 1000.tsv
..........................
│   │   └── qualifiers
│   │       ├── 10.tsv
│   │       ├── 1000.tsv
..........................
```

Where `data/chromosome1/change_log/locations/xxx.tsv` contains changes in location that were introduced in revision `xxx`. The file might be empty if no changes where made in that revision. Otherwise it contains:

* Coordinates of removed/added features.
* Changes in coordinates of features that are present in both versions.

The output looks like this:

```tsv
revision	user	date	systematic_id	primary_name	feature_type	added_or_removed	value
8462	vw253	2022-09-27	SPNCRNA.145		ncRNA	removed	239730..240571
8462	vw253	2022-09-27	SPNCRNA.18		ncRNA	removed	complement(3699381..3700010)
8462	vw253	2022-09-27	SPNCRNA.193		ncRNA	removed	2404684..2405924
8462	vw253	2022-09-27	SPNCRNA.884		ncRNA	removed	complement(3084619..3086087)
8462	vw253	2022-09-27	SPNCRNA.951		ncRNA	removed	3951825..3952588
```

To combine all the changes in a single file, you can then run:

```
python create_single_coordinate_changes_file.py>yourfile.tsv
```

`data/chromosome1/change_log/qualifiers/xxx.tsv` contains changes introduced in revision `xxx` to qualifiers of features that existed in revision `xxx` and the previous one. The output looks like this:

```tsv
revision	user	date	systematic_id	primary_name	feature_type	qualifier_type	added_or_removed	value
7940	vw253	2022-01-04	SPAC20G4.09		intron	controlled_curation	added	('term=misc, confirmed intron',)
7940	vw253	2022-01-04	SPAC20G4.09		intron	controlled_curation	removed	('misc, confirmed',)
7940	vw253	2022-01-04	SPAC23D3.16		intron	controlled_curation	added	('term=misc, confirmed intron',)
7940	vw253	2022-01-04	SPAC23D3.16		intron	controlled_curation	removed	('term=misc, confirmed',)
7940	vw253	2022-01-04	SPAC25G10.03		intron	controlled_curation	added	('term=misc, confirmed intron',)
7940	vw253	2022-01-04	SPAC25G10.03		intron	controlled_curation	removed	('term=misc, confirmed',)
```

To combine all the changes in a single file, you can then run:

```
python create_single_qualifier_changes_file.py>yourfile.tsv
```