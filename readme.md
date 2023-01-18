# Genome Changelog

This repository contains scripts for:

* Downloading previous versions of `.contig` files in the svn repository `curation.pombase.org/var/svn-repos/pombe-embl`.
* Summarise the differences between subsequent versions of the contig files, namely:
  * Coordinates of removed/added features.
  * Changes in coordinates of features that are present in both versions.
  * Changes in qualifiers of features that are present in both versions.
* It also contains a script to perform the diff for any two embl files (`two_genomes_diff.py`)
* More PomBase-specific stuff
  * Download versions that pre-date the use of svn from the PomBase FTP server.

## TL;DR; to update diff files ⏩

All the steps below that need to be re-ran are in [update_file.sh](update_file.sh).

```bash
# install dependencies
poetry install

# activate venv
poetry shell

# run this script (The comments explain what it does)
bash update_file.sh

# Commit changes
```

A github action is there, but it won't work since the svn repo is not public.

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

## Using the code to calculate differences between two genomes (non-PomBase use-case)

Most of the repository was made for a pombe genome analysis, in which many diffs of pombe genome were compared. If you are here only to quickly compare two genomes, you can run the script:

```
python two_genomes_diff.py --new_genome data/chromosome1/8485.contig --old_genome data/chromosome1/8338.contig --output_locations_file 'a.tsv' --output_qualifiers_file 'b.tsv'
```

Arguments:

* `--new_genome` and `--old_genome`: files to be compared, in embl format.
* `--revision_string`: a string with 3 space-separated revision-related values (revision number, user, date). If not provided, it is not printed.
* `--output_locations_file`: the file where the diff in locations will be stored.
* `--output_qualifiers_file`: the file where the diff in qualifiers will be stored.

Before using this, read the next section.

### Feature unique identifiers ⚠️

In order to assess whether a feature has changed or not, the script [two_genomes_diff.py](two_genomes_diff.py) needs to use some unique identifier of a feature to compare it between
two genome versions. For PomBase, that is the feature qualifier `\systematic_id`. Old genomes did not have this qualifier, and used `\gene` instead. However, the `\gene` qualifier
does not have to be unique, and there are some revisions were the usage of both is mixed.
This is handled for the PomBase case in the function `read_pombe_genome` in [genome_functions.py](genome_functions.py).
The data used for this was generated with the pre-svn pombe genomes and with the script [find_missing_synonyms.sh](find_missing_synonyms.sh) and is in the directory
[valid_ids_data](valid_ids_data).

The same unique identifier can refer to multiple features sometimes, such as introns. In that case, we report the ones that are removed or added.

### Resolving conflicts for a feature with multiple unique identifiers

Sometimes, there have been features that were associated with a single 

## Getting the data (PomBase)

> **WARNING:** Downloading all revisions and generating the full diffs will require ~100GB of space.

No matter what you do, you need to download the synonyms list:
```
bash get_valid_ids.sh
```

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

After this, you can download all versions each contig file where changes were made by running (note that you can run each chromosome in parallel to speed up):

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
python create_single_coordinate_changes_file.py --output_file the_file.tsv
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
python create_single_qualifier_changes_file.py --output_file the_file.tsv
```

### Listing changes on main features

> **NOTE** You need to run this first:
> ```
> bash get_data_gene_changes_comments_and_pmids.sh
> ```
>

`get_info_from_changes.py` generates [only_modified_coordinates.tsv](only_modified_coordinates.tsv), a table listing only the changes in gene locations (addition and removal in the same revision) for the main types of features (`CDS`,`ncRNA`,`snRNA`,`repeat_region`,`rRNA`,`tRNA`,`snoRNA`,`misc_RNA`).

```
python get_info_from_changes.py --input_files all_coordinate_changes_file.tsv pre_svn_coordinate_changes_file.tsv --output_modified_coordinates only_modified_coordinates.tsv
```

This can be particularly useful for alleles that refer to previous gene coordinates. This is used in the repository https://github.com/pombase/allele_qc.

 contains a table where only coordinate modifications are shown (modifications as in changes to gene coordinates, rather than new additions or removals of features). It combines the data from the svn server and the pre-svn data.

### Associating publications with changes in gene features.

Pombase-specific. Links changes in genome coordinates in the file `only_modified_coordinates.tsv` to either:

* Comments from PomBase website. The file from PomBase [curation repository]() lists some (not all) of the gene feature changes, and the reasons that led to them (a publication, personal communication, etc.).
* Changes in `dbxref` that occurred in the same revision as a change recorded in `only_modified_coordinates.tsv`.

To run this:
```bash
# Download qualifier changes from latest release
bash get_data_gene_changes_comments_and_pmids.sh

# Make the associations
python associate_comments_with_genome_changes.py --output only_modified_coordinates_with_comments.tsv
```

### Listing revisions where genome sequence changed

Output in the style of [genome_changes_svn](genome_changes_svn.tsv).

```bash
python get_revisions_where_genome_sequence_changes.py --data data/* --output_file output.tsv
```

### Delete analysis data

```
rm pre_svn_data/*/change_log/*/*.tsv
rm data/*/change_log/*/*.tsv
```

## Pre-svn data

Some of the contig files pre-date the use of SVN, to download them and calculate the differences, they are in the ftp server of PomBase: https://www.pombase.org/data/genome_sequence_and_features/artemis_files/OLD/. The full list of those that pre-date svn are in the file ![pre_svn_folder_list.tsv]([pre_svn_folder_list.tsv]). The output files are attached in the release.

```bash
# Download files from first revision of svn and prepare directory structure (pre_svn_data)
bash prepare_pre_svn_folder.sh

# Download the contig files from ftp site and produce an equivalent to the revisions.txt described above
python get_ftp_site_files.py

# Run the diffs on the pre_svn_data directory
python pombe_svn_diff.py --data_folders pre_svn_data/*

# Combine in single files
python create_single_coordinate_changes_file.py --data_folder pre_svn_data/  --output_file pre_svn_coordinate_changes_file.tsv
python create_single_qualifier_changes_file.py --data_folder pre_svn_data/ --output_file pre_svn_qualifier_changes_file.tsv
```

### Find missing synonyms

Some of the code was used to finding missing synonyms of genes and gather them into "tsv dictionaries" (see files `valid_ids_data/missing_synonyms.tsv` and `valid_ids_data/obsoleted_ids.tsv`). The code is in `find_missing_synonyms.sh` and the python scripts called within. To run it you need to have the pre_svn data as well as the latest version of the pombe genome.

