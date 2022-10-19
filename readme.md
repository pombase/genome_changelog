# Genome Changelog

This repository contains scripts for:

* Downloading all previous versions of `.contig` files in the svn repository `curation.pombase.org/var/svn-repos/pombe-embl`.
* Summarise the differences between subsequent versions of the contig files, namely:
  * Coordinates of removed/added features.
  * Changes in coordinates of features that are present in both versions.
  * Changes in qualifiers of features that are present in both versions.


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

### Getting the different revisions

```bash
# If you haven't, activate the local python environment
poetry shell

# Create the basic folder structure, and download the information about revisions
bash coordinate_changelog.sh
```

This will create a directory structure