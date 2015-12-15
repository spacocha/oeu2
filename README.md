% swo

This is a reword of the old OEU-calling pipeline, now purely in R.

# Requirements
You'll need the R packages:
- ggplot2
- dplyr
- optparse
- reshape2

You can install a package by running `install.packages("foo")` in R.


# Historical note
I re-implemented the original pipeline in `oeu-old.R` as a way to double-check
that the implementations were equivalent. I discovered some bugs that way.

That file is a combination of preprocessing and OEU calling.


# Main executables
The two R executables `oeu.R` and `make-plots.R` can call OTUs and make
plots of OTUs in each OEU.

Both scripts use optparse, so you can call, say, `oeu.R --help` to get
reminders about the order of arguments.

## `oeu.R`
This script takes an OTU table that:
- is normalized (i.e., it shows relative abundance in sample), 
- is preprocessed (i.e., replicated have been pooled, bad OTUs/samples have been removed),
- has OTU names in the first column (preferably under "OTU\_ID"), and
- is tab-separated.

This script produces a "groups file", which is a tab-separated file with
two columns: OTU ID and group number.

## `make-plots.R`
This script takes
- an OTU table, the same kind that fed to `oeu.R`, and
- a groups file, the same kind produced by `oeu.R`.

The script will produce a one plot per OEU in a directory you choose.


# 08-aug and 13-aug
These two directories have data from the two timepoints. They are each
preprocessed using some custom scripts (since the header names and number
of replicates is different for both).


# Analysis: `changen`
There's an older script here that I had used to test what happens when
OEUs are called with different numbers of candidate clusters `k`. I plan
to update this into something more script-y that calls `oeu.R`.
