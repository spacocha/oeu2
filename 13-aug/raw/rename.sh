#!/bin/bash

# might have to trim final tab in first line

out=otu_rename.txt

tr '\n' '\t' < headers.txt | sed 's/\t$//' > $out
echo >> $out
sed '1d' raw.txt >> $out
