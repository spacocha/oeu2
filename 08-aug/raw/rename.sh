#!/bin/bash

out=otu_rename.txt

tr '\n' '\t' < headers.txt | sed 's/\t$//' > $out
echo >> $out
sed '1d' raw.txt >> $out
