#!/bin/bash

# might have to trim final tab in first line

out=raw_rename.txt

tr '\n' '\t' < headers.txt > $out
echo >> $out
sed '1d' raw.txt >> $out
