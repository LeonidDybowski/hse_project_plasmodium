#!/bin/bash

while IFS= read -r -u 4 line1 && IFS= read -r -u 5 line2; do
	esearch -db nucleotide -query "$line1" | efetch -format fasta > $line2.fasta
	esearch -db nucleotide -query "$line1" | efetch -format gb > $line2.gb
done 4< files_for_project.txt 5< names_for_files_for_project.txt
