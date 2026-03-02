#!/bin/bash

folder="/Users/carloslaxmolina/Library/Mobile Documents/com~apple~CloudDocs/Lab - Postdoctorado/UniNE - Croll Lab/Holomycota TE Evolution/data/proteins_annot_faa"

cd "$folder"

for file in *.gz; do
    [ -e "$file" ] || continue     # si no hay .gz, no hace nada
    gunzip "$file"                # descomprime y elimina el .gz
done