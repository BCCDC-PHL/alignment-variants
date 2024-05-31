#!/bin/bash

mkdir -p .github/data/refs

curl -o .github/data/refs/NC_000962.3.fa "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=NC_000962.3&db=nucleotide&rettype=fasta"
curl -o .github/data/refs/NC_002973.6.fa "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=NC_002973.6&db=nucleotide&rettype=fasta"
