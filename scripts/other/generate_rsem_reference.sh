#!/bin/bash
# According to instructions at https://deweylab.github.io/RSEM/README.html
# - Build RSEM references using RefSeq, Ensembl, or GENCODE annotations
# - Follow instructions for Ensembl

# Get absolute path of folder where Ensembl genome and matching GTF annotation is stored
OUT_FOLDER="$(realpath temp/RSEM/grch38_ensembl)"
GENOMES="$(realpath downloads/genomes/ensembl)"
cd "$GENOMES"

# Download primary assembly
wget http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Download annotation
wget https://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz
gunzip Homo_sapiens.GRCh38.103.gtf.gz

cd "$OUT_FOLDER"
rsem-prepare-reference --gtf $GENOMES/Homo_sapiens.GRCh38.v103.gtf \
--star -p 8 \
$GENOMES/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
grch38_ensembl
