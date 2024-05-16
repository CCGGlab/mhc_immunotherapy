#!/bin/bash
mkdir temp/RSEM_output/

# Generate STAR reference index
bash 'scripts/other/generate_rsem_reference.sh'

# Activate Conda environment where RSEM is installed
conda activate temp/RSEM/RSEM/
# channels:
#   - conda-forge
#   - bioconda
# dependencies:
#   - rsem=1.3.3
#   - star=2.7.10b
#   - python=3.12

# Analyze 4 samples in parallel
echo 4 > /tmp/jobs_rsem

# Liu
input_folder="downloads/dbGaP/phs000452/Liu/RNA"
output_folder="temp/RSEM_output/Liu_2019"

find "$input_folder" -type d -name 'SRR*' -printf '%f\n'|
parallel \
-j /tmp/jobs_rsem \
-q \
bash 'scripts/other/align_and_count_rsem.sh' '{}' "${input_folder}" "${output_folder}"

# Gide
input_folder="/home/labgroups/ccgg/downloads/ENA/ERP105482"
output_folder="temp/RSEM_output/Gide_2019"

find "$input_folder" -type d -name 'ERR*' -printf '%f\n'|
parallel \
-j /tmp/jobs_rsem \
-q \
bash 'scripts/other/align_and_count_rsem.sh' '{}' "${input_folder}" "${output_folder}"

# Hugo
input_folder="downloads/SRA/SRP070710"
output_folder="temp/RSEM_output/Hugo_2016"

find "$input_folder" -type d -name 'SRR*' -printf '%f\n'|
parallel \
-j /tmp/jobs_rsem \
-q \
bash 'scripts/other/align_and_count_rsem.sh' '{}' "${input_folder}" "${output_folder}"

# Riaz
input_folder="downloads/pub/riaz_2017/full/RNA/"
output_folder="temp/RSEM_output/Riaz_2017"

find "$input_folder" -type d -name 'SRR*' -printf '%f\n'|
parallel \
-j /tmp/jobs_rsem \
-q \
bash 'scripts/other/align_and_count_rsem.sh' '{}' "${input_folder}" "${output_folder}"
