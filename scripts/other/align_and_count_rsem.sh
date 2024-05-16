#!/bin/bash

set -exo pipefail

sample_id="$1"
input_folder="$2"
output_folder="$3"
REFERENCE="temp/RSEM/grch38_ensembl"

quantify_samples() {
  local sample_id="$1"
  local F1="${sample_id}_1.fastq.gz"
  local F2="${sample_id}_2.fastq.gz"

  rsem-calculate-expression --star --paired-end --star-gzipped-read-file -p 8 \
  "$F1" "$F2" \
  "$REFERENCE" \
  "${sample_id}" >& "${sample_id}.rsem.log"
}

ls "${output_folder}" | grep -q "${sample_id}" &&
{ echo File exists; exit 0; }

tempfolder="$(mktemp -d -p "temp/temp_nosnapshot/arne/")"
cleanup() {
    if [ -n "$tempfolder" ]; then
        echo cleanup
        rm -r "$tempfolder"
    fi
}
trap cleanup EXIT

# Perform all analyses inside the temp folder for this sample
cd "$tempfolder"

# Download samples
cp -s "${input_folder}/${sample_id}/"* "./"

# Perform alignment
echo 'Processing' "${sample_id}";
quantify_samples "$sample_id"

# Copy results to remote server
cp "${sample_id}.genes.results" "${sample_id}.isoforms.results" "${output_folder}/"
