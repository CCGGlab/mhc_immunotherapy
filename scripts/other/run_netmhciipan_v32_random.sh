#!/bin/bash

# Split fasta file
mkdir "temp/random_nmp/pep15"
cd "temp/random_nmp/pep15"
split -l 20000 temp/random_15mers.fasta

# Run netMHC
tool='/home/labgroups/ccgg/tools/netMHCIIpan-3.2/netMHCIIpan'
ifile='temp/random_nmp/pep15'
ofolder='temp/random_nmp/netMHCIIpan-3.2_output/'

calc_affinity () {
    allele=$1
    chunk=$2
    output="$ofolder/${allele}_${chunk}.xls"
    if [ ! -e "${output}" ]; then
        "$tool" -a "${allele}" -f "${ifile}/${chunk}" -inptype 0 -length 15 -xls -xlsfile "${output}" > /dev/null
    fi
}

# Expected time: 6 hours per allele
. `which env_parallel.bash`
echo 50 > /tmp/mhcjobs
env_parallel --progress -j /tmp/mhcjobs calc_affinity :::: temp/melanoma_alleles_II.txt <(ls "$ifile")
