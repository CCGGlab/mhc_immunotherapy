#!/bin/bash

# Split fasta file
mkdir "temp/random_nmp/pep9"
cd "temp/random_nmp/pep9"
split -l 20000 temp/random_9mers.fasta

# Run netMHC
tool='/home/labgroups/ccgg/tools/netMHCpan-4.0/netMHCpan'
ifile='temp/random_nmp/pep9'
ofolder='temp/random_nmp/netMHCpan-4.0_output/'

calc_affinity () {
    allele=$2
    chunk=$1
    output="$ofolder/${allele}_${chunk}.xls"
    if [ ! -e "${output}" ]; then
        "$tool" -a "${allele}" -f "${ifile}/${chunk}" -inptype 0 -l 9 -BA -xls -xlsfile "${output}" > /dev/null
    fi
}

. `which env_parallel.bash`
echo 100 > /tmp/mhcjobs
# Update 26/05: order of parameters was switched here compared to
# original (and other scripts)
# Goal: first calculate each peptide for all alleles, then calculate next
# (We can already fit a distribution earlier and look when we reach convergence)
env_parallel --progress -j /tmp/mhcjobs calc_affinity :::: <(ls "$ifile") temp/melanoma_alleles_I.txt
