#!/bin/bash

# This script invokes immcantation scripts that together download IMGT
# sequences and makes igblast databases for BCRs and TCRs for both
# human and mouse

# exit upon any error
set -e

usage () {
    echo -e "Usage: `basename $0` [OPTIONS]"
    echo -e "  -m  makeblastdb directory"
    echo -e "  -c  Conda python3 environment name (or have python3 in path)"
}

MSET=false
CSET=false

# parse args
while getopts "i:m:c:h" option; do
    case "$option" in
    m)  MAKEBLASTDB=${OPTARG}
        MSET=true
        ;;
    c)  CONDA3=${OPTARG}
        CSET=true
        ;;
    h)  usage
        exit
        ;;
    \?) echo "Invalid option: -$OPTARG" >&2
        exit 1
        ;;
    :)  echo "Option -$OPTARG requires an argument" >&2
        exit 1
        ;;
    esac
done

# check for missing args
echo ""
if ! $MSET; then
    echo "Error: missing makeblastdb directory path" >&2
    exit 1
fi

if $CSET; then
    source activate $CONDA3
fi

type python3 >/dev/null 2>&1 || {
    echo >&2 "Error: python3 not found. Aborting.";
    exit 1;
}

#
# Create IMGT databases (format: <species>/<receptor>/*.fasta)
export PATH=$PATH:$PWD/immcantation_scripts:$MAKEBLASTDB  # necessary for internal calls
#bash fetch_imgtdb.sh -o germline_fastas

# check for empty files
if find germline_fastas/ -type f -empty | grep '.'; then
    echo "The above files are empty!"
    echo "Terminating."
    echo ""
    exit 1
fi

# Combine fastas and run makeblastdb
bash imgt2igblast.sh -i $PWD/germline_fastas

echo ""
if grep -q "J-C" database/*_c.fasta; then
    echo ""
    echo "Fasta lines to delete:"
    grep "J-C" database/*_c.fasta
    echo "Terminating."
    echo ""
    exit 1
fi

echo ""
echo "Making BASIC bowtie2 indices:"
#python basic_from_immcantation.py

echo ""
echo "Success. IgBLAST databases created. Exiting."
