#!/usr/bin/env bash

/local10G/dcroote/resources/ncbi-igblast-1.8.0/bin/makeblastdb -dbtype nucl -parse_seqids -in $1
