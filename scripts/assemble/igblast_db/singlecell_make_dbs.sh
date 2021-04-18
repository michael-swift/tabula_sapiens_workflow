#!/bin/bash

./make_dbs.sh -m ~/resources/ncbi-blast-2.7.1+/bin/ -c py3.5 2>&1 | tee out.log
