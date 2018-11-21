#!/bin/bash

module load plink2

FILE=$1
BASE=`basename $FILE`

LD=$2
P=$3
SUMMARY_STATS=$4
SUMMARY_BASE=`basename $SUMMARY_STATS | perl -pi -e 's/\.gz//g'`

cd /home/armartin/pgc_scz_asia/imp_collect/prs_score/armartin

#1. Clump SNPs, considering LD
plink --bfile $LD \
--clump $SUMMARY_STATS \
--clump-field P \
--clump-snp-field SNP \
--clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-kb 500 \
--out $SUMMARY_BASE
