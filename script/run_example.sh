#!/bin/bash

ODIR=example/cis_eqtl_mtx
mkdir -p $ODIR

# bin/fastQTL --vcf example/genotypes.vcf.gz \
#             --bed example/phenotypes.bed.gz \
#             --region 22:17000000-18000000 \
#             --out example/nominals.default.txt.gz \
#             --mtx $ODIR

echo ""
echo ""
echo "CHECK OUTPUT MATRICES:"
echo ""
for I in 0 1 2 3 4
do
    FNAME=example/cis_eqtl_mtx/$I.txt
    wc -l $FNAME 
    head -n3 $FNAME | cut -f1-3
    awk '{print NF}' $FNAME | sort -nu | tail -n 1
done
