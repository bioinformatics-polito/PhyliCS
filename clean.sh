#!/bin/bash

ROOT="$(pwd)"

rm -rf $ROOT/bin/

rm -rf scripts cli genomes

cd $ROOT/phylics/local/src/ginkgo/scripts
rm -f testBED status interval binUnsorted CNVcaller

cd $ROOT/phylics/local/src/ginkgo/genomes/scripts
rm -f GC bin bounds findCentromeres fixed match_genes_to_bins simReads


