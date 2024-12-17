#! /bin/bash

mkdir data/repliseq_274ko/bedgraph/
for bigwigpath in data/repliseq_274ko/bigwig/*; do
    bigwig=$(basename ${bigwigpath})
    bin_external/bigWigToBedGraph $bigwigpath data/repliseq_274ko/bedgraph/$bigwig.bedgraph
done