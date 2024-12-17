#! /bin/bash

# entire signal, on all cell fractions.
wget -P data https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562WaveSignalRep1.bigWig

mkdir data/repliseq_k562
mkdir data/repliseq_k562/bigwig

# a single bigwig per cell fraction, will have to see what the percentage corresponds to.
# G1
wget -P data/repliseq_k562/bigwig https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562G1PctSignalRep1.bigWig
# S1
wget -P data/repliseq_k562/bigwig https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562S1PctSignalRep1.bigWig
# S2
wget -P data/repliseq_k562/bigwig https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562S2PctSignalRep1.bigWig
# S3
wget -P data/repliseq_k562/bigwig https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562S3PctSignalRep1.bigWig
# S4
wget -P data/repliseq_k562/bigwig https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562S4PctSignalRep1.bigWig
# G2
wget -P data/repliseq_k562/bigwig https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562G2PctSignalRep1.bigWig


# transforming the outputs from bigwig to bedgraph, which can be human-read
mkdir data/repliseq_k562/bedgraph/
for bigwigpath in data/repliseq_k562/bigwig/*; do
    bigwig=$(basename ${bigwigpath})
    bin_external/bigWigToBedGraph $bigwigpath data/repliseq_k562/bedgraph/$bigwig.bedgraph
done