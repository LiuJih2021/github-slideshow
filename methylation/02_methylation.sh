#!/usr/bin/env bash
#bedtools v2.30.0
# Time-stamp: <2021年 07月 01日 星期四 11:38:48 CST liujihong>

usage() {
    echo "usage: bash 02_methylation.gene.sh <METHYLATIONGFF> <SAMPLE> "
    echo "where: <METHYLATIONGFF> is a specific DNA methylation file"
    echo "       <SAMPLE> is the strain&DNA methylation name"
    echo "This script does the job of fetching the list of variables to pass to"
    echo "the script a1_BTPMI.sh"
}

# Minimal argument checking

if [ $# -lt 1 ]; then
    usage
    exit
fi
# Set variable for input file
SAMPLE=$2
echo "Start"
date
#get methylatioin site bed
cp $1 tmp.gff
cp $SAMPLE/*.gff genome.gff
cp $SAMPLE/*.fna genome.fna
cp $SAMPLE/*.tsv genome.tsv
sa=`grep ">" genome.fna |sed 's/^.//'`

echo "check genome name $sa"
#sed -i '1,3d' tmp.gff
awk '{print "'$sa'""\t"$4"\t"$4}' tmp.gff > motif.bed
echo "head check methylation bed"
head motif.bed
echo "generate forward and reverse strand"
#positive bed
grep '[[:space:]]+\+' genome.gff  > positive.gtf
grep '[[:space:]]-\+' genome.gff  > negative.gtf
Rscript ./script/get.tss.r
echo "integrate regulaation region file"
awk '{print "'$sa'""\t"$1"\t"$2}' tss.po.bed > tss.bed
awk '{print "'$sa'""\t"$1"\t"$2}' tss.ne.bed >>tss.bed
echo "do a bedtools merge on ${SAMPLE} and regulation region"
bedtools intersect -a tss.bed -b motif.bed -wa > tss.methylation.bed
bedtools intersect -a tss.bed -b motif.bed -wa -wb > tss.methylation.me.bed
echo "generate a fatsa file means this regulation methylated"
sort -n tss.methylation.bed | uniq > tss.methylation.sorted.bed
bedtools getfasta -fi genome.fna  -bed tss.methylation.sorted.bed -fo > $SAMPLE.fimo.fasta
echo "find methylated genes"
cp tss.methylation.me.bed ./tmp.bed
awk '{print $2+100}' tmp.bed >tss.bed
awk '{print $3-100}' tmp.bed >>tss.bed

awk '{print $4"\t"$5}' genome.gff > test.bed
grep -wf tss.bed test.bed > over.txt
grep -wf over.txt genome.gff |cut -f9 |cut -f1 -d";"|sed -e 's/ID\=//g' > locus.tmp
grep -wf locus.tmp genome.tsv |cut -f4,7 > gene.rcc.list.txt
grep -v "hypothetical" gene.rcc.list.txt > RCC.list.txt
rm tss.bed
rm tss.po.bed
rm tss.ne.bed
rm tmp.gff
