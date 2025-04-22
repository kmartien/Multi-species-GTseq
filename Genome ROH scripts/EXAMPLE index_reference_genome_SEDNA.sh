#! /bin/bash

# Sedna has 20 cores per node, 
# low-memory nodes: 96 Gb RAM, 20 compute cores/node (28 nodes) (node names = 1 - 28)
# high-memory nodes: n=3: 1.5 Tb RAM, 24 cores each; local scratch storage

# this is the header for SEDNA
#SBATCH --job-name=index_ref
#SBATCH -e index_ref_%j.e.txt
#SBATCH -o index_ref_%j.log
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  # (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -c 20    # <number of cores to ask for>
#SBATCH --mem=80G   # job memory request increase to 32G after testing
#SBATCH -t 24:00:00   # walltime in mins or mins:secs or hrs:mins:secs. Increase to 48hr after testing
# SBATCH -D /scratch/pmorin/temp		# <folder to change to when starting the job>

# not used:
# SBATCH -p himem  # this would request the hi-memory nodes, if >96G memory needed.
# SBATCH --tasks-per-node=8
################
# example usage:
# cd ~/Ref_genomes/Psin/mPhoSin1.pri_ref_genome
# SCRIPT=~/scripts/misc/index_reference_genome_SEDNA.sh
# FASTA=GCF_008692025.1_mPhoSin1.pri_genomic.fna
# sbatch ${SCRIPT} ${FASTA}
################
# module load java # not needed on SEDNA
module load aligners/bwa-mem2/2.2.1
module load bio/samtools/1.11
module load bio/picard/2.23.9
module load bio/bedtools/2.29.2 
module load R/4.0.3

TEMP_DIR=/scratch/pmorin/temp

FASTA=${1}
echo "Preparing index files for ${FASTA}"

echo ; echo "Generating BWA index files"
bwa-mem2 index ${FASTA} 

echo ; echo "Generating Samtools index"
samtools faidx ${FASTA}

# This requires that the fasta file be decompressed (gunzip)
echo ; echo "Generating sequence dictionary"
java -Xmx8G -jar -Djava.io.tmpdir=${TEMP_DIR} \
${PICARD} CreateSequenceDictionary REFERENCE=${FASTA} OUTPUT=${FASTA%.f*}.dict 

################################################################
# Get scaffold lengths for those that are over 1MB
################################################################

awk '$2 > 1000000 {print $1"\t"$2}' ${1}.fai > ${1}_1MB_scaffold.lengths.txt

################################################################
# Get sliding windows of 1,000,000 bp
################################################################

bedtools makewindows -g ${1}_1MB_scaffold.lengths.txt -w 1000000 | awk '$3 ~ "000000" {print$1":"$2"-"$3}' > ${1}_1MB_windows.txt

bedtools makewindows -g ${1}_1MB_scaffold.lengths.txt -w 100000 | awk '$3 ~ "00000" {print$1":"$2"-"$3}' > ${1}_100kb_windows.txt

################################################################
