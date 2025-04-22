#!/bin/bash

# this is the header for SLURM
#SBATCH --job-name=Het_angsd    ## job name
#SBATCH -e Het_angsd_%j.e.txt    ## error message name
#SBATCH -o Het_angsd.log.%j.out  ## log file output name
#SBATCH --mail-user=karen.martien@noaa.gov
#SBATCH --mail-type=ALL 
#SBATCH -c 5 # <number of cores to ask for>
#SBATCH --mem=20G
#SBATCH -t 24:00:00   ## walltime in mins or mins:secs or hrs:mins:secs. 
#SBATCH --array=1-21 #corresponds to chromosomes/scaffolds in ref genome)
#SBATCH --ntasks=1  
#SBATCH -D /scratch/kmartien/temp		## <folder to change to when starting the job>

### Sliding window heterozygosity based on Angsd variant detection
##########################################################################################
# set parameters and folders
##########################################################################################
module load bio/samtools/1.11
module load bio/angsd/0.933 
module load bio/bedtools/2.29.2 
module load R/4.0.3

# R script for summarizing data from ANGSD
Het_plot=/home/kmartien/scripts/misc/ROH_angsd_plot_chromosomes.R
set -eux

###########################################
# variables for bam and reference files
speciesID=Pcra													# species
#
### NEED VALUES FOR THE NEXT THREE LINES ###########
COV=94															# mean depth of coverage
BAMDIR=/home/pmorin/projects/Ddel/ERR11040185_map2ref			# BAM Directory
BAMFILE=Ddel_ERR11040185_dedup_noRepeats.bam
#############################################
#
REFDIR=../home/khernandez/Genomes # Reference directory
REF=GCF_039906515.1_mPseCra1.hap1_genomic.fna						# Reference genome

OUTDIR=${BAMDIR}/heterozygosity 
TEMPDIR=/scratch/kmartien/temp
LOG=${OUTDIR}/log_${now}
mkdir -p ${BAMDIR}/heterozygosity 
mkdir -p ${OUTDIR}/log_${now}

# get Bam file based on unique part of file name, and extract portions of the name for output files.
short_ID=`echo $BAMFILE | cut -f1,2,3 -d "_"`
sample=${short_ID}

###########################################
# variables (these are needed to make the output file names unique, needed for running 
# multiple runs simultaneously)
MEANDEPTH=${COV} # has to be an integer for mindepth/maxdepth calc below.
	let MINDEPTH=MEANDEPTH/3  # 1/3x average coverage
	let MAXDEPTH=MEANDEPTH*2  # 2x average coverage

genome=RM # RM for repeat-masked, FULL for complete genome  
winsize=1000  #100 (kb) or 1000 (kb=1MB)   
MBQ=20  # minimum base quality filter (Foote et al. 2021 KW ROH pub = 20)
MAPQ=30  # minimum map quality filter (Foote et al. 2021 KW ROH pub = 30)
THREADS=5
now=$(date +'%d.%m.%Y_%H.%M')

RefName=`echo $REF | cut -f1,2 -d "_"`

OUTFILE=${speciesID}_${sample}_${genome}_${winsize}_map2${RefName}_${now}

WINDOWSLIST=${REF}_1MB_windows.txt    
# WINDOWSLIST=${REF}_100kb_windows.txt   
#
#### HOW DO I GENERATE THIS FILE? ############
SCAFFOLDLIST=${REF}_1MB_scaffold.lengths.txt # see README for how to generate list. This list is the same regardless of window length (just lengths of scaffolds)

################################################################
# Get working scaffold based on array number
NUM=$(printf %02d ${SLURM_ARRAY_TASK_ID})

# generate windows list for each scaffold
CHR=$(head -n ${NUM} ${REFDIR}/${SCAFFOLDLIST} | tail -n 1)
CHR=$(echo ${CHR} | awk -F " " '{ print $1 }')
CHR=$(echo ${CHR} | awk -F " " '{$0=$0":"}{print}')
echo ${CHR}

grep ${CHR} ${REFDIR}/${WINDOWSLIST} > ${OUTDIR}/${REF}_${NUM}_windows.txt

################################################################
# Run angsd
################################################################
# .saf files are put into temp directory (scratch) to speed up processing and reduce clutter in output directory. See http://www.popgen.dk/angsd/index.php/Input#Genotype_Likelihood_Files for explanation of options.
# -doSaf 1: Calculate the Site allele frequency likelihood based on individual genotype likelihoods assuming HWE
# If you don't have the ancestral state, you can instead estimate the folded SFS. This is done by supplying the -anc with the reference genome and applying -fold 1 to realSFS.
# baq = base alignment quality. 0: No BAQ calcualtion. 1:normal BAQ (same as default in SAMtools). 2:extended BAQ (same as default in SAMtools).

angsd -GL 2 -setMinDepth ${MINDEPTH} -setMaxDepth ${MAXDEPTH} -minmapq ${MAPQ} -minq ${MBQ} -uniqueonly 1 -remove_bads 1 -only_proper_pairs 1 -C 50 -docounts 1 -i ${BAMDIR}/${BAMFILE} -ref ${REFDIR}/${REF} -P ${THREADS} -out ${TEMPDIR}/${OUTFILE}_${NUM} -doSaf 1 -anc ${REFDIR}/${REF} -rf ${OUTDIR}/${REF}_${NUM}_windows.txt -baq 2 -fold 1 &> ${LOG}/angsd.${NUM}.log

################################################################
# Run realSFS
################################################################
while read -r line; do realSFS -r $line ${TEMPDIR}/${OUTFILE}_${NUM}.saf.idx  -P ${THREADS} -tole 1e-8 2> log >> ${OUTDIR}/${OUTFILE}_combined_sfs_${NUM}.txt; done < ${OUTDIR}/${REF}_${NUM}_windows.txt
# -P = number of threads

################################################################
# Calculate heterozygosity (ends up being the value in column 4)
################################################################
awk '{print $1,$2,$3=$1+$2,$4=$2/$3}' ${OUTDIR}/${OUTFILE}_combined_sfs_${NUM}.txt >> ${OUTDIR}/${OUTFILE}_ROH_${NUM}_summary.het

########### Run in R script to plot and summarize outfile_ROH_summary.het #############
cd ${OUTDIR}
Rscript ${Het_plot} ${speciesID} ${sample} ${genome} ${winsize}
