#! /bin/bash

#SBATCH --job-name=PSMC
#SBATCH -e PSMC_%j.e.txt
#SBATCH -o PSMC_%j.log 
#SBATCH --mail-user= #email address
#SBATCH --mail-type=ALL  # (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -c 20 
#SBATCH --mem=80G
#SBATCH -t 3-0   # walltime in hrs:mins:secs or days-hours. 
#SBATCH -D /scratch/pmorin/temp # <folder to change to when starting the job>

##############################################################
# programs
module load bio/samtools/1.11
module load bio/bcftools/1.11
module load tools/gnuplot/5.4.1
PSMC=~/programs/psmc
set -eux

##############################################################
# variables for bam and reference files
BAMDIR=/home/pmorin/projects/Ddel/ERR11040185_map2ref			# BAM directory
BAMFILE=Ddel_ERR11040185_dedup_noRepeats.bam
REFDIR=/home/pmorin/Ref_genomes/Ddel/GCA_949987515.1_mDelDel1.1 # reference directory
REF=GCA_949987515.1_mDelDel1.1_genomic.fna						# reference genome
OUTDIR=${BAMDIR}/PSMC

GEN=14.8 											# generation time
COV=94												# average depth of coverage (integer)
THREADS=50 											# run 50 bootstraps in parallel

# Extract portions of the BAM filename (separated by "_") for output files.
ID=`echo $BAMFILE | cut -f1,2 -d "_"`

DIPLOID=${ID}_diploid.fq.gz

mkdir -p ${OUTDIR}

##############################################################
# PSMC parameters
TEMP_DIR=/scratch/pmorin/temp
MUT=4.90E-10 #mutation rate (µ/site/yr)
MUTRATE=$(echo $GEN $MUT | awk '{ printf "%.12f", $1*$2 }') # mutation rate (µ/site/yr) multiplied by generation time. 
t=15	# default from Li and Durbin is 15
PSMC_INT="4+25*2+4+6" # default from Li and Durbin is 4+25*2+4+6

# diploid genome parameters
MEANDEPTH=${COV} # has to be an integer for mindepth/maxdepth calc below.
let MIN=MEANDEPTH/3  # 1/3x average coverage
let MAX=MEANDEPTH*2  # 2x average coverage

#########################################################################
#########################################################################
### Create a diploid consensus fastq file from the bam file. This may take ≥24hr
samtools mpileup -uf ${REFDIR}/${REF} ${BAMDIR}/${BAMFILE} | bcftools call -c --threads ${THREADS} | vcfutils.pl vcf2fq -d${MIN} -D${MAX} > ${OUTDIR}/${DIPLOID}

#########################################################################
### PSMC
# generate psmcfa file from diploid genome file
${PSMC}/utils/fq2psmcfa -q20 ${OUTDIR}/${DIPLOID} > ${OUTDIR}/${ID}_diploid.psmcfa

# generate the psmc file using the default settings for humans (-N25, -t20, -r5).
${PSMC}/psmc -N25 -t${t} -r5 -p ${PSMC_INT} -o ${OUTDIR}/${ID}_${MUT}_t${t}.psmc ${OUTDIR}/${ID}_diploid.psmcfa

# Make psmc plots and adapt the scaling using this psmc file.
nice ${PSMC}/utils/psmc_plot.pl -u ${MUTRATE} -g ${GEN} -RM ${ID}"_"${MUT} ${OUTDIR}/${ID}_${MUT}_t${t}_psmc.out ${OUTDIR}/${ID}_${MUT}_t${t}.psmc

######################
### bootstrap PSMC

# split the PSMC file
${PSMC}/utils/splitfa ${OUTDIR}/${ID}_diploid.psmcfa > ${OUTDIR}/split_${ID}_diploid.psmcfa

# PSMC bootstrap, multithread (written to temp directory)
seq 100 | xargs -P ${THREADS} -i ${PSMC}/psmc -N25 -t${t} -r5 -b -p ${PSMC_INT} -o ${TEMP_DIR}/${ID}_round-{}.psmc ${OUTDIR}/split_${ID}_diploid.psmcfa | sh

# merge original.psmc round-*.psmc > merged.psmc
# within folder containing original sample psmc and bootstrap files
cat ${OUTDIR}/${ID}_${MUT}_t${t}.psmc ${TEMP_DIR}/${ID}_round-*.psmc > ${OUTDIR}/merged_${ID}_t${t}_boot.psmc

# And then plot it: (- -> Id placed on plot).
${PSMC}/utils/psmc_plot.pl -u ${MUTRATE} -g ${GEN} -RM "" ${OUTDIR}/merged_${ID}_t${t}_boot.out ${OUTDIR}/merged_${ID}_t${t}_boot.psmc

cp ${BAMDIR}/*PSMCboot_SEDNA.sh ${OUTDIR}
########################################################################################
# for PSMC parameter info, see https://github.com/lh3/psmc
# the `-p' option specifies that there are 64 atomic time intervals and 28 (=1+25+1+1)
# free interval parameters. The first parameter spans the first 4 atomic time
# intervals, each of the next 25 parameters spans 2 intervals, the 27th spans 4
# intervals and the last parameter spans the last 6 time intervals. The `-p' and
# `-t' options are manually chosen such that after 20 rounds of iterations, at
# least ~10 recombinations are inferred to occur in the intervals each parameter
# spans.

### use the custom script provide in the utils folder with PSMC to convert the fastq file to a pseudo-fasta file. It basically codes regions in 100bp bins as being heterozygote or homozygote.

