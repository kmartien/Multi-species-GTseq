#!/bin/bash

# this is the header for SEDNA
#SBATCH --job-name=repeatmask   ## job name
#SBATCH -e repeatmask_%j.e.txt    ## error message name
#SBATCH -o repeatmask.log.%j.out  ## log file output name
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  ## (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -c 20    ## <number of cores to ask for> (cpu cores per task; 20/node)
#SBATCH -p medmem
# SBATCH --mem=180G #(max 96-192/node)
#SBATCH -t 20-0   ## walltime in mins or mins:secs or hrs:mins:secs, or days-hours. 
#SBATCH --ntasks=1                   ## Run a single task (from Morgan's script)
#SBATCH -D /scratch/pmorin/temp		## <folder to change to when starting the job>

##########################################################################################
# run repeatmasker on genome assembly
# run using: 
# script=/home/pmorin/scripts/misc/cetacea_RepeatMasker4.1.6_Dfam3.8_SEDNA.sh
# sbatch ${script} /home/pmorin/Ref_genomes/Ddel/GCA_949987515.1_mDelDel1.1


module load bio/repeatmasker/4.1.6
source /opt/bioinformatics/venv/repeatmasker-4.1.6/bin/activate

REFDIR=${1}

REF=*genomic.fna
OUTDIR=${REFDIR}/RepeatMask4.1.6_cetacea
SP=cetacea
THREADS=10 # should be half the number of cores requested
now=$(date +'%d.%m.%Y_%H.%M')

cd ${REFDIR}
mkdir -p ${OUTDIR}
# cp *RepeatMasker4.1.6* ${OUTDIR}
#########################################################################
# save copy of repeats for the species (or species group) indicated.
famdb.py -i /opt/bioinformatics/bio/repeatmasker/repeatmasker-4.1.6/Libraries/famdb lineage ${SP} -d > ${OUTDIR}/${SP}_Dfam3.8_list_${now}.txt

echo $REF

RepeatMasker -species ${SP} -pa ${THREADS} -s -dir ${OUTDIR} ${REFDIR}/${REF}
# rm generates 2x as many threads as -pa, so set to 10 when using 20 cores

# -pa(rallel) [number]
#         The number of processors to use in parallel (only works for batch
#         files or sequences over 50 kb)
# 
#     -s  Slow search; 0-5% more sensitive, 2-3 times slower than default
# 
#     -q  Quick search; 5-10% less sensitive, 2-5 times faster than default
# 
#     -qq Rush job; about 10% less sensitive, 4->10 times faster than default
#         (quick searches are fine under most circumstances) repeat options

# species repeat database options:
#mammals
#odontoceti
#cetacea

# https://github.com/rmhubley/RepeatMasker/blob/master/repeatmasker.help
#  (copied from help)
# Interspersed repeats mostly are copies of transposable elements in
# different states of erosion. Thus, dependent on the time of activity
# of the source transposable element, interspersed repeats generally are
# specific to a (clade of) species, and different redatabase
# (http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html). In
# principal, all unique clade names occurring in this database can be
# used. Examples are:
# 
# -species "sus scrofa"
# -species chimpanzee
# -species arabidopsis
# -species canidae
# -species mammals
# 
# Capitalization is ignored, multiple words need to bound by apostrophes.
# 
# RepeatMasker builds one or more repeat consensus files the first time
# a species/group has been chosen, or when a new database has been
# downloaded. These will be written in a subdirectory of the Libraries
# directory named after the date of the repeat database version and the
# Latin name of the clade. For example, "-species monocotyledons"
# creates the file
# "..../RepeatMasker/Libraries/20040616/liliopsida/specieslib". 
# Currently, only for mammalian species multiple files are created,
# bearing names like "shortcutlib" and "longlib", which the queries are
# compared to sequentially.
