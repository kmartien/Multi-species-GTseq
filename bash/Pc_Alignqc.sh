#!/bin/bash

#SBATCH --job-name=Pc_Alignqc
#SBATCH -e Pc_Alignqc_%j.e.txt
#SBATCH -o Pc_Alignqc_%j.log
#SBATCH --mail-user=keith.hernandez@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task 10
#SBATCH -t 3-00:00:00
#SBATCH -D /home/khernandez
#SBATCH --array=1-192%30

#load modules
module load bio/samtools
module load bio/bamtools
module load bio/picard
module load bio/bamutil

#Define variables
FILE=/home/khernandez/align_array.txt
IDS=$(cat ${FILE})

#Create job index to run through array and file names
for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	fq_r1=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID}==${job_index} ]]; then
	break
fi
done

sample_id=$(echo $fq_r1 | sed 's!^.*/!!')
sample_id=${sample_id%%_*}

#convert sam to bam
samtools view -bS -F 4 /home/khernandez/bamtools/Pc_${sample_id}.sam > /home/khernandez/bamtools/Pc_${sample_id}.bam

#sort bam file
samtools view -h /home/khernandez/bamtools/Pc_${sample_id}.bam | samtools view -buS | samtools sort -o /home/khernandez/bamtools/Pc_${sample_id}_sorted.bam

#Use Picard to remove duplicates
java -jar ${PICARD} MarkDuplicates I=/home/khernandez/bamtools/Pc_${sample_id}_sorted.bam O=/home/khernandez/bamtools/Pc_${sample_id}_sorted_dedup.bam M=/home/khernandez/bamtools/Pc_${sample_id}_dups.log VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

#Clip overlapping read pairs
bam clipOverlap --in /home/khernandez/bamtools/Pc_${sample_id}_sorted_dedup.bam --out /home/khernandez/bamtools/Pc_${sample_id}_sorted_dedup_clipped.bam --stats

#Calculate depth of coverage per sample and save to a compressed file for genotype likelihood calculations
samtools depth -aa /home/khernandez/bamtools/Pc_${sample_id}_sorted_dedup_clipped.bam | cut -f 3 | gzip > /home/khernandez/bamtools/Pc_${sample_id}.depth.gz

#Index the aligned bam files
samtools index /home/khernandez/bamtools/Pc_${sample_id}_sorted_dedup_clipped.bam
