#!/bin/bash
#SBATCH --ntasks=1                                  # Number of tasks (usually set to 1 for serial jobs)
#SBATCH --cpus-per-task=4                           # Number of CPU cores per task
#SBATCH --mem=1G                                    # Memory per node (e.g., 16GB)
#SBATCH --partition=short                           # Partition name
#SBATCH --time=03:00:00                             # Time

ACCESSION=$1
echo `date`": "${ACCESSION}" Begins"

cd /n/data1/hms/dbmi/farhat/jiazheng/eqtl/rna_analysis
module load fastqc/0.12.1
module load gcc/6.2.0
module load bwa/0.7.17
module load samtools/1.15.1
export PATH=$PATH:/home/jim744/subread-2.0.6-source/bin/

cp ../rnaseq/${ACCESSION}* `pwd`

echo `date`": "${ACCESSION}" fastqc"
fastqc ${ACCESSION}_*.fastq

echo `date`": "${ACCESSION}" bwa"
bwa mem -t 4 ref/h37rv.fna ${ACCESSION}_1.fastq ${ACCESSION}_2.fastq > ${ACCESSION}.sam

echo `date`": "${ACCESSION}" samtools"
samtools view -bS -@ 4 ${ACCESSION}.sam > ${ACCESSION}.bam
samtools sort -@ 4 ${ACCESSION}.bam -o ${ACCESSION}_sorted.bam
samtools index -@ 4 ${ACCESSION}_sorted.bam

echo `date`": "${ACCESSION}" featureCounts"
featureCounts -T 4 -p -t gene -G ref/h37rv.fna -a ref/h37rv.gtf -o ${ACCESSION}_counts.txt ${ACCESSION}_sorted.bam

mkdir align_results/$ACCESSION
mv $ACCESSION* align_results/$ACCESSION
echo `date`": "${ACCESSION}" Ends"
