#!/bin/bash
#SBATCH --ntasks=1

while getopts i:a:c:r:m: flag
do
    case "${flag}" in
        i) INPUT_PREFIX=${OPTARG};;
        a) ACCESSION=${OPTARG};;
        c) NUM_CORE=${OPTARG};;
        r) REFERENCE_GENOME_FOLDER=${OPTARG};;
        m) MEM=${OPTARG};;
    esac
done

echo `date`": "${ACCESSION}" Begins"

module load gcc/6.2.0
module load bwa/0.7.17
module load samtools/1.15.1
module load java/jdk-21.0.2
module load bcftools

cd /n/data1/hms/dbmi/farhat/jiazheng/eqtl/dna_analysis
mkdir ${ACCESSION}
cd ${ACCESSION}
bwa mem -t $NUM_CORE ${REFERENCE_GENOME_FOLDER}/h37rv.fna ${INPUT_PREFIX}/${ACCESSION}_1.fastq > ${ACCESSION}.sam
samtools view -bS -@ $NUM_CORE ${ACCESSION}.sam > ${ACCESSION}.bam
samtools sort -@ $NUM_CORE ${ACCESSION}.bam -o ${ACCESSION}_sorted.bam
samtools index -@ $NUM_CORE ${ACCESSION}_sorted.bam
samtools depth -@ $NUM_CORE -a ${ACCESSION}_sorted.bam > ${ACCESSION}_coverage.txt
cp /home/jim744/standalone/pilon-1.24.jar `pwd`
java -Xmx${MEM} -jar pilon-1.24.jar --genome ${REFERENCE_GENOME_FOLDER}/h37rv.fna --unpaired ${ACCESSION}_sorted.bam --output $ACCESSION --vcf
bcftools view -i 'ALT!="."' ${ACCESSION}.vcf -o ${ACCESSION}_filtered.vcf
rm pilon-1.24.jar
echo `date`": "${ACCESSION}" Ends"

