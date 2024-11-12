#!/bin/bash
#SBATCH --job-name=DOWNLOAD
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=short

while getopts i:o: flag
do
    case "${flag}" in
        i) INPUT_FILE=${OPTARG};;
        o) OUTPUT_FOLDER=${OPTARG};;
    esac
done

echo `date`": START"
cd $OUTPUT_FOLDER
module load sratoolkit/2.10.7
for ACCESSION in `cat $INPUT_FILE`
do
prefetch $ACCESSION
fastq-dump --split-files $ACCESSION
rm -rf $ACCESSION
done
echo `date`": DONE"