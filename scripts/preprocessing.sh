#!/usr/bin/env bash
#### Trim low quality bases using trimmomatic
#### inputs are: 1) sample ID and 2) the data directory 3) output directory
#### & 4) number of threads
#### Written by NIBRT: colin.clarke@nibrt.ie 12-2019
if (($# == 0)); then
        echo "Usage:"
        echo "-s = sample ID"
        echo "-i = input directory"
        echo "-o = output directory"
        echo "-p = num threads"
        exit 2
fi

while getopts s:i:o:p: option
  do
    case "${option}"
      in
      s) SAMPLE_ID=${OPTARG};;
      i) IN_DIR=${OPTARG};;
      o) OUT_DIR=${OPTARG};;
      p) NUM_THREADS=${OPTARG};;
    esac
done

mkdir -p $OUT_DIR
mkdir -p $OUT_DIR/paired $OUT_DIR/unpaired

java -jar  /home/alan/KX576660/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
     -threads $NUM_THREADS \
     $IN_DIR/"$SAMPLE_ID"_L001_R1_001.fastq.gz  $IN_DIR/"$SAMPLE_ID"_L001_R2_001.fastq.gz  \
     $OUT_DIR/paired/"$SAMPLE_ID"_L001_R1_001.fastq.gz  $OUT_DIR/unpaired/"$SAMPLE_ID"_L001_R1_001.fastq.gz  \
     $OUT_DIR/paired/"$SAMPLE_ID"_L001_R2_001.fastq.gz  $OUT_DIR/unpaired/"$SAMPLE_ID"_L001_R2_001.fastq.gz  \
     ILLUMINACLIP:IlluminaDNAPrep.fa:2:30:10 \
     LEADING:3 \
     TRAILING:3 \
     MINLEN:50 \
     -trimlog $OUT_DIR/"$SAMPLE_ID".trimmomatic.log