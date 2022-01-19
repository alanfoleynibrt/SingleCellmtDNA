#!/usr/bin/env bash

samtools faidx reference_genome/cgr_reference.fasta
samtools faidx reference_genome/cgr_reference.fasta cgr_mtdna:8000-16283 > reference_genome/first_segment.fasta
samtools faidx reference_genome/cgr_reference.fasta cgr_mtdna:1-7999 | grep -v ">" > reference_genome/second_segment.fasta

cat reference_genome/first_segment.fasta reference_genome/second_segment.fasta > reference_genome/shifted_reference_sequence.fasta

cat reference_genome/shifted_reference_sequence.fasta | awk '{if (substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' > reference_genome/joinedlineoutput.fasta

mv reference_genome/joinedlineoutput.fasta reference_genome/shift_cgr_reference.fasta
grep -v ">" reference_genome/shift_cgr_reference.fasta | wc | awk '{print $3-$1}'

#bwa index -a is -p cgridx shift_cgr_reference.fasta

# samtools fasta index
#samtools faidx shift_cgr_reference.fasta

#sequence dictonary
#java -jar ~/picard-tools-1.119/CreateSequenceDictionary.jar R=shift_cgr_reference.fasta O=shift_cgr_reference.dict


#check shifted length
#grep -v ">" reference_genome/first_segment.fasta | wc | awk '{print $3-$1}'
#grep -v ">" reference_genome/second_segment.fasta | wc | awk '{print $3-$1}'


rm reference_genome/first_segment.fasta reference_genome/second_segment.fasta
