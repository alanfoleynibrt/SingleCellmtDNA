
---
### "Single CHO cell mtDNA pipeline"

author: "Alan Foley, adapted from Colin Clarke's code"

date: "01/02/2022"

output: html_document

---

## Dependencies
snpEff, VarScan, GenomeAnalysisTK, picard tools, lofreq_star-2.1.2


## Trim adaptor/index sequences
```bash
mkdir data/raw
```

```bash
cat data/sample_info.txt | while read sample; do
./scripts/preprocessing.sh -s $sample  -i data/raw -o data/preprocessed -p 32
done
```


## Prepare reference indexes
```bash
samtools faidx reference_genome/unshifted/KX576660.fasta
samtools faidx reference_genome/shifted/KX576660.fasta

java -jar  /home/alan/bin/picard.jar CreateSequenceDictionary \
R=reference_genome/unshifted/KX576660.fasta \
O=reference_genome/unshifted/KX576660.dict

java -jar  /home/alan/bin/picard.jar CreateSequenceDictionary \
R=reference_genome/shifted/KX576660.fasta \
O=reference_genome/shifted/KX576660.dict
```

## Mapping

## create the mtDNA reference indexes for unshifted & shifted mtDNA sequences
```bash
mkdir -p bwa_index/unshifted/
bwa index -a is -p  bwa_index/unshifted/cgridx \
  reference_genome/unshifted/KX576660.fasta
  
mkdir bwa_index/shifted
bwa index -a is -p  bwa_index/shifted/cgridx \
  reference_genome/shifted/KX576660.fasta
```

## mtDNA genome mapping
```bash
for mode in unshifted shifted; do
  mkdir -p bwa_mapping/$mode
  cat data/sample_info.txt | while read sample; do
    # map
    bwa mem -M -t 32 \
    bwa_index/$mode/cgridx \
    data/preprocessed/paired/"$sample"_L001_R1_001.fastq.gz \
    data/preprocessed/paired/"$sample"_L001_R2_001.fastq.gz \
    > bwa_mapping/$mode/"$sample".sam;

    # sort
    java -jar  /home/alan/bin/picard.jar SortSam \
    INPUT=bwa_mapping/$mode/"$sample".sam \
    OUTPUT=bwa_mapping/$mode/"$sample".unfiltered.bam \
    SORT_ORDER=coordinate;

    # filter
    samtools view -bq 20 \
    bwa_mapping/$mode/"$sample".unfiltered.bam \
    > bwa_mapping/$mode/"$sample".bam

    # remove SAM
    rm bwa_mapping/$mode/"$sample".sam # remove the sam file
  done
done
```

### Mark duplicates
```bash
for mode in unshifted shifted; do
  mkdir bwa_mapping/$mode/duplicate_marked bwa_mapping/$mode/added_read_group
  cat data/sample_info.txt | while read sample; do
    # Mark duplicates
    java -jar  /home/alan/bin/picard.jar MarkDuplicates \
    INPUT=bwa_mapping/$mode/"$sample".bam \
    OUTPUT=bwa_mapping/$mode/duplicate_marked/"$sample".bam \
    METRICS_FILE=/dev/null

    # Add read groups
    java -jar  /home/alan/bin/picard.jar AddOrReplaceReadGroups \
    I=bwa_mapping/$mode/duplicate_marked/"$sample".bam \
    O=bwa_mapping/$mode/added_read_group/"$sample".bam \
    RGSM="$sample" \
    RGLB=mtdna_seq \
    RGPL=illumina \
    RGPU=none \
    VALIDATION_STRINGENCY=LENIENT;
  done
done
```

## Indel realignment
```bash
for mode in unshifted shifted; do
  mkdir -p variant_calling/$mode/indel_realigned/intervals
  mkdir -p variant_calling/$mode/indel_realigned/bam_files
  cat data/sample_info.txt | while read sample; do
    java -jar  /home/alan/bin/picard.jar BuildBamIndex \
    INPUT=bwa_mapping/$mode/added_read_group/$sample.bam \

    sudo docker run --rm -v /home/alan/KX576660:/home broadinstitute/gatk3:3.8-0 \
    java -jar GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -I ../home/bwa_mapping/$mode/added_read_group/"$sample".bam \
    -R ../home/reference_genome/$mode/KX576660.fasta \
    -o ../home/variant_calling/$mode/indel_realigned/intervals/$sample.IndelRealigner.intervals

    sudo docker run --rm -v /home/alan/KX576660:/home broadinstitute/gatk3:3.8-0 \
    java -jar  /usr/GenomeAnalysisTK.jar \
    -T IndelRealigner \
      -I ../home/bwa_mapping/$mode/added_read_group/"$sample".bam \
      -R ../home/reference_genome/$mode/KX576660.fasta \
      --maxReadsForRealignment 15000000 \
      --maxReadsInMemory 150000000 \
      -targetIntervals ../home/variant_calling/$mode/indel_realigned/intervals/$sample.IndelRealigner.intervals \
    -o  ../home/variant_calling/$mode/indel_realigned/bam_files/$sample.bam
  done
done
```

## Base recalibration
```bash
for mode in unshifted shifted; do
mkdir -p variant_calling/$mode/base_recalibration/vcf
mkdir variant_calling/$mode/base_recalibration/bam_files

  cat data/sample_info.txt | while read sample; do
    # initial snp calling with lofreq
    sudo ./lofreq_star-2.1.2/bin/lofreq call-parallel \
    --pp-threads 12 \
    -f reference_genome/$mode/KX576660.fasta \
    variant_calling/$mode/indel_realigned/bam_files/$sample.bam \
    -o variant_calling/$mode/base_recalibration/vcf/$sample.raw.lofreq.vcf

    sudo ./lofreq_star-2.1.2/bin/lofreq filter \
    -i variant_calling/$mode/base_recalibration/vcf/$sample.raw.lofreq.vcf \
    -o variant_calling/$mode/base_recalibration/vcf/$sample.base.recal.lofreq

    # apply "ground truth" SNPs
    sudo docker run --rm -v /home/alan/KX576660:/home broadinstitute/gatk3:3.8-0 \
      java -jar  /usr/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R ../home/reference_genome/$mode/KX576660.fasta \
    -I ../home/variant_calling/$mode/indel_realigned/bam_files/$sample.bam \
    -knownSites ../home/variant_calling/$mode/base_recalibration/vcf/$sample.base.recal.lofreq \
    -o ../home/variant_calling/$mode/base_recalibration/vcf/"$sample".recal_data.table

    # assess the effect of base recalibration
    sudo docker run --rm -v /home/alan/KX576660:/home broadinstitute/gatk3:3.8-0 \
    java -jar  /usr/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R ../home/reference_genome/$mode/KX576660.fasta \
    -I ../home/variant_calling/$mode/indel_realigned/bam_files/$sample.bam \
    -knownSites ../home/variant_calling/$mode/base_recalibration/vcf/$sample.base.recal.lofreq \
    -BQSR ../home/variant_calling/$mode/base_recalibration/vcf/$sample.recal_data.table \
    -o ../home/variant_calling/$mode/base_recalibration/vcf/$sample.post.recal.data.table

    sudo docker run --rm -v /home/alan/KX576660:/home broadinstitute/gatk3:3.8-0 \
    java -jar  /usr/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R ../home/reference_genome/$mode/KX576660.fasta \
    -I ../home/variant_calling/$mode/indel_realigned/bam_files/$sample.bam \
    -BQSR ../home/variant_calling/$mode/base_recalibration/vcf/$sample.recal_data.table \
    -o  ../home/variant_calling/$mode/base_recalibration/bam_files/$sample.bam
  done
done
```

# Variant calling
## Lofreq
```bash
for mode in unshifted shifted; do
mkdir -p variant_calling/$mode/snps/lofreq
mkdir -p variant_calling/$mode/indels/lofreq

  cat data/sample_info.txt | while read sample; do
      sudo ./lofreq_star-2.1.2/bin/lofreq \
      call-parallel \
       --pp-threads 24 \
       variant_calling/$mode/base_recalibration/bam_files/$sample.bam \
       -f reference_genome/$mode/KX576660.fasta \
       -o variant_calling/$mode/snps/lofreq/$sample.lofreq

       sudo ./lofreq_star-2.1.2/bin/lofreq call --call-indels --only-indels \
        variant_calling/$mode/base_recalibration/bam_files/$sample.bam \
       -f reference_genome/$mode/KX576660.fasta \
       -o variant_calling/$mode/indels/lofreq/$sample.lofreq

       sudo ./lofreq_star-2.1.2/bin/lofreq filter \
       -i variant_calling/$mode/snps/lofreq/$sample.lofreq \
       -o variant_calling/$mode/snps/lofreq/$sample.filter.lofreq \
       -a 0.02 \
       -v 500

       sudo ./lofreq_star-2.1.2/bin/lofreq filter \
       -i variant_calling/$mode/indels/lofreq/$sample.lofreq \
       -o variant_calling/$mode/indels/lofreq/$sample.filter.lofreq \
       -a 0.02 \
       -v 500
  done
done
```
## Varscan
Identify variants using the Varscan algorithm
```bash
for mode in unshifted shifted; do
mkdir -p variant_calling/$mode/snps/varscan
mkdir -p variant_calling/$mode/indels/varscan
  cat data/sample_info.txt | while read sample; do
    samtools mpileup -B  -d 10000000 \
    -f reference_genome/$mode/KX576660.fasta \
    variant_calling/$mode/base_recalibration/bam_files/$sample.bam | \
    java -jar VarScan.v2.3.9.jar \
    mpileup2snp \
    --min-avg-qual 20 \
    --min-coverage 500 \
    --min-var-freq 0.02 \
    --p-value 0.05 \
    --output-vcf 1 > variant_calling/$mode/snps/varscan/$sample.filter.varscan

  samtools mpileup -B -d 10000000 \
  -f reference_genome/$mode/KX576660.fasta \
  variant_calling/$mode/base_recalibration/bam_files/$sample.bam | \
  java -jar VarScan.v2.3.9.jar \
  mpileup2indel \
  --min-avg-qual 20 \
  --min-coverage  500 \
  --min-var-freq 0.02 \
  --p-value 0.05 \
  --output-vcf 1 > variant_calling/$mode/indels/varscan/$sample.filter.varscan  
  done
done
```
## compare the result of varscan and lofreq
Identify the variants identified by both algorithms
```bash
for mode in unshifted shifted; do
mkdir -p variant_calling/$mode/calling_overlap/snp/common_sites
mkdir -p variant_calling/$mode/calling_overlap/indels/common_sites
  cat data/sample_info.txt | while read sample; do
  # SNPs
  vcftools --vcf variant_calling/$mode/snps/lofreq/$sample.filter.lofreq \
  --diff variant_calling/$mode/snps/varscan/$sample.filter.varscan \
  --diff-site \
  --out variant_calling/$mode/calling_overlap/snp/common_sites/$sample.snps
  # INDELs
  vcftools --vcf variant_calling/$mode/indels/lofreq/$sample.filter.lofreq \
  --diff variant_calling/$mode/indels/varscan/$sample.filter.varscan \
  --diff-site \
  --out variant_calling/$mode/calling_overlap/indels/common_sites/$sample.indels
  done
done
```

# Create VCF with SNPs/indels identified by both algorithms

```bash
for mode in unshifted shifted; do
mkdir -p variant_calling/$mode/calling_overlap/snp/common_vcf/snp_positions
mkdir -p variant_calling/$mode/calling_overlap/indels/common_vcf/indel_positions

  # modify the lofreq vcf for annotation of variants using snpEff
  cat data/sample_info.txt | while read sample; do
  awk '{print $3}' variant_calling/$mode/calling_overlap/snp/common_sites/$sample.snps.diff.sites_in_files >  \
  variant_calling/$mode/calling_overlap/snp/common_vcf/snp_positions/current.snp.positions

  awk '{print $3}' variant_calling/$mode/calling_overlap/indels/common_sites/$sample.indels.diff.sites_in_files > \
  variant_calling/$mode/calling_overlap/indels/common_vcf/indel_positions/current.indel.positions

  if [ "$mode" = "unshifted" ]; then
    awk 'BEGIN{while((getline<"variant_calling/unshifted/calling_overlap/snp/common_vcf/snp_positions/current.snp.positions")>0)l[$1]=1}/^#/||l[$2]' variant_calling/$mode/snps/lofreq/$sample.filter.lofreq > variant_calling/$mode/calling_overlap/snp/common_vcf/$sample.common.vcf
  else
  awk 'BEGIN{while((getline<"variant_calling/shifted/calling_overlap/snp/common_vcf/snp_positions/current.snp.positions")>0)l[$1]=1}/^#/||l[$2]' variant_calling/$mode/snps/lofreq/$sample.filter.lofreq > variant_calling/$mode/calling_overlap/snp/common_vcf/$sample.common.vcf
  fi

  if [ "$mode" = "unshifted" ]; then
  awk 'BEGIN{while((getline<"variant_calling/unshifted/calling_overlap/indels/common_vcf/indel_positions/current.indel.positions")>0)l[$1]=1}/^#/||l[$2]' variant_calling/$mode/indels/lofreq/$sample.filter.lofreq > variant_calling/$mode/calling_overlap/indels/common_vcf/$sample.common.vcf
  else
  awk 'BEGIN{while((getline<"variant_calling/shifted/calling_overlap/indels/common_vcf/indel_positions/current.indel.positions")>0)l[$1]=1}/^#/||l[$2]' variant_calling/$mode/indels/lofreq/$sample.filter.lofreq > variant_calling/$mode/calling_overlap/indels/common_vcf/$sample.common.vcf
  fi

  done
  rm -r variant_calling/$mode/calling_overlap/snp/common_vcf/snp_positions
  rm -r variant_calling/$mode/calling_overlap/indels/common_vcf/indel_positions
done
```

# Reorientate the variant positions in the shifted sequence
```bash
mkdir variant_calling/shifted/calling_overlap/snp/common_modified_vcf
mkdir variant_calling/shifted/calling_overlap/indels/common_modified_vcf
vcf_path=variant_calling/shifted/calling_overlap/

cat data/sample_info.txt | while read sample; do
echo $sample
grep "#" $vcf_path/snp/common_vcf/$sample.common.vcf > \
$vcf_path/snp/common_modified_vcf/current_vcf_header

grep -v "#" $vcf_path/snp/common_vcf/$sample.common.vcf >  \
$vcf_path/snp/common_modified_vcf/current_variant_information

awk '{if (($2 > 3999 && $2  < 8283) || ($2 > 8283 && $2< 12284)) print $0; }' \
$vcf_path/snp/common_modified_vcf/current_variant_information > \
$vcf_path/snp/common_modified_vcf/inrange.positions

awk '{ if ($2 > 3999 && $2 < 8283 ) print $2+7999; else if ($2 > 8283 && $2< 12284) print $2-8284;}'   \
$vcf_path/snp/common_modified_vcf/inrange.positions > \
$vcf_path/snp/common_modified_vcf/replacement_column
#> reorientated_vcf/snps/"$sample_name".replacement_column

paste $vcf_path/snp/common_modified_vcf/inrange.positions \
$vcf_path/snp/common_modified_vcf/replacement_column |\
awk -v OFS='\t' '{ print $1, $9, $3, $4, $5, $6, $7, $8}' > \
$vcf_path/snp/common_modified_vcf/$sample.vcf

cat $vcf_path/snp/common_modified_vcf/current_vcf_header \
$vcf_path/snp/common_modified_vcf/$sample.vcf > \
$vcf_path/snp/common_modified_vcf/$sample.common.vcf

sed s/:8000-16283// $vcf_path/snp/common_modified_vcf/$sample.common.vcf > \
$vcf_path/snp/common_modified_vcf/$sample.annotated.common.vcf

rm $vcf_path/snp/common_modified_vcf/$sample.common.vcf
rm $vcf_path/snp/common_modified_vcf/$sample.vcf
rm $vcf_path/snp/common_modified_vcf/replacement_column
rm $vcf_path/snp/common_modified_vcf/current_vcf_header; \
rm $vcf_path/snp/common_modified_vcf/current_variant_information; \
rm $vcf_path/snp/common_modified_vcf/inrange.positions;

# indels
grep "#" $vcf_path/indels/common_vcf/$sample.common.vcf > \
$vcf_path/indels/common_modified_vcf/current_vcf_header

grep -v "#" $vcf_path/indels/common_vcf/$sample.common.vcf >  \
$vcf_path/indels/common_modified_vcf/current_variant_information

awk '{if (($2 > 3999 && $2 < 8283) || ($2 > 8283 && $2< 12284)) print $0; }' \
$vcf_path/indels/common_modified_vcf/current_variant_information > \
$vcf_path/indels/common_modified_vcf/inrange.positions

awk '{ if ($2 > 3999 && $2 < 8283 ) print $2+7999; else if ($2 > 8283 && $2< 12284) print $2-8284;}' \
$vcf_path/indels/common_modified_vcf/inrange.positions > \
$vcf_path/indels/common_modified_vcf/replacement_column
#> reorientated_vcf/snps/$sample.replacement_column

paste $vcf_path/indels/common_modified_vcf/inrange.positions \
$vcf_path/indels/common_modified_vcf/replacement_column | \
awk -v OFS='\t' '{ print $1, $9, $3, $4, $5, $6, $7, $8}' > \
$vcf_path/indels/common_modified_vcf/$sample.vcf

cat $vcf_path/indels/common_modified_vcf/current_vcf_header \
$vcf_path/indels/common_modified_vcf/$sample.vcf > \
$vcf_path/indels/common_modified_vcf/$sample.common.vcf

sed s/:8000-16283// $vcf_path/indels/common_modified_vcf/$sample.common.vcf > \
$vcf_path/indels/common_modified_vcf/$sample.annotated.common.vcf

rm $vcf_path/indels/common_modified_vcf/$sample.common.vcf
rm $vcf_path/indels/common_modified_vcf/$sample.vcf
rm $vcf_path/indels/common_modified_vcf/replacement_column
rm $vcf_path/indels/common_modified_vcf/current_vcf_header; \
rm $vcf_path/indels/common_modified_vcf/current_variant_information; \
rm $vcf_path/indels/common_modified_vcf/inrange.positions;
done
```
# annotate the potential variant effects
```bash
for mode in unshifted shifted; do
  mkdir -p results_output/$mode/snps/snpEff_vcf
  mkdir -p results_output/$mode/snps/variant_table
  mkdir -p results_output/$mode/indel/snpEff_vcf
  mkdir -p results_output/$mode/indel/variant_table
  mkdir -p results_output/$mode/all_variants/

  cat data/sample_info.txt | while read sample; do
  touch results_output/$mode/snps/variant_table/$sample.summary
  touch results_output/$mode/indel/variant_table/$sample.summary
    if [ "$mode" = "unshifted" ]; then
      java -jar "/home/alan/KX576660/snpEff/snpEff.jar" ann \
      -no-downstream -no-upstream -no-utr -geneId \
      KX576660 \
      variant_calling/$mode/calling_overlap/snp/common_vcf/$sample.common.vcf  > \
      results_output/$mode/snps/snpEff_vcf/$sample.vcf

      java -jar "/home/alan/KX576660/snpEff/snpEff.jar" ann \
      -no-downstream -no-upstream -no-utr -geneId \
      KX576660 \
      variant_calling/$mode/calling_overlap/indels/common_vcf/$sample.common.vcf  > \
      results_output/$mode/indel/snpEff_vcf/$sample.vcf
    else
      vcf-sort variant_calling/$mode/calling_overlap/snp/common_modified_vcf/$sample.annotated.common.vcf > \
      variant_calling/$mode/calling_overlap/snp/common_modified_vcf/$sample.annotated.common.sorted.vcf
      java -jar "/home/alan/KX576660/snpEff/snpEff.jar" ann \
      -no-downstream -no-upstream -no-utr -geneId \
      KX576660 \
      variant_calling/$mode/calling_overlap/snp/common_modified_vcf/$sample.annotated.common.sorted.vcf  > \
      results_output/$mode/snps/snpEff_vcf/$sample.vcf

      vcf-sort variant_calling/$mode/calling_overlap/indels/common_modified_vcf/$sample.annotated.common.vcf > \
      variant_calling/$mode/calling_overlap/indels/common_modified_vcf/$sample.annotated.common.sorted.vcf
      java -jar "/home/alan/KX576660/snpEff/snpEff.jar" ann \
      -no-downstream -no-upstream -no-utr -geneId \
      KX576660 \
      variant_calling/$mode/calling_overlap/indels/common_modified_vcf/$sample.annotated.common.sorted.vcf  > \
      results_output/$mode/indel/snpEff_vcf/$sample.vcf
    fi

    sudo docker run --rm -v /home/alan/KX576660:/home broadinstitute/gatk3:3.8-0 \
    java -jar /usr/GenomeAnalysisTK.jar \
    -R ../home/reference_genome/$mode/KX576660.fasta \
    -T VariantsToTable \
    -V ../home/results_output/$mode/snps/snpEff_vcf/$sample.vcf \
    -F POS -F REF -F ALT -F DP -F AF -F DP4 -F SB -F ANN \
    -o ../home/results_output/$mode/snps/variant_table/$sample.summary

    sudo docker run --rm -v /home/alan/KX576660:/home broadinstitute/gatk3:3.8-0 \
    java -jar /usr/GenomeAnalysisTK.jar \
    -R ../home/reference_genome/$mode/KX576660.fasta \
    -T VariantsToTable \
    -V ../home/results_output/$mode/indel/snpEff_vcf/$sample.vcf \
    -F POS -F REF -F ALT -F DP -F AF -F DP4 -F SB -F ANN \
    -o ../home/results_output/$mode/indel/variant_table/$sample.summary

    grep -v "#" results_output/$mode/snps/snpEff_vcf/$sample.vcf | \
    cat results_output/$mode/indel/snpEff_vcf/$sample.vcf /dev/stdin | \
    vcf-sort > results_output/$mode/all_variants/$sample.vcf

    cat results_output/$mode/snps/variant_table/$sample.summary \
    results_output/$mode/indel/variant_table/$sample.summary > \
    results_output/$mode/all_variants/$sample.summary
  done
done
```

# Combine variants from the shifted and unshifted mapping
```bash

mkdir -p results_output/final_results/vcf
mkdir results_output/final_results/table
cat data/sample_info.txt | while read sample; do
  grep -v "#" results_output/shifted/all_variants/$sample.vcf | \
  cat results_output/unshifted/all_variants/$sample.vcf /dev/stdin | \
  vcf-sort > results_output/final_results/vcf/$sample.vcf

  echo $sample
  grep -v "#" results_output/final_results/vcf/$sample.vcf | wc -l

  cat results_output/unshifted/all_variants/$sample.summary \
      results_output/shifted/all_variants/$sample.summary > \
      results_output/final_results/table/$sample.table
done
```

# Mapping Statistics
```bash
mkdir -p plotting_data/mapping_statistics/mapping_rates
mkdir  -p plotting_data/mapping_statistics/coverage/perbase_coverage/unshifted
mkdir  -p plotting_data/mapping_statistics/coverage/perbase_coverage/shifted


cat data/sample_info.txt | while read sample; do
  samtools depth variant_calling/unshifted/base_recalibration/bam_files/$sample.bam | \
  awk '{sum+=$3} END { print "'$sample'",sum/NR}' >> plotting_data/mapping_statistics/coverage/sample_depth

  samtools depth variant_calling/unshifted/base_recalibration/bam_files/$sample.bam > \
  plotting_data/mapping_statistics/coverage/perbase_coverage/unshifted/$sample.perbase.coverage

  samtools depth variant_calling/shifted/base_recalibration/bam_files/$sample.bam > \
  plotting_data/mapping_statistics/coverage/perbase_coverage/shifted/$sample.perbase.coverage

  bamtools stats -in variant_calling/shifted/base_recalibration/bam_files/$sample.bam > \
  plotting_data/mapping_statistics/mapping_rates/$sample.map.rates

  total_reads=$(bamtools stats -in bwa_mapping/unshifted/"$sample".unfiltered.bam  | grep 'Total reads' | awk '{print $3}')

  mapped_reads=$(bamtools stats -in variant_calling/unshifted/base_recalibration/bam_files/$sample.bam | grep 'Mapped reads' | awk '{print $3}')

  duplicates=$(bamtools stats -in variant_calling/unshifted/base_recalibration/bam_files/$sample.bam  | grep 'Duplicates' | awk '{print $2}')

  echo -e $sample'\t'$total_reads'\t'$mapped_reads'\t'$duplicates >> plotting_data/mapping_statistics/mapping_rates.txt
done
```
