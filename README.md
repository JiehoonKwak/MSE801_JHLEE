# 2024 MSE801 - Demo for WGS

# Pre-requisite
- download the data from the link provided in [HERE](docs/setup.md)
- install packages from `requirements.txt`
```bash
conda install --file requirements.txt
```

# Case 1 : WGS
## 1. Process Analyze-ready bam file
### 1. Create index files
- create `.fai` file
```bash
gunzip hg38.fa.gz
samtools faidx hg38.fa
```
- create `.dict` file
```bash
gatk CreateSequenceDictionary R=hg38.fa O=hg38.dict
```

### 2. Quality Control
Perform quality control of fastq files with `fastqc`
```bash
mkdir qc
fastqc -t 16 -o qc N724.fastq.gz
```

### 3. Alignment to reference genome
map to reference using `BWA-MEM`
```bash
bwa mem -t 16 -R "@RG\tID:sample\tPL:ILLUMINA\tSM:sample" hg38 read1 read2 > paired.sam
```
\
Then, check output using `samtools`
```bash
samtools flagstat
```

### 4. Mark duplicates
Mark duplicates (Deduplicate) and Sort using `GATK4`
```bash
gatk MarkDuplicatesSpark -I paired.sam -O out.bam
```

### 5. Base quality recalibration
1. Build Model
```bash
gatk BaseRecalibrator -I out.bam -R hg38 --known-sites vcf -O recal_data.table
```

2. Adjust base quality scores
```bash
gatk ApplyBQSR -I out.bam -R hg38 --bqsr-recal-file recal_data.table -O recal.bam
```

### 6. Collect alignment metrics
Collect alignment metrics and Insert size metrics using `GATK4`
```bash
gatk CollectAlignmentSummaryMetrics R=hg38 I=recal.bam O=alignment_metrics.txt
gatk CollectInsertSizeMetrics I=recal.bam O=insert_size_metrics.txt H=insert_size_metrics.txt HISTOGRAM_FILE=insert_size_histogram.pdf
```

## 2. Variant Calling
### 1. HaplotypeCaller
Call variants using `GATK4` (This takes long time)
```bash
gatk HaplotypeCaller -R hg38 -I recal.bam -O raw_variants.vcf
```
\
extract SNP & indels seperately
```bash

```

## Case 2 : Amplicon/Panel Sequencing
- alingment to reference genome
- `dada2`

