# 2024 MSE801 - Demo for WGS
This is real-world WGS whole workflow, for 2024 MSE801 KAIST  
- Reference : [GATK4 Best Practices](https://gatk.broadinstitute.org/hc/en-us), [Biostar Handbook](https://www.biostarhandbook.com/)  


## Before we begin!
- download the data from the link provided in [HERE](docs/setup.md)
- install packages from `requirements.txt` in your own, new conda environment
```bash
conda create -n YOUR_ENV -y
conda activate YOUR_ENV
conda install --file requirements.txt
```

# Case 1 : WGS
## Part 1. Process Analyze-ready bam file
GATK4 : [Data pre-processing for variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery)

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

- create `.bwt` file
It takes about 40minutes
```bash
bwa index hg38.fa
```

### 2. Quality Control
Perform quality control of fastq files with `fastqc`
- `-t` : number of threads
```bash
mkdir qc

for file in raw/*.fastq.gz; do 
    fastqc -t 4 -o qc "$file"
done
```

Inspect the quality control results with `multiqc`
```bash
multiqc qc -o . -n raw_multiqc_report.html

# host http server to view the report
python -m http.server 7222
https://<server.ip>:7222/raw_multiqc_report.html

# or, open the file directly
firefox multiqc_report.html
```

(**Optionally**), preprocess the fastq files with `fastp`. But, in this case we are not going to run this code since QC metrics are already good and bwa supports soft clipping.
```bash
mkdir trimmed
fastp -i raw/SRR7138438_WES_of_homo_sapiens_subventricular_zone_of_brain_tumor_patient_1.fastq.gz -I raw/SRR7138438_WES_of_homo_sapiens_subventricular_zone_of_brain_tumor_patient_2.fastq.gz -o trimmed/SRR7138438_WES_of_homo_sapiens_subventricular_zone_of_brain_tumor_patient_1.fastq.gz -O trimmed/SRR7138438_WES_of_homo_sapiens_subventricular_zone_of_brain_tumor_patient_2.fastq.gz
```
  
(**Optional**) Then, re-run `fastqc` and `multiqc` for the trimmed files
```bash
mkdir reqc
for file in trimmed/*.fastq.gz; do # run this code in demo INSTEAD!!
    fastqc -t 4 -o reqc "$file"
done
multiqc reqc -o . -n trimmed_multiqc_report.html
```

### 3. Alignment to reference genome
**Tip** : run with `tmux` or `screen`  

map to reference using `BWA-MEM`
```bash
mkdir out
bwa mem -Ma -t 16 \
-R "@RG\tID:Blood\tSM:Blood\tPL:ILLUMINA\tLB:ILLUMINA" \
ref/hg38.fa \
raw/SRR7138437_WES_of_homo_sapiens_brain_of_brain_tumor_patient_1.fastq.gz \
raw/SRR7138440_WES_of_homo_sapiens_blood_of_brain_tumor_patient_2.fastq.gz > out/Blood.sam

bwa mem -Ma -t 42 \
-R "@RG\tID:Tumor\tSM:Tumor\tPL:ILLUMINA\tLB:ILLUMINA" \
ref/hg38.fa \
raw/SRR7138437_WES_of_homo_sapiens_brain_of_brain_tumor_patient_1.fastq.gz \
raw/SRR7138437_WES_of_homo_sapiens_brain_of_brain_tumor_patient_2.fastq.gz > out/Tumor.sam

# You can use for-loops
for R1 in raw/*_1.fastq.gz; do
    R2=${R1/_1.fastq.gz/_2.fastq.gz}
    SAMPLE_NAME=$(basename "$R1" _1.fastq.gz)
    OPTIONS="@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tPL:ILLUMINA\tLB:ILLUMINA"

    bwa mem -Ma -t 40 -R "$OPTIONS" ref/hg38.fa "$R1" "$R2" > "out/${SAMPLE_NAME}.sam"
done
```
  
Then, check output using `samtools`
```bash
samtools view Blood.sam | less
samtools flagstat Blood.sam
```

### 4. Mark duplicates
Mark duplicates (Deduplicate) and Sorting alignment file using `GATK4:: MarkDuplicatesSpark`

```bash
gatk MarkDuplicatesSpark -I out/Blood.sam -O out/Blood.bam
gatk MarkDuplicatesSpark -I out/Tumor.sam -O out/Tumor.bam
gatk MarkDuplicatesSpark -I out/SVZ.sam -O out/SVZ.bam
```

Then, again check with `samtools` as above

### 5. Base quality recalibration
adjust base quality based on machine learning, using known variants
1. Build Model
```bash
gatk BaseRecalibrator -I out/Blood.bam -R ref/hg38.fa --known-sites ref/Homo_sapiens_assembly38.dbsnp138.vcf -O out/recal_data_blood.table
```

2. Adjust base quality scores based on model
```bash
gatk ApplyBQSR -I out/Blood.bam -R ref/hg38.fa --bqsr-recal-file out/recal_data_blood.table -O out/Blood_bqsr.bam
```
  
Now, analysis-ready bam file is ready for downstream analysis.  
But before we begin, let's check some metrices.  

### 6. Collect alignment metrics
Collect alignment metrics and Insert size metrics using `GATK4`
```bash
gatk CollectAlignmentSummaryMetrics R=ref/hg38.fa I=out/Blood_bqsr.bam O=out/alignment_metrics_blood.txt
gatk CollectInsertSizeMetrics I=out/Blood_bqsr.bam O=out/insert_metrics_blood.txt H=insert_size_metrics.txt HISTOGRAM_FILE=out/insert_histogram_blod.pdf
```

## Part 2. Variant Calling & Downstream Analysis
Use Analysis-ready bam file from Part 1  
GATK : [Somatic short variant discovery(CNVs)](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels)

### 1. Call candidate variants
Call variants using `GATK4::Mutect2` (This takes long time)
```bash
gatk Mutec2 -R hg38 -I normal.bam -I tumor.bam -tumor tumor -normal normal -O somatic.vcf \
--germline-resource af-only-gnomad.hg38.vcf.gz --panel-of-normals 1000g_pon.hg38.vcf.gz \
--f1r2-tar-gz f1r2.tar.gz
```
### 2. Calculate Contamination
calculate cross-sample contamination using `GATK4::CalculateContamination`  
1. First get pileup summaries
- normal tissue
```bash
gatk GetPileupSummaries -I recal.bam -V af-only-gnomad.hg38.vcf.gz -L exome_calling_regions.v1.1.interval_list -O pileups.table
```

- tumor tissue
```bash
gatk GetPileupSummaries -I recal.bam -V af-only-gnomad.hg38.vcf.gz -L exome_calling_regions.v1.1.interval_list -O pileups.table
```
  
2. Then calculate contamination
- calculate contamination
```bash
gatk CalculateContamination -I tumor_pileups.table -matched normal_pileups.table -O contamination.table
```

### 3. Learn Orientation Bias Artifacts
Trim read-order artifacts
```bash
gatk LearnReadOrientationModel -I f1r2.tar.gz -O read-orientation-model.tar.gz
```

### 4. Filter Variants
```bash
gatk FilterMutectCalls -V somatic.vcf -R hg38 --contamination-table contamination.table --tumor-segmentation read-orientation-model.tar.gz -O filtered.vcf
```

### 5. Annotate variants
```bash
gatk Funcotator -R hg38 -V filtered.vcf -O funcotated.vcf
```

## Case 2 : Amplicon/Panel Sequencing
- alingment to reference genome
- `dada2`

