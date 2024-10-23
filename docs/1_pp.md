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
This takes long time
```bash
bwa index hg38.fa
```

### 2. Quality Control
Perform quality control of fastq files with `fastqc`
- `-t` : number of threads
```bash
mkdir qc

for file in raw/*.fastq.gz; do # run this code in demo INSTEAD!!
    fastqc -t 4 -o qc "$file"
done
```

Inspect the quality control results with `multiqc`
```bash
multiqc qc -o . -n raw_multiqc_report.html

# host http server to view the report
python -m http.server 7222
https://localhost:7222/raw_multiqc_report.html

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
raw/subsampled_SRR7138440_WES_of_homo_sapiens_blood_of_brain_tumor_patient_2.fastq.gz \
raw/subsampled_SRR7138440_WES_of_homo_sapiens_blood_of_brain_tumor_patient_2.fastq.gz > out/Blood.sam

bwa mem -Ma -t 16 \
-R "@RG\tID:Tumor\tSM:Tumor\tPL:ILLUMINA\tLB:ILLUMINA" \
ref/hg38.fa \
raw/subsampled_SRR7138437_WES_of_homo_sapiens_brain_of_brain_tumor_patient_1.fastq.gz \
raw/subsampled_SRR7138437_WES_of_homo_sapiens_brain_of_brain_tumor_patient_2.fastq.gz > out/Tumor.sam

bwa mem -Ma -t 16 \
-R "@RG\tID:Svz\tSM:Svz\tPL:ILLUMINA\tLB:ILLUMINA" \
ref/hg38.fa \
raw/subsampled_SRR7138438_WES_of_homo_sapiens_subventricular_zone_of_brain_tumor_patient_1.fastq.gz \
raw/subsampled_SRR7138438_WES_of_homo_sapiens_subventricular_zone_of_brain_tumor_patient_2.fastq.gz > out/Svz.sam


# You can use for-loops
for R1 in raw/*_1.fastq.gz; do
    R2=${R1/_1.fastq.gz/_2.fastq.gz}
    SAMPLE_NAME=$(basename "$R1" _1.fastq.gz)
    OPTIONS="@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tPL:ILLUMINA\tLB:ILLUMINA"

    bwa mem -Ma -t 16 -R "$OPTIONS" ref/hg38.fa "$R1" "$R2" > "out/${SAMPLE_NAME}.sam"
done
```
  
Let's rename files for make it clean (run if you used for loop)
```bash
cd out

mv subsampled_SRR7138437_WES_of_homo_sapiens_brain_of_brain_tumor_patient.sam Tumor.sam
mv subsampled_SRR7138438_WES_of_homo_sapiens_subventricular_zone_of_brain_tumor_patient.sam Svz.sam
mv subsampled_SRR7138440_WES_of_homo_sapiens_blood_of_brain_tumor_patient.sam Blood.sam
```
  
Then, check output using `samtools`
```bash
samtools view Blood.sam | less
samtools flagstat Blood.sam
```

### 4. Mark duplicates
Mark duplicates (Deduplicate) and Sorting alignment file using `GATK4:: MarkDuplicatesSpark`

```bash
cd ..
gatk MarkDuplicatesSpark -I out/Blood.sam -O out/Blood.bam
gatk MarkDuplicatesSpark -I out/Tumor.sam -O out/Tumor.bam
gatk MarkDuplicatesSpark -I out/Svz.sam -O out/Svz.bam

# or you can do this way
for sample in Blood Tumor Svz; do
  gatk MarkDuplicatesSpark -I out/${sample}.sam -O out/${sample}.bam
done

ls out
```

Then, again check with `samtools` as above

### 5. Base quality recalibration
adjust base quality based on machine learning, using known variants
1. Build Model
```bash
gatk BaseRecalibrator -I out/Blood.bam -R ref/hg38.fa --known-sites ref/Homo_sapiens_assembly38.dbsnp138.vcf -O out/recal_data_blood.table

gatk BaseRecalibrator -I out/Tumor.bam -R ref/hg38.fa --known-sites ref/Homo_sapiens_assembly38.dbsnp138.vcf -O out/recal_data_tumor.table

gatk BaseRecalibrator -I out/Svz.bam -R ref/hg38.fa --known-sites ref/Homo_sapiens_assembly38.dbsnp138.vcf -O out/recal_data_svz.table

# you can do this way
for sample in Blood Tumor Svz; do
  gatk BaseRecalibrator -I out/${sample}.bam -R ref/hg38.fa --known-sites ref/Homo_sapiens_assembly38.dbsnp138.vcf -O out/recal_data_${sample}.table
done

```

2. Adjust base quality scores based on model
```bash
gatk ApplyBQSR -I out/Blood.bam -R ref/hg38.fa --bqsr-recal-file out/recal_data_blood.table -O out/Blood_bqsr.bam

gatk ApplyBQSR -I out/Tumor.bam -R ref/hg38.fa --bqsr-recal-file out/recal_data_tumor.table -O out/Tumor_bqsr.bam

gatk ApplyBQSR -I out/Svz.bam -R ref/hg38.fa --bqsr-recal-file out/recal_data_svz.table -O out/Svz_bqsr.bam

# you can do this way
for sample in Blood Tumor Svz; do
  gatk ApplyBQSR -I out/${sample}.bam -R ref/hg38.fa --bqsr-recal-file out/recal_data_${sample}.table -O out/${sample}_bqsr.bam
done
```
  
Now, analysis-ready bam file is ready for downstream analysis.  
But before we begin, let's check some metrices.  

### 6. (Optional) Collect alignment metrics
Collect alignment metrics and Insert size metrics using `GATK4` (This takes long)
```bash
gatk CollectAlignmentSummaryMetrics -R ref/hg38.fa -I out/Blood_bqsr.bam -O out/alignment_metrics_blood.txt
gatk CollectAlignmentSummaryMetrics -R ref/hg38.fa -I out/Tumor_bqsr.bam -O out/alignment_metrics_tumor.txt
gatk CollectAlignmentSummaryMetrics -R ref/hg38.fa -I out/Svz_bqsr.bam -O out/alignment_metrics_svz.txt

gatk CollectInsertSizeMetrics -I out/Blood_bqsr.bam -O out/insert_metrics_blood.txt -H out/insert_histogram_blod.pdf
gatk CollectInsertSizeMetrics -I out/Tumor_bqsr.bam -O out/insert_metrics_tumor.txt -H out/insert_histogram_tumor.pdf
gatk CollectInsertSizeMetrics -I out/Svz_bqsr.bam -O out/insert_metrics_svz.txt -H out/insert_histogram_svz.pdf

# you can do this
for sample in Blood Tumor Svz; do
  gatk CollectAlignmentSummaryMetrics -R ref/hg38.fa -I out/${sample}_bqsr.bam -O out/alignment_metrics_${sample}.txt
  
  gatk CollectInsertSizeMetrics -I out/${sample}_bqsr.bam -O out/insert_metrics_${sample}.txt -H out/insert_histogram_${sample}.pdf
done
```
