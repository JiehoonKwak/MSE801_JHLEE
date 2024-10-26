## Part 2. Variant Calling & Downstream Analysis
Use Analysis-ready bam file from Part 1  
GATK : [Somatic short variant discovery(CNVs)](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels)

### 1. Call candidate variants
Call variants using `GATK4::Mutect2` (This takes long time)
- `Matched Noraml` : normal samples from same individual, to eliminate germline variants
- `Panel of Normal` : unrelated "normal" samples, to eliminate common/recurring technical artifacts
```bash
# run for Tumor
gatk Mutect2 -R ref/hg38.fa \
    -I out/Tumor_bqsr.bam \
    -I out/Blood_bqsr.bam \
    -normal Blood \
    -O out/Tumor.vcf.gz \
    --germline-resource ref/af-only-gnomad.hg38.vcf.gz \
    --panel-of-normals ref/1000g_pon.hg38.vcf.gz \
    --f1r2-tar-gz out/Tumor_f1r2.tar.gz \ 
    --native-pair-hmm-threads 36

# run for SVZ
gatk Mutect2 -R ref/hg38.fa \
    -I out/Svz_bqsr.bam \
    -I out/Blood_bqsr.bam \
    -normal Blood \
    -O out/Svz.vcf.gz \
    --germline-resource ref/af-only-gnomad.hg38.vcf.gz \
    --panel-of-normals ref/1000g_pon.hg38.vcf.gz \
    --f1r2-tar-gz out/Svz_f1r2.tar.gz \
    --native-pair-hmm-threads 36

# cf, you can run multiple samples at once
# gatk Mutect2 -R ref/hg38.fa \
#     -I out/Tumor_bqsr.bam \
#     -I out/Svz_bqsr.bam \
#     -I out/Blood_bqsr.bam \
#     -normal Blood 
#     -O out/somatic.vcf.gz \
# --germline-resource ref/af-only-gnomad.hg38.vcf.gz --panel-of-normals ref/1000g_pon.hg38.vcf.gz \
# --f1r2-tar-gz out/f1r2.tar.gz
```
  
You can find sample name in the header of the bam file using `samtools`
```bash
samtools view -H out/Blood_bqsr.bam | grep '@RG'
```

### 2. Calculate Contamination
calculate cross-sample contamination using `GATK4::CalculateContamination`  
1. First get pileup summaries
- normal tissue
```bash
gatk GetPileupSummaries -I out/Blood_bqsr.bam -V ref/af-only-gnomad.hg38.vcf.gz -L ref/exome_calling_regions.v1.1.interval_list -O out/Blood_pileups.table
```

- tumor tissue
```bash
gatk GetPileupSummaries -I out/Tumor_bqsr.bam -V ref/af-only-gnomad.hg38.vcf.gz -L ref/exome_calling_regions.v1.1.interval_list -O out/Tumor_pileups.table

gatk GetPileupSummaries -I out/Svz_bqsr.bam -V ref/af-only-gnomad.hg38.vcf.gz -L ref/exome_calling_regions.v1.1.interval_list -O out/Svz_pileups.table
```

for loops
```bash
for sample in Blood Tumor Svz; do
  gatk GetPileupSummaries -I out/${sample}_bqsr.bam -V ref/af-only-gnomad.hg38.vcf.gz -L ref/exome_calling_regions.v1.1.interval_list -O out/${sample}_pileups.table
done
```
  
2. Then calculate contamination
calculate contamination
```bash
gatk CalculateContamination -I out/Tumor_pileups.table -matched out/Blood_pileups.table -O out/Tumor_contamination.table
gatk CalculateContamination -I out/Svz_pileups.table -matched out/Blood_pileups.table -O out/Svz_contamination.table
```

### 3. Learn Orientation Bias Artifacts
Trim read-order artifacts, using f1r2.tar.gz files
```bash
gatk LearnReadOrientationModel -I out/Tumor_f1r2.tar.gz -O out/Tumor_read-orientation-model.tar.gz
gatk LearnReadOrientationModel -I out/Svz_f1r2.tar.gz -O out/Svz_read-orientation-model.tar.gz
```

### 4. Filter Variants
Filter based on contamination and orientation bias
```bash
gatk FilterMutectCalls -V out/Tumor.vcf.gz -R ref/hg38.fa --contamination-table out/Tumor_contamination.table --ob-priors out/Tumor_read-orientation-model.tar.gz -O out/Tumor_filtered.vcf

gatk FilterMutectCalls -V out/Svz.vcf.gz -R ref/hg38.fa --contamination-table out/Svz_contamination.table --ob-priors out/Svz_read-orientation-model.tar.gz -O out/Svz_filtered.vcf

# or you can do this way
for sample in Tumor Svz; do
  gatk FilterMutectCalls -V out/${sample}.vcf.gz -R ref/hg38.fa --contamination-table out/${sample}_contamination.table --ob-priors out/${sample}_read-orientation-model.tar.gz -O out/${sample}_filtered.vcf
done
```

### 5. Annotate variants
```bash
gatk Funcotator -R ref/hg38.fa -V out/Tumor_filtered.vcf -O out/Tumor_funcotated.vcf
gatk Funcotator -R ref/hg38.fa -V out/Svz_filtered.vcf -O out/Svz_funcotated.vcf

# or you can do this way
for sample in Tumor Svz; do
  gatk Funcotator -R ref/hg38.fa -V out/${sample}_filtered.vcf -O out/${sample}_funcotated.vcf
done

for sample in Tumor Svz; do
    gatk Funcotator \
        -V out/${sample}_filtered.vcf \
        -R ref/hg38.fa \
        --ref-version hg38 \
        --data-sources-path ref/funcotator_dataSources.v1.8.hg38.20230908s \
        -O out/${sample}_funcotated.vcf \
        --output-file-format VCF
done
```

### 6. Convert to Table format
```bash
gatk VariantsToTable \
    -V out/Tumor_funcotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
    -O out/Tumor.table

gatk VariantsToTable \
    -V out/Svz_funcotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
    -O out/Svz.table

# or you can do this way
for sample in Tumor Svz; do
    gatk VariantsToTable \
        -V out/${sample}_funcotated.vcf \
        -F CHROM -F POS -F TYPE -F AC -F AN -F DP -F AF -F FUNCOTATION \
        -O out/${sample}_output_snps.table
done
```