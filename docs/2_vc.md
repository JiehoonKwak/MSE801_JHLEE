## Part 2. Variant Calling & Downstream Analysis
Use Analysis-ready bam file from Part 1  
GATK : [Somatic short variant discovery(CNVs)](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels)  

- since we subsetted original fastq files and only aligned to Chr5, the following codes wold not work. So, we should parse vcf and interval files with following codes (I already created this files in `ref` directory,it's not necessary for real-world data)
```bash
# bcftools view -r chr5 ref/af-only-gnomad.hg38.vcf.gz -Oz -o ref/af-only-gnomad.chr5.vcf.gz
# tabix -p vcf ref/af-only-gnomad.chr5.vcf.gz  
# grep -P "^@|^chr5\s" ref/exome_calling_regions.v1.1.interval_list > ref/exome_calling_regions.chr5.interval_list

# bcftools view -r chr5 ref/af-only-gnomad.hg38.vcf.gz -Oz -o ref/af-only-gnomad.chr5.hg38.vcf.gz
# tabix -p vcf ref/af-only-gnomad.chr5.hg38.vcf.gz

# bcftools view -r chr5 ref/1000g_pon.hg38.vcf.gz -Oz -o ref/1000g_pon.chr5.hg38.vcf.gz
# tabix -p vcf ref/1000g_pon.chr5.hg38.vcf.gz
```
Use `af-only-gnomad.hg38.vcf.gz` and `exome_calling_regions.v1.1.interval_list` in real-world data

### 1. Call candidate variants
Call variants using `GATK4::Mutect2` (This takes **VERY long** time)
- `Matched Noraml` : normal samples from same individual, to eliminate germline variants
- `Panel of Normal` : unrelated "normal" samples, to eliminate common/recurring technical artifacts
- remove `-L` parameter for real-world data
```bash
# run for Tumor
# real	5m48.879s + 5m57.886s
gatk Mutect2 -R ref/hg38.fa \
    -I out/Tumor_bqsr.bam \
    -I out/Blood_bqsr.bam \
    -normal Blood \
    -O out/Tumor.vcf.gz \
    -L chr5 \
    --germline-resource ref/af-only-gnomad.hg38.vcf.gz \
    --panel-of-normals ref/1000g_pon.hg38.vcf.gz \
    --f1r2-tar-gz out/Tumor_f1r2.tar.gz \
    --native-pair-hmm-threads 4

# run for SVZ
gatk Mutect2 -R ref/hg38.fa \
    -I out/Svz_bqsr.bam \
    -I out/Blood_bqsr.bam \
    -normal Blood \
    -O out/Svz.vcf.gz \
    -L chr5 \
    --germline-resource ref/af-only-gnomad.hg38.vcf.gz \
    --panel-of-normals ref/1000g_pon.hg38.vcf.gz \
    --f1r2-tar-gz out/Svz_f1r2.tar.gz \
    --native-pair-hmm-threads 4

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

```bash
1. First get pileup summaries  
```bash
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
# real 1m20.753s
for sample in Blood Tumor Svz; do
  gatk GetPileupSummaries -I out/${sample}_bqsr.bam -V ref/af-only-gnomad.hg38.vcf.gz -L ref/exome_calling_regions.v1.1.interval_list -O out/${sample}_pileups.table
done
```
  
2. Then calculate contamination
calculate contamination
```bash
# real	0m12.980s
gatk CalculateContamination -I out/Tumor_pileups.table -matched out/Blood_pileups.table -O out/Tumor_contamination.table
gatk CalculateContamination -I out/Svz_pileups.table -matched out/Blood_pileups.table -O out/Svz_contamination.table
```
### 3. Learn Orientation Bias Artifacts
Trim read-order artifacts, using f1r2.tar.gz files
```bash
# real	0m35.544s
gatk LearnReadOrientationModel -I out/Tumor_f1r2.tar.gz -O out/Tumor_read-orientation-model.tar.gz
gatk LearnReadOrientationModel -I out/Svz_f1r2.tar.gz -O out/Svz_read-orientation-model.tar.gz
```

### 4. Filter Variants
Filter based on contamination and orientation bias
```bash
# real	0m13.377s
gatk FilterMutectCalls -V out/Tumor.vcf.gz -R ref/hg38.fa --contamination-table out/Tumor_contamination.table --ob-priors out/Tumor_read-orientation-model.tar.gz --stats out/Tumor.vcf.gz.stats -O out/Tumor_filtered.vcf

gatk FilterMutectCalls -V out/Svz.vcf.gz -R ref/hg38.fa --contamination-table out/Svz_contamination.table --ob-priors out/Svz_read-orientation-model.tar.gz --stats out/Svz.vcf.gz.stats -O out/Svz_filtered.vcf

# or you can do this way
for sample in Tumor Svz; do
  gatk FilterMutectCalls -V out/${sample}.vcf.gz -R ref/hg38.fa --contamination-table out/${sample}_contamination.table --ob-priors out/${sample}_read-orientation-model.tar.gz --stats out/${sample}.vcf.gz.stats -O out/${sample}_filtered.vcf
done
```

### 5. Annotate variants
you can create output of either `maf` or `vcf` format
```bash
# real	0m21.662s (MAF funcotation)
# VCF format
gatk Funcotator -R ref/hg38.fa --ref-version hg38 --data-sources-path ref/funcotator_dataSources.v1.8.hg38.20230908s -V out/Tumor_filtered.vcf -O out/Tumor_funcotated.vcf --output-file-format VCF 
gatk Funcotator -R ref/hg38.fa --ref-version hg38 --data-sources-path ref/funcotator_dataSources.v1.8.hg38.20230908s -V out/Svz_filtered.vcf -O out/Svz_funcotated.vcf --output-file-format VCF

# or you can do this way
for sample in Tumor Svz; do
    gatk Funcotator \
        -V out/${sample}_filtered.vcf \
        -R ref/hg38.fa \
        --ref-version hg38 \
        --data-sources-path ref/funcotator_dataSources.v1.8.hg38.20230908s \
        -O out/${sample}_funcotated.vcf \
        --output-file-format VCF
done

# MAF format
gatk Funcotator -R ref/hg38.fa --ref-version hg38 --data-sources-path ref/funcotator_dataSources.v1.8.hg38.20230908s -V out/Tumor_filtered.vcf -O out/Tumor_funcotated.maf --output-file-format MAF
gatk Funcotator -R ref/hg38.fa --ref-version hg38 --data-sources-path ref/funcotator_dataSources.v1.8.hg38.20230908s -V out/Svz_filtered.vcf -O out/Svz_funcotated.maf --output-file-format MAF

# or you can do this way
for sample in Tumor Svz; do
    gatk Funcotator \
        -V out/${sample}_filtered.vcf \
        -R ref/hg38.fa \
        --ref-version hg38 \
        --data-sources-path ref/funcotator_dataSources.v1.8.hg38.20230908s \
        -O out/${sample}_funcotated.maf \
        --output-file-format MAF
done
```

### (Optional) 6. Convert to Table format
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
        -F CHROM -F POS -F TYPE -F REF -F ALT \
        -F AC -F AN -F DP -F AF -GF DP -GF AD -F FUNCOTATION \
        -O out/${sample}_output_snps.table
done
```
