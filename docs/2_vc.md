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
