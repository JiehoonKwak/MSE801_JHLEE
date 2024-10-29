# Data Preparation
- **This is data prepartion process for demo workflow. If you want to follow whole workflow with real-world dataset, refer to [this](whole_setup.md)**
- **TL;DR** : run following code
```bash
bash <(curl -s https://raw.githubusercontent.com/JiehoonKwak/MSE801_JHLEE/main/download_demo.sh)
```

## How datas were created
- original fastq files were too large, so i tired to reverse-create fastq.gz files that were only aligned to Chr5 : [subsample.sh](subsample.sh), but it leads to decreased depth and quality of the data. So, I decided to random sample 1% of the original fastq files.
```bash
# random subsampling do not effectively call variants
ls *.fastq.gz | parallel -j 6 'seqtk sample -s777 {} 0.01 | gzip > sample/subsampled_{}'
```

1. Download sequence of Chr5 (Reference genome)
- [GRCh38's Chr5](https://hgrdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr6.fa.gz)
```bash
curl -L https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr5.fa.gz -o hg38.fa.gz
```

2. Download Base Quality Score Recalibration (BQSR)
- BQSR : [Link](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR)
```bash
curl -L https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf -o Homo_sapiens_assembly38.dbsnp138.vcf
curl -L https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx -o Homo_sapiens_assembly38.dbsnp138.vcf.idx
```

3. Download Mutect2 resources
- gnomeAD
```bash
curl -L https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz -o af-only-gnomad.hg38.vcf.gz
curl -L https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi -o af-only-gnomad.hg38.vcf.gz.tbi
```
 
- PoN
```bash
curl -L https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz -o 1000g_pon.hg38.vcf.gz
curl -L https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi -o 1000g_pon.hg38.vcf.gz.tbi
```
- Interval files
```bash
curl -L https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/exome_calling_regions.v1.1.interval_list -o exome_calling_regions.v1.1.interval_list
```

- Funcotator source file
```bash
gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download --hg38
```