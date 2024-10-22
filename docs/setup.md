# Data Preparation

## Raw sequence reads & reference genome download
1. Download sequence reads using `SRA explorer`
- Original paper : [Human glioblastoma arises from subventricular zone cells with low-level driver mutations](https://www.nature.com/articles/s41586-018-0389-3)
    - SRA ID : SRP145073
- `SRA explorer` : [Link](https://sra-explorer.info/)
```bash
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR713/000/SRR7138440/SRR7138440_1.fastq.gz -o SRR7138440_WES_of_homo_sapiens_blood_of_brain_tumor_patient_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR713/000/SRR7138440/SRR7138440_2.fastq.gz -o SRR7138440_WES_of_homo_sapiens_blood_of_brain_tumor_patient_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR713/009/SRR7138439/SRR7138439_1.fastq.gz -o SRR7138439_WES_of_homo_sapiens_tumor_of_brain_tumor_patient_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR713/009/SRR7138439/SRR7138439_2.fastq.gz -o SRR7138439_WES_of_homo_sapiens_tumor_of_brain_tumor_patient_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR713/008/SRR7138438/SRR7138438_1.fastq.gz -o SRR7138438_WES_of_homo_sapiens_subventricular_zone_of_brain_tumor_patient_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR713/008/SRR7138438/SRR7138438_2.fastq.gz -o SRR7138438_WES_of_homo_sapiens_subventricular_zone_of_brain_tumor_patient_2.fastq.gz
```

2. Download reference genome
- [GRCh38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)
```bash
curl -L https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -o hg38.fa.gz
```

3. Download Base Quality Score Recalibration (BQSR)
- BQSR : [Link](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR)
```bash
curl -L https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf -o Homo_sapiens_assembly38.dbsnp138.vcf
curl -L https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx -o Homo_sapiens_assembly38.dbsnp138.vcf.idx
```

4. Download Mutect2 resources
- gnomeAD
```bash
curl -L https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz -o af-only-gnomad.hg38.vcf.gz
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi -o af-only-gnomad.hg38.vcf.gz.tbi
```
\
- PoN
```bash
curl -L https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz -o 1000g_pon.hg38.vcf.gz
curl -L https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi -o 1000g_pon.hg38.vcf.gz.tbi
```