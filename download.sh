#!/bin/bash
TARGET_DIR="$HOME/demo/ref"

TARGET_DIR="$HOME/demo/ref"

if [ ! -d "$TARGET_DIR" ]; then
  mkdir -p "$TARGET_DIR"
fi

download_if_not_exists() {
  local url=$1
  local output=$2
  if [ ! -f "$output" ]; then
    curl -L "$url" -o "$output"
  else
    echo "File $output already exists. Skipping download."
  fi
}

echo "Downloading files to $TARGET_DIR, it may take more than 3hours..."

# Download sequence reads
download_if_not_exists "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR713/000/SRR7138440/SRR7138440_1.fastq.gz" "$TARGET_DIR/SRR7138440_WES_of_homo_sapiens_blood_of_brain_tumor_patient_1.fastq.gz"
download_if_not_exists "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR713/000/SRR7138440/SRR7138440_2.fastq.gz" "$TARGET_DIR/SRR7138440_WES_of_homo_sapiens_blood_of_brain_tumor_patient_2.fastq.gz"
download_if_not_exists "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR713/009/SRR7138439/SRR7138439_1.fastq.gz" "$TARGET_DIR/SRR7138439_WES_of_homo_sapiens_tumor_of_brain_tumor_patient_1.fastq.gz"
download_if_not_exists "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR713/009/SRR7138439/SRR7138439_2.fastq.gz" "$TARGET_DIR/SRR7138439_WES_of_homo_sapiens_tumor_of_brain_tumor_patient_2.fastq.gz"
download_if_not_exists "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR713/008/SRR7138438/SRR7138438_1.fastq.gz" "$TARGET_DIR/SRR7138438_WES_of_homo_sapiens_subventricular_zone_of_brain_tumor_patient_1.fastq.gz"
download_if_not_exists "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR713/008/SRR7138438/SRR7138438_2.fastq.gz" "$TARGET_DIR/SRR7138438_WES_of_homo_sapiens_subventricular_zone_of_brain_tumor_patient_2.fastq.gz"

# Download reference genome
download_if_not_exists "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz" "$TARGET_DIR/hg38.fa.gz"

# Download Base Quality Score Recalibration (BQSR) files
download_if_not_exists "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf" "$TARGET_DIR/Homo_sapiens_assembly38.dbsnp138.vcf"
download_if_not_exists "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx" "$TARGET_DIR/Homo_sapiens_assembly38.dbsnp138.vcf.idx"

# Download Mutect2 resources
download_if_not_exists "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz" "$TARGET_DIR/af-only-gnomad.hg38.vcf.gz"
download_if_not_exists "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi" "$TARGET_DIR/af-only-gnomad.hg38.vcf.gz.tbi"
download_if_not_exists "https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz" "$TARGET_DIR/1000g_pon.hg38.vcf.gz"
download_if_not_exists "https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi" "$TARGET_DIR/1000g_pon.hg38.vcf.gz.tbi"
download_if_not_exists "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/exome_calling_regions.v1.1.interval_list" "$TARGET_DIR/exome_calling_regions.v1.1.interval_list"

echo "All files have been downloaded to $TARGET_DIR"