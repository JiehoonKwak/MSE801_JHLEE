#!/bin/bash
set -x
TARGET_DIR="$PWD/raw"

mkdir -p "$TARGET_DIR"

download_if_not_exists() {
  local url=$1
  local output=$2
  if [ ! -f "$output" ]; then
    curl -L "$url" -o "$output"
  else
    echo "File $output already exists. Skipping download."
  fi
}


### Following sections were prepared for the demo ###

# 1. Download files in raw/ from GitHub directory
GITHUB_RAW_URL="https://raw.githubusercontent.com/JiehoonKwak/MSE801_JHLEE/main/raw"
download_if_not_exists "$GITHUB_RAW_URL/subsampled_SRR7138440_WES_of_homo_sapiens_blood_of_brain_tumor_patient_1.fastq.gz" "$TARGET_DIR/subsampled_SRR7138440_WES_of_homo_sapiens_blood_of_brain_tumor_patient_1.fastq.gz"
download_if_not_exists "$GITHUB_RAW_URL/subsampled_SRR7138440_WES_of_homo_sapiens_blood_of_brain_tumor_patient_2.fastq.gz" "$TARGET_DIR/subsampled_SRR7138440_WES_of_homo_sapiens_blood_of_brain_tumor_patient_2.fastq.gz"
download_if_not_exists "$GITHUB_RAW_URL/subsampled_SRR7138435_WES_of_homo_sapiens_tumor_of_brain_tumor_patient_1.fastq.gz" "$TARGET_DIR/subsampled_SRR7138435_WES_of_homo_sapiens_tumor_of_brain_tumor_patient_1.fastq.gz"
download_if_not_exists "$GITHUB_RAW_URL/subsampled_SRR7138435_WES_of_homo_sapiens_tumor_of_brain_tumor_patient_2.fastq.gz" "$TARGET_DIR/subsampled_SRR7138435_WES_of_homo_sapiens_tumor_of_brain_tumor_patient_2.fastq.gz"
download_if_not_exists "$GITHUB_RAW_URL/subsampled_SRR7138438_WES_of_homo_sapiens_subventricular_zone_of_brain_tumor_patient_1.fastq.gz" "$TARGET_DIR/subsampled_SRR7138438_WES_of_homo_sapiens_subventricular_zone_of_brain_tumor_patient_1.fastq.gz"
download_if_not_exists "$GITHUB_RAW_URL/subsampled_SRR7138438_WES_of_homo_sapiens_subventricular_zone_of_brain_tumor_patient_2.fastq.gz" "$TARGET_DIR/subsampled_SRR7138438_WES_of_homo_sapiens_subventricular_zone_of_brain_tumor_patient_2.fastq.gz"

# 2. Download sequence of Chr5 
# download_if_not_exists "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr5.fa.gz" "$TARGET_DIR/hg38.fa.gz"


# 3. Download Base Quality Score Recalibration (BQSR) files
# download_if_not_exists "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf" "$TARGET_DIR/Homo_sapiens_assembly38.dbsnp138.vcf"
# download_if_not_exists "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx" "$TARGET_DIR/Homo_sapiens_assembly38.dbsnp138.vcf.idx"

# 4. Download Mutect2 resources
# download_if_not_exists "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz" "$TARGET_DIR/af-only-gnomad.hg38.vcf.gz"
# download_if_not_exists "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi" "$TARGET_DIR/af-only-gnomad.hg38.vcf.gz.tbi"
# download_if_not_exists "https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz" "$TARGET_DIR/1000g_pon.hg38.vcf.gz"
# download_if_not_exists "https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi" "$TARGET_DIR/1000g_pon.hg38.vcf.gz.tbi"
# download_if_not_exists "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/exome_calling_regions.v1.1.interval_list" "$TARGET_DIR/exome_calling_regions.v1.1.interval_list"



echo "All files have been downloaded to $TARGET_DIR"