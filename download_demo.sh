#!/bin/bash
set -x
TARGET_DIR="$PWD/raw"
REPO="JiehoonKwak/MSE801_JHLEE"

mkdir -p "$TARGET_DIR"

# Check if gh CLI is installed
if ! command -v gh &> /dev/null
then
    echo "gh CLI could not be found, re-check Step 1."
    exit 1
fi

if [ ! -d "$REPO" ]; then
  gh repo clone "$REPO" -- --filter=blob:none --sparse
  cd MSE801_JHLEE
else
  cd MSE801_JHLEE
  git pull
fi

git sparse-checkout set raw
git lfs pull

cp raw/* "$TARGET_DIR/"
echo "Files have been downloaded to $TARGET_DIR"


### Following sections were prepared for the demo ###

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