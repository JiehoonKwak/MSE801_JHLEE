#!/usr/bin/env bash
set -euo pipefail

WORKING_DIR="$(pwd)"
TARGET_DIR="${1:-$WORKING_DIR/raw}"
REF_DIR="${2:-$WORKING_DIR/ref}"
REPO_URL="https://github.com/JiehoonKwak/MSE801_JHLEE.git"

log() {
  echo "[MSE801] $*"
}

usage() {
  cat <<'EOF'
Usage:
  bash download_demo.sh [RAW_DIR] [REF_DIR]

Arguments:
  RAW_DIR    Path to store raw FASTQ files (default: <cwd>/raw)
  REF_DIR    Path to store reference files (default: <cwd>/ref)
EOF
}

download_if_not_exists() {
  local url=$1
  local out=$2
  local fallback_urls=()
  shift 2
  fallback_urls=("$@")

  mkdir -p "$(dirname "$out")"

  if [[ -f "$out" ]]; then
    log "Skip (exists): $out"
    return
  fi

  for u in "$url" "${fallback_urls[@]}"; do
    local temp_out="${out}.part"
    log "Downloading: $out"
    rm -f "$temp_out"
    if curl -fL --retry 5 --retry-delay 3 --retry-all-errors --max-time 1200 --output "$temp_out" "$u"; then
      mv "$temp_out" "$out"
      return
    fi
    rm -f "$temp_out"
    log "Download failed from $u"
  done

  log "ERROR: cannot download $out. Stop."
  return 1
}

find_raw_commit() {
  local repo_dir=$1
  local commit
  commit="HEAD"

  while true; do
    if git -C "$repo_dir" ls-tree -r --name-only "$commit" | grep -q '^raw/'; then
      echo "$commit"
      return
    fi
    commit=$(git -C "$repo_dir" rev-parse "${commit}^") || {
      return 1
    }
  done
}

prepare_demo_fastq() {
  local tmp_dir=$1
  log "Prepare raw FASTQ from repository history (filtered Chr5 dataset)"

  mkdir -p "$tmp_dir"
  git clone --filter=blob:none --no-checkout "$REPO_URL" "$tmp_dir/repo"
  local raw_commit
  raw_commit=$(find_raw_commit "$tmp_dir/repo")
  if [[ -z "$raw_commit" ]]; then
    log "ERROR: cannot find a commit that contains raw/ files."
    return 1
  fi

  cd "$tmp_dir/repo"
  git sparse-checkout init --cone
  git sparse-checkout set raw
  git checkout "$raw_commit"

  if command -v git-lfs >/dev/null 2>&1; then
    git lfs pull
  fi

  mkdir -p "$TARGET_DIR"
  cp -a raw/* "$TARGET_DIR/"
  cd - >/dev/null
}

prepare_reference_resources() {
  mkdir -p "$REF_DIR"

  log "Download reference genome and Mutect2 resources"
  download_if_not_exists \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr5.fa.gz" \
    "$REF_DIR/hg38.fa.gz"
  if [[ ! -f "$REF_DIR/hg38.fa" ]]; then
    log "Extract reference: $REF_DIR/hg38.fa"
    gzip -dc "$REF_DIR/hg38.fa.gz" > "$REF_DIR/hg38.fa"
  fi

  download_if_not_exists \
    "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz" \
    "$REF_DIR/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
  download_if_not_exists \
    "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi" \
    "$REF_DIR/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi"

  download_if_not_exists \
    "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz" \
    "$REF_DIR/af-only-gnomad.hg38.vcf.gz"
  download_if_not_exists \
    "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi" \
    "$REF_DIR/af-only-gnomad.hg38.vcf.gz.tbi"

  download_if_not_exists \
    "https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz" \
    "$REF_DIR/1000g_pon.hg38.vcf.gz"
  download_if_not_exists \
    "https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi" \
    "$REF_DIR/1000g_pon.hg38.vcf.gz.tbi"

  download_if_not_exists \
    "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/exome_calling_regions.v1.1.interval_list" \
    "$REF_DIR/exome_calling_regions.v1.1.interval_list"
  if [[ ! -f "$REF_DIR/exome_calling_regions.chr5.interval_list" ]]; then
    awk '/^@/ || $1 == "chr5" { print }' \
      "$REF_DIR/exome_calling_regions.v1.1.interval_list" > "$REF_DIR/exome_calling_regions.chr5.interval_list"
  fi

  if command -v bcftools >/dev/null 2>&1 && command -v tabix >/dev/null 2>&1; then
    if [[ ! -f "$REF_DIR/af-only-gnomad.chr5.hg38.vcf.gz" || ! -f "$REF_DIR/af-only-gnomad.chr5.hg38.vcf.gz.tbi" ]]; then
      bcftools view -r chr5 "$REF_DIR/af-only-gnomad.hg38.vcf.gz" -Oz -o "$REF_DIR/af-only-gnomad.chr5.hg38.vcf.gz"
      tabix -p vcf "$REF_DIR/af-only-gnomad.chr5.hg38.vcf.gz"
    fi

    if [[ ! -f "$REF_DIR/1000g_pon.chr5.hg38.vcf.gz" || ! -f "$REF_DIR/1000g_pon.chr5.hg38.vcf.gz.tbi" ]]; then
      bcftools view -r chr5 "$REF_DIR/1000g_pon.hg38.vcf.gz" -Oz -o "$REF_DIR/1000g_pon.chr5.hg38.vcf.gz"
      tabix -p vcf "$REF_DIR/1000g_pon.chr5.hg38.vcf.gz"
    fi
  else
    log "bcftools/tabix not found. chr5-filtered VCF files were not created."
    log "Install bcftools and tabix and rerun script if you need chr5-filtered VCFs."
  fi
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

mkdir -p "$TARGET_DIR" "$REF_DIR"

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

prepare_demo_fastq "$tmpdir"
prepare_reference_resources

log "Done."
log "  FASTQ : $TARGET_DIR"
log "  REF   : $REF_DIR"
