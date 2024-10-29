#!/bin/bash
# FAIL : subsetting fastq from bam not works...

# process_bam() {
#     bam_file="$1"
#     base_name=$(basename "$bam_file" .bam)

#     case "$base_name" in
#         "Blood")
#             srr_id="SRR7138440"
#             sample_desc="blood_of_brain_tumor"
#             ;;
#         "Svz")
#             srr_id="SRR7138438"
#             sample_desc="subventricular_zone_of_brain_tumor"
#             ;;
#         "Tumor")
#             srr_id="SRR7138435"
#             sample_desc="tumor_of_brain_tumor"
#             ;;
#         *)
#             echo "Unknown sample type: $bam_file"
#             return
#             ;;
#     esac

#     samtools view -b "$bam_file" "chr5" > "raw/${base_name}_chr5.bam"

#     samtools sort -o "raw/${base_name}_chr5_sorted.bam" "raw/${base_name}_chr5.bam"
#     samtools index "raw/${base_name}_chr5_sorted.bam"

#     samtools fastq -1 "raw/${base_name}_R1.fastq" \
#                    -2 "raw/${base_name}_R2.fastq" \
#                    -0 /dev/null -s /dev/null -n "raw/${base_name}_chr5_sorted.bam"

#     gzip "raw/${base_name}_R1.fastq"
#     gzip "raw/${base_name}_R2.fastq"

#     # Rename to mach the original file names
#     mv "raw/${base_name}_R1.fastq.gz" "raw/subsampled_${srr_id}_WES_of_homo_sapiens_${sample_desc}_patient_1.fastq.gz"
#     mv "raw/${base_name}_R2.fastq.gz" "raw/subsampled_${srr_id}_WES_of_homo_sapiens_${sample_desc}_patient_2.fastq.gz"

#     # Remove intermediate files
#     rm "raw/${base_name}_chr5.bam" "raw/${base_name}_chr5_sorted.bam" "raw/${base_name}_chr5_sorted.bam.bai"
# }

# export -f process_bam
# parallel -j 3 process_bam ::: raw/*.bam