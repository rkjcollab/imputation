#!/bin/bash

set -e
set -u

while getopts p:d: opt; do
   case "${opt}" in
      p) plink_prefix=${OPTARG};;
      d) dir=${OPTARG};;
      \?) echo "Invalid option -$OPTARG" >&2
      exit 1;;
   esac
done


# Make file_list for input into bcftools concat
ls ${dir}/*_clean.vcf.gz | sort -V > "${dir}/imputed_cleaned_file_list.txt"

# Concatenate VCFs
bcftools concat \
   --file-list "${dir}/imputed_cleaned_file_list.txt" \
   --output-type z \
   --output "${dir}/chr_all_concat.vcf.gz"

# VCF sample IDs are FID_IID as written by PLINK 1.9 --recode vcf. Make
# temp PSAM with correct IDs instead of using delim which errors when
# IDs have internal delim
awk 'BEGIN{OFS="\t"; print "#IID","SEX"} {print $1"_"$2, $5}' \
   "${plink_prefix}.fam" > "${dir}/tmp_concat.psam"
awk 'BEGIN{OFS="\t"} {print "0", $1"_"$2, $1, $2}' \
    "${plink_prefix}.fam" > "${dir}/tmp_concat_ids.txt"

# Convert to PLINK, adding sex info from pre_qc files
plink2 --vcf "${dir}/chr_all_concat.vcf.gz" \
   --make-pgen --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 10000 \
   --psam "${dir}/tmp_concat.psam" \
   --update-ids "${dir}/tmp_concat_ids.txt" \
   --out "${dir}/chr_all_concat"

# Clean up
rm -f ${dir}/tmp_*
