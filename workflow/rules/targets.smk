rule all:
    input:
        f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/chr_all_concat.pvar",
        f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/chr_all_concat.psam",
        f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/chr_all_concat.pgen"
