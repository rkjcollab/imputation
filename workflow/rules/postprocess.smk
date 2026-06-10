rule download_results:
    threads: 16
    output:
        [f"{out_dir}/imputed/chr_{c}.zip" for c in chr_noY]
    log:
        f"{out_dir}/imputed/download_results.log"
    params:
        script=Path(code_dir, "scripts/download_results.sh")
    shell:
        """
        bash {params.script} \
            -i {imp} \
            -c {code_dir} \
            -o {out_dir}/imputed \
            -j {imp_job_id} \
            > {log} 2>&1
        """


# Different from the other rules, this script in this rule runs once for each chr
rule unzip_results:
    input:
        [f"{out_dir}/imputed/chr{c}.dose.vcf.gz" for c in chr_noY]


rule unzip_results_helper:
    input:
        f"{out_dir}/imputed/chr_{{chr}}.zip"
    output:
        f"{out_dir}/imputed/chr{{chr}}.dose.vcf.gz"
    log:
        f"{out_dir}/imputed/chr{{chr}}_unzip_results.log"
    params:
        script=Path(code_dir, "scripts/unzip_results.sh")
    shell:
        """
        bash {params.script} \
            -d {out_dir}/imputed \
            -p "{zip_pw}" \
            -c {wildcards.chr} \
            > {log} 2>&1
        """


# Different from the other rules, this script in this rule runs once for each chr
rule filter_info_and_vcf_files:
    input:
        [f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/chr{c}_clean.vcf.gz" for c in chr_noY]


rule filter_info_and_vcf_files_helper:
    input:
        f"{out_dir}/imputed/chr{{chr}}.dose.vcf.gz"
    output:
        f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/chr{{chr}}_clean.vcf.gz"
    log:
        f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/chr{{chr}}_filter_info_and_vcf_files.log"
    params:
        script=f'{code_dir}/scripts/filter_info_and_vcf_files{"_bcftools" if opt == "all" else ""}.sh'
    shell:
        """
        bash {params.script} \
            -n {wildcards.chr} \
            -r {rsq} \
            -m {maf} \
            -p {out_dir}/pre_qc/pre_qc \
            -d {out_dir} \
            -o {opt} \
           > {log} 2>&1
        """


rule concat_convert_to_plink:
    input:
        [f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/chr{c}_clean.vcf.gz" for c in chr_noY]
    output:
        f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/chr_all_concat.pvar",
        f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/chr_all_concat.psam",
        f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/chr_all_concat.pgen"
    log:
        f"{out_dir}/imputed_clean_maf{maf}_rsq{rsq}/concat_convert_to_plink.log"
    params:
        script=Path(code_dir, "scripts/concat_convert_to_plink.sh")
    shell:
        """
        bash {params.script} \
            -p {out_dir}/pre_qc/pre_qc \
            -d {out_dir}/imputed_clean_maf{maf}_rsq{rsq} \
            > {log} 2>&1
        """
