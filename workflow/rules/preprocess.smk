rule create_initial_input:
    input:
        f"{plink_dir}/{plink_prefix_name}.bed",
        f"{plink_dir}/{plink_prefix_name}.bim",
        f"{plink_dir}/{plink_prefix_name}.fam"
    output:
        [f"{out_dir}/pre_qc/chr{c}_pre_qc.vcf.gz" for c in chr]
    log:
        f"{out_dir}/pre_qc/create_initial_input.log"
    params:
        script=Path(code_dir, "scripts/create_initial_input.sh"),
        summary=f"{out_dir}/pre_qc/summary_pre_qc.txt"
    run:
        cmd = (
            f"bash {params.script} "
            f"-p {plink_dir}/{plink_prefix_name} "
            f"-o {out_dir}/pre_qc "
            f"-c {code_dir} "
            f"-n \"{chr_str}\" "
            f"-b {orig_build} "
            f"-t {to_build} "
            f"-s \"{sexcheck_str}\" "
            f"-l {params.summary} "
        )
        # If id_list given for filtering
        if id_list:
            id_list_name=Path(id_list).name
            if use_cont:
                id_list_dir: str = config["id_list_dir_cont"]
            else:
                id_list_dir = Path(id_list).parent
            cmd += f" -k {id_list_dir}/{id_list_name}"
        # If id_list_hwe given for HWE calculation
        if id_list_hwe:
            id_list_hwe_name=Path(id_list_hwe).name
            if use_cont:
                id_list_hwe_dir: str = config["id_list_hwe_dir_cont"]
            else:
                id_list_hwe_dir = Path(id_list_hwe).parent
            cmd += f" -h {id_list_hwe_dir}/{id_list_hwe_name}"
        # Redirect logs
        cmd += f" > {log} 2>&1"
        # Run
        shell(cmd)


# Need "" around chr_str to keep all as single input
rule submit_initial_input:
    input:
        [f"{out_dir}/pre_qc/chr{c}_pre_qc.vcf.gz" for c in chr]
    output:
        log_final=f"{out_dir}/pre_qc/submit_initial_input.log"
    params:
        script=Path(code_dir, "scripts/submit.py"),
        log_tmp=f"{out_dir}/pre_qc/tmp_submit_initial_input.log"
    shell:
        """
        python {params.script} \
            --dir {out_dir}/pre_qc \
            --chr "{chr_noY_str}" \
            --imp {imp} \
            --build {to_build} \
            --mode "qconly" \
            --imp-name {imp_name} \
            > {params.log_tmp} 2>&1

        mv {params.log_tmp} {output.log_final}
        """


rule download_qc:
    output:
        log_final=f"{out_dir}/pre_qc/download_qc.log",
        snps_excl=f"{out_dir}/pre_qc/snps-excluded.txt"
    params:
        script=Path(code_dir, "scripts/download_results.sh"),
        log_tmp=f"{out_dir}/pre_qc/tmp_download_qc.log"
    shell:
        """
        bash {params.script} \
            -i {imp} \
            -c {code_dir} \
            -o {out_dir}/pre_qc \
            -j {imp_job_id} \
            > {params.log_tmp} 2>&1

        mv {params.log_tmp} {output.log_final}
        """


rule fix_strands:
    input:
        f"{out_dir}/pre_qc/snps-excluded.txt",
        f"{out_dir}/pre_qc/download_qc.log"
    output:
        [f"{out_dir}/post_qc/chr{c}_post_qc.vcf.gz" for c in chr_noY]
    log:
        f"{out_dir}/post_qc/fix_strands.log"
    params:
        script=Path(code_dir, "scripts/fix_strands.sh")
    shell:
        """
        bash {params.script} \
            -o {out_dir} \
            -c {code_dir} \
            -n "{chr_noY_str}" \
            -t {to_build} \
            -i {imp} \
            > {log} 2>&1
        """


rule submit_fix_strands:
    input:
        [f"{out_dir}/post_qc/chr{c}_post_qc.vcf.gz" for c in chr_noY]
    output:
        log_final=f"{out_dir}/post_qc/submit_fix_strands.log"
    params:
        script=Path(code_dir, "scripts/submit.py"),
        log_tmp=f"{out_dir}/post_qc/tmp_submit_fix_strands.log"
    shell:
        """
        python {params.script} \
            --dir {out_dir}/post_qc \
            --chr "{chr_noY_str}" \
            --imp {imp} \
            --build {to_build} \
            --mode imputation \
            --rsq-filt {imp_rsq_filt} \
            --imp-name {imp_name} \
            > {params.log_tmp} 2>&1

        mv {params.log_tmp} {output.log_final}
        """
