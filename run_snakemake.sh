#!/bin/bash

set -e
set -u

# Help message
usage() {
    echo "Usage: $0 -s <step> -f <config> [-d <dry-run>] [-c <cores>] [-u <unlock>] [-h]"
    echo "Options:"
    echo "  -s  Step to run (required): submit_initial_input, submit_fix_strands, or filter_info_and_vcf_files"
    echo "  -f  Path to config file"
    echo "  -d  Run snakemake --dry-run"
    echo "  -c  Number of cores to use (default: 6)"
    echo "  -u  Unlock snakemake"
    echo "  -h  Show this help message"
    exit 1
}

# Run help message if no args given
if [ "$#" -eq 0 ]; then
    usage
fi

# Default args
dry_run=0
n_cores=6
unlock=0

# Parse args
while getopts "s:f:c:duh" opt; do
  case $opt in
    s) step="$OPTARG" ;;
    f) config="$OPTARG" ;;
    d) dry_run=1 ;;
    c) n_cores="$OPTARG" ;;
    u) unlock=1 ;;
    h) usage ;;
    \?) usage ;;
  esac
done

# Set dry-run flag
if [ "$dry_run" -eq 1 ]; then
    dry_flag="--dry-run"
else
    dry_flag=""
fi

# Set unlock flag
if [ "$unlock" -eq 1 ]; then
    unlock_flag="--unlock"
else
    unlock_flag=""
fi

# Get info about config file path passed in
config_name=$(basename "$config")
config_path=$(dirname "$config")

# Get values set in config file
    # Note that grep only works on simple, single-level .yml file. Using "yq"
    # would be better, but requires tool to be available outside of container
plink_prefix=$(sed 's/#.*//' "$config" | grep "^plink_prefix:" | cut -d':' -f2 | tr -d ' "')
plink_dir=$(dirname "$plink_prefix")

id_list=$(sed 's/#.*//' "$config" | grep "^id_list:" | cut -d':' -f2 | tr -d ' "')
id_list_dir=$(dirname "$id_list")

id_list_hwe=$(sed 's/#.*//' "$config" | grep "^id_list_hwe:" | cut -d':' -f2 | tr -d ' "')
id_list_hwe_dir=$(dirname "$id_list_hwe")

out_dir=$(sed 's/#.*//' "$config" | grep "^out_dir:" | cut -d':' -f2 | tr -d ' "')
repo=$(sed 's/#.*//' "$config" | grep "^repo:" | cut -d':' -f2 | tr -d ' "')
use_cont=$(sed 's/#.*//' "$config" | grep "^use_cont:" | cut -d':' -f2 | tr -d ' "')
cont_engine=$(sed 's/#.*//' "$config" | grep "^cont_engine:" | cut -d':' -f2 | tr -d ' "')
cont_engine=${cont_engine:-apptainer}
docker_image="ghcr.io/rkjcollab/imputation:1.1.3"  # image matching this release

# Container paths - fixed here rather than in the config file, and passed to
# snakemake with --config since the Snakefile needs them inside the container
plink_dir_cont="/input_data"
id_list_dir_cont="/id_list"
id_list_hwe_dir_cont="/id_list_hwe"
out_dir_cont="/output_data"
repo_cont="/repo"
cont_config=(--config \
    "plink_dir_cont=${plink_dir_cont}" \
    "id_list_dir_cont=${id_list_dir_cont}" \
    "id_list_hwe_dir_cont=${id_list_hwe_dir_cont}" \
    "out_dir_cont=${out_dir_cont}")

if [ "$use_cont" = "false" ]; then
    # Run snakemake on local machine using provided conda environment
    snakemake --rerun-triggers mtime --snakefile ${repo}/Snakefile \
        --configfile ${config_path}/${config_name} \
        --cores "$n_cores" "$step" $dry_flag $unlock_flag

elif [ "$use_cont" = "true" ] && [ "$cont_engine" = "docker" ]; then
    # Run snakemake in Docker container. Docker requires absolute host paths,
    # so resolve them first, and only mount optional ID lists if given.
    abs_path() { (cd "$1" && pwd); }
    repo_abs=$(abs_path "$repo")
    plink_dir_abs=$(abs_path "$plink_dir")
    config_path_abs=$(abs_path "$config_path")
    out_dir_abs=$(abs_path "$out_dir")

    id_list_mounts=()
    if [ -n "$id_list" ]; then
        id_list_mounts+=(-v "$(abs_path "$id_list_dir"):${id_list_dir_cont}")
    fi
    if [ -n "$id_list_hwe" ]; then
        id_list_mounts+=(-v "$(abs_path "$id_list_hwe_dir"):${id_list_hwe_dir_cont}")
    fi

    # Run as host user so output is not owned by root, and mount the working
    # dir as apptainer does so .snakemake is written to the same place
    docker run --rm \
        --platform linux/amd64 \
        --user "$(id -u):$(id -g)" \
        --env HOME=/tmp \
        --env USER="$(id -un)" \
        -v "${repo_abs}:${repo_cont}" \
        -v "${plink_dir_abs}:${plink_dir_cont}" \
        ${id_list_mounts[@]+"${id_list_mounts[@]}"} \
        -v "${config_path_abs}:/proj_repo" \
        -v "${out_dir_abs}:${out_dir_cont}" \
        -v "${PWD}:${PWD}" \
        -w "${PWD}" \
        "$docker_image" \
        snakemake --rerun-triggers mtime --snakefile ${repo_cont}/Snakefile \
            --configfile /proj_repo/${config_name} "${cont_config[@]}" \
            --cores "$n_cores" "$step" $dry_flag $unlock_flag

elif [ "$use_cont" = "true" ] && [ "$cont_engine" = "apptainer" ]; then
    # Run snakemake in apptainer container (default)
    apptainer exec \
        --writable-tmpfs \
        --bind "${repo}:${repo_cont}" \
        --bind "${plink_dir}:${plink_dir_cont}" \
        --bind "${id_list_dir}:${id_list_dir_cont}" \
        --bind "${id_list_hwe_dir}:${id_list_hwe_dir_cont}" \
        --bind "${config_path}:/proj_repo" \
        --bind "${out_dir}:${out_dir_cont}" \
        ${repo}/envs/env_imputation.sif \
        snakemake --rerun-triggers mtime --snakefile ${repo_cont}/Snakefile \
            --configfile /proj_repo/${config_name} "${cont_config[@]}" \
            --cores "$n_cores" "$step" $dry_flag $unlock_flag

else
    echo "Config file use_cont must be either true or false, and cont_engine must be either apptainer or docker"
    exit 1
fi
