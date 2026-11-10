## **imputation**

*PLEASE NOTE that this repo is still under active development and testing.*

Snakemake-run pipeline for imputing autosomal genetic array data with hg19 or hg38
coordinates to the TOPMed r3 reference panel, the Michigian 1000 genomes phase 3
version 5 reference panel, or the Michigan HLA four-digit multi-ethnic v1 or v2 panels.

Adapted from Michelle Daya's topmed_imputation pipeline.

## Setup

### Apptainer

The apptainer image is built from the Docker image matching this release
(`ghcr.io/rkjcollab/imputation:1.1.0`). From the top level of this repository run:

```
apptainer build envs/env_imputation.sif docker://ghcr.io/rkjcollab/imputation:1.1.0
```

### Docker

To use Docker instead of apptainer, set `cont_engine: "docker"` in the config file.
The image matching this release (`ghcr.io/rkjcollab/imputation:1.1.0`) is pulled
from GHCR. The image is built for linux/amd64, so it runs under emulation (slower)
on Apple Silicon Macs. If needed, the image can instead be built from the top level
of this repository with:

```
docker buildx build --platform linux/amd64 -t ghcr.io/rkjcollab/imputation:1.1.0 - < envs/Dockerfile
```

### Conda Environment

To run without a container (`use_cont: false`), or for development, build the
mamba environment on a linux-64, osx-64, or osx-arm64 operating system with:

```
mamba env create -f envs/env_imputation.yml
```

Tool versions match the container, and the environment also includes
development tools (lintr, styler, and pre-commit).

## Input Files

* Genetic array data in PLINK1.9 file format.

## How to Run

### Step 1

Copy "config.yml" to the project repo want to run imputation for, and update paths and
settings specific to the data:

```
### TODO: Update pipeline settings
chr: [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]  # list of chromosome numbers with format [1,10,22]
orig_build: "hg19"  #  either "hg19" or "hg38"
to_build: "hg38"  # either "hg19" or "hg38"
imp: "topmed"  # should be 'topmed', 'mich_hla_v1', 'mich_hla_v2', 'mich_1kg_p3_v5', 'mich_hrc'
imp_rsq_filt: "0"  # DEFAULT: "0", filter to apply on imputation server (0, 0.001, 0.1, 0.2, or 0.3)
imp_name: "job-name"  # job name for imputation server
imp_job_id: ""  # should be "job-####-##-###"
zip_pw: ""  # add from email when imputed job finishes
opt: "gt"  # DEFAULT: "gt", option for "gt" only or "all" dosage information (much slower)
use_cont: true  # DEFAULT: true, true to use container, false to use local machine
cont_engine: "apptainer"  # DEFAULT: "apptainer", container engine to use, either "apptainer" or "docker"
sexcheck: [0.8, 0.2]  # DEFAULT: [0.8, 0.2]
maf: "0"  # DEFAULT: "0", MAF filter applied after imputation
rsq: "0.3"  # DEFAULT: "0.3", Rsq filter applied to imputed variants after imputation

### TODO: Update host paths outside container
# PLINK1.9 file prefix with path - must be single file, but can contain one or more chr
plink_prefix: "/path/to/plink_prefix"
# OPTIONAL: ID list (FID IID) want to run imputation for, if none use ""
id_list: ""
# OPTIONAL: ID list to calculate HWE filtering in (likely controls), if none use ""
id_list_hwe: ""
# Top-level directory for all output - must contain the pre-QC PLINK1.9 file at
# out_dir/pre_qc/pre_qc
out_dir: "/path/to/output"
# Repo with snakemake pipeline
repo: "/path/to/imputation"
```

Note that some variables (imp_job_id & zip_pw) will need to be added to the config file as
the pipeline is run.

Note that imp_rsq_filt refers to a filter based on imputation quality Rsq that would be
applied by the imputation server. To preserve information about all variants the server
attempted to impute, it is best to set this to 0. If necessary, due to large file size,
a filter can be applied on the server by setting this option. However, although the server
generally returns TYPED input variants, it does almost attempt to impute them and applying
the server-side Rsq filter appears to removed TYPED variants that imputed with an Rsq
less than the filter as of 10/29/25.

Note that for TOPMed imputation, the 'population' option is set to 'all', which means allele
frequencies will be compared between the input data and the TOPMed panel to generate a QC
report. For Michigan imputation, the 'population' option is set to 'off', since (at least
for 1000G whole genome imputation) there is no good option for allele frequency comparison
for mixed populations. Either way, imputation results should not be affected.

If `use_cont` is set to false, activate the conda environment before running the
pipeline with `mamba activate imputation`.

### Step 2

Use the "run_snakemake.sh" script to run the snakefile through initial data preparation and
submission to the imputation server for QC only. This will create and upload the initial
input files to the imputation server and panel selected in the config file.

```
bash run_snakemake.sh \
    -s submit_initital_input \
    -f /path/to/project/config.yml

```

Note that you can add "-d" to run a snakemake --dry-run, use "-c" to set the number of
cores used (the default is 6), or use "-h" to show the run script help message.

This will create and upload the initial input files to the imputation server and panel
selected in the config file.

The pre-imputation QC includes liftover (if necessary), removal of strand ambiguous SNPs,
updating all variant IDs to chr:pos:ref:alt, filtering by PLINK2 --maf (1e-6), --geno (0.05),
and --hwe (1e-20 for chr6 MHC region, 1e-6 for all other regions/chromosomes). A summary of
each step is included in the log file.

### Step 3

Once the QC automatically submitted has run, you will receive an email. The pipeline here
will download the results for you, modify the input based on the server QC, and resubmit
to the server for imputation. However, you must log onto the imputation server, get the
job ID, and add it to the project-specific config.yml first:

```
### TODO: Update pipeline settings
# ...
imp_job_id: "job-1234-56-789"  # should be "job-####-##-###"
zip_pw: ""  # add from email when imputed job finishes

```

Only add the imp_job_id for now. You will add the zip_pw later.

Use the "run_snakemake.sh" to run the pipeline through the "submit_fix_strands" step:

```
bash run_snakemake.sh \
    -s submit_fix_strands \
    -f /path/to/project/config.yml

```

This step flips strands of variants identified as such in the snps-exlcuded.txt file - this
will produce the final post QC VCF files for imputation.

*TO NOTE: if using Michigan imputation server, may need to move files out of*
*subfolders manually.*

### Step 4

Once the imputation has run, you will receive an email. The pipeline here
will download the results for you, filter them, and create PLINK files. However,
you must log onto the imputation server, get the updated job ID, and add it and
the zip password that was emailed to you to the project-specific config.yml first:

```
### TODO: Update pipeline settings
# ...
imp_job_id: "job-2345-67-891"  # should be "job-####-##-###"
zip_pw: "abc123zyx987"  # add from email when imputed job finishes

```

Use the "run_snakemake.sh" to run the pipeline through the "filter_info_and_vcf_files" step:

```
bash run_snakemake.sh \
    -s filter_info_and_vcf_files \
    -f /path/to/project/config.yml

```


## Development Notes

### CI/CD

TODO: add detail for this. CI/CD on GitHub will check styler and linter.
Recommend using pre commit hook locally first always:

```
mamba activate imputation  # pre-commit is included in the conda environment
pre-commit install

```

### TODOs

* update below!
* add instructions for getting TOPMed API key and/or manual submit?
* add instructions for key.py folder with API keys
* add checks that config values are only the allowable ones?
* add more documentation to python functions
* think could remove a config arg so that to_build is auto set based on chosen imputation
* not sure how to better handle relative paths within the pipeline directory?
* add documentation about geno in all individuals followed by hwe in controls -
    may lead to warning since subsetting increases variant missingness, but
    think should be okay to ignore: https://groups.google.com/g/plink2-users/c/xjmPNfN0Swc/m/79_CAtMRBwAJ

* discuss sex chromosome options, note that SEX column in PLINK file must be filled out
* awknowledge that better than duplicate removal would be to consider fixing mismapped
    SNPs using reference panel... potential future step
* add note that create initial input should be run twice so that plot with sex checks
    can be examined and thresholds can be adjusted
* automate getting job-id from log files
* add note the pre-filtering, sex-checks done for all chr present in given file, but
    SNP dedup, HWE and imputation prep only done for requested chr to prepare
* add note that chrY can be prepared but not imputed
* split into two snakefiles - one pre- and one post-QC (would help with post QC paths
    being different)
* need to automate removal of imputed files that don't want to keep
* need to confirm automated download works for Michigan imputation server
* Update names for pre_qc (final input) vs. pre_qc_hg38 (after liftover but before
    strang ambig removal)
