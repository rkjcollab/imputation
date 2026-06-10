from pathlib import Path

# Setup -----------------------------------------------------------------------------------------

### Get variables from config file
# Required pipeline settings
chr: list[str] = config["chr"]
orig_build: str = config["orig_build"]
to_build: str = config["to_build"]
imp: str = config["imp"]
imp_name: str = config["imp_name"]
zip_pw: str = config["zip_pw"]
imp_job_id: str = config["imp_job_id"]

# Required pipeline settings with defaults
imp_rsq_filt: str = config.get("imp_rsq_filt", "0")
opt: str = config.get("opt", "gt")
use_cont: bool = config.get("use_cont", True)
sexcheck: list[str] = config.get("sexcheck", ["0.8", "0.2"])

# Required host paths outside container
plink_prefix: str = config["plink_prefix"]
plink_prefix_name = Path(plink_prefix).name

# Optional host paths outside container
id_list: str = config.get("id_list", None)
id_list_hwe: str = config.get("id_list_hwe", None)

if use_cont:
    # Container paths
    plink_dir: str = config["plink_dir_cont"]
    out_dir: str = config["out_dir_cont"]
else:
    plink_dir = Path(plink_prefix).parent
    out_dir: str = config["out_dir"]

# Modify arrays for bash script
chr_str = " ".join(map(str, chr))
chr_noY = [c for c in chr if c != "Y"]
chr_noY_str = " ".join(map(str, chr_noY))
sexcheck_str = " ".join(map(str, sexcheck))

### Other prep
# Set default values currently not controlled by config file
maf = "0"
rsq = "0.3"

# Get dir of pipeline
code_dir = workflow.current_basedir

# Rules -----------------------------------------------------------------------------------------

include: "workflow/rules/targets.smk"
include: "workflow/rules/preprocess.smk"
include: "workflow/rules/postprocess.smk"
