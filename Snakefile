from collections import defaultdict
from snakemake.utils import validate
from workflow_utils import preprocess_config


configfile: "default.yaml"

include "workflow/mapping.smk"
include "workflow/calling.smk"
