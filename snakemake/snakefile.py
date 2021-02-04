# Bash script to run on EBI cluster:

# conda activate snakemake
# snakemake \
#   --jobs 5000 \
#   --latency-wait 100 \
#   --cluster-config mikk_genome/code/snakemake/config/cluster.json \
#   --cluster 'bsub -g /snakemake_bgenie -J {cluster.name} -n {cluster.n} -M {cluster.memory} -o {cluster.output} -e {cluster.error}' \
#   --keep-going \
#   --rerun-incomplete \
#   --use-conda \
#   -s pilot_paper/snmk/snakefile.py \
#   -p

# Import functions and packages

from os.path import join
import pandas as pd

# Load config file and provide content as config object

configfile: "mikk_genome/snakemake/config/config.yaml"

# Load samples to process

SAMPLES = pd.read_csv(config["samples_file"], comment="#", skip_blank_lines=True, index_col=0)

# Globals

#ASSAYS = ["open_field", "novel_object"]
#QUADRANTS = ["q1", "q2", "q3", "q4"]

# Rules

rule all:
    input:
        expand("emfs/ensembl_102/50_fish.epo.{sample}.emf.gz",
            sample = SAMPLES.index)

rule download:
    input:
        join(config["input_dir"], "50_fish.epo.{sample}.emf.gz")
    output:
        "emfs/ensembl_102/50_fish.epo.{sample}.emf.gz"
    shell:
        "wget -P emfs/ensembl_102 {input}"
