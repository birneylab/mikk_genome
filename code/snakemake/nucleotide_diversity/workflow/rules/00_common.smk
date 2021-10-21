######################
# Libraries
######################

import pandas as pd
import numpy as np
import os

######################
# Config file
######################

configfile: "code/snakemake/nucleotide_diversity/config/config.yaml"

######################
# Samples
######################

MIKK_SAMPLES = pd.read_csv(config["samples_line_only"],
                           header = None)[0].tolist()

WILD_SAMPLES = pd.read_csv(config["kiyosu_samples_only"],
                           header = None)[0].tolist()