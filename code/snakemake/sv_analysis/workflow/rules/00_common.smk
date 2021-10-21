######################
# Libraries
######################

import pandas as pd
import numpy as np
import os

######################
# Config file
######################

configfile: "code/snakemake/sv_analysis/config/config.yaml"

SAMPLES = config["samples"]