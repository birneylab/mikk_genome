####################
# EBI codon cluster
####################

ssh codon
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
bsub -M 20000 -Is bash
cd /hps/software/users/birney/ian/repos/mikk_genome
conda activate snakemake_6.7.0
snakemake \
  --jobs 5000 \
  --latency-wait 100 \
  --cluster-config code/snakemake/nucleotide_diversity/config/cluster.yaml \
  --cluster 'bsub -g /snakemake_bgenie -J {cluster.name} -q {cluster.queue} -n {cluster.n} -M {cluster.memory} -o {cluster.outfile}' \
  --keep-going \
  --rerun-incomplete \
  --use-conda \
  --use-singularity \
  -s code/snakemake/nucleotide_diversity/workflow/Snakefile \
  -p

####################
# RStudio Server
####################

# Set container path
CONT=/hps/software/users/birney/ian/containers/mikk_genome/nucleotide_diversity/R_4.1.0.sif

# Build container
singularity build --remote \
    $CONT \
    code/snakemake/nucleotide_diversity/workflow/envs/R_4.1.0.def

# Run RStudio Server
ssh proxy-codon
bsub -M 20000 -Is bash
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp 
singularity shell --bind /hps/software/users/birney/ian/rstudio_db:/var/lib/rstudio-server \
                  --bind /hps/software/users/birney/ian/tmp:/tmp \
                  --bind /hps/software/users/birney/ian/run:/run \
                  $CONT
# Then run rserver, setting path of config file containing library path
rserver --rsession-config-file /hps/software/users/birney/ian/repos/mikk_genome/code/snakemake/nucleotide_diversity/workflow/envs/rsession.conf

ssh -L 8787:hl-codon-12-01:8787 proxy-codon
