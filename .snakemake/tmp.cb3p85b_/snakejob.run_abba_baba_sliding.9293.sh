#!/bin/sh
# properties = {"type": "single", "rule": "run_abba_baba_sliding", "local": false, "input": ["../introgression/genos/final_no-sibs/14.txt", "data/no_sibs.pop.txt"], "output": ["data/introgression/abba_sliding/melastigma_hdrr_14.txt"], "wildcards": {"p1": "melastigma", "p2": "hdrr", "chr": "14"}, "params": {"p1": "melastigma", "p2": "hdrr", "window_length": "25000", "minimum_sites": "250"}, "log": [], "threads": 1, "resources": {}, "jobid": 9293, "cluster": {"memory": "20000", "n": "1", "name": "run_abba_baba_sliding.chr=14,p1=melastigma,p2=hdrr", "output": "../log/run_abba_baba_sliding_chr=14,p1=melastigma,p2=hdrr.out", "error": "../log/run_abba_baba_sliding_chr=14,p1=melastigma,p2=hdrr.err"}}
 cd /hps/research1/birney/users/ian/mikk_paper/mikk_genome && \
PATH='/nfs/research1/birney/users/brettell/anaconda3/envs/snakemake/bin':$PATH /nfs/research1/birney/users/brettell/anaconda3/envs/snakemake/bin/python3.9 \
-m snakemake data/introgression/abba_sliding/melastigma_hdrr_14.txt --snakefile /hps/research1/birney/users/ian/mikk_paper/mikk_genome/code/snakemake/introgression/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.cb3p85b_ ../introgression/genos/final_no-sibs/14.txt data/no_sibs.pop.txt --latency-wait 100 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules run_abba_baba_sliding --nocolor --notemp --no-hooks --nolock \
--mode 2  --use-conda  --use-singularity  && touch /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.cb3p85b_/9293.jobfinished || (touch /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.cb3p85b_/9293.jobfailed; exit 1)

