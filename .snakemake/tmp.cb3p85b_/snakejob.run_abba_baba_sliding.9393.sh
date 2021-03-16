#!/bin/sh
# properties = {"type": "single", "rule": "run_abba_baba_sliding", "local": false, "input": ["../introgression/genos/final_no-sibs/8.txt", "data/no_sibs.pop.txt"], "output": ["data/introgression/abba_sliding/javanicus_hsok_8.txt"], "wildcards": {"p1": "javanicus", "p2": "hsok", "chr": "8"}, "params": {"p1": "javanicus", "p2": "hsok", "window_length": "25000", "minimum_sites": "250"}, "log": [], "threads": 1, "resources": {}, "jobid": 9393, "cluster": {"memory": "20000", "n": "1", "name": "run_abba_baba_sliding.chr=8,p1=javanicus,p2=hsok", "output": "../log/run_abba_baba_sliding_chr=8,p1=javanicus,p2=hsok.out", "error": "../log/run_abba_baba_sliding_chr=8,p1=javanicus,p2=hsok.err"}}
 cd /hps/research1/birney/users/ian/mikk_paper/mikk_genome && \
PATH='/nfs/research1/birney/users/brettell/anaconda3/envs/snakemake/bin':$PATH /nfs/research1/birney/users/brettell/anaconda3/envs/snakemake/bin/python3.9 \
-m snakemake data/introgression/abba_sliding/javanicus_hsok_8.txt --snakefile /hps/research1/birney/users/ian/mikk_paper/mikk_genome/code/snakemake/introgression/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.cb3p85b_ ../introgression/genos/final_no-sibs/8.txt data/no_sibs.pop.txt --latency-wait 100 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules run_abba_baba_sliding --nocolor --notemp --no-hooks --nolock \
--mode 2  --use-conda  --use-singularity  && touch /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.cb3p85b_/9393.jobfinished || (touch /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.cb3p85b_/9393.jobfailed; exit 1)

