#!/bin/sh
# properties = {"type": "single", "rule": "emf_2_geno_remove_sibs", "local": false, "input": ["../introgression/genos/final_full-run/21.txt", "data/20200227_panel_lines_excluded.txt"], "output": ["../introgression/genos/final_no-sibs/21.txt"], "wildcards": {"chr": "21"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 9267, "cluster": {"memory": "20000", "n": "1", "name": "emf_2_geno_remove_sibs.chr=21", "output": "../log/emf_2_geno_remove_sibs_chr=21.out", "error": "../log/emf_2_geno_remove_sibs_chr=21.err"}}
 cd /hps/research1/birney/users/ian/mikk_paper/mikk_genome && \
PATH='/nfs/research1/birney/users/brettell/anaconda3/envs/snakemake/bin':$PATH /nfs/research1/birney/users/brettell/anaconda3/envs/snakemake/bin/python3.9 \
-m snakemake ../introgression/genos/final_no-sibs/21.txt --snakefile /hps/research1/birney/users/ian/mikk_paper/mikk_genome/code/snakemake/introgression/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.cb3p85b_ ../introgression/genos/final_full-run/21.txt data/20200227_panel_lines_excluded.txt --latency-wait 100 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules emf_2_geno_remove_sibs --nocolor --notemp --no-hooks --nolock \
--mode 2  --use-conda  --use-singularity  && touch /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.cb3p85b_/9267.jobfinished || (touch /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.cb3p85b_/9267.jobfailed; exit 1)

