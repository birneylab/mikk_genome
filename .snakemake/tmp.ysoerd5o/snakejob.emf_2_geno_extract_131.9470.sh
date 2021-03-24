#!/bin/sh
# properties = {"type": "single", "rule": "emf_2_geno_extract_131", "local": false, "input": ["../introgression/genos/final_full-run/19.txt", "data/sv_analysis/20210323_full-run_except_131-2.txt"], "output": ["../introgression/genos/final_131-2/19.txt"], "wildcards": {"line_131_sib": "131-2", "chr": "19"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 9470, "cluster": {"memory": "20000", "n": "1", "name": "emf_2_geno_extract_131.chr=19,line_131_sib=131-2", "output": "../log/emf_2_geno_extract_131_chr=19,line_131_sib=131-2.out", "error": "../log/emf_2_geno_extract_131_chr=19,line_131_sib=131-2.err"}}
 cd /hps/research1/birney/users/ian/mikk_paper/mikk_genome && \
PATH='/nfs/research1/birney/users/brettell/anaconda3/envs/snakemake/bin':$PATH /nfs/research1/birney/users/brettell/anaconda3/envs/snakemake/bin/python3.9 \
-m snakemake ../introgression/genos/final_131-2/19.txt --snakefile /hps/research1/birney/users/ian/mikk_paper/mikk_genome/code/snakemake/introgression/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.ysoerd5o ../introgression/genos/final_full-run/19.txt data/sv_analysis/20210323_full-run_except_131-2.txt --latency-wait 100 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules emf_2_geno_extract_131 --nocolor --notemp --no-hooks --nolock \
--mode 2  --use-conda  --use-singularity  && touch /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.ysoerd5o/9470.jobfinished || (touch /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.ysoerd5o/9470.jobfailed; exit 1)

