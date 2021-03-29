#!/bin/sh
# properties = {"type": "single", "rule": "kaga_merge_mikk_vcf", "local": false, "input": ["../introgression/vcfs/merged/kaga.vcf.gz", "../introgression/vcfs/full-run_line-ids.vcf.gz"], "output": ["../introgression/vcfs/full-run_line-ids_with-kaga.vcf.gz"], "wildcards": {}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 103, "cluster": {"memory": "5000", "n": "1", "name": "kaga_merge_mikk_vcf.", "output": "../log/kaga_merge_mikk_vcf_.out", "error": "../log/kaga_merge_mikk_vcf_.err"}}
 cd /hps/research1/birney/users/ian/mikk_paper/mikk_genome && \
PATH='/nfs/research1/birney/users/brettell/anaconda3/envs/snakemake/bin':$PATH /nfs/research1/birney/users/brettell/anaconda3/envs/snakemake/bin/python3.9 \
-m snakemake ../introgression/vcfs/full-run_line-ids_with-kaga.vcf.gz --snakefile /hps/research1/birney/users/ian/mikk_paper/mikk_genome/code/snakemake/introgression/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.kykj44vn ../introgression/vcfs/merged/kaga.vcf.gz ../introgression/vcfs/full-run_line-ids.vcf.gz --latency-wait 100 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules kaga_merge_mikk_vcf --nocolor --notemp --no-hooks --nolock \
--mode 2  --use-conda  --use-singularity  && touch /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.kykj44vn/103.jobfinished || (touch /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.kykj44vn/103.jobfailed; exit 1)

