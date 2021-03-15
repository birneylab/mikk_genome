#!/bin/sh
# properties = {"type": "single", "rule": "split_mikk_vcf", "local": false, "input": ["../introgression/vcfs/ill_no-sibs_no-missing.vcf.gz", "data/introgression/20200914_abbababa_mikk_1.txt"], "output": ["../introgression/vcfs/ill_no-sibs_no-missing_bi-snps_with-af_1.vcf.gz"], "wildcards": {"group": "1"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 9034, "cluster": {"memory": "10000", "n": "1", "name": "split_mikk_vcf.group=1", "output": "../log/split_mikk_vcf_group=1.out", "error": "../log/split_mikk_vcf_group=1.err"}}
 cd /hps/research1/birney/users/ian/mikk_paper/mikk_genome && \
PATH='/nfs/research1/birney/users/brettell/anaconda3/envs/snakemake/bin':$PATH /nfs/research1/birney/users/brettell/anaconda3/envs/snakemake/bin/python3.9 \
-m snakemake ../introgression/vcfs/ill_no-sibs_no-missing_bi-snps_with-af_1.vcf.gz --snakefile /hps/research1/birney/users/ian/mikk_paper/mikk_genome/code/snakemake/introgression/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.lgsughb2 ../introgression/vcfs/ill_no-sibs_no-missing.vcf.gz data/introgression/20200914_abbababa_mikk_1.txt --latency-wait 100 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules split_mikk_vcf --nocolor --notemp --no-hooks --nolock \
--mode 2  --use-conda  --use-singularity  && touch /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.lgsughb2/9034.jobfinished || (touch /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.lgsughb2/9034.jobfailed; exit 1)

