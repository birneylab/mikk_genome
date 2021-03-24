#!/bin/sh
# properties = {"type": "single", "rule": "add_kaga_read_groups", "local": false, "input": ["../introgression/sams/kaga.sam"], "output": ["../introgression/bams/kaga_with_read_groups.bam"], "wildcards": {}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 9954, "cluster": {"memory": "20000", "n": "16", "name": "add_kaga_read_groups.", "output": "../log/add_kaga_read_groups_.out", "error": "../log/add_kaga_read_groups_.err"}}
 cd /hps/research1/birney/users/ian/mikk_paper/mikk_genome && \
PATH='/nfs/research1/birney/users/brettell/anaconda3/envs/snakemake/bin':$PATH /nfs/research1/birney/users/brettell/anaconda3/envs/snakemake/bin/python3.9 \
-m snakemake ../introgression/bams/kaga_with_read_groups.bam --snakefile /hps/research1/birney/users/ian/mikk_paper/mikk_genome/code/snakemake/introgression/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.ysoerd5o ../introgression/sams/kaga.sam --latency-wait 100 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules add_kaga_read_groups --nocolor --notemp --no-hooks --nolock \
--mode 2  --use-conda  --use-singularity  && touch /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.ysoerd5o/9954.jobfinished || (touch /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.ysoerd5o/9954.jobfailed; exit 1)

