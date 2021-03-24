#!/bin/sh
# properties = {"type": "single", "rule": "align_kaga_to_hdrr", "local": false, "input": ["../introgression/fastas/DRR002225_1.fastq.gz", "../introgression/fastas/DRR002225_2.fastq.gz"], "output": ["../introgression/sams/DRR002225.sam"], "wildcards": {"kaga_sample": "DRR002225"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 9, "cluster": {"memory": "20000", "n": "16", "name": "align_kaga_to_hdrr.kaga_sample=DRR002225", "output": "../log/align_kaga_to_hdrr_kaga_sample=DRR002225.out", "error": "../log/align_kaga_to_hdrr_kaga_sample=DRR002225.err"}}
 cd /hps/research1/birney/users/ian/mikk_paper/mikk_genome && \
PATH='/nfs/research1/birney/users/brettell/anaconda3/envs/snakemake/bin':$PATH /nfs/research1/birney/users/brettell/anaconda3/envs/snakemake/bin/python3.9 \
-m snakemake ../introgression/sams/DRR002225.sam --snakefile /hps/research1/birney/users/ian/mikk_paper/mikk_genome/code/snakemake/introgression/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.0zwmzoyt ../introgression/fastas/DRR002225_1.fastq.gz ../introgression/fastas/DRR002225_2.fastq.gz --latency-wait 100 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules align_kaga_to_hdrr --nocolor --notemp --no-hooks --nolock \
--mode 2  --use-conda  --use-singularity  && touch /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.0zwmzoyt/9.jobfinished || (touch /hps/research1/birney/users/ian/mikk_paper/mikk_genome/.snakemake/tmp.0zwmzoyt/9.jobfailed; exit 1)

