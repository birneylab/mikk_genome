#!/bin/bash

set -euo pipefail

bsub \
      	-J NanoSnake_RNA_illumina_medaka \
        -eo bsub_err.txt \
        -oo bsub_out.txt \
        -n 3 \
	-M 5000 \
	NanoSnake RNA_illumina \
		-g ftp://ftp.ensembl.org/pub/release-98/fasta/oryzias_latipes/dna/Oryzias_latipes.ASM223467v1.dna_sm.toplevel.fa.gz \
	        -t ftp://ftp.ensembl.org/pub/release-98/fasta/oryzias_latipes/cdna/Oryzias_latipes.ASM223467v1.cdna.all.fa.gz \
		-a ftp://ftp.ensembl.org/pub/release-98/gff3/oryzias_latipes/Oryzias_latipes.ASM223467v1.98.gff3.gz \
		-s all_sample_sheet.tsv \
		--cluster_config cluster_config.yaml \
		--restart_times 1 \
		--conda_prefix /nfs/software/birney/adrien/miniconda3/envs/NanoSnake_wrappers
