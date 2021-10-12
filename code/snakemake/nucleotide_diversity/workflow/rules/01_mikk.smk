## Filter VCF for non-sibling MIKK lines (N = 63)
rule keep_no_sibs:
    input:
        vcf = config["mikk_vcf"],
        samples_cram = config["samples_cram"]
    output:
        os.path.join(config["working_dir"], "vcfs/mikk_no-sibs_cram-id.vcf.gz")
    log:
        os.path.join(config["working_dir"], "logs/keep_no_sibs.log")
    container:
        config["bcftools"]
    shell:
        """
        bcftools view \
            --output-file {output[0]} \
            --samples-file {input.samples_cram} \
            {input.vcf} \
                2> {log}
        """

## Replace CRAM IDs with MIKK line IDs
rule rehead:
    input:
        vcf = os.path.join(config["working_dir"], "vcfs/mikk_no-sibs_cram-id.vcf.gz"),
        samples = config["samples"]
    output:
        os.path.join(config["working_dir"], "vcfs/mikk_no-sibs_line-id.vcf")
    log:
        os.path.join(config["working_dir"], "logs/rehead.log")
    container:
        config["bcftools"]
    shell:
        """
        bcftools reheader \
            --output {output[0]} \
            --samples {input.samples} \
            {input.vcf} \
                2> {log}
        """

## Compress
rule compress:
    input:
        os.path.join(config["working_dir"], "vcfs/mikk_no-sibs_line-id.vcf")
    output:
        os.path.join(config["working_dir"], "vcfs/mikk_no-sibs_line-id.vcf.gz")
    log:
        os.path.join(config["working_dir"], "logs/compress.log")
    container:
        config["bcftools"]
    shell:
        """
        bcftools view \
            --output-type z \
            --output-file {output[0]} \
            {input[0]} \
                2> {log}
        """

## Index
rule index:
    input:
        os.path.join(config["working_dir"], "vcfs/mikk_no-sibs_line-id.vcf.gz")
    output:
        os.path.join(config["working_dir"], "vcfs/mikk_no-sibs_line-id.vcf.gz.tbi")
    log:
        os.path.join(config["working_dir"], "logs/index.log")
    container:
        config["bcftools"]
    shell:
        """
        bcftools index \
            --tbi \
            {input[0]} \
                2> {log}
        """

## Caculate nucleotide divergence
rule nucleotide_divergence:
    input:
        vcf = os.path.join(config["working_dir"], "vcfs/mikk_no-sibs_line-id.vcf.gz"),
        index = os.path.join(config["working_dir"], "vcfs/mikk_no-sibs_line-id.vcf.gz.tbi")
    output:
        os.path.join(config["lts_dir"], "nucleotide_divergence/mikk/all/{window_size}.windowed.pi")
    log:
        os.path.join(config["working_dir"], "logs/nucleotide_divergence/{window_size}.log")
    params:
        window_size = lambda wildcards: wildcards.window_size,
        prefix = lambda wildcards: os.path.join(config["lts_dir"], "nucleotide_divergence", wildcards.window_size)
    container:
        config["vcftools"]
    shell:
        """
        vcftools \
            --gzvcf {input.vcf} \
            --window-pi {params.window_size} \
            --out {params.prefix} \
                2> {log}
        """

# Take random sample of 7 MIKK lines
rule filter_random_mikk:
    input:
        vcf = os.path.join(config["working_dir"], "vcfs/mikk_no-sibs_line-id.vcf.gz"),
        index = os.path.join(config["working_dir"], "vcfs/mikk_no-sibs_line-id.vcf.gz.tbi"),
        samples_file = config["random_7_mikk"]
    output:
        os.path.join(config["working_dir"], "vcfs/mikk_no-sibs_line-id_random-7.vcf.gz")
    log:
        os.path.join(config["working_dir"], "logs/filter_random_mikk.log")
    container:
        config["bcftools"]
    shell:
        """
        bcftools view \
            --output-file {output[0]} \
            --samples-file {input.samples_file} \
            {input.vcf} \
                2> {log}
        """

rule index_random_7:
    input:
        os.path.join(config["working_dir"], "vcfs/mikk_no-sibs_line-id_random-7.vcf.gz")
    output:
        os.path.join(config["working_dir"], "vcfs/mikk_no-sibs_line-id_random-7.vcf.gz.tbi")
    log:
        os.path.join(config["working_dir"], "logs/index_random_7.log")
    container:
        config["bcftools"]
    shell:
        """
        bcftools index \
            --tbi \
            {input[0]} \
                2> {log}
        """

rule pi_random_7:
    input:
        vcf = os.path.join(config["working_dir"], "vcfs/mikk_no-sibs_line-id.vcf.gz"),
        index = os.path.join(config["working_dir"], "vcfs/mikk_no-sibs_line-id.vcf.gz.tbi")
    output:
        os.path.join(config["lts_dir"], "nucleotide_divergence/mikk/random/500.windowed.pi")
    log:
        os.path.join(config["working_dir"], "logs/pi_random_7.log")
    params:
        prefix = lambda wildcards: os.path.join(config["lts_dir"], "nucleotide_divergence/mikk/random/500")
    container:
        config["vcftools"]
    shell:
        """
        vcftools \
            --gzvcf {input.vcf} \
            --keep {config[random_7_mikk]} \
            --window-pi 500000 \
            --out {params.prefix} \
                2> {log}
        """

# Pull mapping quality
rule mapping_quality:
    input:
        vcf = os.path.join(config["working_dir"], "vcfs/mikk_no-sibs_line-id.vcf.gz"),
        index = os.path.join(config["working_dir"], "vcfs/mikk_no-sibs_line-id.vcf.gz.tbi")
    output:
        os.path.join(config["lts_dir"], "mapping_quality/mapping_quality.csv")
    log:
        os.path.join(config["working_dir"], "logs/mapping_quality.log")
    container:
        config["bcftools"]
    shell:
        """
        bcftools query \
            --format '%CHROM,%POS,%INFO/MQ\\n' \
            --output {output[0]} \
            {input.vcf} \
                2> {log}
        """

# Split VCF per-sample
rule split_mikk:
    input:
        vcf = os.path.join(config["working_dir"], "vcfs/mikk_no-sibs_line-id.vcf.gz"),
        index = os.path.join(config["working_dir"], "vcfs/mikk_no-sibs_line-id.vcf.gz.tbi")
    output:
        expand(os.path.join(config["working_dir"], "vcfs/per_sample/mikk/{mikk_sample}.vcf.gz"),
                mikk_sample = MIKK_SAMPLES
        )
    log:
        os.path.join(config["working_dir"], "logs/split.log")
    params:
        out_dir = lambda wildcards: os.path.join(config["working_dir"], "vcfs/per_sample/mikk")
    container:
        config["bcftools"]
    shell:
        """
        bcftools +split \
            {input.vcf} \
            --output {params.out_dir} \
            --output-type z \
                2> {log}
        """

rule index_per_sample_mikk:
    input:
        os.path.join(config["working_dir"], "vcfs/per_sample/mikk/{mikk_sample}.vcf.gz")
    output:
        os.path.join(config["working_dir"], "vcfs/per_sample/mikk/{mikk_sample}.vcf.gz.tbi")
    log:
        os.path.join(config["working_dir"], "logs/index_per_sample_mikk/{mikk_sample}.log")
    container:
        config["bcftools"]
    shell:
        """
        bcftools index \
            --tbi \
            {input[0]} \
                2> {log}
        """

rule nd_per_sample_mikk:
    input:
        vcf = os.path.join(config["working_dir"], "vcfs/per_sample/mikk/{mikk_sample}.vcf.gz"),
    output:
        os.path.join(config["lts_dir"], "nucleotide_divergence/mikk/per_sample/{mikk_sample}.windowed.pi")
    log:
        os.path.join(config["working_dir"], "logs/nd_per_sample_mikk/{mikk_sample}.log")
    params:
        window_size = 500000,
        prefix = lambda wildcards, output: os.path.join(os.path.dirname(output[0]), wildcards.mikk_sample)
    container:
        config["vcftools"]
    shell:
        """
        vcftools \
            --gzvcf {input.vcf} \
            --window-pi {params.window_size} \
            --out {params.prefix} \
                2> {log}
        """

# Ewan: On the calculation of Pi, do you have a minimum count per allele (how many reads) or Quality of base calls if count of 1?
# Get various stats for single MIKK sample to determine cutoffs
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##  Info on GQ: https://gatk.broadinstitute.org/hc/en-us/articles/360035531692
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
rule get_test_stats:
    input:
        vcf = os.path.join(config["working_dir"], "vcfs/per_sample/mikk/10_1.vcf.gz"),
    output:
        os.path.join(config["working_dir"], "qc_stats/10_1.csv")
    log:
        os.path.join(config["working_dir"], "logs/get_test_stats.log")
    container:
        config["bcftools"]
    shell:
        """
        bcftools query \
            --format '%CHROM,%POS,%INFO/DP,%INFO/MQ,%INFO/QD[,%GQ]\\n' \
            --output {output[0]} \
            {input.vcf} \
                2> {log}
        """

rule filter_mikk:
    input:
        vcf = os.path.join(config["working_dir"], "vcfs/per_sample/mikk/{mikk_sample}.vcf.gz"),
        index = os.path.join(config["working_dir"], "vcfs/per_sample/mikk/{mikk_sample}.vcf.gz.tbi"),
    output:
        os.path.join(config["working_dir"], "vcfs/per_sample/filtered/mikk/{mikk_sample}.vcf.gz"),
    log:
        os.path.join(config["working_dir"], "logs/filter_mikk/{mikk_sample}.log")
    params:
        min_MQ = config["min_MQ"],
        min_DP = config["min_DP"],
        min_GQ = config["min_GQ"]
    container:
        config["bcftools"]
    shell:
        """
        bcftools view \
            --genotype ^miss \
            --include 'FORMAT/DP >= {params.min_DP} & MQ >= {params.min_MQ} & FORMAT/GQ >= {params.min_GQ}' \
            --output {output[0]} \
            --output-type z \
            {input.vcf} \
                2> {log}
        """

rule nd_filtered_mikk:
    input:
        vcf = os.path.join(config["working_dir"], "vcfs/per_sample/filtered/mikk/{mikk_sample}.vcf.gz"),
    output:
        os.path.join(config["lts_dir"], "nucleotide_divergence/mikk/filtered/{mikk_sample}.windowed.pi")
    log:
        os.path.join(config["working_dir"], "logs/nd_filtered_mikk/{mikk_sample}.log")
    params:
        window_size = 500000,
        prefix = lambda wildcards, output: os.path.join(os.path.dirname(output[0]), wildcards.mikk_sample)
    container:
        config["vcftools"]
    shell:
        """
        vcftools \
            --gzvcf {input.vcf} \
            --window-pi {params.window_size} \
            --out {params.prefix} \
                2> {log}
        """