## Filter VCF for wild Kiyosu samples (N = 7)
rule keep_no_sibs_kw:
    input:
        vcf = config["mikk_vcf"],
        samples_cram = config["kiyosu_samples_cram"]
    output:
        os.path.join(config["working_dir"], "vcfs/wild-kiyosu_cram-id.vcf.gz")
    log:
        os.path.join(config["working_dir"], "logs/keep_no_sibs_kw.log")
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
rule rehead_kw:
    input:
        vcf = os.path.join(config["working_dir"], "vcfs/wild-kiyosu_cram-id.vcf.gz"),
        samples = config["kiyosu_samples"]
    output:
        os.path.join(config["working_dir"], "vcfs/wild-kiyosu_line-id.vcf")
    log:
        os.path.join(config["working_dir"], "logs/rehead_kw.log")
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
rule compress_kw:
    input:
        os.path.join(config["working_dir"], "vcfs/wild-kiyosu_line-id.vcf")
    output:
        os.path.join(config["working_dir"], "vcfs/wild-kiyosu_line-id.vcf.gz")
    log:
        os.path.join(config["working_dir"], "logs/compress_kw.log")
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
rule index_kw:
    input:
        os.path.join(config["working_dir"], "vcfs/wild-kiyosu_line-id.vcf.gz")
    output:
        os.path.join(config["working_dir"], "vcfs/wild-kiyosu_line-id.vcf.gz.tbi")
    log:
        os.path.join(config["working_dir"], "logs/index_kw.log")
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
rule nucleotide_divergence_kw:
    input:
        vcf = os.path.join(config["working_dir"], "vcfs/wild-kiyosu_line-id.vcf.gz"),
        index = os.path.join(config["working_dir"], "vcfs/wild-kiyosu_line-id.vcf.gz.tbi")
    output:
        os.path.join(config["lts_dir"], "nucleotide_divergence/wild_kiyosu/{window_size}.windowed.pi")
    log:
        os.path.join(config["working_dir"], "logs/nucleotide_divergence_kw/{window_size}.log")
    params:
        window_size = lambda wildcards: wildcards.window_size,
        prefix = lambda wildcards: os.path.join(config["lts_dir"], "nucleotide_divergence/wild_kiyosu", wildcards.window_size)
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