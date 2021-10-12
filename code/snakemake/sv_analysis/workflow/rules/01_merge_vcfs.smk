# Move polished VCFs from EBI cluster to Codon
rule copy_polished_vcfs:
    input:
        os.path.join(config["polished_on_ebi"], "{sample}.vcf")
    output:
        os.path.join(config["working_dir"], "vcfs/indiv/{sample}.vcf")
    log:
        os.path.join(config["working_dir"], "logs/copy_polished_vcfs/{sample}.log")
    shell:
        """
        cp {input[0]} {output[0]} \
            2> {log}
        """

rule make_new_headers:
    input:
        os.path.join(config["working_dir"], "vcfs/indiv/{sample}.vcf")
    output:
        os.path.join(config["working_dir"], "vcfs/new_headers/{sample}.txt")
    log:
        os.path.join(config["working_dir"], "logs/make_new_headers/{sample}.log")
    container:
        config["bcftools"]
    shell:
        """
        cp {config[new_header]} {output[0]} ;
        bcftools view --header {input[0]} | grep -v "^##" >> {output[0]}
        """

rule rehead_vcfs:
    input:
        vcf = os.path.join(config["working_dir"], "vcfs/indiv/{sample}.vcf"),
        header = os.path.join(config["working_dir"], "vcfs/new_headers/{sample}.txt")
    output:
        os.path.join(config["working_dir"], "vcfs/indiv_rehead/{sample}.vcf")
    log:
        os.path.join(config["working_dir"], "logs/rehead_vcfs/{sample}.log")
    container:
        config["bcftools"]
    shell:
        """
        bcftools reheader \
            --header {input.header} \
            --output {output[0]} \
            {input.vcf} \
                2> {log}
        """

rule compress_vcfs:
    input:
        os.path.join(config["working_dir"], "vcfs/indiv_rehead/{sample}.vcf")
    output:
        os.path.join(config["working_dir"], "vcfs/indiv_rehead/{sample}.vcf.gz")
    log:
        os.path.join(config["working_dir"], "logs/compress_vcfs/{sample}.log")
    container:
        config["bcftools"]
    shell:
        """
        bcftools view \
            --output {output[0]} \
            --output-type z \
            {input} \
                2> {log}
        """

rule sort_vcfs:
    input:
        os.path.join(config["working_dir"], "vcfs/indiv_rehead/{sample}.vcf.gz")
    output:
        os.path.join(config["working_dir"], "vcfs/indiv_sorted/{sample}.vcf.gz")
    log:
        os.path.join(config["working_dir"], "logs/index_vcfs/{sample}.log")
    container:
        config["bcftools"]
    shell:
        """
        bcftools sort \
            --output {output[0]} \
            --output-type z \
            {input[0]} \
                2> {log}
        """

rule index_vcfs:
    input:
        os.path.join(config["working_dir"], "vcfs/indiv_sorted/{sample}.vcf.gz")
    output:
        os.path.join(config["working_dir"], "vcfs/indiv_sorted/{sample}.vcf.gz.csi")
    log:
        os.path.join(config["working_dir"], "logs/index_vcfs/{sample}.log")
    container:
        config["bcftools"]
    shell:
        """
        bcftools index \
            --csi \
            {input[0]} \
                2> {log}
        """

rule merge_vcfs:
    input:
        vcfs = expand(os.path.join(config["working_dir"], "vcfs/indiv_sorted/{sample}.vcf.gz"),
                        sample = config["samples"]),
        indexes = expand(os.path.join(config["working_dir"], "vcfs/indiv_sorted/{sample}.vcf.gz.csi"),
                        sample = config["samples"])
    output:
        os.path.join(config["working_dir"], "vcfs/merged/all.vcf")
    log:
        os.path.join(config["working_dir"], "logs/merge_vcfs.log")
    container:
        config["bcftools"]
    shell:
        """
        bcftools merge \
            --output {output[0]} \
            --output-type z \
            {input} \
                2> {log}
        """

