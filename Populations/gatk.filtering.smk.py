rule split_SNPs:
    input:
        "03.vcfs/{sample}.vcf.gz",
    output:
        temp("04.filtering/snps/{sample}.snp.vcf.gz")
    resources:
        tmpdir="/BioData/tmp",
        mem_gb=64
    threads: 2
    run:
        if "Norm_depth" in config["Other_paras"]:
            if get_fraction(wildcards.sample) <= 1:
                shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\" \
                SelectVariants  -V {input} -select-type SNP -O {output} --tmp-dir {resources.tmpdir}")
            else:
                shell("touch {output}")
        else:
            shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\" \
                SelectVariants  -V {input} -select-type SNP -O {output} --tmp-dir {resources.tmpdir}")

rule split_Indels:
    input:
        "03.vcfs/{sample}.vcf.gz",
    output:
        temp("04.filtering/indels/{sample}.indel.vcf.gz")
    resources:
        tmpdir="/BioData/tmp",
        mem_gb=64
    threads: 2
    run:
        if "Norm_depth" in config["Other_paras"]:
            if get_fraction(wildcards.sample) <= 1:
                shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\" \
                SelectVariants  -V {input} -select-type INDEL -select-type MIXED -O {output} --tmp-dir {resources.tmpdir}")
            else:
                shell("touch {output}")
        else:
            shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\" \
                SelectVariants  -V {input} -select-type INDEL -select-type MIXED -O {output} --tmp-dir {resources.tmpdir}")

rule SNP_filtering:
    input:
        rules.split_SNPs.output
    output:
        "04.filtering/snps/{sample}.snp.filtered.vcf.gz"
    resources:
        tmpdir="/BioData/tmp",
        mem_gb=64
    threads: 2
    run:
        if "Norm_depth" in config["Other_paras"]:
            if get_fraction(wildcards.sample) <= 1:
                shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\" \
                VariantFiltration -V {input} -O {output} -filter 'QD <2.0' --filter-name 'QD2' --tmp-dir {resources.tmpdir}\
                -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'MQ < 40.0' --filter-name 'MQ40'\
                -filter 'FS > 60.0' --filter-name 'FS60' -filter 'MQRankSum < -20.0' --filter-name 'MQRankSum-20'\
                -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8' -filter 'SOR > 3.0' --filter-name 'SOR3'")
            else:
                shell("touch {output}")
        else:
            shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\" \
                VariantFiltration -V {input} -O {output} -filter 'QD <2.0' --filter-name 'QD2' --tmp-dir {resources.tmpdir}\
                -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'MQ < 40.0' --filter-name 'MQ40'\
                -filter 'FS > 60.0' --filter-name 'FS60' -filter 'MQRankSum < -20.0' --filter-name 'MQRankSum-20'\
                -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8' -filter 'SOR > 3.0' --filter-name 'SOR3'")

rule Indel_filtering:
    input:
        rules.split_Indels.output
    output:
        "04.filtering/indels/{sample}.indel.filtered.vcf.gz"
    resources:
        tmpdir="/BioData/tmp",
        mem_gb=64
    threads: 2
    run:
        if "Norm_depth" in config["Other_paras"]:
            if get_fraction(wildcards.sample) <= 1:
                shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\" \
                VariantFiltration -V {input} -O {output} -filter 'QD < 2.0' --filter-name 'QD2' \
                -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'FS > 200.0' --filter-name 'FS200' \
                -filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' --tmp-dir {resources.tmpdir}")
            else:
                shell("touch {output}")
        else:
            shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\" \
                VariantFiltration -V {input} -O {output} -filter 'QD < 2.0' --filter-name 'QD2' \
                -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'FS > 200.0' --filter-name 'FS200' \
                -filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' --tmp-dir {resources.tmpdir}")

rule merge_filtered_VCF:
    input:
        snp=rules.SNP_filtering.output,
        indel=rules.Indel_filtering.output
    output:
        temp("04.filtering/merged/{sample}.sorted.vcf.gz")
    resources:
        tmpdir="/BioData/tmp",
        mem_gb=64
    threads: 2
    run:
        if "Norm_depth" in config["Other_paras"]:
            if get_fraction(wildcards.sample) <= 1:
                shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\" \
                SortVcf -I {input.snp} -I {input.indel} -O {output}")
            else:
                shell("touch {output}")
        else:
            shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\" \
                SortVcf -I {input.snp} -I {input.indel} -O {output}")

rule filtered_removed_VCF:
    input:
        rules.merge_filtered_VCF.output
    output:
        "04.filtering/merged/{sample}.filtered.vcf.gz"
    resources:
        tmpdir="/BioData/tmp",
        mem_gb=64
    threads: 2
    run:
        if "Norm_depth" in config["Other_paras"]:
            if get_fraction(wildcards.sample) <= 1:
                shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\" \
                SelectVariants -V {input} -O {output} --exclude-filtered true --tmp-dir {resources.tmpdir}")
            else:
                shell("touch {output}")
        else:
            shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\" \
                SelectVariants -V {input} -O {output} --exclude-filtered true --tmp-dir {resources.tmpdir}")


