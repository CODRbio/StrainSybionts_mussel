rule get_index:
    input:
        "00.data/ref.fa"
    output:
        "00.data/ref.fa.fai"
    shell:
        "set +eu "
        " && PS1=dummy"
        " && . $(conda info --base)/etc/profile.d/conda.sh"
        " && conda activate population_gatk;"
        "samtools faidx {input}"


rule bam_index:
    input:
        "00.data/ref.fa"
    output:
        "00.data/ref.fa.img"
    resources:
        tmpdir="/BioData/tmp",
        mem_gb=100
    shell:
        "set +eu "
        " && PS1=dummy"
        " && . $(conda info --base)/etc/profile.d/conda.sh"
        " && conda activate population_gatk;"
        "/nfs_genome/anaconda/envs/population_gatk/bin/gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=4\" BwaMemIndexImageCreator -I {input}"

rule dict_building:
    input:
        "00.data/ref.fa"
    output:
        "00.data/ref.dict"
    shell:
        "set +eu "
        " && PS1=dummy"
        " && . $(conda info --base)/etc/profile.d/conda.sh"
        " && conda activate population_gatk;"
        "/nfs_genome/anaconda/envs/population_gatk/bin/gatk4 CreateSequenceDictionary R={input}"


rule bwa_map_deduplication:
    message:
        "mapping reads to the reference"
    input:
        ref="00.data/ref.fa",
        reads="00.data/{sample}.unaligned.bam",
        fai="00.data/ref.fa.fai",
        img="00.data/ref.fa.img",
        dict="00.data/ref.dict"
    output:
        temp("01.mapped_reads/{sample}.derep.bam"),
        temp("01.mapped_reads/{sample}.derep.bam.sbi")
    priority:
        9
    resources:
        tmpdir="/BioData/tmp",
        mem_gb=800
    threads: 64
    shell:
        "set +eu "
        " && PS1=dummy"
        " && . $(conda info --base)/etc/profile.d/conda.sh"
        " && conda activate population_gatk;"
        "/nfs_genome/anaconda/envs/population_gatk/bin/gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\" \
                BwaAndMarkDuplicatesPipelineSpark \
                -- --spark-runner LOCAL --spark-master local[*] \
                --conf spark.executor.cores=8 --conf spark.executor.memory=90g --conf spark.num.executors=8 \
                -R {input.ref} --tmp-dir {resources.tmpdir} -I {input.reads} -O {output[0]}"

rule sort_bam:
    message:
        "Making the bam sorted"
    input:
        rules.bwa_map_deduplication.output[0]
    output:
        temp("01.mapped_reads/{sample}.sorted.bam"),
    resources:
        tmpdir="/BioData/tmp",
        mem_gb=800
    threads: 60
    shell:
        "set +eu "
        " && PS1=dummy"
        " && . $(conda info --base)/etc/profile.d/conda.sh"
        " && conda activate population_gatk;"
        "/nfs_genome/anaconda/envs/population_gatk/bin/gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\" \
                SortSamSpark \
                -- --spark-runner LOCAL --spark-master local[*] \
                --conf spark.executor.cores=6 --conf spark.executor.memory=90g --conf spark.num.executors=10 \
                -I {input[0]} -O {output[0]} --tmp-dir {resources.tmpdir}"

rule indel_refinement_pre:
    message:
        "Refining the bams in the intdel areas"
    input:
        ref="00.data/ref.fa",
        bam=rules.sort_bam.output[0]
    threads: 4
    output:
        temp("01.mapped_reads/{sample}.intervals")
    shell:
        "export _JAVA_OPTIONS=-Xmx128G\n"
        "gatk3 -T RealignerTargetCreator -R {input.ref} -I {input.bam} --minReadsAtLocus 4 -o {output}"

rule indel_refinement:
    input:
        ref="00.data/ref.fa",
        bam=rules.sort_bam.output[0],
        interval=rules.indel_refinement_pre.output
    threads: 4
    output:
        "01.mapped_reads/{sample}.realigned.bam"
    shell:
        "export _JAVA_OPTIONS=-Xmx128G\n"
        "gatk3 -T IndelRealigner -R {input.ref} -I {input.bam} -targetIntervals {input.interval} -LOD 3.0 -o {output}"

rule sample_coverage:
    input:
        rules.indel_refinement.output
        #"mapped_reads/{sample}.realigned.bam"
    threads: 2
    output:
        "01.mapped_reads/{sample}.stats"
    shell:
        """samtools depth {input} | awk '{{sum+=$3;cnt++}} END {{print "{wildcards.sample}""	"sum/cnt}}' > {output}"""
rule average_coverage:
    input:
        expand("01.mapped_reads/{sample}.stats", sample=config["samples"])
    output:
        "01.mapped_reads/average_cov.stats"
    script:
        "scripts/cov_mean.py"
