rule bwa_map_deduplication:
    message:
        "mapping reads to the reference"
    input:
        ref="data/ref.fa",
        reads="data/sample/{sample}.unaligned.bam"
    output:
        temp("mapped_reads/{sample}.derep.bam")
    shell:
        "gatk4 --java-options \"-Xmx690G -XX:+UseParallelGC -XX:ParallelGCThreads=4”\" \
                BwaAndMarkDuplicatesPipelineSpark \
                -- --spark-runner LOCAL --spark-master local[*] \
                --conf spark.executor.cores=4 --conf spark.executor.memory=27g --conf spark.num.executors=25 \
                -I {input.ref} -R {input.reads} -O {output}"

rule sort_bam
    message:
        "Making the bam sorted"
    input:
        rules.bwa_map_deduplication.output
    output:
        temp("mapped_reads/{sample}.sorted.bam")
    shell:
        "gatk4 --java-options \"-Xmx690G -XX:+UseParallelGC -XX:ParallelGCThreads=4”\" \
                SortSamSpark \
                -- --spark-runner LOCAL --spark-master local[*] \
                --conf spark.executor.cores=4 --conf spark.executor.memory=27g --conf spark.num.executors=25 \
                -I {input} -O {output}"

rule indel_refinement_pre
    message:
        "Refining the bams in the intdel areas"
    input:
        ref="data/ref.fa",
        bam=rules.sort_bam.output
    output:
        temp("mapped_reads/{sample}.intervals")
    shell:
        "gatk3 -T RealignerTargetCreator -R {input.ref} -I {input.bam} --minReadsAtLocus 4 -o {output}"

rule indel_refinement
    input:
        ref="data/ref.fa",
        bam=rules.sort_bam.output,
        interval=rules.indel_refinement_pre.output
    output:
        "mapped_reads/{sample}.realigned.bam"
    shell:
        "gatk3 -T IndelRealigner -R {input.ref} -I {input.bam} -targetIntervals {input.interval} -LOD 3.0 -o {output}"

rule all
    input:
        interval=rules.indel_refinement_pre.output
