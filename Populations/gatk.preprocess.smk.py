rule pre_propossing: 
    input:
        fq1 = "00.data/{sample}.fq1.gz",
        fq2 = "00.data/{sample}.fq2.gz",
        ref = "00.data/ref.fa"
    output:
        bam="00.data/{sample}.unaligned.bam",
    conda:
        "/nfs_genome/anaconda/envs/population_gatk"
    resources:
        tmpdir="/BioData/tmp",
        mem_gb=64
    threads: 5
    shell:
        "set +eu "
        " && PS1=dummy"
        " && . $(conda info --base)/etc/profile.d/conda.sh"
        " && conda activate population_gatk;"
        "/nfs_genome/anaconda/envs/population_gatk/bin/gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\" FastqToSam F1={input.fq1} F2={input.fq2} O={output.bam} SM={wildcards.sample} RG=rg0013 PL=ILLUMINA"
        

rule interval_gen:
    input:
        ref = "00.data/ref.fa"
    output:
        lis="00.data/interval.list"
    run:
        shell("grep \">\" ref.fa | perl -ne 'chomp; my @info=split(/\s+/,); print \"$info[0]\n\"' | sed \"s/>//\" >{output.lis}")
