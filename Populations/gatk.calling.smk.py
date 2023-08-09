MIN_REF_LENGTH=int(config["Other_paras"]["Minimum_contig_length"]) if "Minimum_contig_length" in config["Other_paras"] else 1000
NUMBER_OF_HC_INTERVALS= int(config["Other_paras"]["Intervals"]) if "Intervals" in config["Other_paras"] else 1

def get_fraction(sample):
    if "Norm_depth" in config["Other_paras"]:
        cov_obj =  float(config["Other_paras"]["Norm_depth"])
        files = "01.mapped_reads/"+sample+".stats"
        if os.path.exists(files):
            df = pd.read_table(files,header=None)
            cov_sample = float(df.iloc[0,1])
            fraction = cov_obj / cov_sample
            return(fraction)
        else:
            return(0)
    else:
        return(1)

selected = []
filtered = []

for item in samples:
    if get_fraction(item)>0 and get_fraction(item)<=1:
        selected.append(item)
    else:
        filtered.append(item)

def pick_samples(samples):
    picked = []
    for item in samples:
        if get_fraction(item) <= 1:
            picked.append(item)
    return(picked)

def reads_select():
    if "Norm_depth" in config["Other_paras"]: 
        return("01.mapped_reads/{sample}.subsampled.bam")
    else:
        return("01.mapped_reads/{sample}.realigned.bam")

def aggregate_input(wildcards):
    intervalDir = checkpoints.scatter_intervals.get(**wildcards).output[0]
    return expand("02.gvcfs/{{sample}}.interval_{interval}.g.vcf.gz", 
                 interval=glob_wildcards(intervalDir + "/interval_{match}.list").match) 

def aggregate_vcf_input(wildcards):
    intervalDir = checkpoints.scatter_intervals.get(**wildcards).output[0]
    return expand("03.vcfs/{{sample}}.interval_{interval}.s.vcf.gz",
                 interval=glob_wildcards(intervalDir + "/interval_{match,\d+}.list").match)

if "Norm_depth" in config["Other_paras"]:
    rule subsample:
        input:
            bam="01.mapped_reads/{sample}.realigned.bam",
            stats="01.mapped_reads/{sample}.stats"
        output:
            "01.mapped_reads/{sample}.subsampled.bam",
            "01.mapped_reads/{sample}.subsampled.bam.bai"
        run:
            fraction = get_fraction(wildcards.sample)
            if fraction <= 1:
                shell("samtools view -@ 20 -s {fraction} -b {input.bam} -o {output[0]}")
                shell("samtools index -@ 20 {output[0]}")
                print(f"{wildcards.sample} has been subsampled at {fraction}")
                shell("echo {wildcards.sample} >>01.mapped_reads/samples.txt")
            else:
                shell("touch {output[0]}")
                shell("touch {output[1]}")
                print(f"Work around {wildcards.sample}")
                shell("echo {wildcards.sample} >>01.mapped_reads/remove.txt")

checkpoint scatter_intervals:
    input:
        genome=RefGenome
    output:
        intervals=directory("interval_lists")
    params:
        cutoff=MIN_REF_LENGTH,
        intervals=NUMBER_OF_HC_INTERVALS
    resources:
        mem_gb=1 # will only use a few Mb; float not accepted here

    shell:
        "mkdir interval_lists\n"
        "scripts/createIntervalLists.py {input} interval_lists {params.cutoff} {params.intervals} "


rule haplotypeCaller:
    # the HaplotypeCaller uses several threads, causing a peak load at
    # the beginning and then falling back to approx. 2 threads for the 
    # rest of the run time
    input:
        bam=reads_select(),
        genome=RefGenome,
        intervals="interval_lists/{interval}.list"
    output:
        temp("02.gvcfs/{sample}.{interval}.g.vcf.gz")
    priority:
        8
    params:
        ploidy = config["Other_paras"]["Ploidy"],
        # no option (= empty string) or append path to option
        pairHMM= "" if not NUMBER_OF_HC_HMM_THREADS else "--native-pair-hmm-threads " +
                      str(NUMBER_OF_HC_HMM_THREADS),
    threads:
        #2 # can be influenced through the option --native-pair-hmm-threads
        10
    resources:
        tmpdir="/BioData/tmp",
        mem_gb=96
    run:
        if "Norm_depth" in config["Other_paras"]:
            if get_fraction(wildcards.sample) <= 1:
                shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\" HaplotypeCaller -ERC GVCF \
                       -I {input.bam} -R {input.genome} --annotate-with-num-discovered-alleles true \
                       -ploidy {params.ploidy} -O {output} -L {input.intervals} {params.pairHMM} --tmp-dir {resources.tmpdir}")
            else:
                shell("touch {output}")
        else:
            shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\" HaplotypeCaller -ERC GVCF \
                -I {input.bam} -R {input.genome} --annotate-with-num-discovered-alleles true \
                -ploidy {params.ploidy} -O {output} -L {input.intervals} {params.pairHMM} --tmp-dir {resources.tmpdir}")

rule run_genotyping:
    input:
        ref=RefGenome,
        vcfs="02.gvcfs/{sample}.interval_{int}.g.vcf.gz"
    output:
        temp("03.vcfs/{sample}.interval_{int}.s.vcf.gz")
    params:
        ploidy = config["Other_paras"]["Ploidy"]
    threads: 5
    resources:
        tmpdir="/BioData/tmp",
        mem_gb=128
    run:
        if "Norm_depth" in config["Other_paras"]:
            if get_fraction(wildcards.sample) <= 1:
                shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=4\" GenotypeGVCFs \
                      --annotate-with-num-discovered-alleles true -R {input.ref} -V {input.vcfs} -O {output} -stand-call-conf 30 -ploidy {params.ploidy}")
            else:
                shell("touch {output}")
        else:
            shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=4\" GenotypeGVCFs \
                --annotate-with-num-discovered-alleles true -R {input.ref} -V {input.vcfs} -O {output} -stand-call-conf 30 -ploidy {params.ploidy}")

rule gather_gvcf_intervals:
    input:
        inGVCFs=aggregate_input
    output:
        "02.gvcfs/{sample}.gvcf.gz"
    priority:
        5
    threads: 4
    params:
        options=lambda wildcards, input: " ".join(
                    [ "-I "+x for x in sorted(
                            input.inGVCFs,
                            key=lambda key: [ int(s) if s.isdigit() else s for s in re.split('([0-9]+)', key) ]
                        )
                    ]
                ),
        intervalCount=lambda wildcards, input: len(input.inGVCFs),
        # no option (= empty string) or append path to option
    resources:
        tmpdir="/BioData/tmp",
        mem_gb=64
    run:
        if "Norm_depth" in config["Other_paras"]:
            if get_fraction(wildcards.sample) <= 1:
                if params.intervalCount >1:
                    shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\"\
                            GatherVcfs -O {output} {params.options} --TMP_DIR {resources.tmpdir}")
                    shell("tabix -p vcf {output}")
                else:
                    shell("ln {input.inGVCFs[0]} {output}")
                    shell("ln {input.inGVCFs[0]}.tbi {output}.tbi ")
            else:
                shell("touch {output}")
        else:
            if params.intervalCount >1:
                shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\"\
                        GatherVcfs -O {output} {params.options} --TMP_DIR {resources.tmpdir}")
                shell("tabix -p vcf {output}")
            else:
                shell("ln -s {input.inGVCFs[0]} {output}")
                shell("ln -s {input.inGVCFs[0]}.tbi {output}.tbi ")

rule gather_vcf_intervals:
    input:
        inGVCFs=aggregate_vcf_input
    output:
        "03.vcfs/{sample}.vcf.gz"
    priority:
        5
    params:
        options=lambda wildcards, input: " ".join(
                    [ "-I "+x for x in sorted(
                            input.inGVCFs,
                            key=lambda key: [ int(s) if s.isdigit() else s for s in re.split('([0-9]+)', key) ]
                        )
                    ]
                ),
        intervalCount=lambda wildcards, input: len(input.inGVCFs),
    resources:
        tmpdir="/BioData/tmp",
        mem_gb=64
    threads: 2
    run:
        if "Norm_depth" in config["Other_paras"]:
            if get_fraction(wildcards.sample) <= 1:
                if params.intervalCount >1:
                    shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\"\
                            GatherVcfs -O {output} {params.options} --TMP_DIR {resources.tmpdir}")
                    shell("tabix -p vcf {output}")
                else:
                    shell("ln {input.inGVCFs[0]} {output}")
                    shell("ln {input.inGVCFs[0]}.tbi {output}.tbi ")
            else:
                shell("touch {output}")
        else:
            if params.intervalCount >1:
                shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=2\"\
                        GatherVcfs -O {output} {params.options} --TMP_DIR {resources.tmpdir}")
                shell("tabix -p vcf {output}")
            else:
                shell("ln {input.inGVCFs[0]} {output}")
                shell("ln {input.inGVCFs[0]}.tbi {output}.tbi ")

rule combine_calls:
    input:
        ref=RefGenome,
        interval = "00.data/interval.list",
        gvcfs= expand("02.gvcfs/{sample}.gvcf.gz", sample =  (lambda : samples if "Norm_depth" not in config["Other_paras"] else pick_samples(samples))())
    output:
        directory("10.gvcfDB")
    resources:
        tmpdir="/BioData/tmp",
        mem_gb=512
    threads: 2
    params:
        gvcfs= expand("-V 02.gvcfs/{sample}.gvcf.gz", sample =  (lambda : samples if "Norm_depth" not in config["Other_paras"] else pick_samples(samples))())
    run:
        shell("gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=4\" GenomicsDBImport {params.gvcfs} -L {input.interval} --genomicsdb-workspace-path {output} --batch-size 50 --tmp-dir {resources.tmpdir} --reader-threads 20 --max-num-intervals-to-import-in-parallel 15")


rule genotype_variants:
    input:
        vcfs = "10.gvcfDB",
        ref="00.data/ref.fa"
    params:
        ploidy = config["Other_paras"]["Ploidy"]
    resources:
        tmpdir="/BioData/tmp",
        mem_gb=512
    threads: 2
    output:
        "03.vcfs/all.genotyped.vcf.gz"
    shell:
        "gatk4 --java-options \"-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads=4\" GenotypeGVCFs --tmp-dir {resources.tmpdir} --annotate-with-num-discovered-alleles true -R {input.ref} -V gendb://{input.vcfs} -O {output} -stand-call-conf 30 -ploidy {params.ploidy}"

