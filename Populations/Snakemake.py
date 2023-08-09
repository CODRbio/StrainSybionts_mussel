configfile: "config.yaml"
samples=list(config["samples"].keys())
RefGenome="00.data/ref.fa"
NUMBER_OF_HC_HMM_THREADS = 4
RefLength=config["Ref_length_kb"]

import pandas as pd

rule all:
    input:
        "00.data/ref.fa.fai",
        "00.data/ref.fa.img",
        "00.data/ref.dict",
        "01.mapped_reads/average_cov.stats",
        expand("01.mapped_reads/{sample}.stats", sample=samples),
        "03.vcfs/all.genotyped.vcf.gz",
        expand("03.vcfs/{sample}.vcf.gz", sample=samples),
        expand("04.filtering/merged/{sample}.filtered.vcf.gz", sample=samples),
        "04.filtering/merged/all.genotyped.filtered.vcf.gz",
        "05.statistics/all.snv.tsv"
        "05.statistics/PCA.snp.pdf",
        "05.statistics/PCA.all.pdf",
        "05.statistics/PCoA.all.pdf",
        "05.statistics/PCoA.snp.pdf",
        "05.statistics/allele_all_frequency.csv",
        "05.statistics/allele_snp_frequency.csv",
        "05.statistics/snp_density.csv",
        "05.statistics/variation_density.csv",
        expand("05.statistics/snv_count/{sample}.snv.count.csv",sample=samples),
        "05.statistics/FstHeatmap.pdf",
        "05.statistics/all_pn_ps_genes.tsv",
        "05.statistics/all_pn_ps_genome.tsv",
        "06.desman/desmanInput_sel_var.csv",
        "06.desman/res/Dev.csv",
        "06.desman/res/Dev.pdf",
        "06.desman/res/best_genotype.csv",
        lambda wildcards: expand("06.desman/desmanResN_{g}_{r}",
                                 g= config["desman_nr_genomes_range"],
                                 r=config["desman_nr_variable_positions"])



include: "gatk.preprocess.smk.py"        
include: "gatk.bammapping.smk.py"
include: "gatk.calling.smk.py"
include: "gatk.filtering.smk.py"
include: "gatk.statistics.py" 

def count_samples(samples):
    picked = 0
    for item in samples:
        if get_fraction(item) <= 1:
            picked = picked+1
    return(picked)

counts=count_samples(samples)
ends=3+4*counts

rule Desman_raw:
    input:
        cnt="05.statistics/sample.list_Fst_pos.txt"
    output:
        tmp="06.desman/desman_tmp.csv",
        all="06.desman/desman_all.csv"
    run:
        import pandas as pd
        shell("cp {input} {output.tmp}")
        shell("sed -i \"1s/\.\([ATGC]\)\(\s\)/-\\1\\2/g\" {output.tmp}") 
        print(ends)
        df = pd.read_table(output.tmp , header=0 ,sep='\t')
        df2 = df.iloc[:,[0,1]+list(range(3,ends))]
        df2.to_csv(output.all,index=False)

rule Desman_core_data:
    input:
        core="00.data/Core_all.tsv",
        allD=rules.Desman_raw.output.all
    output:
        csv="06.desman/desman_core.csv"
    run:
        ind = pd.read_table(input.core,header=None)[0].to_list()
        df = pd.read_csv(input.allD,header=0,index_col=0)
        ind_use = list(set(ind) & set(df.index))  
        df_sel = df.loc[ind_use]
        print(df_sel)
        df_sel.to_csv(output.csv,index=True)

rule Desman_filter:
    input:
        coreV=rules.Desman_core_data.output.csv
    output:
        expand("06.desman/desmanInput_{output_file_type}",
                output_file_type = ["sel_var.csv", "p_df.csv", "q_df.csv", "r_df.csv", "tran_df.csv", "log.txt"])
    shell:
        "/nfs_genome/anaconda/envs/desman/bin/Variant_Filter.py --output_stub 06.desman/desmanInput_ {input}"

rule Desman_runner:
    input:
        sel_var = "06.desman/desmanInput_sel_var.csv",
        tran_df = "06.desman/desmanInput_tran_df.csv"
    output:
        desman_results = "06.desman/desmanResN_{g}_{r}/log_file.txt",
        fit = "06.desman/desmanResN_{g}_{r}/fit.txt",
        desman_log = "06.desman/desmanResN_{g}_{r}.log",
        dir = directory("06.desman/desmanResN_{g}_{r}")
    resources:
        tmpdir="/BioData/tmp",
        mem_gb=600
    threads: 60
    priority: 50
    run:
        output_dir = os.path.dirname(output.desman_results)
        shell("""
            set +u; conda activate desman; set -u;
            /nfs_genome/anaconda/envs/desman/bin/desman {input.sel_var} -e {input.tran_df} -o {output_dir} -r 1000 -i 100 -g {wildcards.g} -s {wildcards.r} > {output.desman_log} 2>&1
        """)

rule run_all_desman:
    input: 
        expand("06.desman/desmanResN_{g}_{r}/fit.txt",
                                 g= config["desman_nr_genomes_range"],
                                 r=config["desman_nr_variable_positions"])
    output: 
        "06.desman/res/Dev.csv"
    resources:
        tmpdir="/BioData/tmp",
        mem_gb=100
    threads: 10
    priority: 50
    shell:
        """cat <(echo 'H,G,LP,Dev') <(cat {input} | cut -d',' -f2-) > {output}"""

rule desman_fit_pdf:
    input:
        rules.run_all_desman.output
    output: "06.desman/res/Dev.pdf"
    shell:
        "/nfs_genome/anaconda/envs/desman/bin/PlotDev.R -l {input} -o {output}"

rule find_best:
    input:
        expand("06.desman/desmanResN_{g}_{r}",g = config["desman_nr_genomes_range"], r = config["desman_nr_variable_positions"])
    output:
        "06.desman/res/best_genotype.csv"
    run:
        shell("/nfs_genome/anaconda/envs/desman/bin/resolvenhap.py 06.desman/desmanResN >{output}")
