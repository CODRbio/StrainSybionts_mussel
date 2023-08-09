if "Group" in config["Other_paras"]:
    rule PCA_snp_drawing:
        input:
            "04.filtering/snps/all.genotyped.snp.filtered.vcf.gz",
            config["Other_paras"]["Group"]
        output:
            "05.statistics/PCA.snp.pdf"
        script:
            "scripts/variationPCA_group.R"


    rule PCA_merge_drawing:
        input:
            "04.filtering/merged/all.genotyped.filtered.vcf.gz",
            config["Other_paras"]["Group"]
        output:
            "05.statistics/PCA.all.pdf"
        script:
            "scripts/variationPCA_group.R"

    rule PCoA_merge_drawing:
        input:
            "04.filtering/merged/all.genotyped.filtered.vcf.gz",
            config["Other_paras"]["Group"]
        output:
            "05.statistics/PCoA.all.pdf"
        script:
            "scripts/variationPCoA_group.R"

    rule PCoA_snp_drawing:
        input:
            "04.filtering/snps/all.genotyped.snp.filtered.vcf.gz",
            config["Other_paras"]["Group"]
        output:
            "05.statistics/PCoA.snp.pdf"
        script:
            "scripts/variationPCoA_group.R"
else:
    rule PCA_snp_drawing:
        input:
            "04.filtering/snps/all.genotyped.snp.filtered.vcf.gz"
        output:
            "05.statistics/PCA.snp.pdf"
        script:
            "scripts/variationPCA.R"


    rule PCA_merge_drawing:
        input:
            "04.filtering/merged/all.genotyped.filtered.vcf.gz"
        output:
            "05.statistics/PCA.all.pdf"
        script:
            "scripts/variationPCA.R"

    rule PCoA_merge_drawing:
        input:
            "04.filtering/merged/all.genotyped.filtered.vcf.gz"
        output:
            "05.statistics/PCoA.all.pdf"
        script:
            "scripts/variationPCoA.R"

    rule PCoA_snp_drawing:
        input:
            "04.filtering/snps/all.genotyped.snp.filtered.vcf.gz"
        output:
            "05.statistics/PCoA.snp.pdf"
        script:
            "scripts/variationPCoA.R"

rule allele_merged_pre:
    input:
        "04.filtering/merged/all.genotyped.filtered.vcf.gz"
    output:
        temp("05.statistics/all_genotype.tsv")
    shell:
        "gatk4 VariantsToTable -V {input} -F CHROM -F POS -F TYPE -GF AD -GF DP -O {output}"

rule allele_merged_frequency:
    input:
        "05.statistics/all_genotype.tsv"
    params:
        len=RefLength
    output:
        "05.statistics/allele_all_frequency.csv",
        "05.statistics/variation_density.csv"
    shell:
        "scripts/Allele_fre.py {params.len} {input} {output[0]} {output[1]}"


rule allele_snp_pre:
    input:
        "04.filtering/snps/all.genotyped.snp.filtered.vcf.gz"
    output:
        temp("05.statistics/snp_genotype.tsv")
    shell:
        "gatk4 VariantsToTable -V {input} -F CHROM -F POS -F TYPE -GF AD -GF DP -O {output}"

rule allele_snp_frequency:
    input:
        "05.statistics/snp_genotype.tsv"
    params:
        len=RefLength
    output:
        "05.statistics/allele_snp_frequency.csv",
        "05.statistics/snp_density.csv"
    shell:
        "scripts/Allele_fre.py {params.len} {input} {output[0]} {output[1]}"

rule vcf2count:
    input:
        "04.filtering/snps/all.genotyped.snp.filtered.vcf.gz"
    output:
        "05.statistics/all.snv.tsv"
    shell:
        "gatk4 VariantsToTable -V {input} -F CHROM -F POS -F REF -F ALT -GF AD -O {output}"

rule countTable:
    input:
        rules.vcf2count.output
    output:
        expand("05.statistics/snv_count/{sample}.snv.count.csv",sample=selected),
        temp("05.statistics/sample.list")
    run:
        shell("scripts/vcf_counts.py {input}")
        with open("05.statistics/sample.list",'w') as out:
            for item in selected:
                out.write(f"{item}\t05.statistics/snv_count/{item}.snv.count.csv\n")

rule countTable2:
    output:
        expand("05.statistics/snv_count/{sample}.snv.count.csv",sample=filtered)
    run:
        for item in filtered:
            shell("touch 05.statistics/snv_count/{item}.snv.count.csv")

rule Fst_gen:
    input:
        "05.statistics/sample.list"
    output:
        "05.statistics/sample.list_Fst_pos.txt"
    shell:
        "scripts/structure.r 05.statistics/sample.list"

rule Fst_whole:
    input:
        list="05.statistics/sample.list_Fst_pos.txt",
        ref="00.data/ref.fa"
    output:
        "05.statistics/Fst_fst.mat",
        "05.statistics/Fst_fst.txt",
        "05.statistics/Fst_pi.txt"
    shell:
        "scripts/genome_wise_caculation.py {input.list} {input.ref} Fst"

rule Fst_pdf:
    input:
        "05.statistics/Fst_fst.mat"
    output:
        "05.statistics/FstHeatmap.pdf"
    shell:
        "scripts/heatmap_draw.R -m {input} -R -C -o {output}"

rule pn_ps:
    input:
        ref="00.data/ref.fa",
        tab="05.statistics/all.snv.tsv"
    output:
        "05.statistics/pn_ps_genes.tsv",
        "05.statistics/pn_ps_genome.tsv"
    shell:
        "scripts/pn_ps_caculation.py {input.tab} {input.ref}" 
