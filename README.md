# StrainSybionts_mussel
 Strain analysis of symbionts in the deep-sea mussel, v2.0
Strain-level assembly of Symbionts in mussels from vents and cold seeps in the Northwest Pacific Ocean

This repository contains scripts used in the assembly of MAGs to uncover the within-species diversity of methane-oxidizing (MOX) symbionts in Gigantidasplatifrons platifrons. The efficiency of the pipeline can be assessed using the scripts located in the mockTest_scripts folder. Additionally, these scripts facilitate analysis on positive selection, clade-specific gene recovery, and differential gene distribution across specific loci in the genomes. These scripts are provided to assist readers of the manuscript titled "Mosaic environment-driven evolution of the deep-sea mussel Gigantidas platifrons bacterial endosymbiont" in replicating the analysis. We have not verified their applicability to other datasets, and it's important to note that they are not intended as tools for broader use.

For most of the scripts, -h option following each command will give the general introduction on their usages.

The modified binning scripts can be found in the 'metawrap-modules' directory. While many of these scripts originate from the official metawrap website, please note that the binning module is an exception. This module has undergone significant modifications, and many of the binner packages are not included in this repository. If required, you'll need to download them separately. The script RPKG_fromBam_minimap2.py and merge_DEPTH.py, which parse the coverage information, are available in this folder. Ensure they're placed in a directory that's included in your system's PATH.

Scripts for phylogenetic tree construction using both single-copy and multiple-copy orthologs can be found in the phylogeny folder. This includes astral_tre.pl and adj_AstrialTree.py. To further construct the cds trees, follow these steps: Align and trim the corresponding cds using alignCDSonFaaALN.py, then concatenate them with concanateOrtho.py. Employ fasttree or raxML to construct trees as needed.

The mockTest_scripts folder contains scripts tailored for assessing the performance of the reads-binning pipeline. This includes pacbio_reads_stimulator.py for creating mock PacBio reads, illu_reads_stimulator.py for generating mock Illumina reads aimed at coverage recovery, and assemble_mock_flye.py to aid in batch assembly. Additionally, colinearity_analysis.py is provided to visualize the colinearity between the assembled genome and reference genomes from which the reads were derived.

The Population folder houses modified Snakemake pipelines intended for population genetic analysis and DESMAN strain recovery. Use get_strain_fasta.py to extract the fasta files for individual strains for subsequent analysis. Customize the pipeline to your datasets by adjusting the config.yaml file. Before initiating the pipeline, ensure that both DESMAN and GATK4 are properly installed.

Scripts designed for statistical analysis of PanGenomes are located in the Pan_statistics directory. Among these, get_specialOrtholog.py identifies clade-specific genes from orthologs; selection_paml_group.py executes PAML analysis to pinpoint genes under selection for individual clades; bilateral_gene_cruise.py is used to compare gene density around a specified locus. Additionally, pan_genome_catalog.py, COG_summarize.py, and compare_core_specific.py are employed to determine core genes and analyze differences in their COG distribution.

The Evolutionary_rates folder houses scripts that compare the evolutionary rates of clades. Within this directory, dnds_paml_calculation.py is employed to retrieve the dN/dS values for genes that meet the filtration criteria (filtered by dnds_filtering.py) as well as for all concatenated genes. Meanwhile, dnds_permutation.py calculates the dN/dS values for permuted gene blocks, each comprised of 20 randomly selected genes.

The Other folder contains additional scripts such as GOenrich.r and KEGGenrich.R, designed for GO and KEGG enrichment analyses respectively. Additionally, the FastANI_heatmap.py script is included for visualizing FastANI results in heatmap format.
