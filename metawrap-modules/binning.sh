#!/usr/bin/env bash

##############################################################################################################################################################
#
# This script is meant to be run on the outputs of assembly.sh pipeline to split the assembly contigs into metagenomic bins.
# Ideally it should take in the assembly file of all of your samples, followed by the reads of all the samples that went into the assembly.
# The more samples, the better the binning. 
#
# The script uses metaBAT2 and/or CONCOCT and/or MaxBin2 to bin the contigs. MetaBAT2 is the defualt due to its speed and great performance,
# but all these binners have their advantages and disadvantages, so it recomended to run the bin_refinement module to QC the bins, get the 
# best bins of all of each method, and to reassembly and refine the final bins. 
#
# Author of pipeline: German Uritskiy. I do not clain any authorship of the many programs this pipeline uses.
# For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu.
# The pipeline were modified by wangminxiao@qdio.ac.cn
#
##############################################################################################################################################################


help_message () {
	echo ""
	echo "Usage: metaWRAP binning [options] -a assembly.fa -o output_dir readsA_1.fastq | readsA.fq1.gz readsA_2.fastq | readsA.fq2.gz ... [readsX.fq1.gz readsX.fq2.gz]"
	echo "Note1: Make sure to provide all your separately replicate read files, not the joined file."
	echo "Note2: You may provide single end or interleaved reads as well with the use of the correct option"
	echo "Note3: If the output already has the .bam alignments files from previous runs, the module will skip re-aligning the reads"
	echo ""
	echo "Options:"
	echo ""
	echo "	-a STR          metagenomic assembly file"
	echo "	-o STR          output directory"
	echo "	-t INT          number of threads (default=1)"
	echo "	-m INT		amount of RAM available (default=4)"
	echo "	-M INT          how to deal with multiple hitting, 0 discarding, 1 pick as 1, 2 thresholding, 3 proportion, default 1"
	echo "	-l INT		minimum contig length to bin (default=1000bp). Note: metaBAT will default to 1500bp minimum"
	echo "	-k INT		kmer length for concoct composition calculation"
	echo "	-P INT		num of parellization (default 2)"
	echo "	-p INT STR	num of pacbio to align, list of files"
	echo "	-T STR		temporary dictionary to use (default=\\tmp)"
	echo "	-c INT		Groop Core minimum contig length (default=2000)"
	echo ""
        echo "	--mapOnly       map reads without bin operations"
	echo "	--metabat2      bin contigs with metaBAT2"
	echo "	--metabat1	bin contigs with the original metaBAT"
	echo "	--maxbin2	bin contigs with MaxBin2"
	echo "	--concoct	bin contigs with CONCOCT"
	echo "	--bmc3c		bin contigs with BMC3C"
	echo "	--vamb		bin contigs with VAMP"
	echo "	--binsanity     bin contigs with BinSanity"
	echo "	--solidbin      bin contigs with SolidBin"
	echo "	--groopm	bin contigs with GroopM"
	echo "	--semibin	bin contigs with SemiBin"
	echo "	--strain	map reads considering overlap strains"
	echo "	--universal	use universal marker genes instead of bacterial markers in MaxBin2 (improves Archaea binning)"
	echo "	--run-checkm	immediately run CheckM on the bin results (requires 40GB+ of memory)"
	echo "	--single-end	non-paired reads mode (provide *.fastq or *.*gz files)"
	echo "	--interleaved	the input read files contain interleaved paired-end reads"
	echo "";}

comm () { ${SOFT}/print_comment.py "$1" "-"; }
error () { ${SOFT}/print_comment.py "$1" "*"; exit 1; }
warning () { ${SOFT}/print_comment.py "$1" "*"; }
announcement () { ${SOFT}/print_comment.py "$1" "#"; }
run_checkm () {
	comm "Running CheckM on ${1} bins"

	# determine --pplacer_threads count. It is either the max thread count or RAM/4, whichever is higher
	ram_max=$(($mem / 40))
	if (( $ram_max < $threads )); then
		p_threads=$ram_max
	else
		p_threads=$threads
	fi
	comm "There is $mem RAM and $threads threads available, and each pplacer thread uses <40GB, so I will use $p_threads threads for pplacer"

	mkdir ${1}.tmp
	checkm lineage_wf -x fa ${1} ${1}.checkm -t $threads --tmpdir ${1}.tmp --pplacer_threads $p_threads
	if [[ ! -s ${1}.checkm/storage/bin_stats_ext.tsv ]]; then error "Something went wrong with running CheckM. Exiting..."; fi
	rm -r ${1}.tmp
	${SOFT}/summarize_checkm.py ${1}.checkm/storage/bin_stats_ext.tsv ${1##*/}\
	| (read -r; printf "%s\n" "$REPLY"; sort -rn -k2) > ${1}.stats
	if [[ $? -ne 0 ]]; then error "Cannot make checkm summary file. Exiting."; fi
	num=$(cat ${1}.stats | awk '{if ($2>=70 && $2<=100 && $3>=0 && $3<=10) print $1 }' | wc -l)
	comm "There are $num 'good' bins found in ${1}! (>70% completion and <10% contamination)"
}


########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################


# setting scripts and databases from config file (should be in same folder as main script)
config_file=$(which config-metawrap)
source $config_file

# default params
threads=1; mem=4; kmer=4; len=1000; out=false; ASSEMBLY=false; tmpD="/BioData/tmp"; core_length=2000;
# long options defaults
mapOnly=false;metabat1=false; metabat2=false; binsanity=false; maxbin2=false; concoct=false; binsanity=false; bmc3c=false; vamb=false; solidbin=false;
checkm=false; read_type=paired
markers=107
paralization=2
multiple=1
# load in params
OPTS=`getopt -o ht:m:M:o:a:l:k:T:p:P:c --long help,semibin,metabat1,groopm,strain,mapOnly,metabat2,binsanity,maxbin2,concoct,vamb,bmc3c,solidbin,run-checkm,single-end,universal,interleaved -- "$@"`
if [ $? -ne 0 ]; then help_message; exit 1; fi
strain=false
# loop through input params
while true; do
        case "$1" in
                -t) threads=$2; shift 2;;
		-c) core_length=$2; shift 2;;
		-T) tmpD=$2; shift 2;;
		-m) mem=$2; shift 2;;
		-M) multiple=$2; shift 2;;
                -o) out=$2; shift 2;;
                -a) ASSEMBLY=$2; shift 2;;
		-l) len=$2; shift 2;;
		-p) num_file=$2; for((i=3;i<=$[num_file+2];i++)); do  pac_file="$pac_file ${!i}" ; done; shift $[num_file+2];; 
		-P) paralization=$2; shift 2;;
		-k) kmer=$2; shift 2;;
                -h | --help) help_message; exit 1; shift 1;;
		--groopm) groopm=true; shift 1;;
		--metabat2) metabat2=true; shift 1;;
		--metabat1) metabat1=true; shift 1;;
		--mapOnly) mapOnly=true; shift 1;;
		--maxbin2) maxbin2=true; shift 1;;
		--concoct) concoct=true; shift 1;;
		--bmc3c) bmc3c=true; shift 1;;
		--vamb)  vamb=true; shift 1;;
		--binsanity) binsanity=true; shift 1;;
		--solidbin) solidbin=true; shift 1;;
		--semibin) semibin=true; shift 1;;
		--run-checkm) checkm=true; shift 1;;
		--single-end) read_type=single; shift 1;;
		--interleaved) read_type=interleaved; shift 1;;
		--strain) strain=true; shift 1;;
		--universal) markers=40; shift 1;;
                --) help_message; exit 1; shift; break ;;
                *) break;;
        esac
done
#echo "Hello:\t$[num_file+2]\t$pac_file"
########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################
# Make sure at least one binning method was chosen
if [ $maxbin2 = false ] && [ $vamb = false ] && [ $mapOnly = flase ] && [ $metabat2 = false ] && [ $bmc3c = false ] && [ $concoct = false ] && [ $binsanity = false ]; then
	help_message
	error "You must select at least one binning method or mapOnly"
fi

if [ $strain = true ]; then
	map_max=250
else
	map_max=5
fi

# check if all parameters are entered
if [ $out = false ] || [ $ASSEMBLY = false ] ; then 
	comm "Non-optional parameters -a and/or -o were not entered"
	help_message; exit 1
fi

#check if the assembly file exists
if [ ! -s $ASSEMBLY ]; then error "$ASSEMBLY does not exist. Exiting..."; fi

comm "Entered read type: $read_type"

if [ $read_type = paired ]; then
	# check for at least one pair of read *gz files:
	F="no"; R="no"
	for num in "$@"; do
		if [[ $num == *".fq1.gz" || $num == *"_1.fastq" ]]; then F="yes"; fi
		if [[ $num == *".fq2.gz" || $num == *"_2.fastq" ]]; then R="yes"; fi
	done
	if [ $F = "no" ] || [ $R = "no" ]; then
		comm "Unable to find proper *gz or fastq read pairs in the format *fq1.gz and *fq2.gz or in the format of *_1.fastq and *_2.fastq"
		help_message; exit 1
	fi
else
	# check for at least one *gz read
	F="no"
	for num in "$@"; do
		if [[ $num == *".fq.gz" || $num == *".fastq" ]]; then F="yes"; fi
	done
	if [ $F = "no" ]; then
		comm "Unable to find read files in format *.*gz (for single-end or interleaved reads)"
		help_message; exit 1
	fi
fi

if [ $read_type = paired ]; then
	#determine number of *gz read files provided:
	num_of_F_read_files=$(for I in "$@"; do echo $I | grep -E 'fq1.gz|_1.fastq'; done | wc -l)
	num_of_R_read_files=$(for I in "$@"; do echo $I | grep -E 'fq2.gz|_2.fastq'; done | wc -l)

	comm "$num_of_F_read_files forward and $num_of_R_read_files reverse read files detected"
	if [ ! $num_of_F_read_files == $num_of_R_read_files ]; then error "Number of F and R reads must be the same!"; fi
fi

if [ $len -lt 1500 ]; then
	metabat_len=1500
else
	metabat_len=$len
fi


# Checks for correctly configures meta-scripts folder
if [ ! -s $SOFT/sort_contigs.py ]; then
	error "The folder $SOFT doesnt exist. Please make sure config.sh is in the same filder as the mains scripts and all the paths in the config.sh file are correct"
fi


########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################


########################################################################################################
########################         ALIGNING READS TO MAKE COVERAGE FILES          ########################
########################################################################################################
announcement "ALIGNING READS TO MAKE COVERAGE FILES"

# setting up the output folder
if [ ! -d $out ]; then mkdir $out;
else 
	echo "Warning: $out already exists."
	rm -r ${out}/*checkm
fi

if [ ! -d ${out}/work_files ]; then mkdir ${out}/work_files; fi

if [ -f ${out}/work_files/assembly.fa ]; then
	comm "Looks like the assembly file is already coppied. Skipping..."
else
	comm "making copy of assembly file $ASSEMBLY"
	cp $ASSEMBLY ${out}/work_files/assembly.fa
fi

tmp=${ASSEMBLY##*/}
sample=${tmp%.*}

# Index the assembly
if [ -f ${out}/work_files/ref.mmi ]; then
	comm "Looks like there is a index of the assembly already. Skipping..."
else
	comm "Indexing assembly file"
	minimap2 -I 40g -t $threads -d ${out}/work_files/ref.mmi ${out}/work_files/assembly.fa
	if [[ $? -ne 0 ]] ; then error "Something went wrong with indexing the assembly. Exiting."; fi
fi

# If there are several pairs of reads passed, they are processed sepperately

if [ num_file>0 ]; then
	if [ -f ${out}/work_files/ref.mmi ]; then
        	comm "Looks like there is a index of the assembly already. Skipping..."
	else
        	comm "Indexing assembly file for pacbio reads"
		minimap2 -I 40g -t $[ threads/2 ] -d ${out}/work_files/ref.mmi ${out}/work_files/assembly.fa	
		if [[ $? -ne 0 ]] ; then error "Something went wrong with indexing the assembly. Exiting."; fi
	fi
	pac_file_adj=$(echo $pac_file | sed "s/^\s\+//" | sed "s/\s\+$//")
	Pac_arr=($pac_file_adj)
	for id in $(seq 0 $[num_file-1])
	do
		#echo "hello "${Pac_arr[$id]}
		tmp=`basename ${Pac_arr[$id]}`;
		BASE=${tmp%.*}
		if [[ ! -f ${out}/work_files/${BASE}_pacmapped.sorted.bam ]]; then
		comm "Aligning ${Pac_arr[$id]} to assembly"
		minimap2 -t $threads -N ${map_max} --secondary yes -ax map-hifi $out/work_files/ref.mmi ${Pac_arr[$id]} -o ${out}/work_files/${BASE}_pacmapped.sam
		RPKG_fromBam_minimap2.py -s -r 1 -f ${out}/work_files/${BASE}_pacmapped.sam -m $multiple -o ${out}/work_files/${BASE}_RPKG -i 0.98 -c 0.75;	
		sleep 1
		#--sort -m 40G $out/work_files/ref.mmi -j $threads ${Pac_arr[$id]} $out/work_files/${BASE}_pacmapped.bam
		else comm "skipping aligning $BASE reads to assembly because $out/work_files/${BASE}_pacmapped.bam already exists"
		fi
	done	
fi
rm -f ${out}/work_files/mapping.sh
alias mv="mv -f"
for num in "$@"; do
	# paired end reads
	if [ $read_type = paired ]; then
		if [[ $num == *"fq1.gz"* || $num == *"_1.fastq"* ]]; then 
			reads_1=$num
			if [[ $num =~ 'gz' ]];then reads_2=${num%.fq*.gz}.fq2.gz; fi
			if [[ $num =~ 'fastq' ]];then reads_2=${num%_1.fastq}_2.fastq; fi
			if [ ! -s $reads_1 ]; then error "$reads_1 does not exist. Exiting..."; fi
			if [ ! -s $reads_2 ]; then error "$reads_2 does not exist. Exiting..."; fi
			if [[ $reads_1 =~ 'fastq' ]];
			then
			tmp=${reads_1##*/}
			sample=${tmp%_*}
			fi
			if [[ $reads_1 =~ 'gz' ]];
			then
			tmp=${reads_1##*/}
			sample=${tmp%.fq*.gz}
			fi
			if [[ ! -f ${out}/work_files/${sample}_DEPTH.tsv ]]; then
				comm "Aligning $reads_1 and $reads_2 back to assembly"
				out2=`basename $out`
				if [[ $strain = true ]]; then
					echo "minimap2 -ax sr --secondary yes -t $threads -N ${map_max} ${out}/work_files/ref.mmi $reads_1 $reads_2 | awk '{if (\$1 ~ /^@/ || \$12~/NM:i:/) {split(\$12,a,\":\"); if(a[3]<=3 || \$1 ~ /^@/) print \$0}}' > ${tmpD}/${out2}_${sample}_tmp.sam; RPKG_fromBam_minimap2.py -s -r 1 -S -f ${tmpD}/${out2}_${sample}_tmp.sam -m $multiple -M 80000 -o ${out}/work_files/${sample}_RPKG -i 0.98 -c 0.9; wait" >>${out}/work_files/mapping.sh
				else
					echo "minimap2 -ax sr --secondary yes -t $threads -N ${map_max} ${out}/work_files/ref.mmi $reads_1 $reads_2 -o ${tmpD}/${out2}_${sample}_tmp.sam; RPKG_fromBam_minimap2.py -s -r 1 -f ${tmpD}/${out2}_${sample}_tmp.sam -m $multiple -M 500000 -o ${out}/work_files/${sample}_RPKG -i 0.98 -c 0.9" >>${out}/work_files/mapping.sh
					#echo "minimap2 -ax sr -t $threads ${out}/work_files/ref.mmi $reads_1 $reads_2 -o ${tmpD}/${out2}_${sample}_tmp.sam; RPKG_fromBam_minimap2.py -s -r 1 -f ${tmpD}/${out2}_${sample}_tmp.sam -m $multiple -M 150000 -o ${out}/work_files/${sample}_RPKG -i 0.98 -c 0.9" >>${out}/work_files/mapping.sh
				fi
			fi
		fi
	# single end or interleaved reads
	else
		if [[ $num == *".fq.gz"* || $num == *".fastq"* ]]; then
			reads=$num
			if [ ! -s $reads ]; then error "$reads does not exist. Exiting..."; fi
			tmp=${reads##*/}
			sample=${tmp%.*}
			if [[ ! -f ${out}/work_files/${sample}.bam ]]; then
				comm "Aligning $reads back to assembly, and sorting the alignment"
				if [ $read_type = single ]; then
					bwa mem -t $threads ${out}/work_files/assembly.fa $reads > ${out}/work_files/${sample}.sam
					if [[ $? -ne 0 ]]; then error "Something went wrong with aligning the reads to the assembly!"; fi
				elif [ $read_type = interleaved ]; then
					bwa mem -p -t $threads ${out}/work_files/assembly.fa $reads > ${out}/work_files/${sample}.sam
					if [[ $? -ne 0 ]]; then error "Something went wrong with aligning the reads to the assembly!"; fi
				else
					error "something from with the read_type (=$read_type)"
				fi
				
				comm "Sorting the $sample alignment file"
				samtools sort -T ${out}/work_files/tmp-samtools -@ $threads -O BAM -o ${out}/work_files/${sample}.bam ${out}/work_files/${sample}.sam && rm -rf ${out}/work_files/${sample}.sam 
				if [[ $? -ne 0 ]]; then error "Something went wrong with sorting the alignments. Exiging..."; fi
			else
				comm "skipping aligning $sample reads to assembly because ${out}/work_files/${sample}.bam already exists."
			fi
		fi
	fi
done

if [ -s "${out}/work_files/mapping.sh" ]; then
	ParaFly -c ${out}/work_files/mapping.sh -CPU ${paralization}
	if [[ $? -ne 0 ]]; then error "Something went wrong with sorting the alignments. Exiging..."; fi 
fi

if [ ! -f "${out}/work_files/mergedRPKG-pseudoVariation.tsv" ]; then 	
	merge_RPKG.py -f ${out}/work_files -o ${out}/work_files/mergedRPKG -v -s ${out}/work_files/assembly.fa
fi

if [ ! -f "${out}/work_files/mergedDEPTH-pseudoVariation.tsv" ]; then
        merge_DEPTH.py  -f ${out}/work_files -o ${out}/work_files/mergedDEPTH -v -s ${out}/work_files/assembly.fa
fi


alias mv='mv -i'

if [ $metabat2 = true ]  && [ ! -f "${out}/metabat2.done" ]; then
	########################################################################################################
	########################                   RUNNING METABAT2                     ########################
	########################################################################################################
	announcement "RUNNING METABAT2"

	comm "Starting binning with metaBAT2..."
	metabat2 -i ${out}/work_files/assembly.fa -a ${out}/work_files/mergedDEPTH-pseudoVariation.tsv \
	 -o ${out}/metabat2_bins/bin -m $metabat_len -t $threads --unbinned
	if [[ $? -ne 0 ]]; then error "Something went wrong with running MetaBAT2. Exiting"; fi
	comm "metaBAT2 finished successfully, and found $(ls -l ${out}/metabat2_bins | grep .fa | wc -l) bins!"

	if [ $checkm = true ]; then
		run_checkm ${out}/metabat2_bins
	fi
	touch ${out}/metabat2.done
fi


if [ $metabat1 = true ]  && [ ! -f "${out}/metabat1.done" ]; then
        ########################################################################################################
        ########################                   RUNNING METABAT1                     ########################
        ########################################################################################################
        announcement "RUNNING METABAT1"

        comm "making contig depth file..."
        jgi_summarize_bam_contig_depths --outputDepth ${out}/work_files/metabat_depth.txt ${out}/work_files/*.bam
        if [[ $? -ne 0 ]]; then error "Something went wrong with making contig depth file. Exiting."; fi

        comm "Starting binning with metaBAT1..."
        metabat1 -i ${out}/work_files/assembly.fa -a ${out}/work_files/metabat_depth.txt\
         -o ${out}/metabat1_bins/bin -m $metabat_len --minContigByCorr $len -t $threads --unbinned --superspecific
        if [[ $? -ne 0 ]]; then error "Something went wrong with running MetaBAT1. Exiting"; fi
        comm "metaBAT1 finished successfully, and found $(ls -l ${out}/metabat1_bins | grep .fa | wc -l) bins!"

        if [ $checkm = true ]; then
                run_checkm ${out}/metabat1_bins
        fi
	touch ${out}/metabat1.done
fi



if [ $maxbin2 = true ]  && [ ! -f "${out}/maxbin2.done" ]; then
        ########################################################################################################
        ########################                   RUNNING MAXBIN2                     ########################
        ########################################################################################################
        announcement "RUNNING MAXBIN2"

	comm "making contig depth file..."

	#calculate total numper of columns
	A=($(head -n 1 ${out}/work_files/mergedDEPTH.tsv)) 
	N=${#A[*]}
	
	# split the contig depth file into multiple files
	comm "split master contig depth file into individual files for maxbin2 input"
	if [ -f ${out}/work_files/mb2_abund_list.txt ]; then rm ${out}/work_files/mb2_abund_list.txt; fi
	for i in $(seq 2 $N); do 
		sample=$(head -n 1 ${out}/work_files/mergedDEPTH.tsv | cut -f $i)
		echo "processing $sample depth file..."
		grep -v totalAvgDepth ${out}/work_files/mergedDEPTH.tsv | cut -f 1,$i > ${out}/work_files/mb2_${sample%.*}.txt
		if [[ $out == /* ]]; then
			echo ${out}/work_files/mb2_${sample%.*}.txt >> ${out}/work_files/mb2_abund_list.txt
		else
			echo $(pwd)/${out}/work_files/mb2_${sample%.*}.txt >> ${out}/work_files/mb2_abund_list.txt
		fi
	done

	run_MaxBin.pl
	if [[ $? -ne 0 ]]; then
		comm "looks like our default perl libraries are not the conda ones. Manually setting perl5 library directory"
        	conda_path=$(which metawrap)
		echo "metawrap path: $conda_path"
		conda_path=${conda_path%/*}
		if [ $(echo -n $conda_path | tail -c 1) = "/" ]; then conda_path=${conda_path%/*}; fi
		conda_path=${conda_path%/*}
		if [ ! -d ${conda_path}/lib/perl5/site_perl/5.22.0 ]; then 
			error "${conda_path}/lib/perl5/site_perl/5.22.0 does not exixt. Cannot set manual path to perl5 libraries. Exiting..."
		fi
	
	        perl_libs=${conda_path}/lib/perl5/site_perl/5.22.0
	        echo "Will use perl5 libraries located in $perl_libs - hopefully they are there..."
		export PERL5LIB="$perl_libs"
	fi

	
	comm "Starting binning with MaxBin2..."
	mkdir ${out}/work_files/maxbin2_out
	run_MaxBin.pl -contig ${out}/work_files/assembly.fa -markerset $markers -thread $threads -min_contig_length $len\
	-out ${out}/work_files/maxbin2_out/bin \
	-abund_list ${out}/work_files/mb2_abund_list.txt

	if [[ $? -ne 0 ]]; then error "Something went wrong with running MaxBin2. Exiting."; fi
	if [[ $(ls ${out}/work_files/maxbin2_out/ | grep bin | grep .fasta | wc -l) -lt 1 ]]; then error "MaxBin2 did not pruduce a single bin. Something went wrong. Exiting."; fi

	mkdir ${out}/maxbin2_bins
	N=0
	for i in $(ls ${out}/work_files/maxbin2_out/ | grep bin | grep .fasta); do
		cp ${out}/work_files/maxbin2_out/$i ${out}/maxbin2_bins/bin.${N}.fa
		N=$((N + 1))
	done
	comm "MaxBin2 finished successfully, and found $(ls -l ${out}/maxbin2_bins | grep .fa | wc -l) bins!"

	if [ $checkm = true ]; then
		run_checkm ${out}/maxbin2_bins
	fi
	touch ${out}/maxbin2.done
fi

	

if [ $concoct = true ]  && [ ! -f "${out}/concoct.done" ]; then
	########################################################################################################
	########################                    RUNNING CONCOCT                     ########################
	########################################################################################################
        announcement "RUNNING CONCOCT"

        comm "Starting binning with CONCOCT..."
        # I have to do some directory changing because CONCOCT dumps all files into current directory...
        home=$(pwd)
        mkdir ${out}/work_files/concoct_out
	/nfs_genome/anaconda/envs/concoct-env/bin/concoct -i 2000 -l $len -t $threads -k $kmer\
		--coverage_file ${out}/work_files/mergedDEPTH.tsv \
		--composition_file ${out}/work_files/assembly.fa \
		-b ${out}/work_files/concoct_out

	if [[ $? -ne 0 ]]; then error "Something went wrong with binning with CONCOCT. Exiting..."; fi


        comm "splitting contigs into bins"
        mkdir ${out}/concoct_bins
	cd ${out}/concoct_bins
        bin_fasta.pl ../work_files/concoct_out/clustering_gt${len}.csv ../work_files/assembly.fa 1 300000
        cd $home
	
	comm "CONCOCT finished successfully, and found $(ls -l ${out}/concoct_bins | grep .fa | wc -l) bins!"
	sed -i "s/\(>\S\+\)\s.*/\1/" ${out}/concoct_bins/*.fa
	if [ $checkm = true ]; then
		run_checkm ${out}/concoct_bins
	fi
	touch ${out}/concoct.done
fi


if [ $bmc3c = true ]  && [ ! -f "${out}/bmc3c.done" ]; then
        ########################################################################################################
        ########################                    RUNNING BMC3C                       ########################
        ########################################################################################################
        announcement "RUNNING BMC3C"
	mydir=`pwd`
	mkdir ${out}/work_files/bmc3c_out
	cat ${out}/work_files/metabat_depth.txt | awk -F "\t" '{a=$1;for(i=4;i<=NF;i+=2)a=a FS $i;print a}' >${out}/work_files/bmc3c_out/cov_input.tsv
	removeemptycontigs.pl ${out}/work_files/bmc3c_out
	prodigal -p meta -i ${out}/work_files/assembly.fa -o  ${out}/work_files/bmc3c_out/contigs.describe -a ${out}/work_files/bmc3c_out/contigs.proteins -d ${out}/work_files/bmc3c_out/contigs.genes
	removeemptygenes.pl ${out}/work_files/bmc3c_out
	codon_usage.py ${out}/work_files/bmc3c_out/contigs_noempty.fasta ${out}/work_files/bmc3c_out/contigs_noempty.genes ${out}/work_files/bmc3c_out/
	count_contig=`grep ">" -c ${out}/work_files/bmc3c_out/contigs_gene_filter.fasta`
	
	fasta_to_features.py ${out}/work_files/bmc3c_out/contigs_gene_filter.fasta $count_contig $kmer ${out}/work_files/bmc3c_out/kmer_tmp.csv
	# generating codon_usage.txt,kmer_tmp.csv,cov_input.tsv
	awk '{for(i=2;i<=NF;i++) if($i >0 && $i < 0.001) $i =0; print $0}' ${out}/work_files/bmc3c_out/cov_input.tsv >${out}/work_files/bmc3c_out/cov_input_rev.tsv;
	for id in $(grep ">" ${out}/work_files/bmc3c_out/contigs_gene_filter.fasta | sed "s/^>//"); do grep -P "^$id\t" ${out}/work_files/bmc3c_out/cov_input_rev.tsv; done >${out}/work_files/bmc3c_out/cov_final.tsv	
	title=`head ${out}/work_files/bmc3c_out/cov_input.tsv -n1`;
	sed -i "1i$title" ${out}/work_files/bmc3c_out/cov_final.tsv
	cp /nfs_genome/BioInformatics/BMC3C/bmc3c.m ${out}/work_files/bmc3c_out
	sed -i 'N;2iPATH="'"$mydir"'/'"$out"'/work_files/bmc3c_out"' ${out}/work_files/bmc3c_out/bmc3c.m
	sed -i "1iname,1,2,3,4,5" ${out}/work_files/bmc3c_out/codon_usage.txt
	cd ${out}/work_files/bmc3c_out/
	matlab -nodesktop -nosplash -r bmc3c
	if [[ $? -ne 0 ]]; then error "Something went wrong with bmc3c. Exiting."; fi
	mkdir ${mydir}/${out}/bmc3c_bins
	cd ${mydir}/${out}/bmc3c_bins
	bin_fasta.pl ../work_files/bmc3c_out/results.csv ../work_files/assembly.fa 1 50000	
	cd $mydir
	comm "BMC3C finished successfully, and found $(ls -l ${odut}/bmc3c_bins | grep .fa | wc -l) bins!"
	if [ $checkm = true ]; then
                run_checkm ${out}/bmc3c_bins
        fi
	touch ${out}/bmc3c.done
fi


if [ $vamb = true ]  && [ ! -f "${out}/vamb.done" ]; then
        ########################################################################################################
        ########################                    RUNNING VAMB                        ########################
        ########################################################################################################
        announcement "RUNNING VAMB"
        mydir=`pwd`
        mkdir ${out}/work_files/vamb_out
	cd ${out}/work_files/vamb_out
	#vamb -o SEP -z 0.95 -s 30 --outdir OUT --fasta FASTA --bamfiles
	lengthFilter.py -l 1500 -s ../assembly.fa -p ../mergedDEPTH-pseudoVariation.tsv
	/nfs_genome/anaconda3/envs/vamb/bin/vamb -m 1000 --outdir ${out}/work_files/vamb_out/out --fasta ${out}/work_files/vamb_out/assembly_1500.fa --jgi ${out}/work_files/vamb_out/L1500_mergedDEPTH-pseudoVariation.tsv
	if [[ $? -ne 0 ]]; then error "Something went wrong with vamb binning. Exiting."; fi
	mkdir ${out}/vamb_bins
        cd ${out}/vamb_bins
        bin_fasta.pl ../work_files/vamb_out/out/clusters.tsv ../work_files/assembly.fa 2 200000
        cd $mydir
        comm "VAMB finished successfully, and found $(ls -l ${out}/vamb_bins | grep .fa | wc -l) bins!"
        if [ $checkm = true ]; then
                run_checkm ${out}/vamb_bins
        fi
	touch ${out}/vamb.done
fi


if [ $binsanity = true ] && [ ! -f "${out}/binsanity.done" ]; then
        ########################################################################################################
        ########################                   RUNNING BinSanity                    ########################
        ########################################################################################################
        announcement "RUNNING BinSanity"
        mydir=`pwd`
        mkdir ${out}/work_files/BinSanity_out
	contig_num=`grep ">" -c ${out}/work_files/assembly.fa`
	if [ $contig_num -le 110000 ]; then
	let thread_checkm=threads/8
	binsoft="/nfs_genome/anaconda3/envs/BinSanity/bin/Binsanity-wf --threads $thread_checkm"
	else
	binsoft="/nfs_genome/anaconda3/envs/BinSanity/bin/Binsanity-lc -C 70 --checkm_threads $thread_checkm --kmean_threads $threads"
	fi
	if [[ $? -ne 0 ]]; then error "Something went wrong with BinSanity. Exiting."; fi
	rm_enter_fas ${out}/work_files/assembly.fa
	mv sim_assembly.fa ${out}/work_files/BinSanity_out
	mkdir ${out}/work_files/BinSanity_out/bam
	ln -sf ${mydir}/${out}/work_files/*.bam ${mydir}/${out}/work_files/BinSanity_out/bam
	grep ">" ${out}/work_files/BinSanity_out/sim_assembly.fa | sed "s/>//" >${out}/work_files/BinSanity_out/contig_id
	cd ${out}/work_files/BinSanity_out/
	/nfs_genome/anaconda3/envs/BinSanity/bin/Binsanity-profile -i $mydir/${out}/work_files/BinSanity_out/sim_assembly.fa -s $mydir/${out}/work_files/BinSanity_out/bam --transform scale -T 50 -c Binsanity.cov
        
	#cd ${out}/work_files/BinSanity_out/	
	$binsoft -c $mydir/${out}/work_files/mergedDEPTH.tsv -f $mydir/${out}/work_files/BinSanity_out/ -p -10 -l sim_assembly.fa -o $dir/${out}/work_files/BinSanity_out/
	cd $mydir
	mv ${out}/work_files/BinSanity_out/BinSanity-Final-bins ${out}
	for id in $(ls ${out}/BinSanity-Final-bins/*.fna | sed "s/\.fna//"); do mv $id.fna $id.fa; done
        comm "Binsanity finished successfully, and found $(ls -l ${out}/BinSanity-Final-bins | grep .fa | wc -l) bins!"
        if [ $checkm = true ]; then
                run_checkm ${out}/BinSanity-Final-bins
        fi
	touch ${out}/binsanity.done
fi

if [ $solidbin = true ] && [ ! -f "${out}/SolidBin.done" ]; then
        ########################################################################################################
        ########################                   RUNNING SolidBin                     ########################
        ########################################################################################################
        announcement "RUNNING SolidBin"
        mydir=`pwd`
        mkdir ${out}/work_files/solidbin_out
	cd ${out}/work_files/
	samtools faidx assembly.fa
	awk -v OFS='\t' {'print $1,$2'} assembly.fa.fai >length.txt
	cd solidbin_out
	bash /nfs_genome/BioInformatics/SolidBin/scripts/run.sh ../assembly.fa 1000 4
	/nfs_genome/anaconda3/envs/solidbin/bin/python /nfs_genome/BioInformatics/SolidBin/SolidBin.py --contig_file ../assembly.fa --thread $threads --coverage_profiles ../mergedDEPTH.tsv --composition_profiles ../kmer_4_f1000.csv --output ./solidbin.tsv --log solidbin.log --use_sfs
	if [[ $? -ne 0 ]]; then error "Something went wrong with SolidBin binning. Exiting."; fi
	cp ./good_bins/ ../../SolidBin_bins -rf
	cd ../../SolidBin_bins
	for id in $(ls *.bin); do mv $id $id.fa; done
	cd $mydir
        comm "SolidBin finished successfully, and found $(ls -l ${out}/SolidBin_bins | grep .fa | wc -l) bins!"
        if [ $checkm = true ]; then
                run_checkm ${out}/SolidBin_bins
        fi
	touch ${out}/SolidBin.done
fi


if [ "$groopm" = true ] && [ ! -f "${out}/GroopM.done" ]; then
        ########################################################################################################
        ########################                   RUNNING GroopM                       ########################
        ########################################################################################################
        announcement "RUNNING GroopM"
        mydir=`pwd`
        mkdir ${out}/work_files/groopm_out
	if [[ `ls ${out}/work_files/primary/*.bam 2> /dev/null | wc -l` -eq 0 ]]; then
		mkdir ${out}/work_files/primary
		for id in $(ls ${out}/work_files/*.bam); do
			base=`basename $id` 
			samtools view -F 256 -@ 20 $id -b -o ${out}/work_files/primary/$base.primary.bam
			if [ ! -f ${out}/work_files/primary/$base.primary.bam.bai ]; then samtools index ${out}/work_files/primary/$base.primary.bam; fi
			if [[ $? -ne 0 ]]; then error "Something went wrong with sorting. Exiting."; fi
		done
	fi
        cd ${out}/work_files/
	mkdir ../GroopM_bins

        cd groopm_out 
        ln ../assembly.fa ./
	bam_groopm=$(ls ../primary/*.bam )
	if [ ! -s "db.gm" ]; then
	/nfs_genome/anaconda/envs/GroopMv0.3/bin/groopm parse db.gm assembly.fa $bam_groopm -t $threads -f; fi
	if [[ $? -ne 0 ]]; then error "Something went wrong with groopm database building. Exiting."; fi	
	/nfs_genome/anaconda/envs/GroopMv0.3/bin/groopm core db.gm -f -c $core_length -b 600000
	if [[ $? -ne 0 ]]; then error "Something went wrong with groopm clustering. Exiting."; fi
	/nfs_genome/anaconda3/envs/GroopMv0.3/bin/groopm refine db.gm -a
	if [[ $? -ne 0 ]]; then error "Something went wrong with groopm refining. Exiting."; fi
	/nfs_genome/anaconda/envs/GroopMv0.3/bin/groopm recruit db.gm -c 1000
	if [[ $? -ne 0 ]]; then error "Something went wrong with groopm recruiting. Exiting."; fi
	/nfs_genome/anaconda/envs/GroopMv0.3/bin/groopm extract db.gm assembly.fa -o ../../GroopM_bins
	if [[ $? -ne 0 ]]; then error "Something went wrong with groopm bin extraction. Exiting."; fi
	cd ../../GroopM_bins
        for id in $(ls *.fna | sed "s/.fna//"); do mv $id.fna $id.fa; done
        cd $mydir
        comm "GroopM finished successfully, and found $(ls -l ${out}/GroopM_bins | grep .fa | wc -l) bins!"
        if [ $checkm = true ]; then
                run_checkm ${out}/GroopM_bins
        fi
        touch ${out}/GroopM.done
fi

if [ "$semibin" = true ] && [ ! -f "${out}/semibin.done" ]; then
        ########################################################################################################
        ########################                   RUNNING SemiBin                      ########################
        ########################################################################################################
        announcement "RUNNING SemiBin"
	mydir=`pwd`
	export OPENBLAS_NUM_THREADS=32
	mkdir ${out}/work_files/semibin
	if [[ `ls ${out}/work_files/primary/*.bam 2> /dev/null | wc -l` -eq 0 ]]; then
                mkdir ${out}/work_files/primary
                for id in $(ls ${out}/work_files/*.bam); do
                        base=`basename $id`
                        samtools view -F 256 $id -@ 20 -b -o ${out}/work_files/primary/$base.primary.bam
                        if [ ! -f ${out}/work_files/primary/$base.primary.bam.bai ]; then samtools index ${out}/work_files/primary/$base.primary.bam; fi
			if [[ $? -ne 0 ]]; then error "Something went wrong with sorting. Exiting."; fi
                done
        fi
	
	bam=$(ls ${out}/work_files/primary/*.bam | grep -v pacmapped)
	/nfs_genome/anaconda/envs/semibin/bin/SemiBin2 single_easy_bin -i ${out}/work_files/assembly.fa -b $bam -o ${out}/work_files/semibin_out
	if [[ $? -ne 0 ]]; then error "Something went wrong with semibin. Exiting."; fi
	mv ${mydir}/${out}/work_files/semibin_out/output_bins/ ${out}/semibin_bins
	cd ${out}/semibin_bins
	for id in $(ls *.gz); do gzip -d $id; done
	cd $mydir
	if [ $checkm = true ]; then
                run_checkm ${out}/semibin_bins
        fi
	touch ${out}/semibin.done
fi

#comm "cleaning up *.sam to save space..." #change to rm after mapping
#rm ${out}/work_files/*sam

########################################################################################################
########################      BINNING PIPELINE SUCCESSFULLY FINISHED!!!         ########################
########################################################################################################
announcement "BINNING PIPELINE SUCCESSFULLY FINISHED!!!"

