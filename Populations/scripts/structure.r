#!/usr/bin/env Rscript

#Loading package



library('seqinr')

# This function reads genes in gff format
get_genes_gff<-function(gff){
    genes=read.table(gff,sep="\t",quote="")
    sub=genes$V3=="CDS"
    genes=genes[sub,]
    genes2=data.frame(Contig=genes$V1,Start=genes$V4,End=genes$V5)
    s1=genes$V9
    s2=unlist(regmatches(s1,gregexpr("ID=.*?;",s1))) #extract gene names
    genes2$Name=substr(s2,4,nchar(s2)-1)
    genes2
}
# This function reads genes in fasta format
get_genes_fasta<-function(fasta){
  genes=read.fasta(fasta,seqtype = c("DNA"), as.string = F, set.att = FALSE)
  gene_lengths <- data.frame(Contig=character(0), Length=numeric(0))
  for (gene in names(genes)){
    gene_lengths <- rbind(gene_lengths, data.frame(Contig=gene, Length=length(genes[[gene]])))
  }
  gene_lengths
}



get_samples<-function(samplefile){

	data=read.table(samplefile)
	samples=as.character(data$V2)
	names(samples)=data$V1
	samples
}

get_snps<-function(samples){

	print("Reading data")
	alldata=list()
	for (s in names(samples)){
		print(s)
	  alldata[[s]]=read.table(samples[s],header=T)
	}

	data=alldata[[1]]
	for (i in seq(2,length(samples))){
		data=cbind(data,alldata[[i]][,c("A","C","G","T")])
	}

	names=c("Contig","Pos","Ref")
	for (s in names(samples)){
		names=c(names,paste(s,c("A","C","G","T"),sep="."))
	}
	names(data)=names
	data
}

add_pi<-function(snps,samples){
  
	ns=dim(snps)[1]
	print ("Pi calculation")
	for (s in samples){
		print(s)
		snps[,paste("pi__",s,sep="")]=rep(0,ns)
		names=paste(s,c("A","C","G","T"),sep=".")
		for (i in seq(1,ns)){
			subd=snps[i,names]
			sd=rowSums(subd)
			dem=1/sd/(sd-1)
			dem[!is.finite(dem)]=0 #do not count the one with 1 or 0 entries
			pi=2*(subd[,1]*subd[,2]*dem+subd[,1]*subd[,3]*dem+subd[,1]*subd[,4]*dem+subd[,2]*subd[,3]*dem+subd[,2]*subd[,4]*dem+subd[,3]*subd[,4]*dem)
			snps[i,paste("pi__",s,sep="")]=pi
		}
	}
	snps
}

add_pi2<-function(snps,samples){

	print("Between-sample pi calculation")
	ns=dim(snps)[1]
	print(ns)
	for (i1 in seq(1,length(samples)-1)){
		s1=samples[i1]
		for (i2 in seq(i1+1,length(samples))){
			s2=samples[i2]
			print(c(s1,s2))
			name=paste("pi",s1,s2,sep="__")
			snps[,name]=rep(0,ns)
			names1=paste(s1,c("A","C","G","T"),sep=".")
			names2=paste(s2,c("A","C","G","T"),sep=".")
			for (i in seq(1,ns)){
				sub1=snps[i,names1]
				sub2=snps[i,names2]
				sum1=rowSums(sub1)
				sum2=rowSums(sub2)
				dem=1/sum1/sum2
				dem[!is.finite(dem)]=0
				pi=sub1[,1]*sub2[,2]*dem+sub1[,1]*sub2[,3]*dem+sub1[,1]*sub2[,4]*dem+sub1[,2]*sub2[,3]*dem+sub1[,2]*sub2[,4]*dem+sub1[,3]*sub2[,4]*dem+sub2[,1]*sub1[,2]*dem+sub2[,1]*sub1[,3]*dem+sub2[,1]*sub1[,4]*dem+sub2[,2]*sub1[,3]*dem+sub2[,2]*sub1[,4]*dem+sub2[,3]*sub1[,4]*dem
				snps[i,name]=pi
			}
		}
	}
	snps
}

add_fst<-function(snps,samples){
	
	print("Position Fst calculation")
	
	nsamp=length(samples)
	names1=paste("pi",samples,sep="__") #sample pi
	n=names(snps)
	n=n[substr(n,1,4)=="pi__"]
	names2=setdiff(n,names1)  #between sample pi
	
	pi1=snps[,names1]
	pi2=snps[,names2]
	pi1a=rowSums(pi1)/nsamp
	pi2a=rowSums(pi2)/(nsamp*(nsamp-1)/2)
	snps$Fst=1-(pi1a/pi2a)
	
	snps
}

print_fst_pos<-function(snps,samplefile){

	name=paste(samplefile,"Fst_pos.txt",sep="_")
	print(paste("Writing",name))
	write.table(snps,name,sep="\t",quote=F,row.names=F)
}


###### GENES #######

#Add pi to genes in gff format
add_pi_genes_gff<-function(snps,genes,samples, globals){
	print("Pi gene calculation")
	total_length <- 0
  	ng=dim(genes)[1]
	genes$SNPs=rep(0,ng)
	for (g in seq(1,ng)){
		start=genes[g,"Start"]
		end=genes[g,"End"]
		sub=snps$Contig==as.character(genes[g,"Contig"])&snps$Pos>=start&snps$Pos<=end
		#print (sub)
		#print (sum(sub))
		genes[g,"SNPs"]=sum(sub)
	}

	for (s in samples){
		name=paste("pi__",s,sep="")	
		genes[,name]=rep(0,ng)
	}

	for (g in seq(1,ng)){
		if(!genes[g,"SNPs"]){next;}
		start=genes[g,"Start"]
		end=genes[g,"End"]
		total_length <- total_length + (end-start+1)
		sub=snps$Contig==as.character(genes[g,"Contig"])&snps$Pos>=start&snps$Pos<=end
		
		for (s in samples){
			name=paste("pi__",s,sep="")
			pi=sum(snps[sub,name])
			globals$total_pi[(globals$sample==s)] <- globals$total_pi[(globals$sample==s)] + pi
			pi=pi/(end-start+1)
			genes[g,name]=pi
		}
	}
	assign('globals',globals,envir=.GlobalEnv)
	assign('total_length',total_length,envir=.GlobalEnv)
	genes
}



add_pi_genes_fasta<-function(snps,genes,samples, globals){
  
  print("Pi gene calculation")
  total_length <- 0
  ng=dim(genes)[1]
  genes$SNPs=rep(0,ng)
  for (g in seq(1,ng)){
    sub=snps$Contig==as.character(genes[g,"Contig"])
    genes[g,"SNPs"]=sum(sub)
    total_length <- total_length + genes[g,"Length"]
  }
  count = 1
  for (s in samples){
    name=paste("pi__",s,sep="")
    genes[,name]=rep(0,ng)
    count <- count + 1
    
  }
  for (g in seq(1,ng)){
    if(!genes[g,"SNPs"]){next;}
    sub=snps$Contig==as.character(genes[g,"Contig"])
    for (s in samples){
      name=paste("pi__",s,sep="")
      pi=sum(snps[sub,name])
      globals$total_pi[(globals$sample==s)] <- globals$total_pi[(globals$sample==s)] + pi
      pi=pi/genes[g,"Length"]
      genes[g,name]=pi
    }
  }
  assign('globals',globals,envir=.GlobalEnv)
  assign('total_length',total_length,envir=.GlobalEnv)
  genes
}

add_pi_genes2_gff<-function(snps,genes,samples, globals){
  new_df <- data.frame(sample=as.character(), total_pi=as.integer())
	print("Between-sample pi gene calculation")
	ng=dim(genes)[1]
	for (i1 in seq(1,length(samples)-1)){
	  s1=samples[i1]
		for (i2 in seq(i1+1,length(samples))){
		  s2=samples[i2]
			print(c(s1,s2))
			name=paste("pi",s1,s2,sep="__")
			df_name = paste(s1,s2,sep="__")
			genes[,name]=rep(0,ng)
			new_df2 <- data.frame(sample=as.character(df_name), total_pi = as.integer(0))
			new_df <- rbind(new_df,new_df2)
		}
	}

	for (g in seq(1,ng)){
		if(!genes[g,"SNPs"]){next;}
		start=genes[g,"Start"]
		end=genes[g,"End"]
		sub=snps$Contig==as.character(genes[g,"Contig"])&snps$Pos>=start&snps$Pos<=end
		
		for (i1 in seq(1,length(samples)-1)){
			s1=samples[i1]
			for (i2 in seq(i1+1,length(samples))){
				s2=samples[i2]
				name=paste("pi",s1,s2,sep="__")
				df_name = paste(s1,s2,sep="__")
				pi=sum(snps[sub,name])
				new_df$total_pi[(new_df$sample==df_name)] <- new_df$total_pi[(new_df$sample==df_name)] + pi
				pi=pi/(end-start+1)
				genes[g,name]=pi
			}
		}
	}
	globals <- rbind(globals, new_df)
	globals$global_pi <- globals$total_pi/total_length
	assign('globals',globals,envir=.GlobalEnv)
	genes
}


add_pi_genes2_fasta<-function(snps,genes,samples, globals){
  new_df <- data.frame(sample=as.character(), total_pi=as.integer())
  print("Between-sample pi gene calculation")
  ng=dim(genes)[1]
  for (i1 in seq(1,length(samples)-1)){
    s1=samples[i1]
    for (i2 in seq(i1+1,length(samples))){
      s2=samples[i2]
      print(c(s1,s2))
      name=paste("pi",s1,s2,sep="__")
      df_name = paste(s1,s2,sep="__")
      genes[,name]=rep(0,ng)
      new_df2 <- data.frame(sample=as.character(df_name), total_pi = as.integer(0))
      new_df <- rbind(new_df,new_df2)
    }
  }
  
  for (g in seq(1,ng)){
    if(!genes[g,"SNPs"]){next;}

    sub=snps$Contig==as.character(genes[g,"Contig"])
    
    for (i1 in seq(1,length(samples)-1)){
      s1=samples[i1]
      for (i2 in seq(i1+1,length(samples))){
        s2=samples[i2]
        name=paste("pi",s1,s2,sep="__")
        df_name = paste(s1,s2,sep="__")
        pi=sum(snps[sub,name])
        new_df$total_pi[(new_df$sample==df_name)] <- new_df$total_pi[(new_df$sample==df_name)] + pi
        pi=pi/genes[g,"Length"]
        genes[g,name]=pi
      }
    }
  }
  globals <- rbind(globals, new_df)
  globals$global_pi <- globals$total_pi/total_length
  assign('globals',globals,envir=.GlobalEnv)
  genes
}


add_fst_genes<-function(genes,samples, globals){
  globals$total_fst=rep(0,dim(globals)[1])

	print("Gene Fst calculation")
	
	nsamp=length(samples)
	names1=paste("pi",samples,sep="__")
	n=names(genes)
	n=n[substr(n,1,4)=="pi__"]
	print (n)
	print (names1)
	names2=setdiff(n,names1)  #between sample pi
	pi1=genes[,names1]
	pi2=genes[,names2]
	pi1a=rowSums(pi1)/nsamp
	pi2a=rowSums(pi2)/(nsamp*(nsamp-1)/2)
	genes$Fst=1-(pi1a/pi2a)

	for (i1 in seq(1,length(samples)-1)){
	  s1=samples[i1]
		for (i2 in seq(i1+1,length(samples))){
			s2=samples[i2]  
			name=paste("Fst",s1,s2,sep="__")
			df_name = paste(s1,s2,sep="__")
			global_fst <- 1-(globals$global_pi[(globals$sample==s1)]+globals$global_pi[(globals$sample==s2)])/2/globals$global_pi[(globals$sample==paste(s1,s2,sep="__"))]
			globals$total_fst[(globals$sample==paste(s1,s2,sep="__"))] <- global_fst
			genes[,name]=1-(genes[,paste("pi",s1,sep="__")]+genes[,paste("pi",s2,sep="__")])/2/genes[,paste("pi",s1,s2,sep="__")]
		}
	}
	assign('globals',globals,envir=.GlobalEnv)
	genes
}

print_fst_genes<-function(genes,samplefile){

	name=paste(samplefile,"Fst.txt",sep="_")
	print(paste("Writing",name))
	write.table(genes,name,sep="\t",quote=F,row.names=F)
}


  
  

  
#### Arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("Usage: structure.r <samples>\nPlease, report bugs to Anne Kupczok (akupczok@ifam.uni-kiel.de)\n", call.=FALSE)
}

samplefile=args[1]
coding_seq=args[2]


samples=get_samples(samplefile)
snps=get_snps(samples)
snps=add_pi(snps,names(samples))
snps=add_pi2(snps,names(samples))
snps=add_fst(snps,names(samples))
print_fst_pos(snps,samplefile)
globals <- data.frame(sample = names(samples), total_pi=rep(0, length(samples)))
if (length(args)>1) {
  if (grepl('.gff', coding_seq)) {
	  genes=get_genes_gff(coding_seq)
	  genes=add_pi_genes_gff(snps,genes,names(samples),globals)
	  genes=add_pi_genes2_gff(snps,genes,names(samples),globals)
  } else if (grepl('.fasta', coding_seq)){
    genes=get_genes_fasta(coding_seq)
    genes=add_pi_genes_fasta(snps,genes,names(samples),globals)
    genes=add_pi_genes2_fasta(snps,genes,names(samples), globals)
  }
  genes=add_fst_genes(genes,names(samples), globals)
  print (globals)
  print_fst_genes(genes,samplefile)
}
