# Author: Raunak Shrestha
# Last Updated: 11 June 2016
# runGSEAfisher.R
# GeneSet Enrichmet Analysis based on one-tail Fisher's Exact Test

### LoadLibraries -----------------------------------------------------------
require("stringr")

### SetDirectories ----------------------------------------------------------
dir.wrk <- getwd()
dir.geneset <- file.path(dir.wrk, "genesets")
dir.refset <- file.path(dir.wrk, "backgroundset")
dir.output <- file.path(dir.wrk, "enrichment")

dir.create(dir.output, showWarnings=FALSE)
	
### Parse GMT Genesets -------------------------------------------------------
parseGMT <- function(gmt.name, dir.gmt=dir.geneset){
	cat(paste(Sys.time()), "START PARSE GMT ...","\n", sep="\t")
	
	dir.batch <- file.path(dir.gmt, gmt.name)
	
	file.genesetdb <- list.files(path=dir.batch, pattern=".gmt", full.names = T)
	file.names <- list.files(path=dir.batch, pattern=".gmt", full.names = F)

	for(j in 1:length(file.genesetdb)){
		### Load GeneSet Data
		no_col <- max(count.fields(file.genesetdb[j],  sep="\t"))
		dat.genesetdb <- read.delim(file.genesetdb[j], header=FALSE, fill=TRUE, col.names=1:no_col, as.is=paste("X",1:no_col, sep=""))		
	
		dat.genesets <- data.frame(Category=dat.genesetdb[,1])
		for(i in 1:nrow(dat.genesetdb)){
			elements <- as.character(dat.genesetdb[i,3:ncol(dat.genesetdb)])
			elements <- paste(elements[which(elements != "")], collapse=",")
  
			dat.genesets$Genesets[i] <- elements
			dat.genesets$Description[i] <- dat.genesetdb[i,2]
		} 
	
		file.output <- file.path(dir.batch, gsub(".gmt",".txt", file.names[j]))
		write.table(dat.genesets, file.output, sep="\t", row.names=F, col.names=T, quote=F)
		
		cat("FILE GENERATED:", file.output, "\n", sep="\t")
	}
	cat(paste(Sys.time()), "END PARSE GMT ...","\n", sep=" ")
}


#### FISHER EXACT TEST
get.fisher.exact.test <- function(dat.genesets, genes.queryset, genes.refset, ct){  
	dat.result <- data.frame(Category=as.character(dat.genesets$Category))
	for(i in 1:nrow(dat.genesets)){
		genes.genesets <- sort(unique(unlist(str_split(dat.genesets$Genesets[i], ","))),decreasing=F)
		genes.interest <- intersect(genes.queryset, genes.genesets)
  
		yy <- length(intersect(genes.genesets, genes.interest))
		yn <- length(intersect(genes.genesets, setdiff(genes.refset, genes.interest)))
		ny <- length(intersect(setdiff(genes.refset, genes.genesets), genes.interest))
		nn <- length(intersect(setdiff(genes.refset,genes.interest), setdiff(genes.refset, genes.genesets)))
  
		fisherRes <- fisher.test(rbind(c(yy,yn),c(ny,nn)), alternative="greater")
		dat.result$pvalue[i] <- fisherRes$p.value
		dat.result$fdr[i] <- 0

		dat.result$overlap.percent[i] <- round(length(genes.interest)/length(genes.genesets) * 100, digit=2)
		dat.result$overlap.genes[i] <- paste(genes.interest, collapse=",")
	}

	# Multiple-test correction
	dat.result$fdr <- p.adjust(dat.result$pvalue, method="BH")
	
	dat.result$Description <- dat.genesets$Description
	dat.result <- dat.result[order(dat.result$pvalue, decreasing=F),]
	dat.result <- subset(dat.result, dat.result$fdr <= ct)
  
	dat.result$pvalue <- format(dat.result$pvalue, scientific = T, digits = 4)
	dat.result$fdr <- format(dat.result$fdr, scientific = T, digits = 4)
	
	return(dat.result)
}  


### GeneSet Enrichment Test -------------------------------------------------------
get.enrichment.test <- function(file.query, batch, file.background, dir.gmt, threshold, dir.out=dir.output){
	cat(paste(Sys.time()), "START ENRICHMENT ANALYSIS ...","\n", sep="\t")
	
	#Set Directories
	dir.batch <- file.path(dir.out, batch) 
	dir.create(dir.batch, showWarnings=FALSE)
	
	#Get geneset files
	file.genesetdb <- list.files(path=dir.gmt, pattern=".txt", full.names = T)
	file.names <- list.files(path=dir.gmt, pattern=".txt", full.names = F)

	#Load ReferenceGeneset and QueryGeneset
	genes.refset <- read.table(file.background, header=F, as.is="V1")$V1
	genes.queryset <- read.table(file.query,header=F, as.is="V1")$V1

	for(j in 1:length(file.genesetdb)){
		#LOAD GENESETS
		dat.genesets <- read.delim(file.genesetdb[j], header=T, stringsAsFactor=F)
	
		#FISHER EXACT TEST
		dat.result <- get.fisher.exact.test(dat.genesets, genes.queryset, genes.refset, ct=threshold)
	
		if(nrow(dat.result) != 0){ 
			#GENERATE OUTPUT
			file.output <- file.path(dir.batch,  paste("enrichment",file.names[j],sep="_"))
			write.table(dat.result, file.output, sep="\t", row.names=F, col.names=T, quote=F)
			cat("FILE GENERATED:", file.output, "\n", sep=" ")
		}
	}

	cat(paste(Sys.time()), "END ENRICHMENT ANALYSIS ...","\n", sep="\t")
}
