# Author: Raunak Shrestha
# Last Updated: 01 June 2019
# runGSEAfisher.R
# GeneSet Enrichmet Analysis based on one-tail Fisher's Exact Test (Hypergeometric test)

### LOAD LIBRARIES ---
require("stringr")

### DEFINE PATH ---
dir.wrk <- getwd()
dir.geneset <- file.path(dir.wrk, "genesets")
dir.output <- file.path(dir.wrk, "enrichment")

### CREAT OUTPUT FOLDER ---
dir.create(dir.output, showWarnings=FALSE)
	
### FUNCTION: PARSE GMT GENESETS ---
parseGMT <- function(dir.gmt){
	cat(paste(Sys.time()), "PARSE GMT: START ...","\n", sep="\t")
	
	# GET GENESET FILES ---
	file.genesetdb <- list.files(path=dir.gmt, pattern=".gmt", full.names = TRUE)
	file.names <- list.files(path=dir.gmt, pattern=".gmt", full.names = FALSE)

	# PARSE EACH FILE AT A TIME ---
	for(j in 1:length(file.genesetdb)){
	     
		# LOAD GENESET DATA ---
		no_col <- max(count.fields(file.genesetdb[j],  sep="\t"))
		dat.genesetdb <- read.delim(file.genesetdb[j], header=FALSE, fill=TRUE, col.names=1:no_col, as.is=paste("X",1:no_col, sep=""))		
	
		# PARSE EACH GENESETS ENTRY ---
		dat.genesets <- data.frame(Category=dat.genesetdb[,1])
		for(i in 1:nrow(dat.genesetdb)){
			elements <- as.character(dat.genesetdb[i,3:ncol(dat.genesetdb)])
			elements <- paste(elements[which(elements != "")], collapse=":")
  
			dat.genesets$Genesets[i] <- elements
			dat.genesets$Description[i] <- dat.genesetdb[i,2]
		} 
		
	     # WRITE OUTPUT ---
		file.output <- file.path(dir.gmt, gsub(".gmt", ".tsv", file.names[j]))
		write.table(dat.genesets, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
		
		cat("FILE GENERATED:", file.output, "\n", sep="\t")
	}
	
	cat(paste(Sys.time()), "PARSE GMT: END ...","\n", sep=" ")
}


### FUNCTION: FISHER EXACT TEST ---
get.fisher.exact.test <- function(dat.genesets, genes.queryset, genes.refset, ct){  
	
     # CREATE DATAFRAME ---
     dat.result <- data.frame(Category=dat.genesets$Category)
	
     # COMPUTE FISHER EXACT TEST FOR EACH GENESET ---
	for(i in 1:nrow(dat.genesets)){
		genes.genesets <- sort(unique(unlist(str_split(dat.genesets$Genesets[i], ":"))), decreasing=FALSE)
		genes.interest <- intersect(genes.queryset, genes.genesets)
  
		yy <- length(intersect(genes.genesets, genes.interest))
		yn <- length(intersect(genes.genesets, setdiff(genes.refset, genes.interest)))
		ny <- length(intersect(setdiff(genes.refset, genes.genesets), genes.interest))
		nn <- length(intersect(setdiff(genes.refset,genes.interest), setdiff(genes.refset, genes.genesets)))
  
		fisherRes <- fisher.test(rbind(c(yy,yn),c(ny,nn)), alternative="greater")
		
		dat.result$pvalue[i] <- fisherRes$p.value
		dat.result$fdr[i] <- NA

		dat.result$overlap.percent[i] <- round(length(genes.interest)/length(genes.genesets) * 100, digit=2)
		dat.result$overlap.genes[i] <- paste(genes.interest, collapse=":")
	}

	# MULTIPLE TEST CORRECTION ---
	dat.result$fdr <- p.adjust(dat.result$pvalue, method="BH")

     # ADD DESCRIPTION LINK ---
	dat.result$Description <- dat.genesets$Description
	
	# RE-ARRANGE DATA ---
	dat.result <- dat.result[order(dat.result$pvalue, decreasing=FALSE),]
	
	# CUT-OFF BY P-VALUE ---
	dat.result <- subset(dat.result, dat.result$fdr <= ct)
  
	# FORMAT P-VALUES ---
	dat.result$pvalue <- format(dat.result$pvalue, scientific = TRUE, digits = 4)
	dat.result$fdr <- format(dat.result$fdr, scientific = TRUE, digits = 4)
	
	return(dat.result)
}  


### FUNCTION: GENESET ENRICHMENT TEST ---
get.enrichment.test <- function(file.query, batch, file.background, dir.gmt, threshold, dir.out=dir.output){
	cat(paste(Sys.time()), "START ENRICHMENT ANALYSIS ...","\n", sep="\t")
	
	# DEFINE PATH ---
	dir.batch <- file.path(dir.out, batch) 
	dir.create(dir.batch, showWarnings=FALSE)
	
	# GET GENESETS FILES ---
	file.genesetdb <- list.files(path=dir.gmt, pattern=".tsv", full.names = TRUE)
	file.names <- list.files(path=dir.gmt, pattern=".tsv", full.names = FALSE)

	# LOAD REFERENCE GENESET AND QUERY GENESET ---
	genes.refset <- read.table(file.background, header=FALSE, stringsAsFactors=FALSE)$V1
	genes.queryset <- read.table(file.query, header=FALSE, stringsAsFactors=FALSE)$V1

	# COMPUTE ENRICHMENT ---
	for(j in 1:length(file.genesetdb)){
		
	     # LOAD GENESETS ---
		dat.genesets <- read.delim(file.genesetdb[j], header=TRUE, stringsAsFactor=FALSE)
	
		# FISHER EXACT TEST ---
		dat.result <- get.fisher.exact.test(dat.genesets, genes.queryset, genes.refset, ct=threshold)
	
		if(nrow(dat.result) != 0){ 
			
		     # GENERATE OUTPUT ---
			file.output <- file.path(dir.batch,  paste("enrichment", file.names[j], sep="_"))
			write.table(dat.result, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
			
			cat("FILE GENERATED:", file.output, "\n", sep=" ")
		}
	}

	cat(paste(Sys.time()), "END ENRICHMENT ANALYSIS ...","\n", sep="\t")
}
