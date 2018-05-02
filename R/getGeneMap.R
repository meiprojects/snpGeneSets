getGeneMap<-function(gene, isGeneID=TRUE,GRCh=37) {
	
	if(missing(gene)) stop("Need to specify the gene value")	
	
	if(GRCh==37) {
		genedb=system.file("extdata", "gene105.db", package="snpGeneSets")		
	} else if (GRCh==38) {
		genedb=system.file("extdata", "gene106.db", package="snpGeneSets")
	} else {
		stop("The reference genome build, GRCh, has to be 37 or 38")
	}
	if(!file.exists(genedb)) stop("genedb does not exists")
	
	require(RSQLite)
	
	if(isGeneID) {
		if(!(is.numeric(gene)|is.integer(gene))) 
			stop("When isGeneID=TRUE, the gene has to be numeric or integer vector of entrez Gene ID")	
		
		gene_map<-sqlQuery(paste("select * from genemaps INDEXED BY idx_id where", 
						paste("gene_id", "IN (",paste(shQuote(gene), collapse=","),")",sep=" "), sep=" "),
				genedb)
		
		other=gene[!(gene %in% gene_map$gene_id)]
	} else {
		if(!is.character(gene)) 
			stop("When isGeneID=FALSE, the gene has to be a character vector of HGNC gene symbol")	
		
		gene=toupper(gene)
		gene_map<-sqlQuery(paste("select * from genemaps INDEXED BY idx_name where", 
						paste("gene_name", "IN (",paste(shQuote(gene), collapse=","),")",sep=" "), sep=" "),
				genedb)
		
		other=gene[!(gene %in% gene_map$gene_name)]
	}	
	return(list(gene_map=gene_map, other=other))		
}
