getSNPMap<-function(snps,GRCh=37) {
	if(missing(snps)) stop("snps cannot be empty")
	if(!is.character(snps)) stop("snps has to be a character vector")
	
	if(GRCh==37) {
		snpdb=system.file("extdata", "snp138.db", package="snpGeneSets")		
	} else if (GRCh==38) {
		snpdb=system.file("extdata", "snp142.db", package="snpGeneSets")
	} else {
		stop("The reference genome build, GRCh, has to be 37 or 38")
	}	
	
	snps<-tolower(snps)	
	rsid_map<-sqlQuery(paste("select * from snpmaps INDEXED BY idx_snp where", 
					paste("snp", "IN (",paste(shQuote(snps), collapse=","),")",sep=" "), sep=" "),
			snpdb)
	other=snps[!(snps %in% rsid_map$snp)]
	return(list(rsid_map=rsid_map, other=other))	
}
