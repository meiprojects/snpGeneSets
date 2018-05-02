snp2Gene<-function(snps, up=2000, down=2000, GRCh=37) {
	snpmaps=getSNPMap(snps, GRCh=GRCh)
	snpDF=snpmaps$rsid_map
	other=snpmaps$other
	
	if(nrow(snpDF)==0) return(list(map=NULL, other=other))	
	
	if(GRCh==37) {
		genedb=system.file("extdata", "gene105.db", package="snpGeneSets")		
	} else if (GRCh==38) {
		genedb=system.file("extdata", "gene106.db", package="snpGeneSets")
	} else {
		stop("The reference genome build, GRCh, has to be 37 or 38")
	}
	if(!file.exists(genedb)) stop("genedb does not exist!")
	
	require(RSQLite)
	driver <- dbDriver("SQLite")
	connect <- dbConnect(driver, dbname=genedb);
	closeup <- function(){
		dbDisconnect(connect)
	}
	
	rst=apply(snpDF,1, function(isnp) {
				chr=isnp["chr"]; pos=as.integer(isnp["pos"]); 
				sql=paste("select gene_id from genemaps INDEXED BY idx_pos where chr=", 
						shQuote(chr),"AND start<=", pos+up, "AND end>=", pos-down, sep=" ")
				igene=tryCatch(dbGetQuery(connect, sql), finally=closeup)
				if(nrow(igene)==0) {
					return(NULL)
				} else {
					cbind(snp=t(isnp["snp"]), igene, stringsAsFactors=F)
				}				
			})	
	if(!is.null(rst))
		rst=as.data.frame(do.call("rbind",rst[!sapply(rst, is.null)]))	
	dbDisconnect(connect)
	return(list(map=rst, other=c(other, snpDF$snp[!(snpDF$snp %in% rst$snp)])))
}
