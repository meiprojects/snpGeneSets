getSetGenes<-function(setType=2) {
	
	setdb=system.file("extdata", "msigdb.db", package="snpGeneSets")	
	if(!file.exists(setdb)) 
		stop("MSigDB database is not existing, please use msigdb_build function to install msigdb gene sets first!\n")
	
	if(!((is.integer(setType)|is.numeric(setType))&(length(setType)==1)&(setType>=0 & setType<=20))) 
		stop("setType has to be a number from 0 to 20")
	setInfo=getSetType();
	return(switch((setType>=1&setType<=20)+1, 
					sqlQuery("select distinct gene_id from geneInfo",  setdb)[,1],
					sqlQuery(paste("select distinct gene_id from geneInfo where", with(setInfo, symbol[id==setType]), "=1"),  setdb)[,1]
			))
}
