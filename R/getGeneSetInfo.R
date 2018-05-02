getGeneSetInfo<-function(setID) {
	
	setdb=system.file("extdata", "msigdb.db", package="snpGeneSets")	
	if(!file.exists(setdb)) 
		stop("MSigDB database is not existing, please use msigdb_build function to install msigdb gene sets first!\n")
	
	if(missing(setID)) stop("Need to specify setID value")
	if(!(is.numeric(setID)|is.integer(setID))) 
		stop("setID has to be a numeric or integer value")
	if(!file.exists(setdb)) stop("setdb does not exist!")
	require(RSQLite)
	set_rst=sqlQuery(paste("select * from geneSets INDEXED BY idx_pid where pid=",setID), setdb)
	geneid=strsplit(set_rst$geneid,":")[[1]]
	geneid=as.integer(geneid[!is.na(geneid)&geneid!=""])
	set_type=as.logical(set_rst[getSetType()[-1,2]])
	names(set_type)=names(set_rst[getSetType()[-1,2]])
	list(setID=set_rst$pid, set_name=set_rst$pathway, set_link=set_rst$link, set_type=set_type, set_geneid=geneid)	
}
