setID2setType<-function(setID) {
	
	setdb=system.file("extdata", "msigdb.db", package="snpGeneSets")	
	if(!file.exists(setdb)) 
		stop("MSigDB database is not existing, please use msigdb_build function to install msigdb gene sets first!\n")
	
	if(missing(setID)) stop("Need to specify setType value")
	if(!((is.integer(setID)|is.numeric(setID))&(length(setID)==1)))  
		stop("setID has to be a integer or numeric value")
	if(!file.exists(setdb)) stop("setdb does not exist")
	
	tids=sqlQuery(paste("select * from geneSets where pid=", setID), setdb)[,6:24]	
	tids=which(tids==1)
	if(length(tids)>0) tids=c(0,tids)	
	return(tids)	
}
