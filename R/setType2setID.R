setType2setID<-function(setType, setSize=NULL) {
	
	setdb=system.file("extdata", "msigdb.db", package="snpGeneSets")	
	if(!file.exists(setdb)) 
		stop("MSigDB database is not existing, please use msigdb_build function to install msigdb gene sets first!\n")
	
	if(missing(setType)) stop("Need to specify setType value")
	if(!((is.integer(setType)|is.numeric(setType))&(length(setType)==1)))  
		stop("setType has to be a number from 0 to 20")
	if(setType>20 | setType<0) setType=0
	if(!file.exists(setdb)) stop("setdb does not exist")
	
	setInfo=getSetType(); setSym=with(setInfo, symbol[id==setType])	
	
	pids=NA
	if(setType==0) {
		pids=sqlQuery("select distinct pid from geneSets",setdb)[,1]
		if(!is.null(setSize)) {
			pids=sqlQuery(paste("select distinct pid from geneSets where", 
							paste("size", "IN (",paste(setSize, collapse=","),")",sep=" "), sep=" "),
					setdb)[,1]
		}
	} else {
		pids=sqlQuery(paste("select distinct pid from geneSets where", setSym, "=", 1, sep=" "),setdb)[,1]
		if(!is.null(setSize)) {
			pids=sqlQuery(paste("select pid from geneSets where", setSym, "=", 1,
							paste("AND size", "IN (",paste(setSize, collapse=","),")",sep=" "), sep=" "),
					setdb)[,1]
		}		
	}
	
	return(pids)	
}
