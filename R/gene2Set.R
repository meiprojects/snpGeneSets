gene2Set<-function(geneID, setType=2) {
	
	setdb=system.file("extdata", "msigdb.db", package="snpGeneSets")	
	if(!file.exists(setdb)) 
		stop("MSigDB database is not existing, please use msigdb_build function to install msigdb gene sets first!\n")
	
	if(missing(geneID)) stop("Need to specify a geneID value")	
	
	if(!((is.integer(geneID)|is.numeric(geneID))&(length(geneID)==1))) 
		stop("The geneID has to be a integer or numeric Entrez gene ID")
	if(!((is.integer(setType)|is.numeric(setType))&(length(setType)==1)))  
		stop("setType has to be a number from 0 to 20")
	if(setType>20 | setType<0) setType=0
	if(!file.exists(setdb)) stop("setdb does not exist")
	
	require(RSQLite)
	
	setInfo=getSetType(); setSym=with(setInfo, symbol[id==setType])		
	
	if(setType==0) {
		setRst=sqlQuery(paste("select pid from geneSets where geneid LIKE", shQuote(paste("%:",geneID,":%",sep="")), sep=" "),setdb)[,1]		
	} else {
		setRst=sqlQuery(paste("select pid from geneSets where", setSym, "=", 1,  "AND geneid LIKE", shQuote(paste("%:",geneID,":%",sep="")), sep=" "),setdb)[,1]		
	}
	return(setRst)
}
