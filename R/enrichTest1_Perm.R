enrichTest1_Perm<-function(genesSize, hitSize=0, setType=2,times=100, seed=1) {
	
	setdb=system.file("extdata", "msigdb.db", package="snpGeneSets")	
	if(!file.exists(setdb)) 
		stop("MSigDB database is not existing, please use msigdb_build function to install msigdb gene sets first!\n")
	
	if(missing(genesSize)) stop("Need to specify genesSize value")
	
	if(!((is.integer(genesSize)|is.numeric(genesSize)) & (genesSize>=hitSize) )) 
		stop(paste("genesSize has to be a integer or numeric value >=",hitSize))	
	if(!((is.integer(setType)|is.numeric(setType))&(length(setType)==1))) 
		stop("setType has to be a number from 0 to 20")
	if(setType>20 | setType<0) setType=0
	
	setInfo=getSetType(); setSym=with(setInfo, symbol[id==setType])	
	
	require(RSQLite)
	
	if(setType==0){
		pid_geneid=sqlQuery("select distinct pid,  geneid from geneSets",  setdb)	
	} else {
		pid_geneid=sqlQuery(paste("select distinct pid, geneid from geneSets where", setSym, "=1"),  setdb)
	}		
	pid_geneid_list<- sapply(pid_geneid$geneid, function(x) as.integer(strsplit(x,":")[[1]][-1]))
	names(pid_geneid_list)=pid_geneid$pid
	
	
	setGenes=getSetGenes(setType)
	rst=matrix(NA, nrow=nrow(pid_geneid), ncol=times, dimnames = list(pid_geneid$pid, 1:times))	
	
	n_balls=length(setGenes)
	n_whites=genesSize	
	i_drawn=sapply(pid_geneid_list, length) #	
	
	set.seed(seed)
	for (i_time in 1:times) {
		cat("Permutation", i_time, "\n")
		gene_sample<-sample(setGenes, n_whites) #
		i_drawnWhites=sapply(pid_geneid_list, function(x) {					
					sum(gene_sample %in% x)
				}) #
		rst[, i_time]=phyper(i_drawnWhites-1, n_whites, 
				n_balls-n_whites, i_drawn, lower.tail = FALSE)
		if(hitSize>0) rst[i_drawnWhites<hitSize,i_time]=NA		
		#gc()
	}	
	return(rst)
}
