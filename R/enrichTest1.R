enrichTest1<-function(genes, hitSize=0, setType=2) {
	
	setdb=system.file("extdata", "msigdb.db", package="snpGeneSets")	
	if(!file.exists(setdb)) 
		stop("MSigDB database is not existing, please use msigdb_build function to install msigdb gene sets first!\n")
	
	if(missing(genes)) stop("Need to specify genes value")
	if(!file.exists(setdb)) stop("setdb does not exist")
	if(!((is.integer(genes)|is.numeric(genes))&(length(genes)>=hitSize))) 
		stop(paste("genes has to be a integer or numeric value >=",hitSize))	
	if(!((is.integer(setType)|is.numeric(setType))&(length(setType)==1))) 
		stop("setType has to be a number from 0 to 20")
	if(setType>20 | setType<0) setType=0
	setInfo=getSetType(); setSym=with(setInfo, symbol[id==setType])	
	
	require(RSQLite)
	
	if(setType==0){
		pid_geneid=sqlQuery("select distinct pid,size,  geneid from geneSets",  setdb)	
	} else {
		pid_geneid=sqlQuery(paste("select distinct pid, size, geneid from geneSets where", setSym, "=1"),  setdb)
	}
	pid_geneid_list<- sapply(pid_geneid$geneid, function(x) as.integer(strsplit(x,":")[[1]][-1]))
	names(pid_geneid_list)=pid_geneid$pid
	
	setGenes=getSetGenes(setType)	
	n_balls=length(setGenes)
	existing=genes[genes %in% setGenes]
	n_whites=length(existing)
	#n_drawn = length(existing)		
	
	pid_drawn=sapply(pid_geneid_list, length)
	pid_drawnWhites=sapply(pid_geneid_list, function(x) sum(x %in% existing))
	p0=n_whites/n_balls
	pid_stat=pid_drawnWhites/pid_drawn-p0
	pid_sd=sqrt(p0*(1-p0)/pid_drawn)
	pid_p=phyper(pid_drawnWhites-1, n_whites, n_balls-n_whites, pid_drawn, lower.tail = FALSE)
	na_id= (pid_drawn==0) | (pid_drawnWhites<hitSize)
	if(sum(na_id)>0) {
		pid_stat[na_id]<-NA; pid_sd[na_id]<-NA; pid_p[na_id]<-NA;
	}
	#pid_geneid$genes<-pid_drawn
	pid_geneid$genesSize<-pid_drawnWhites
	pid_geneid$effect<-pid_stat
	pid_geneid$sd<-pid_sd
	pid_geneid$pval<-pid_p	
	setInfo=setInfo[setInfo$id==setType,]; rownames(setInfo)=NULL
	
	return(list(enrich_test=pid_geneid[,-3], useGenes=existing, nGenes=n_balls, 
					nTopGenes=n_whites, setTypeInfo =as.list(setInfo)))
}	
