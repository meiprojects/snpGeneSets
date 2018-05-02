enrichTest2<-function(geneDF, cut_score=0.05, hitSize=0, setType=2) {
	
	setdb=system.file("extdata", "msigdb.db", package="snpGeneSets")	
	if(!file.exists(setdb)) 
		stop("MSigDB database is not existing, please use msigdb_build function to install msigdb gene sets first!\n")
	
	if(missing(geneDF)) stop("Need to specify the geneDF value")
	if(!(class(geneDF) %in% c("data.frame","matrix","array" ))) 
		stop("geneDF has to be a data frame or matrix or array")
	if(class(geneDF)!="data.frame") geneDF<-as.data.frame(geneDF)
	if (!all(c("gene_id", "score") %in% colnames(geneDF))) stop("the geneDF has to contain column names of 'gene_id' and 'score'")
	
	
	if(!((is.numeric(geneDF[,"gene_id"])|is.integer(geneDF[,"gene_id"])) & (is.numeric(geneDF[,"score"]))))
		stop("The 'gene_id' column has to be integer or numeric vector of Entrez gene ID and the 'score' column has to be numeric value of p value")	
	if(!((is.integer(setType)|is.numeric(setType))&(length(setType)==1))) 
		stop("setType has to be a number from 0 to 20")
	if(setType>20 | setType<0) setType=0				
	
	setGenes=getSetGenes(setType)	
	geneDF=geneDF[geneDF[,"gene_id"] %in% setGenes,]
	geneDF=geneDF[!is.na(geneDF[,"score"]),]	
	geneDF[,"score"]=uscore(geneDF[,"score"])	
	
	n_balls=nrow(geneDF)
	sigGenes=geneDF$gene_id[geneDF$score<=cut_score]
	n_whites=length(sigGenes)
	
	setInfo=getSetType(); setSym=with(setInfo, symbol[id==setType])
	if(setType==0){
		pid_geneid=sqlQuery("select distinct pid, size, geneid from geneSets",  setdb)	
	} else {
		pid_geneid=sqlQuery(paste("select distinct pid, size, geneid from geneSets where", setSym, "=1"),  setdb)
	}	
	pid_geneid_list<- sapply(pid_geneid$geneid, function(x) as.integer(strsplit(x,":")[[1]][-1]))
	names(pid_geneid_list)=pid_geneid$pid
	
	pid_drawn=sapply(pid_geneid_list, function(x) sum(x %in% geneDF$gene_id))
	pid_drawnWhites=sapply(pid_geneid_list, function(x) sum(x %in% sigGenes))
	p0=n_whites/n_balls
	pid_stat=pid_drawnWhites/pid_drawn-p0
	pid_sd=sqrt(p0*(1-p0)/pid_drawn)
	pid_p=phyper(pid_drawnWhites-1, n_whites, n_balls-n_whites, pid_drawn, lower.tail = FALSE)
	na_id= (pid_drawn==0) | (pid_drawnWhites<hitSize) | (n_whites==0)
	#na_id=(pid_drawn<=hitSize)
	if(sum(na_id)>0) {
		pid_stat[na_id]<-NA; pid_sd[na_id]<-NA; pid_p[na_id]<-NA;
	}
	pid_geneid$genes<-pid_drawn
	pid_geneid$sigGenes<-pid_drawnWhites
	pid_geneid$effect<-pid_stat
	pid_geneid$sd<-pid_sd
	pid_geneid$pval<-pid_p
	
	setInfo=setInfo[setInfo$id==setType,]; rownames(setInfo)=NULL
	return(list(enrich_test=pid_geneid[,-3], useGenes=geneDF[,"gene_id"], 
					nGenes=n_balls, nSigGenes=n_whites, setTypeInfo =as.list(setInfo)))
}
