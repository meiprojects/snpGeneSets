aligator<-function(snpGeneP,setType=2, 
		Nsample=5000, Btimes=1000, pcut=0.05, seed=1) {
	
	setdb=system.file("extdata", "msigdb.db", package="snpGeneSets")	
	if(!file.exists(setdb)) 
		stop("MSigDB database is not existing, please use msigdb_build function to install msigdb gene sets first!\n")
	
	if(missing(snpGeneP)) stop("Need to specify the snpGeneP value")
	if(class(snpGeneP)!="data.frame") stop("snpGeneP has to be a data frame")
	if (!all(c("snp", "gene_id", "p") %in% colnames(snpGeneP))) 
		stop("the geneDF has to contain column names of 'snp', 'gene_id' and 'p'")
	colT=(class(snpGeneP$snp)=="character")&
			(class(snpGeneP$gene_id)=="integer"|class(snpGeneP$gene_id)=="numeric")&
			(class(snpGeneP$p)=="numeric")
	if(!colT) 
		stop("the 'snp', 'gene_id' and 'p' of snpGeneP should be 'character', 'integer'/'numeric' and 'numeric'")
	
	#sapply(colnames(snpGeneP), function(x) class(snpGeneP[[x]]))==c("character", "integer", "")
			
	setGenes=getSetGenes(setType)	
	snpGeneP=snpGeneP[snpGeneP[,"gene_id"] %in% setGenes,]
	snpGeneP=snpGeneP[!is.na(snpGeneP[,"p"]),]
		
	sigid <- which(snpGeneP$p <= pcut)
	sig_gene<-unique(snpGeneP$gene_id[sigid])
	sig_gene_N<-length(sig_gene) #168	
		
	setInfo=getSetType(); setSym=with(setInfo, symbol[id==setType])
	if(setType==0){
		pid_geneid=sqlQuery("select distinct pid,  geneid from geneSets",  setdb)	
	} else {
		pid_geneid=sqlQuery(paste("select distinct pid, geneid from geneSets where", setSym, "=1"),  setdb)
	}	
	pid_geneid_list<- sapply(pid_geneid$geneid, function(x) as.integer(strsplit(x,":")[[1]][-1]))	
	names(pid_geneid_list)=pid_geneid$pid
	
	path_count<-sapply(pid_geneid_list, function(x) sum(x %in% sig_gene))
		
	set.seed(seed)
	cat("generate random", Nsample, "samples:\t")
	t1=Sys.time()
	sig_gene_list <- lapply(1:Nsample, function(x) sampleGenebySNP(snpGeneP[,c("snp","gene_id")],sig_gene_N))
	t2=Sys.time()
	cat(as.numeric(t2-t1,units="mins"),"minutes\n")
		
	sim_count<-matrix(0, nrow=length(pid_geneid_list), ncol=length(sig_gene_list))	
	rownames(sim_count)=names(pid_geneid_list)
	path_N=length(pid_geneid_list)
	for (i in 1:path_N) {
		cat("analyzing pathway i",i,"/",path_N, "in" , Nsample, "random samples:\t")
		t1=Sys.time()		
		path<-pid_geneid_list[[i]]
		sim_count[i,]<-sapply(sig_gene_list, function(x) sum(path %in% x))
		t2=Sys.time()
		cat(as.numeric(t2-t1,units="mins"),"minutes\n")
	}
	
	#cat("***Calculating pathway-specific unadjusted p-value***\n")
	path_unadj_p<-apply(cbind(path_count, sim_count), 1, function(x) {
				mean(x[-1]>=x[1])
			})	
	
	#cat("***Calculating  permutation-adjusted p-value in random samples***\n")
	perm_table<-NULL
	sim_count_N=ncol(sim_count)
	boot_n=sim_count_N #the size of every bootstrap sample
	for (i in 1:Btimes) {
		cat("Bootstrap", i, "/", Btimes, ":")
		t1=Sys.time()
		obi=sample(sim_count_N,1)
		booti<-sample((1:sim_count_N)[-obi],boot_n, replace=TRUE )
		pi<-apply(sim_count,1, function(x) {
					mean(x[booti]>=x[obi])
				})
		perm_table<-c(perm_table, min(pi))
		t2=Sys.time()
		cat(as.numeric(t2-t1,units="mins"),"minutes\n")
	}
	path_adj_p<-sapply(path_unadj_p, function(x) mean(perm_table<=x))
	
	rst=data.frame(pid=as.integer(names(pid_geneid_list)), p=path_unadj_p, adj_p=path_adj_p)
	rownames(rst)=NULL	
	return(rst)	
}


