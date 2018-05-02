enrichTest2_Perm<-function(geneDF=NULL, cut_score=0.05, hitSize=0, 
		setType=2, times=1000, seed=1, snpWeightDF=NULL) {
	
	setdb=system.file("extdata", "msigdb.db", package="snpGeneSets")	
	if(!file.exists(setdb)) 
		stop("MSigDB database is not existing, please use msigdb_build function to install msigdb gene sets first!\n")
	
	if(!((is.integer(setType)|is.numeric(setType))&(length(setType)==1))) 
		stop("setType has to be a number from 0 to 20")
	if(setType>20 | setType<0) setType=0
	
	if(!is.null(geneDF)) {
		if(is.data.frame(geneDF)) {
			if (!all(c("gene_id", "score") %in% colnames(geneDF))) stop("the geneDF has to contain column names of 'gene_id' and 'score'")
		} else {
			stop("the geneDF has to be a data frame")
		}
		if(!(is.numeric(geneDF$gene_id)|is.integer(geneDF$gene_id))) {
			stop("the gene_id column of geneDF has to be a numeric or integer vector")
		}
		if(!is.numeric(geneDF$score)) {
			stop("the score column of geneDF has to be a numeric")
		}
	}
	if(!is.null(snpWeightDF)) {
		if(class(snpWeightDF)!="data.frame") 
			stop("snpWeightDF has to be a data frame")
		if (!all(c("snp", "gene_id") %in% colnames(snpWeightDF))) 
			stop("the snpWeightDF has to contain column names of 'snp' and 'gene_id'")		
		if(!is.character(snpWeightDF$snp))
			stop("the snp column of snpWeightDF has to be character vector")
		if(!(is.numeric(snpWeightDF$gene_id)|is.integer(snpWeightDF$gene_id))) {
			stop("the gene_id column of snpWeightDF has to be a numeric or integer vector")
		}
	}
	
	set.seed(seed)		
	setInfo=getSetType(); setSym=with(setInfo, symbol[id==setType])		
	if(setType==0){
		pid_geneid=sqlQuery("select distinct pid,  geneid from geneSets",  setdb)	
	} else {
		pid_geneid=sqlQuery(paste("select distinct pid, geneid from geneSets where", setSym, "=1"),  setdb)
	}	
	pid_geneid_list<- sapply(pid_geneid$geneid, function(x) as.integer(strsplit(x,":")[[1]][-1]))	
	names(pid_geneid_list)=pid_geneid$pid
	
	setGenes=getSetGenes(setType)
	if(!is.null(geneDF)) {
		geneDF=geneDF[!is.na(geneDF[,"score"]),]
		geneDF=geneDF[geneDF$gene_id %in% setGenes,]
		geneDF[,"score"]=uscore(geneDF[,"score"])
		n_whites=sum(geneDF$score<=cut_score)
		n_balls=nrow(geneDF)
		geneids=geneDF$gene_id 
		i_drawn=sapply(pid_geneid_list, function(x) sum(x %in% geneids))
		if(!is.null(snpWeightDF)) {
			snpWeightDF=merge(data.frame(gene_id=geneDF$gene_id),snpWeightDF,all.x=TRUE,all.y=FALSE)[,c("snp","gene_id")]
		}
	} else {
		n_balls=length(setGenes)
		n_whites=round(n_balls*cut_score)
		geneids=setGenes
		i_drawn=sapply(pid_geneid_list, length)
		if(!is.null(snpWeightDF)) {
			snpWeightDF=merge(data.frame(gene_id=setGenes),snpWeightDF,all.x=TRUE,all.y=FALSE)[,c("snp","gene_id")]
		}
	}
	i_drawn[i_drawn==0]=NA
	
	rst=matrix(NA, nrow=nrow(pid_geneid), ncol=times, dimnames = list(pid_geneid$pid, 1:times))		
	for (i_time in 1:times) {
		cat("Permutation", i_time, "\n")
		if(!is.null(snpWeightDF)) {
			i_genes=sampleGenebySNP(snpWeightDF, n_whites)
		} else {
			i_genes=sample(geneids, n_whites)
		}		
		i_drawnWhites=sapply(pid_geneid_list, function(x) sum(x %in%  i_genes))
		rst[, i_time]=phyper(i_drawnWhites-1, n_whites, 
				n_balls-n_whites, i_drawn, lower.tail = FALSE)
		if(hitSize>0) rst[i_drawnWhites<hitSize,i_time]=NA
	}	
	return(rst)
}
