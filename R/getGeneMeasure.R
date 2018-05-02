getGeneMeasure<-function(snpGeneP) {
	
	if(missing(snpGeneP)) stop("Need to specify the snpGeneP value")
	
	if(class(snpGeneP) !="data.frame") 
		stop("snpGeneP has to be a data frame")
	
	if (!all(c("snp", "gene_id", "p") %in% colnames(snpGeneP))) 
		stop("the snpGeneP has to contain column names of 'snp', 'gene_id' and 'p'")
	
	if(class(snpGeneP$snp)!="character") 
		stop("the snpGeneP$snp has to be a character vector")
	
	if(!any(class(snpGeneP$gene_id) %in% c("numeric","integer"))) 
		stop("the snpGeneP$gene_id has to be a integer or numeric vector")
	
	if(class(snpGeneP$p)!="numeric") 
		stop("the snpGeneP$p has to be a numeric vector")
	
	snpGeneP=snpGeneP[!is.na(snpGeneP$p),]
	
	
	snpGeneP<-snpGeneP[order(snpGeneP$gene_id),]
	gids<-as.character(unique(snpGeneP$gene_id))
	
	
	minp<-by(snpGeneP[,"p"], snpGeneP[,"gene_id"], function(x) min(x))	
	sndp<-by(snpGeneP[,"p"], snpGeneP[,"gene_id"], 
			function(x) {
				if(length(x)==1) x else	sort(x)[2]
			} )
	simp<-by(snpGeneP[,"p"], snpGeneP[,"gene_id"], 
			function(x) min(sort(x)*length(x)/(1:length(x))))
	fishp<-by(snpGeneP[,"p"], snpGeneP[,"gene_id"], 
			function(x) {
				df=2*length(x) 
				ch2=(-2*sum(log(x))) 
				pchisq(ch2, df=df, lower.tail=FALSE)
			})
	
	geneMeasure<-data.frame(gene_id=as.integer(gids), 
			minp=as.numeric(minp[gids]), sndp=as.numeric(sndp[gids]), 
			simp=as.numeric(simp[gids]), fishp=as.numeric(fishp[gids]))
	
#	if(us) {
#		for (imeasure in c("minp","sndp","simp","fishp")) {
#			geneMeasure[,imeasure]<-uscore(geneMeasure[,imeasure])
#		}
#	}	
	
	return(geneMeasure)	
}
