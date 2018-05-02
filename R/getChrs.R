getChrs <-
function() {return(c(1:22,"X","Y","M"))}

sqlQuery <-
function(sql, dbname){
	require(RSQLite)
	driver <- dbDriver("SQLite")
	connect <- dbConnect(driver, dbname=dbname);
	closeup <- function() dbDisconnect(connect)
	dd <- tryCatch(dbGetQuery(connect, sql), finally=closeup)
	dbDisconnect(connect)
	return(dd)
}

sampleGenebySNP<-function(snpGene, size) {	
	if(missing(snpGene)) stop("Need to specify the snpGene value")
	if(!(class(snpGene) %in% c("data.frame","matrix","array" ))) 
		stop("snpGene has to be a data frame or matrix or array")
	if(class(snpGene)!="data.frame") 
		snpGene<-as.data.frame(snpGene, stringsAsFactors = FALSE)	
	if (!all(c("gene_id", "snp") %in% colnames(snpGene))) 
		stop("the snpGene has to contain column names of 'snp' and 'gene_id'")	
	if(!(is.integer(snpGene$gene_id)|is.numeric(snpGene$gene_id))) {
		stop("The gene_id of snpGene has to be a integer or numeric Entrez gene ID")
	}	
	
	smpGenes=sample(snpGene$gene_id)		
	smp_sig_gene=unique(smpGenes[1:size])
	
	smp_sig_gene_N=length(smp_sig_gene)
	n0=size
	
	while(smp_sig_gene_N<size) {
		nd=size-smp_sig_gene_N
		n1=n0+nd
		n0=n0+1		
		smp_sig_gene=unique(c(smp_sig_gene, smpGenes[n0:n1]))
		smp_sig_gene_N=length(smp_sig_gene)
		n0=n1			
	}		
	return(smp_sig_gene)
}
