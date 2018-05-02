msigdb_build <-function (gmt_dir) {
	
	if(missing(gmt_dir)) stop("gmt_dir cannot be missing")
	
	setdb=paste(system.file("extdata", package="snpGeneSets"),"/msigdb.db",sep="")
	
	if(file.exists(setdb)) {
		file.remove(setdb)
	}
	
	owd<-getwd()
	setwd(gmt_dir)	
	
	path.symbol.file="msigdb.v6.1.symbols.gmt"
	path.entrez.file="msigdb.v6.1.entrez.gmt"
	
	path.symbol=readLines(path.symbol.file)
	path.entrez=readLines(path.entrez.file)	
	
	###extract gene id and gene symbol used in msigdb pathway
	genes=NULL
	path.counts=length(path.entrez)
	for (p in 1:path.counts) {
		cat("total=", path.counts, ":\tpath=", p,"\n")
		p.s=path.symbol[p]; p.e=path.entrez[p]
		p.s.vector=strsplit(p.s, split="\t")[[1]]
		p.e.vector=strsplit(p.e, split="\t")[[1]]
		if(length(p.s.vector)!=length(p.e.vector)) {
			stop("path=", p, " is wrong\n")
			next
		}
		genes=unique(c(genes, paste(p.e.vector[-(1:2)], p.s.vector[-(1:2)], sep=":")))
		gc()
	}
	ng=length(genes)
	gene_id=integer(ng); gene_sym=character(ng)
	for (ng_i in 1:ng) {
		cat("total",ng,":",ng_i,"\n")
		ng_i.tmp=strsplit(genes[ng_i],":")[[1]]
		gene_id[ng_i]=as.integer(ng_i.tmp[1])
		gene_sym[ng_i]=ng_i.tmp[2]
	} 
	msigdb_id_sym=data.frame(gene_id,gene_sym, stringsAsFactors = FALSE)
	
	path.table=NULL
	path.counts=length(path.entrez)
	
	for (p in 1:path.counts) {
		cat("total=", path.counts, ":\tpath=", p,"\n")
		p.e=path.entrez[p]
		p.e.vector=strsplit(p.e, split="\t")[[1]]
		path.table=rbind(path.table, c(p, p.e.vector[1:2], length(p.e.vector[-(1:2)]), paste(":", paste(p.e.vector[-(1:2)],collapse=":"),":",sep="")))
	}
	colnames(path.table)=c("pid", "pathway", "link", "size", "geneid")
	path.table=data.frame(path.table, stringsAsFactors = FALSE)
	path.table$pid=as.integer(path.table$pid)
	path.table$size=as.integer(path.table$size)	
	
	#build category
	setTypes<-c("c1","c2","c2_cgp","c2_cp","c2_biocarta","c2_kegg","c2_reactome","c3","c3_mir",
			"c3_tft","c4","c4_cgn","c4_cm","c5","c5_bp","c5_cc","c5_mf","c6","c7", "hallmark")
	
	gmtfiles<-c("c1.all.v6.1.entrez.gmt", "c2.all.v6.1.entrez.gmt","c2.cgp.v6.1.entrez.gmt",
			"c2.cp.v6.1.entrez.gmt", "c2.cp.biocarta.v6.1.entrez.gmt","c2.cp.kegg.v6.1.entrez.gmt",
			"c2.cp.reactome.v6.1.entrez.gmt","c3.all.v6.1.entrez.gmt","c3.mir.v6.1.entrez.gmt",
			"c3.tft.v6.1.entrez.gmt","c4.all.v6.1.entrez.gmt","c4.cgn.v6.1.entrez.gmt",
			"c4.cm.v6.1.entrez.gmt","c5.all.v6.1.entrez.gmt","c5.bp.v6.1.entrez.gmt", 
			"c5.cc.v6.1.entrez.gmt","c5.mf.v6.1.entrez.gmt","c6.all.v6.1.entrez.gmt",
			"c7.all.v6.1.entrez.gmt", "h.all.v6.1.entrez.gmt")
	
	gmtlogic<-file.exists(gmtfiles)
	if(!all(gmtlogic)) {
		cat("!!!These files are not exising: ", paste(gmtfiles[gmtlogic],collapse=","),"!!!\n")
		return(NULL)
	}
	
	for (gmti in 1:length(gmtfiles)) {
		cat("parsing", gmtfiles[gmti],"....\n")
		gmtdata=readLines(gmtfiles[gmti])
		gmtdata=gsub("^(.+?)\\t.+$","\\1", gmtdata)
		path.table[,setTypes[gmti]] <- path.table$pathway %in% gmtdata
		if(sum(path.table[,setTypes[gmti]])!=length(gmtdata)) {
			stop(gmtfiles[gmti], " have at least one pathway not identified!\n")
		}
		set.geneid=paste(path.table[(path.table$pathway %in% gmtdata), "geneid"],collapse="")
		set.geneid=sapply(msigdb_id_sym$gene_id, 
				function(x) grepl(paste(":",x,":",sep=""), set.geneid, fixed = TRUE,useBytes = TRUE))
		
		msigdb_id_sym[,setTypes[gmti]]=set.geneid
		gc()
	}	
	
	##start: check if there is duplicated gene_id (msigdb v6.1)
	if(sum(duplicated(msigdb_id_sym$gene_id))>0) {
		dupgids=msigdb_id_sym[(duplicated(msigdb_id_sym$gene_id)),"gene_id"]
		dupgids_syms=getGeneMap(dupgids)$gene_map[,c("gene_id", "gene_symbol")]
		for(i in 1:nrow(dupgids_syms)) {
			msigdb_id_sym[msigdb_id_sym$gene_id==dupgids_syms[i,"gene_id"],"gene_sym"]=dupgids_syms[i,"gene_symbol"]
		}		
		msigdb_id_sym=msigdb_id_sym[!duplicated(msigdb_id_sym$gene_id),]
		
	}	
	##end
	
	
	require("RSQLite")
	drv=dbDriver("SQLite")
	con=dbConnect(drv, setdb)
	dbWriteTable(con, "geneSets", path.table, row.names=FALSE)
	dbWriteTable(con, "geneInfo", msigdb_id_sym, row.names=FALSE)
	dbExecute(con, 'CREATE UNIQUE INDEX idx_pid on geneSets (pid)')
	dbExecute(con, 'CREATE UNIQUE INDEX idx_geneid on geneInfo (gene_id)')	
	dbDisconnect(con)		
	setwd(owd)	
	cat(setdb, "has been successfully installed!\n")
	return(0)
}
