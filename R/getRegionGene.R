getRegionGene<-function(regionDF, GRCh=37) {
	if(missing(regionDF)) 
		stop("regionDF has to be assigned a single value")
	if(class(regionDF)!="data.frame") stop("regionDF has to be a data.frame")
	if(!all(c("chr", "start") %in% colnames(regionDF))) stop("regionDF has to contain columns of 'chr' and 'start'")
	
	isRegion=all(c("chr", "start", "end") %in% colnames(regionDF))
	if(!is.character(regionDF$chr)) stop("'chr' of regionDF has to be 'character'")
	if(!(is.numeric(regionDF$start)|is.integer(start))) stop("start has to be integer or numeric")
	if(isRegion) {
		if(!(is.numeric(regionDF$end)|is.integer(end))) stop("end has to be integer or numeric")
	}
	
	regionDF$chr=toupper(regionDF$chr)
	if(!all(unique(regionDF$chr) %in% getChrs())) warning("some 'chr' of 'regionDF' do not exist. Please run getChrs() to get the supported chromosome ID")
	
	if(GRCh==37) {
		genedb=system.file("extdata", "gene105.db", package="snpGeneSets")		
	} else if (GRCh==38) {
		genedb=system.file("extdata", "gene106.db", package="snpGeneSets")
	} else {
		stop("The reference genome build, GRCh, has to be 37 or 38")
	}
	if(!file.exists(genedb)) stop("The file of genedb does not exists")
	
	if(isRegion) {	
		geneMap<-by(regionDF, regionDF[,1], function(x) {				
					regionQuery<-paste("(", 
							paste(apply(x, 1, function(xi) {	paste("((start>=", xi["start"], "AND start<=", xi["end"],")",
														"OR (end>=", xi["start"], "AND end<=", xi["end"],")",
														"OR (start<=", xi["start"], "AND end>=", xi["end"],"))",
														sep=" ")
											}), collapse=" OR ")
							,")", sep="")	
					
					sqlQuery(paste("select * from genemaps INDEXED BY idx_pos where chr=", shQuote(unique(x[,1])),"AND", 
									regionQuery),genedb)				
				})
	} else { 
		geneMap<-by(regionDF, regionDF[,1], function(x) {
					regionQuery<-paste("(", 
							paste(apply(x, 1, function(xi) {	paste(
														"(start<=", xi["start"], "AND end>=", xi["start"],")",
														sep=" ")
											}), collapse=" OR ")
							,")", sep="")
					
					sqlQuery(paste("select * from genemaps INDEXED BY idx_pos where chr=", shQuote(unique(x[,1])),"AND", 
									regionQuery),genedb)
					
				})	
		
	}	
	geneMap<-as.data.frame(do.call("rbind",geneMap))
	rownames(geneMap)=NULL	
	return(geneMap)
}
