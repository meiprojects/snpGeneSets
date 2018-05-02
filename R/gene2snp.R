gene2snp<-function(gene, isGeneID=TRUE, up=2000, down=2000, GRCh=37) {
	regions=getGeneMap(gene=gene, isGeneID=isGeneID, GRCh=GRCh)
	if(nrow(regions$gene_map)==0) return(NULL)	
	regions$start=regions$start-up
	regions$end=regions$end+down	
	snps=getRegionSNP(regions$gene_map, GRCh=GRCh)
	return(snps)
}
