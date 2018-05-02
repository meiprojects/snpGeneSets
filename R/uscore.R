uscore<-function(gm) {	
	
	if(!any(class(gm) %in% c("numeric","integer"))) 
		stop("the gm has to be a integer or numeric vector")
	
	pr<-rep(NA, length(gm))
	eid<-!is.na(gm)
	gm<-gm[eid]
	
	pr[eid]<-sapply(gm, function(x) sum(gm<x)+0.5*sum(gm==x))/length(gm)
	return(pr)
}
