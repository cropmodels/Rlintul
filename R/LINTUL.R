# Author: Robert J. Hijmans, r.hijmans@gmail.com
# License GPL3
# Version 0.1  July 2016


#if (!isGeneric("crop<-")) { setGeneric("crop<-", function(x, value) standardGeneric("crop<-")) }	
#if (!isGeneric("soil<-")) { setGeneric("soil<-", function(x, value) standardGeneric("soil<-")) }	
#if (!isGeneric("control<-")) { setGeneric("control<-", function(x, value) standardGeneric("control<-")) }	
#if (!isGeneric("weather<-")) { setGeneric("weather<-", function(x, value) standardGeneric("weather<-")) }	
#if (!isGeneric("run")) { setGeneric("run", function(x, ...) standardGeneric("run")) }	

	

.trim2 <- function(x) return(gsub("^ *|(?<= ) | *$", "", x, perl=TRUE))


.getNumLst <- function(ini) {
	v <- ini[,3]
	vv <- sapply(v, function(i) strsplit(i, ','), USE.NAMES = FALSE)
	vv <- sapply(vv, as.numeric)
	lst <- lapply(vv, function(i) if(length(i) > 1) { matrix(i, ncol=2, byrow=TRUE) } else {i})
	names(lst) <- ini[,2]
	lst
}


