# Author: Robert J. Hijmans, r.hijmans@gmail.com
# License GPL3
# Version 0.1  July 2016



lintul1 <- function(crop, control, weather) {
	m <- Lintul1Model$new()
	if (!missing(crop)) { crop(m) <- crop }
	if (!missing(control)) { control(m) <- control }
	if (!missing(weather)) { weather(m) <- weather }
	return(m)
}


setMethod("run", signature('Rcpp_Lintul1Model'), 
	function(x, ...) {
		x$run()
		#ff <- names(getRefClass("Rcpp_Lintul1Output")$fields())
		out <- x$out
		date <- as.Date(x$control$emergence, origin="1970-01-01") + out$step
		v <- data.frame(date, out$DLV, out$LAI, out$TSUM, out$WRT, out$WLV, out$WLVD, out$WLVG, out$WST, out$WSO)	
		colnames(v) <- c("date", "DLV", "LAI", "TSUM", "WRT", "WLV", "WLVD", "WLVG", "WST", "WSO")
		v$Wtot <- v$WRT + v$WLVD + v$WLVG + v$WST + v$WSO
		v
	}
)


setMethod("crop<-", signature('Rcpp_Lintul1Model', 'list'), 
	function(x, value) {
		parameters <- c("LAIi", "SLA", "Tbase", "RGRL", "Tsum1", "Tsum2", "LAIcr", "RDRSHM", "RUE", "K", "RDRT", "FRTTB", "FLVTB", "FSTTB", "FSOTB")
		nms <- names(value)
		if (!all(parameters %in% nms)) stop(paste("parameters missing:", paste(parameters[!(parameters %in% nms)], collapse=", ")))
		value <- value[parameters]
		nms <- names(value)
		lapply(1:length(value), function(i) eval(parse(text = paste0("x$crop$", nms[i], " <- ", value[i]))))
		return(x)
	}
)


setMethod("weather<-", signature('Rcpp_Lintul1Model', 'list'), 
	function(x, value) {
		parameters <- c("date", "srad", "tmin", "tmax")
		nms <- names(value)
		if (!all(parameters %in% nms)) stop(paste("parameters missing:", paste(parameters[!(parameters %in% nms)], collapse=", ")))
		x$setWeather(value$date, value$tmin, value$tmax, value$srad, value$prec, value$wind, value$vapr)
		return(x)
	}
)


setMethod("weather<-", signature('Rcpp_Lintul1Model', 'list'), 
	function(x, value) {
		parameters <- c("date", "srad", "tmin", "tmax")
		nms <- names(value)
		if (!all(parameters %in% nms)) stop(paste("parameters missing:", paste(parameters[!(parameters %in% nms)], collapse=", ")))
		x$setWeather(value$date, value$tmin, value$tmax, value$srad)
		return(x)
	}
)


setMethod("control<-", signature('Rcpp_Lintul1Model', 'list'), 
	function(x, value) {
		#parameters <- c("emergence", "maxdur", "long_output")
		parameters <- c("emergence", "maxdur")
		nms <- names(value)
		if (!all(parameters %in% nms)) stop(paste("parameters missing:", paste(parameters[!(parameters %in% nms)], collapse=", ")))
		value <- value[parameters]
		nms <- names(value)
		lapply(1:length(value), function(i) eval(parse(text = paste0("x$control$", nms[i], " <- ", value[i]))))
		return(x)
	}
)


lintul1_crop <- function() {
	params <- list(LAIi=0.012, SLA=0.022, Tbase=0, RGRL=0.009, Tsum1=1110, Tsum2=2080, LAIcr=4, RDRSHM=0.03, RUE=3.0, K=0.6)
	RDRT <- matrix(c(-10.,0.03, 10.,0.03, 15.,0.04, 30.,0.09, 50.,0.09), ncol=2, byrow=TRUE)
	FRTTB <- matrix(c(0.,0.50, 110, 0.50, 275 ,0.34, 555, 0.12, 780, 0.07, 1055, 0.03, 1160, 0.02, 1305, 0, 2500, 0), ncol=2, byrow=TRUE)
	FLVTB <- matrix(c(0, 0.33, 110, 0.33, 275, 0.46, 555, 0.44, 780, 0.14, 1055, 0, 2500, 0), ncol=2, byrow=TRUE)
	FSTTB <- matrix(c(0, 0.17, 110, 0.17, 275, 0.20, 555, 0.44, 780, 0.79, 1055.,0.97, 1160, 0, 2500, 0), ncol=2, byrow=TRUE)
	FSOTB <- matrix(c(0, 0, 1055, 0, 1160.,0.98, 1305, 1, 2500, 1), ncol=2, byrow=TRUE)
	crop <- c(params, list(RDRT=RDRT, FRTTB=FRTTB, FLVTB=FLVTB, FSTTB=FSTTB, FSOTB=FSOTB))
	return(crop)
}


setMethod ('show' , 'Rcpp_Lintul1Model', function(object) { utils::str(object) } )	
setMethod ('show' , 'Rcpp_Lintul1Output', function(object) { utils::str(object) } )	
setMethod ('show' , 'Rcpp_Lintul1Crop', function(object) { utils::str(object) } )	
setMethod ('show' , 'Rcpp_Lintul1Control', function(object) { utils::str(object) })	
setMethod ('show' , 'Rcpp_Lintul1Weather', function(object) { utils::str(object) })	


readLIN1output <- function(f) {
	X <- readLines(f)
	H <- unlist(strsplit(X[7], '\t'))
	X <- X[-c(1:8)]
	X = .trim2(X)
	X <- strsplit(X, '\t')
	s = lapply(X, function(i) as.numeric(i))
	ss = do.call(rbind, s)
	colnames(ss) = H
	ss
}

