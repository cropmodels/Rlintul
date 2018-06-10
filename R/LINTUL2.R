# Author: Robert J. Hijmans, r.hijmans@gmail.com
# License GPL3
# Version 0.1  July 2016


lintul2 <- function(crop, soil, control, weather) {
	m <- Lintul2Model$new()
	if (!missing(crop)) { crop(m) <- crop }
	if (!missing(soil)) { soil(m) <- soil }
	if (!missing(control)) { control(m) <- control }
	if (!missing(weather)) { weather(m) <- weather }
	return(m)
}


setMethod("run", signature('Rcpp_Lintul2Model'), 
	function(x, ...) {
		x$run()
		out <- x$out
		date <- as.Date(x$control$emergence, origin="1970-01-01") + out$step
		Wtot <- out$WRT + out$WLVD + out$WLVG + out$WST + out$WSO
		v <- data.frame(date, out$TSUM, out$LAI, out$WRT, out$WLV, out$WLVD, out$WLVG, out$WST, out$WSO, Wtot,
							  out$EVAP, out$TRAN, out$TRANRF, out$WA, out$WC, out$RWA)
		colnames(v) <- c("date", "TSUM", "LAI", "WRT", "WLV", "WLVD", "WLVG", "WST", "WSO", "Wtot",
							"EVAP", "TRAN", "TRANRF", "WA", "WC", "RWA")
		v
	}
)


setMethod("crop<-", signature('Rcpp_Lintul2Model', 'list'), 
	function(x, value) {
		parameters <- c("LAIi", "SLA", "Tbase", "RGRL", "Tsum1", "Tsum2", "LAIcr", "RDRSHM", "RUE", "K", "RDRT", "FRTTB", "FLVTB", "FSTTB", "FSOTB", "ROOTDi", "ROOTDM", "RRDMAX", "TRANCO")
		nms <- names(value)
		if (!all(parameters %in% nms)) stop(paste("parameters missing:", paste(parameters[!(parameters %in% nms)], collapse=", ")))
		value <- value[parameters]
		nms <- names(value)
		lapply(1:length(value), function(i) eval(parse(text = paste0("x$crop$", nms[i], " <- ", value[i]))))
		return(x)
	}
)

setMethod("soil<-", signature('Rcpp_Lintul2Model', 'list'), 
	function(x, value) {
		parameters <- c("WCi", "WCAD", "WCWP", "WCFC", "WCWET", "WCST", "DRATE", "IRRIGF")
		nms <- names(value)
		if (!all(parameters %in% nms)) stop(paste("parameters missing:", paste(parameters[!(parameters %in% nms)], collapse=", ")))
		value <- value[parameters]
		nms <- names(value)
		lapply(1:length(value), function(i) eval(parse(text = paste0("x$soil$", nms[i], " <- ", value[i]))))
		return(x)
	}
)


setMethod("weather<-", signature('Rcpp_Lintul2Model', 'list'), 
	function(x, value) {
		parameters <- c("date", "srad", "tmin", "tmax", "prec", "wind", "vapr")
		nms <- names(value)
		if (!all(parameters %in% nms)) stop(paste("parameters missing:", paste(parameters[!(parameters %in% nms)], collapse=", ")))
		x$setWeather(value$date, value$tmin, value$tmax, value$srad, value$prec, value$wind, value$vapr)
		return(x)
	}
)


setMethod("control<-", signature('Rcpp_Lintul2Model', 'list'), 
	function(x, value) {
		#parameters <- c("emergence", "maxdur", "long_output")
		parameters <- c("start", "emergence", "maxdur")
		nms <- names(value)
		if (!all(parameters %in% nms)) stop(paste("parameters missing:", paste(parameters[!(parameters %in% nms)], collapse=", ")))
		value <- value[parameters]
		nms <- names(value)
		lapply(1:length(value), function(i) eval(parse(text = paste0("x$control$", nms[i], " <- ", value[i]))))
		return(x)
	}
)



lintul2_crop <- function() {
	c(lintul1_crop(), list(ROOTDi = 0.1, ROOTDM = 1.2, RRDMAX = 0.012, TRANCO = 8.))
}

lintul2_soil <- function() {
	list(WCi=0.36, WCAD=0.08, WCWP=0.23, WCFC=0.36, WCWET=0.48, WCST=0.55, DRATE=50, IRRIGF=0)
}

setMethod ('show' , 'Rcpp_Lintul2Model', function(object) { utils::str(object) })	
setMethod ('show' , 'Rcpp_Lintul2Output', function(object) { utils::str(object) })	
setMethod ('show' , 'Rcpp_Lintul2Crop', function(object) { utils::str(object) })	
setMethod ('show' , 'Rcpp_Lintul2Control', function(object) { utils::str(object) })	
setMethod ('show' , 'Rcpp_Lintul2Weather', function(object) { utils::str(object) })	

