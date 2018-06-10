# Author: Robert J. Hijmans, r.hijmans@gmail.com
# License GPL3
# Version 0.1  July 2016


lintul3 <- function(crop, soil, control, weather) {
	m <- Lintul3Model$new()
	if (!missing(crop)) { crop(m) <- crop }
	if (!missing(soil)) { soil(m) <- soil }
	if (!missing(control)) { control(m) <- control }
	if (!missing(weather)) { weather(m) <- weather }
	return(m)
}


setMethod("run", signature('Rcpp_Lintul3Model'), 
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


setMethod("crop<-", signature('Rcpp_Lintul3Model', 'list'), 
	function(x, value) {
		parameters <- c("IDSL", "TSUM1", "TSUM2", "DVSI", "DVSEND", "TDWI", "RGRLAI", "SPA", "DTSMTB", "SLATB", "SSATB", "TBASE", "KDIFTB", "RUETB", "TMPFTB", "KDIFTB", "RUETB", "TMPFTB", "TMNFTB", "COTB", "FRTB", "FLTB", "FSTB", "FOTB", "RDRL", "RDRLTB", "RDRSHM", "RDRNS", "RDRRTB", "RDRSTB", "CFET", "DEPNR", "IAIRDU", "RDI", "RRI", "RDMCR", "DVSDR", "DVSDLT", "DVSNLT", "DVSNT", "TBASEM", "TEFFMX", "TSUMEM", "FNTRT", "FRNX", "FRPX", "FRKX", "LAICR", "LRNR", "LSNR", "LRPR", "LSPR", "LRKR", "LSKR", "NLAI", "NLUE", "NMAXSO", "PMAXSO", "KMAXSO", "NPART", "NFIXF", "NSLA", "RNFLV", "RNFRT", "RNFST", "TCNT", "NMXLV", "RPFLV", "RPFRT", "RPFST", "TCPT", "PMXLV", "RKFLV", "RKFRT", "RKFST", "TCKT", "KMXLV", "PHOTTB", "RDI")
		if (is.null(value$IARDU)) {
			value$IAIRDU = 0
		}
		nms <- names(value)
		if (!all(parameters %in% nms)) stop(paste("parameters missing:", paste(parameters[!(parameters %in% nms)], collapse=", ")))
		value <- value[parameters]
		nms <- names(value)
		lapply(1:length(value), function(i) eval(parse(text = paste0("x$crop$", nms[i], " <- ", value[i]))))
		return(x)
	}
)

setMethod("soil<-", signature('Rcpp_Lintul3Model', 'list'), 
	function(x, value) {
		parameters <- c("SMDRY", "SMW", "SMFC", "SM0", "SMI", "SMLOWI", "RDMSO", "RUNFR", "CFEV", "KSUB", "CRAIRC", "PMINS", "NMINS", "KMINS", "RTPMINS", "RTNMINS", "RTKMINS", "PRFTAB", "NRFTAB", "KRFTAB")
		nms <- names(value)
		if (!all(parameters %in% nms)) stop(paste("parameters missing:", paste(parameters[!(parameters %in% nms)], collapse=", ")))
		value <- value[parameters]
		nms <- names(value)
		lapply(1:length(value), function(i) eval(parse(text = paste0("x$soil$", nms[i], " <- ", value[i]))))
		return(x)
	}
)


setMethod("weather<-", signature('Rcpp_Lintul3Model', 'list'), 
	function(x, value) {
		parameters <- c("date", "srad", "tmin", "tmax", "prec", "wind", "vapr")
		nms <- names(value)
		if (!all(parameters %in% nms)) stop(paste("parameters missing:", paste(parameters[!(parameters %in% nms)], collapse=", ")))
		x$setWeather(value$date, value$tmin, value$tmax, value$srad, value$prec, value$wind, value$vapr)
		return(x)
	}
)


setMethod("control<-", signature('Rcpp_Lintul3Model', 'list'), 
	function(x, value) {
		parameters <- c("start", "emergence","maxdur", "IOPT", "IDPL", "DAYPL", "PL", "IRRI", "IRRTAB", "FERNTAB", "FERPTAB", "FERKTAB")		
		nms <- names(value)
		if (!all(parameters %in% nms)) stop(paste("parameters missing:", paste(parameters[!(parameters %in% nms)], collapse=", ")))
		value <- value[parameters]
		nms <- names(value)
		lapply(1:length(value), function(i) eval(parse(text = paste0("x$control$", nms[i], " <- ", value[i]))))
		return(x)
	}
)


lintul3_control <- function() {
	f <- system.file("lintul/control.ini", package="Rlintul")
	ini <- .readIniFile(f)
	.getNumLst(ini)
}


lintul3_soil <- function(name='p1') {
	f <- list.files(system.file("lintul/soil", package="Rlintul"), full.names=TRUE)
	soils <- gsub('LINTUL_soil_', '', basename(f))
	soils <- gsub('.ini', '', soils)
	if (name %in% soils) {
		i <- which (name == soils)
		lst <- .getNumLst(.readIniFile(f[i]))

		f2 <- system.file("lintul/soil.ini", package="Rlintul")
		lst2 <- .getNumLst(.readIniFile(f2))

		return(c(lst, lst2))
	} else {
		stop(paste('not available. Choose one of:', paste(soils, collapse=', ')))
	}
}



lintul3_crop <- function(name) {
    f <- list.files(system.file("lintul/crop", package="Rlintul"), full.names=TRUE)
    crops <- gsub('lintul3_', '', basename(f))
    crops <- gsub('.ini', '', crops)
    if (name %in% crops) {
		i <- which (name == crops)
		ini <- .readIniFile(f[i])
		lst <- .getNumLst(ini)
		return(lst)
    } else {
		stop(paste('not available. Choose one of:', paste(crops, collapse=', ')))
    }
}



readLIN3output <- function(f) {
	X <- readLines(f)
	H <- .trim2(X[17])
	H <- unlist(strsplit(H, ' '))
	X <- X[-c(1:17)]
	X = .trim2(X)
	X <- strsplit(X, ' ')
	s = lapply(X, function(i) as.numeric(i))
	ss = do.call(rbind, s)

	colnames(ss) = H
	ss
}



setMethod ('show' , 'Rcpp_Lintul3Model', function(object) { utils::str(object) })	
setMethod ('show' , 'Rcpp_Lintul3Output', function(object) { utils::str(object) })	
setMethod ('show' , 'Rcpp_Lintul3Crop', function(object) { utils::str(object) })	
setMethod ('show' , 'Rcpp_Lintul3Control', function(object) { utils::str(object) })	
setMethod ('show' , 'Rcpp_Lintul3Weather', function(object) { utils::str(object) })	
