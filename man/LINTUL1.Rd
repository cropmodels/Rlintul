\name{lintul1}

\alias{readLIN1output}

\alias{lintul1}
\alias{lintul1_crop}
\alias{crop<-,Rcpp_Lintul1Model-method}
\alias{weather<-,Rcpp_Lintul1Model-method}
\alias{control<-,Rcpp_Lintul1Model-method}
\alias{run,Rcpp_Lintul1Model-method}


\title{
LINTUL1 crop growth model
}

\description{
LINTUL1 is a relatively simple crop growth simulation model. It computes potential production.
}

\usage{
lintul1(crop, control, weather)
lintul1_crop()
\S4method{crop}{Rcpp_Lintul1Model}(x)<-value
\S4method{control}{Rcpp_Lintul1Model}(x)<-value
\S4method{weather}{Rcpp_Lintul1Model}(x)<-value
\S4method{run}{Rcpp_Lintul1Model}(x, ...)
}


\arguments{
 \item{x}{Rcpp_Lintul1Model object}
 \item{crop}{list with named crop parameters. See Details. An example is returned by \code{lintul1_crop()}} 
 \item{control}{list with named control parameters} 
 \item{weather}{data.frame with daily weather data} 
 \item{value}{a list with crop or control parameters, or a data.frame with weather data. As above}
 \item{...}{additional arguments. None implemented}
}

\value{
matrix
}

\references{
Spitters, C.J.T. and A.H.C.M. Schapendonk, 1990. Evaluation of breeding strategies for drought tolerance in potato by means of crop growth simulation. Plant and Soil 123: 193-203.
}


\examples{

library(Rlintul)
w <- example_weather()
crop <- lintul1_crop()
contr <- list(emergence=as.Date('2012-04-01'), maxdur=200)

m <- lintul1(crop, contr, w)
# or
m <- lintul1()
crop(m) <- crop
control(m) <- contr
weather(m) <- w

x <- run(m)
head(x)

# or 
m$run()
m


# change parameters
m$crop$RUE <- 3.5
m$control$emergence=as.Date('2012-04-01')
m$run()

out <- m$out
out


date <- as.Date(m$control$emergence, origin="1970-01-01") + out$step
par(mfrow=c(1,2))
plot(date, out$LAI)

prod <- cbind(out$WRT, out$WLV, out$WST, out$WSO)

wtot <- rowSums(prod)
plot(date, wtot, type='l', lwd=2)
lines(date, prod[,2], col='green', lwd=2)
lines(date, prod[,4], col='red', lwd=2)


# test LINTUL1

library(meteor)
library(Rlintul)
f1 <- system.file("lintul/test/1/res.dat", package="Rlintul")
s <- readLIN1output(f1)
f2 <- system.file("extdata/Netherlands_Wageningen.csv", package="meteor")
wth <- read.csv(f2)
wth$date <- as.Date(wth$date)
wth$srad <- wth$srad / 1000

crop <- lintul1_crop()
sdate <- dateFromDoy(58, 1971)
contr <- list(emergence=sdate+2, start = sdate, maxdur=365)
x <- run(lintul1(crop, contr, wth))
x$DOY <- doyFromDate(x$date)

plot(s[,'TIME'], s[,'LAI'], type='l')
points(x[, 'DOY'], x[, 'LAI'])

}

