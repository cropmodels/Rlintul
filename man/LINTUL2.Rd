\name{lintul2}

\alias{lintul2}
\alias{lintul2_crop}
\alias{lintul2_soil}
\alias{crop<-,Rcpp_Lintul2Model-method}
\alias{soil<-,Rcpp_Lintul2Model-method}
\alias{weather<-,Rcpp_Lintul2Model-method}
\alias{control<-,Rcpp_Lintul2Model-method}
\alias{run,Rcpp_Lintul2Model-method}


\title{
LINTUL2 crop growth model
}

\description{
LINTUL is a relatively simple crop growth simulation model
}

\references{
Spitters, C.J.T. and A.H.C.M. Schapendonk, 1990. Evaluation of breeding strategies for drought tolerance in potato by means of crop growth simulation. Plant and Soil 123: 193-203.
}



\usage{
lintul2(crop, soil, control, weather)
lintul2_crop()
lintul2_soil()
\S4method{crop}{Rcpp_Lintul2Model}(x)<-value
\S4method{soil}{Rcpp_Lintul2Model}(x)<-value
\S4method{control}{Rcpp_Lintul2Model}(x)<-value
\S4method{weather}{Rcpp_Lintul2Model}(x)<-value
\S4method{run}{Rcpp_Lintul2Model}(x, ...)
}

\arguments{
 \item{x}{Rcpp_Lintul2Model object}
 \item{crop}{list with named crop parameters. See Details. An example is returned by \code{lintul2_crop()}} 
 \item{soil}{list with named soil parameters. See Details. An example is returned by \code{lintul2_soil()}} 
 \item{control}{list with named control parameters} 
 \item{weather}{data.frame with daily weather data} 
 \item{value}{a list with crop, soil, or control parameters, or a data.frame with weather data. As above}
 \item{...}{additional arguments. None implemented}
}

\value{
matrix
}

\examples{

library(Rlintul)
m <- lintul2()
emerg <- as.Date('2012-04-01')
contr <- list(emergence=emerg, start=emerg-1, maxdur=200)
w <- example_weather()
acrop <- lintul2_crop()
asoil <- lintul2_soil()

crop(m) <- acrop
soil(m) <- asoil
control(m) <- contr
weather(m) <- w

m$run()

m$out

# or do 

m2 <- lintul2(acrop, asoil, contr, w)
r <- run(m2)
tail(r)


# Compare with FORTRAN version results)

f1 <- system.file("lintul/test/2/res.dat", package="Rlintul")
s <- readLIN1output(f1)
f2 <- system.file("extdata/Netherlands_Wageningen.csv", package="meteor")
wth <- read.csv(f2)
wth$date <- as.Date(wth$date)
wth$srad <- wth$srad / 1000

m <- lintul2()
crop(m) <- lintul2_crop()
soil(m) <- lintul2_soil()
sdate <- dateFromDoy(58, 1971)
control(m) <- list(emergence=sdate+2, start = sdate, maxdur=365)
weather(m) <- wth
m$run()

x <- m$out

plot(s[,'TIME'], s[,'LAI'], type='l')
points(m$out$step+58, m$out$LAI)

plot(s[,'TIME'], s[,'WSO'], type='l')
points(m$out$step+58, m$out$WSO)

}

