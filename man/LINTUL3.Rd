\name{lintul3}

\alias{lintul3}
\alias{lintul3_crop}
\alias{lintul3_soil}
\alias{lintul3_control}

\alias{crop<-,Rcpp_Lintul3Model-method}
\alias{soil<-,Rcpp_Lintul3Model-method}
\alias{weather<-,Rcpp_Lintul3Model-method}
\alias{control<-,Rcpp_Lintul3Model-method}
\alias{run,Rcpp_Lintul3Model-method}

\alias{readLIN3output}


\title{
LINTUL3 crop growth model
}

\description{
LINTUL3 simulates the growth of a crop as function of intercepted radiation, temperature and light use efficiency. Soil water (free drainage) and simple nitrogen, phosphorus and potassium balances are simulated and also the effects of water and N, P and K supplies on crop growth. 

LINTUL3 was developed by Joost Wolf at Wageningen University.
}

\usage{
lintul3(crop, soil, control, weather)
lintul3_crop(name)
lintul3_soil(name='p1')
lintul3_control()
\S4method{crop}{Rcpp_Lintul3Model}(x)<-value
\S4method{soil}{Rcpp_Lintul3Model}(x)<-value
\S4method{control}{Rcpp_Lintul3Model}(x)<-value
\S4method{weather}{Rcpp_Lintul3Model}(x)<-value
\S4method{run}{Rcpp_Lintul3Model}(x, ...)
}


\arguments{
 \item{x}{Rcpp_Lintul3Model object}
 \item{crop}{list with named crop parameters. See Details. An example is returned by \code{lintul3_crop()}} 
 \item{soil}{list with named soil parameters. See Details. An example is returned by \code{lintul3_soil()}} 
 \item{control}{list with named control parameters} 
 \item{weather}{data.frame with daily weather data} 
 \item{value}{a list with crop, soil, or control parameters, or a data.frame with weather data. As above}
 \item{...}{additional arguments. None implemented}
 \item{name}{character}
 
}

\value{
matrix
}


\examples{

library(Rlintul)
wth <- example_weather()
crop <- lintul3_crop('POTATO')
soil <- lintul3_soil(name='p1')

contr <- lintul3_control()
contr$emergence <- as.Date('2012-04-01')
contr$start <- contr$emergence - 10
contr$PL <- FALSE
contr$DAYPL <- contr$emergence
contr$IOPT <- 2
# IOPT 1=potential, 2=water limited, 3=water and N limited and 4=water and N, P and K limited

x <- lintul3(crop, soil, contr, wth)


\dontrun{

z <- run(x)

#or 
x$run()
x$out


# test LINTUL3

wf <- system.file("extdata/Netherlands_Wageningen.csv", package="meteor")
wth <- read.csv(wf)
wth$date <- as.Date(wth$date)
wth$srad <- wth$srad / 1000
crop <- lintul3_crop('WWHEAT')
soil <- lintul3_soil('p1')
sdate <- dateFromDoy(1, 1969)
contr <- lintul3_control()
contr$emergence=sdate
contr$start = sdate
contr$long_output=TRUE
contr$PL = FALSE


contr$IOPT = 1
contr.IRRI = 1
f1 <- system.file("lintul/test/3/WAGI1PTS.OUT", package="Rlintul")
s1 <- readLIN3output(f1)
x <- lintul3(crop, wth, soil, contr)
plot(s1[,'DAY'], s1[,'LAI'], type='l')
points(x[, 'step']+1, x[, 'LAI'])

contr$IOPT = 2
contr.IRRI = 1
f2 <- system.file("lintul/test/3/WAGW2PTS.OUT", package="Rlintul")
s2 <- readLIN3output(f2)
x2 <- lintul3(crop, wth, soil, contr)
plot(s2[,'DAY'], s2[,'LAI'], type='l')
points(x2[, 'step']+1, x2[, 'LAI'])

contr$IOPT = 3
f3 <- system.file("lintul/test/3/WAGW3PTS.OUT", package="Rlintul")
s3 <- readLIN3output(f3)
x3 <- lintul3(crop, wth, soil, contr)
plot(s3[,'DAY'], s3[,'LAI'], type='l')
points(x3[, 'step']+1, x3[, 'LAI'])

contr$IOPT = 4
f4 <- system.file("lintul/test/3/WAGW4PTS.OUT", package="Rlintul")
s4 <- readLIN3output(f4)
x4 <- lintul3(crop, wth, soil, contr)
plot(s4[,'DAY'], s4[,'LAI'], type='l')
points(x4[, 'step']+1, x4[, 'LAI'])
}
}


