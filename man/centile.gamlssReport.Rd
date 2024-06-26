% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gamlssReport.R
\name{centile.gamlssReport}
\alias{centile.gamlssReport}
\title{Calculate the centiles}
\usage{
centile.gamlssReport(object, y, newdata, gamlss.prediction = FALSE)
}
\arguments{
\item{object}{the object of class \code{"gamlssReport"} generated by \code{\link{gamlssReport}}.}

\item{y}{the score for which to calculate the centile.}

\item{newdata}{a dataframe containing the Xs for which to evaluate the centile.}

\item{gamlss.prediction}{logical, if \code{TRUE} it uses the spline interpolation approach implemented in \code{predict.gamlss}, otherwise it uses the interpolation of the B-spline basis and the estimated penalized coefficients. Defaults to \code{FALSE}. Note: if set to \code{TRUE}, the \code{\link{gamlssReport}} had to be used with \code{extract.smooth=TRUE}, oterwise the error is returned.}
}
\value{
a vector containing the centiles which are values in (0,100)
}
\description{
Using the object generated by \code{\link{gamlssReport}}, the function calculates the centile for a given y and x
}
\details{
\code{y} can be a vector, in this case \code{newdata} need to be a dataframe with single row, or dataframe with the same no of rows as \code{y}; if \code{y} is a scalar \code{newdata} can have as many rows as desired; see examples.
}
\examples{
library(gamlss)


data(aids)
aids1<-gamlss(y~pb(x,df=4)+qrt,data=aids,family=PO)
m.aids1<-gamlssReport(aids1)
centile.gamlssReport(m.aids1,y=c(14,11,20),newdata=data.frame(x=4,qrt=4))
centile.gamlssReport(m.aids1,y=14,newdata=data.frame(x=c(4,5),qrt= c(4,2)))
centile.gamlssReport(m.aids1,y=c(14,12),newdata=data.frame(x=c(4,5),qrt= c("reference",2) ))
}
\seealso{
\code{\link{gamlssReport}}, \code{\link{predict.gamlssReport}}
}
\author{
Rok Blagus, \email{rok.blagus@mf.uni-lj.si}
}
