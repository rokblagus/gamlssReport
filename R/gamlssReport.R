




#' @title Extract information from the GAMLSS fitted objects
#' @description The function extracts the paramaters of the fitted GAMLSS object.
#'
#'
#' @param object an object of class \code{"gamlss"}.
#' @param extract.smooth logical indicating if the prediction for the smooth terms should be extracted (see details). Defaults to \code{FALSE}.
#' Note: multiple smooth terms per parameter are allowed (these smooth terms had to be obtained using \code{pb()}).
#'
#' @details In the presence of additive terms there are two approaches for making interpolation (for the observed data points the two approaches return the same result): 1) interpolating the values of the B-spline basis for the new datapoints and using the estimated penalized coefficients to obtain predictions and 2) spline interpolation is performed using the estimated smooth terms for the observed datapoints. The latter approach (implemented in \code{\link{predict.gamlss}} using natural splines) requires that the fitted smooth terms and original datapoints for the variable(s) used to form the spline are saved. Setting \code{extract.smooth=TRUE} will save these and enable to exactly replicate the results of \code{\link{predict.gamlss}} when using \code{\link{predict.gamlssReport}}. Note that this will enable sharing the potentially sensible information!
#' @return An object of class \code{"gamlssReport"} for which \code{print}, \code{plot}, and \code{predict} functions are available.
#' A list with components
#' \itemize{
#' \item{\code{family}} {the distribution family used when fitting GAMLSS}
#' \item{\code{params}} {the names of the fitted parameters, i.e. mu, sigma, nu, tau}
#' \item{\code{link}} {the link function for the model used for each parameter (a named list, with names corresponding to the parameters)}
#' \item{\code{coef.beta}} {the linear coefficients for each parameter (a named list, with names corresponding to the parameters)}
#' \item{\code{coef.spline}} {the penalized coefficients for the P-spline for each parameter where pb was used to model non-linear association (a named list, with names corresponding to the parameters).}
#' \item{\code{knots.spline}} {the knots used to form the B-spline basis for the P-spline (a named list, with names corresponding to the parameters)}
#' \item{\code{range.x}} {the range of variable that was used for P-spline (a named list, with names corresponding to the parameters; \code{NULL} for parameters where \code{pb()} was not used)}
#' \item{\code{degree.spline}} {the degree of the spline (a named list, with names corresponding to the parameters; \code{NULL} for parameters where \code{pb()} was not used)}
#' \item{\code{terms}} {details about the fixed effects and nonlinear term}
#' \itemize{
#' \item{\code{fixformula}} {formula for the linear effects (a named list, with names corresponding to the parameters)}
#' \item{\code{splinevar}} {names of the variables that were used for P-splines (a named list, with names corresponding to the parameters; \code{NULL} for parameters where \code{pb()} was not used)}
#' }
#' \item{\code{term.smooth} if \code{extract.smooth=TRUE} the fitted smooth terms for the observed datapoints}
#' }
#' @seealso \code{\link{plot.gamlssReport}}, \code{\link{print.gamlssReport}}, \code{\link{predict.gamlssReport}}, \code{\link{centile.gamlssReport}}, \code{\link{score.gamlssReport}}
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @export
#' @examples
#' library(gamlss)
#'
#'
#' data(aids)
#' aids1<-gamlss(y~pb(x,df=4)+qrt,data=aids,family=PO)
#' m.aids1<-gamlssReport(aids1)
#' m.aids1
#'
#' data(abdom)
#' mod<-gamlss(y~pb(x),sigma.fo=~pb(x),family=BCT, data=abdom, method=mixed(1,20))
#' m.mod<-gamlssReport(mod)
#' m.mod



gamlssReport<-function(object,extract.smooth=FALSE){

  .family<-object$family

  params<-object$parameters
	.terms<-extract_terms(object)

  .link<-coef.beta<-coef.spline<-knots.spline<-range.x<-degree.spline<-vector("list",length=length(params))

  names(.link)<-names(coef.beta)<-names(coef.spline)<-names(knots.spline)<-names(range.x)<-names(degree.spline)<-params



  for (i in params){ #currently assuming only one pb term!, assuming there is no offset!
    #.terms[[i]]<-object[[paste(i,".terms",sep="")]]
    .link[[i]]<-object[[paste(i,".link",sep="")]]
    coef.beta[[i]]<-object[[paste(i,".coefficients",sep="")]]

		if (!is.null(object[[paste(i,".coefSmo",sep="")]])){
			coef.spline[[i]]<-knots.spline[[i]]<-degree.spline[[i]]<-range.x[[i]]<-
					vector("list",length=length(object[[paste(i,".coefSmo",sep="")]]))
			names(coef.spline[[i]])<-names(knots.spline[[i]])<-names(degree.spline[[i]])<-names(range.x[[i]])<-.terms$splinevar[[i]]
			zz=0
			for (ii in .terms$splinevar[[i]]){
			zz=zz+1
    				coef.spline[[i]][[ii]]<-object[[paste(i,".coefSmo",sep="")]][[zz]]$coef
    				knots.spline[[i]][[ii]]<-object[[paste(i,".coefSmo",sep="")]][[zz]]$knots
     				degree.spline[[i]][[ii]]<-3 #not general!

    		      	x<-object$mu.x[, attributes(object$mu.terms)$term.labels[grepl(object$mu.coefSmo[[zz]]$name,attributes(object$mu.terms)$term.labels)] ]#this could be an issue!!
      			range.x[[i]][[ii]]<- c(min(x),max(x))
    				}
			}
  }



  res<-list(family=.family,params=params,link=.link,
       coef.beta=coef.beta,
       coef.spline=coef.spline,knots.spline=knots.spline,range.x=range.x,degree.spline=degree.spline,terms=.terms)

  fun.extract<-function(x){
    what<-x$parameters
    nms<-paste0(what,".s")
    y<-x[nms]
    names(y)<-what
    nms2<-paste0(what,".coefSmo")
    x<-lapply(x[nms2],function(x) lapply(x,function(x){
      Call <- object$call;
      data<-eval(Call$data);
      data[,x$name]
    }  ))
    x<-lapply(x, function(x){
      if (length(x)==0) m<-NULL else {
      nsmooth<-length(x)
      nrows<-length(x[[1]])
      m<-matrix(NA,ncol=nsmooth,nrow=nrows)
      for (ii in 1:length(x)){
         m[,ii]<- x[[ii]]
      }
      }
      m
    })
    names(x)<-what
    list(x=x,y=y)
  }
  if (extract.smooth==TRUE){
    .term.smooth<-fun.extract(object)
    res$term.smooth<-.term.smooth
  }

  class(res)<-"gamlssReport"
res
}





#' @title Print function
#' @description Prints the extracted information that can be used in the publications
#'
#' @param x an object of class \code{"gamlssReport"} obtained as a call to \code{\link{gamlssReport}}.
#' @param digits_range a numeric value determining how many digits to display for range of x; defaults to \code{Inf} (no rounding).
#' @param digits_coefs a numeric value determining how many digits to display for coefficients (fixed and penalized); defaults to \code{Inf} (no rounding).
#' @param digits_knots a numeric value determining how many digits to display for the knots of the B-spline basis; defaults to \code{Inf} (no rounding).
#' @param ... for extra arguments.
#'
#'
#' @seealso \code{\link{gamlssReport}}, \code{\link{plot.gamlssReport}}
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @export
#' @examples
#' library(gamlss)
#'
#'
#' data(aids)
#' aids1<-gamlss(y~pb(x,df=4)+qrt,data=aids,family=PO)
#' m.aids1<-gamlssReport(aids1)
#' print.gamlssReport(m.aids1)
#'
#' data(abdom)
#' mod<-gamlss(y~pb(x),sigma.fo=~pb(x),family=BCT, data=abdom, method=mixed(1,20))
#' m.mod<-gamlssReport(mod)
#' print.gamlssReport(m.mod)

print.gamlssReport<-function(x,digits_range=Inf,digits_coefs=Inf,digits_knots=Inf,...){
  txt1<-"********************************************************\n"
  txt2<-"--------------------------------------------------------\n"

  cat(txt1,"Family: ", paste(x$family,collapse=", "), "; parameters: ",paste(x$params,collapse=", "),"\n",
      txt2,sep="")

  for (i in x$params){
    cat(i," link function: ", x$link[[i]], "\n", i," coefficients:\n"
        ,sep="" )

    cfs<-as.matrix(round(x$coef.beta[[i]],digits_coefs),ncol=1)
    colnames(cfs)<-"Estimate"

    print(cfs)

    if (!is.null(x$terms$splinevar[[i]])){

		splvar<-x$terms$splinevar[[i]]
		for (ii in 1:length(splvar)){
      		cat("\n",i," P-spline details for ",splvar[ii],": range of ",splvar[ii] ," ",paste("from ", round(x$range.x[[i]][[ii]][1],digits_range)," to ",round(x$range.x[[i]][[ii]][2],digits_range) ,sep=""),"; P-spline degree ",x$degree.spline[[i]][[ii]],
          		"\n",i," P-spline knots and coefficients\n",sep="")



      		resi<-rbind(round(x$knots.spline[[i]][[ii]],digits_knots),t(round(x$coef.spline[[i]][[ii]],digits_coefs)))
      		rownames(resi)<-c("knots","coefs")
      		colnames(resi)<-paste("P-spline ",1:length(x$knots.spline[[i]][[ii]]),sep="")

      		print(t(resi))

			}
    }
    cat("\n",txt2,sep="")

  }



}







#' @title Plot centiles
#'
#' @description Function \code{"plot.gamlssReport"} plots the centile curves for the object created by \code{\link{gamlssReport}}.
#'
#'
#' @param x An object of class \code{"gamlssReport"} representing the GAMLSS model summarized by the function \code{\link{gamlssReport}}.
#' @param xname name (a character of length one) of the x variable used on the x-axis of the plot.
#' Note, the name must exactly match the name of the variable that was used when fitting GAMLSS.
#' @param range.x a numeric vector of length 2 specifying the range of x for which to show the centiles.
#' @param gamlss.prediction logical, if \code{TRUE} it uses the spline interpolation approach implemented in \code{predict.gamlss}, otherwise it uses the interpolation of the B-spline basis and the estimated penalized coefficients. Defaults to \code{FALSE}. Note: if set to \code{TRUE}, the \code{\link{gamlssReport}} had to be used with \code{extract.smooth=TRUE}, oterwise the error is returned.
#' @param x.transform a function used to transform the x-axis, see details. Defaults to \code{NULL} which means that no transformation is used.
#' @param centiles a numeric vector with entries in (0,1) containing centiles to be shown. Defaults to \code{c(0.1,0.5,0.9)}.
#' @param newdata the data.frame with a single row containig the values of the other variables that were used when fitting the model, set to \code{NULL} (default) if the model only has a single x. Needs to be nonnull if there are more xs. The names of the variables must exactly match the names used when fitting GAMLSS.
#' @param return.object logical, if TRUE an object that can be used to make a plot is returned instead of the plot. Defaults to \code{FALSE}.
#' @param seq.length length of the sequence for "new x". Defaults to \code{1e5}.
#' @param ... other arguments controling the appearance of the plot.
#'
#' @return if \code{return.object=TRUE} returns a list with components
#' \itemize{
#' \item{\code{x}} {vector of length \code{seq.length} containing xs for which the scores are calculated}
#' \item{\code{y}} {a list of length \code{length(centiles)} containig the scores; each element of the list is a vector of length \code{length(centiles)}}
#' \item{\code{cent}} {vector of length \code{length(centiles)} listing the centiles}
#' }
#' @details When fitting the model, a power transformation of some variable, say xa=x^p, is used. The argument \code{x.transform} can be used to plot the centiles on the x scales (or any other if desired by the user) instead of the default (\code{x.transform=NULL}) scale. E.g. if  xa=x^(1/2) was used, specifying \code{x.transform=function(x) x^2} will plot the centiles on the x scale. Note that the \code{range.x} always refers to the variable that was used when fitting the model. See examples.
#' @seealso \code{\link{gamlssReport}}, \code{\link{print.gamlssReport}}
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @export
#' @examples
#' library(gamlss)
#'
#'
#' data(aids)
#' aids1<-gamlss(y~pb(x,df=4)+qrt,data=aids,family=PO)
#' m.aids1<-gamlssReport(aids1)
#' plot.gamlssReport(x=m.aids1,xname="x",range.x=c(1,45),
#'    x.transform=NULL,centiles=seq(from=0.1,to=0.9,by=0.2),
#'    newdata=data.frame(qrt=c("reference")),xlab="x",ylab="score")
#'
#' aids$nx<-sqrt(aids$x)
#' aids2<-gamlss(y~pb(nx,df=4)+qrt,data=aids,family=PO)
#' m.aids2<-gamlssReport(aids2)
#'
#' ##plot on the fitted scale
#' plot.gamlssReport(x=m.aids2,xname="nx",range.x=sqrt(c(1,45)),
#'    x.transform=NULL,centiles=seq(from=0.1,to=0.9,by=0.2),
#'    newdata=data.frame(qrt=c("reference")),xlab="nx",ylab="score")
#'
#' ##plot on the original scale
#' plot.gamlssReport(x=m.aids2,xname="nx",range.x=sqrt(c(1,45)),
#'    x.transform=function(x) x**2,centiles=seq(from=0.1,to=0.9,by=0.2),
#'    newdata=data.frame(qrt=c("reference")),xlab="x",ylab="score")
#'
#'
#' data(abdom)
#' mod<-gamlss(y~pb(x),sigma.fo=~pb(x),family=BCT, data=abdom, method=mixed(1,20))
#' m.mod<-gamlssReport(mod)
#' plot.gamlssReport(x=m.mod,xname="x",range.x=c(13,40),
#'    x.transform=NULL,centiles=seq(from=0.1,to=0.9,by=0.2),
#'    xlab="x",ylab="score")


plot.gamlssReport<-function(x,xname,range.x,gamlss.prediction=FALSE,x.transform=NULL,centiles=c(0.1,0.5,0.9),newdata=NULL,return.object=FALSE,seq.length=1e5,...){
object<-x
  params<-object$params


  #make predictions for a sequence of x

  x.new<-seq(from=range.x[1], to=range.x[2],length.out=seq.length)

  dfn<-data.frame(v1=x.new)
  names(dfn)<-xname

  if (is.null(newdata)) df<-dfn else {

    df<-cbind(dfn,newdata)

  }



  yf<-vector("list",length=length(centiles))
  ii=0
  for (cet in centiles){
    ii=ii+1
    yf[[ii]]<-score.gamlssReport(object,cet,df,gamlss.prediction)
  }

  if (return.object==FALSE){
    mn<-min(unlist(yf))
    mx<-max(unlist(yf))
    if (!is.null(x.transform)){
      range.x<-x.transform(range.x)
      x.new<-x.transform(x.new)
    }
    plot(1:10,col="white",xlim=range.x,ylim=c(mn,mx),...)
    for (j in 1:length(yf)){
      lines(x.new,yf[[j]],col=j)
    }
    legend("topleft",legend=centiles,lty=1,col=1:length(centiles),bty="n")
  } else {
    res<-list(x=x.new,y=yf,cent=centiles)
    res
  }



}








#' @title Predict function
#' @description Using the object generated by \code{\link{gamlssReport}}, it makes predictions (for the parameters), based on the new data
#'
#' @param object the object of class \code{"gamlssReport"} generated by \code{\link{gamlssReport}}.
#' @param newdata A dataframe containing all the variables (the names must match comepletely) needed to make predictions; see details for factors.
#' @param gamlss.prediction logical, if \code{TRUE} it uses the spline interpolation approach implemented in \code{predict.gamlss}, otherwise it uses the interpolation of the B-spline basis and the estimated penalized coefficients. Defaults to \code{FALSE}. Note: if set to \code{TRUE}, the \code{\link{gamlssReport}} had to be used with \code{extract.smooth=TRUE}, oterwise the error is returned.
#' @param ... for extra arguments.
#'
#'
#' @return a list with components
#' \itemize{
#' \item{\code{lp}} {a named list containing the linear predictors (the names corespond to the fitted parameters, e.g. \code{mu}, \code{sigma}, \code{nu}, and \code{tau})}
#' \item{\code{fv}} {a named list containing the fitted values (the names corespond to the fitted parameters, e.g. \code{mu}, \code{sigma}, \code{nu}, and \code{tau})}
#' }
#' @details For factors enter a value; this value must be an existing level of the factor. For the reference value enter \code{"reference"}.
#'
#' @seealso \code{\link{gamlssReport}}
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @export
#' @examples
#' library(gamlss)
#'
#'
#' data(aids)
#' aids1<-gamlss(y~pb(x,df=4)+qrt,data=aids,family=PO)
#' m.aids1<-gamlssReport(aids1)
#' new_data<-data.frame(x=c(10,15,2,4),qrt=c("reference",3,4,2) )
#' predict.gamlssReport(m.aids1,newdata=new_data   )
#'
#' #multiple smooth terms
#' set.seed(1)
#' aids$z<-runif(nrow(aids))
#' fit<-gamlss(y~pb(z,df=10)+qrt+pb(x,df=5),data=aids,family=PO)
#' m<-gamlssReport(aids1)
#' newm<-data.frame(qrt="3",z=0.25,x=2)
#' predict(m,newm)$lp$mu
#' #gamlss
#' predict(fit,newdata=data.frame(qrt=3,z=0.25,x=2),what="mu",type="link")
#'
#' m2<-gamlssReport(aids1,extract.smooth=TRUE)
#' predict(m2,newm,gamlss.prediction=TRUE)$lp$mu #the same as gamlss



predict.gamlssReport<-function(object,newdata,gamlss.prediction=FALSE,...){
  splinevar<-object$terms$splinevar
  fixformula<-lapply(object$terms$fixformula,as.formula)


  params<-object$params
  list.nl<-list.l<-vector("list",length=length(params))
  names(list.nl)<-names(list.l)<-params

  for (i in params){


    list.l[[i]]<-model.matrix(fixformula[[i]],newdata)

    if (!is.null(object$terms$splinevar[[i]])){
      list.nl[[i]]<-vector("list",length(object$terms$splinevar[[i]]))
	names(list.nl[[i]])<-object$terms$splinevar[[i]]
	for (ii in object$terms$splinevar[[i]]){
	list.nl[[i]][[ii]]<-newdata[,ii]
	}

    }

  }

  make_prediction(object,list.nl,list.l,gamlss.prediction)


}




#' @title Calculate the centiles
#' @description Using the object generated by \code{\link{gamlssReport}}, the function calculates the centile for a given y and x
#'
#' @param object the object of class \code{"gamlssReport"} generated by \code{\link{gamlssReport}}.
#' @param y the score for which to calculate the centile.
#' @param newdata a dataframe containing the Xs for which to evaluate the centile.
#' @param gamlss.prediction logical, if \code{TRUE} it uses the spline interpolation approach implemented in \code{predict.gamlss}, otherwise it uses the interpolation of the B-spline basis and the estimated penalized coefficients. Defaults to \code{FALSE}. Note: if set to \code{TRUE}, the \code{\link{gamlssReport}} had to be used with \code{extract.smooth=TRUE}, oterwise the error is returned.
#' @return a vector containing the centiles
#' @details \code{y} can be a vector, in this case \code{newdata} need to be a dataframe with single row, or dataframe with the same no of rows as \code{y}; if \code{y} is a scalar \code{newdata} can have as many rows as desired; see examples.
#'
#' @seealso \code{\link{gamlssReport}}, \code{\link{predict.gamlssReport}}
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @export
#' @examples
#' library(gamlss)
#'
#'
#' data(aids)
#' aids1<-gamlss(y~pb(x,df=4)+qrt,data=aids,family=PO)
#' m.aids1<-gamlssReport(aids1)
#' centile.gamlssReport(m.aids1,y=c(14,11,20),newdata=data.frame(x=4,qrt=4))
#' centile.gamlssReport(m.aids1,y=14,newdata=data.frame(x=c(4,5),qrt= c(4,2)))
#' centile.gamlssReport(m.aids1,y=c(14,12),newdata=data.frame(x=c(4,5),qrt= c("reference",2) ))


centile.gamlssReport<-function(object,y,newdata,gamlss.prediction=FALSE){

  pred<-predict.gamlssReport(object,newdata,gamlss.prediction)$fv
  family<-object$family[1]

  f<-paste("p",family,sep="")
  fun <- get(f)
  if (length(pred[[1]])==1&length(y)>1) pred<-lapply(pred,function(x,y) rep(x,length(y)),y)

  pred[[names(formals(fun))[1]]]<-y

  do.call(fun,pred)

}








#' @title Calculate the score for given centiles
#' @description Using the object generated by \code{\link{gamlssReport}}, the function calculates the score for a given centile and x
#'
#' @param object the object of class \code{"gamlssReport"} generated by \code{\link{gamlssReport}}.
#' @param centile the centile for which to calculate the score.
#' @param newdata a dataframe containing the Xs for which to evaluate the centile.
#' @param gamlss.prediction logical, if \code{TRUE} it uses the spline interpolation approach implemented in \code{predict.gamlss}, otherwise it uses the interpolation of the B-spline basis and the estimated penalized coefficients. Defaults to \code{FALSE}. Note: if set to \code{TRUE}, the \code{\link{gamlssReport}} had to be used with \code{extract.smooth=TRUE}, oterwise the error is returned.
#' @return a vector containing the scores
#' @details \code{centile} can be a vector, in this case \code{newdata} need to be a dataframe with single row, or dataframe with the same no of rows as \code{centile}; if \code{centile} is a scalar \code{newdata} can have as many rows as desired; see examples.
#'
#' @seealso \code{\link{gamlssReport}}, \code{\link{predict.gamlssReport}}
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @export
#' @examples
#' library(gamlss)
#'
#'
#' data(aids)
#' aids1<-gamlss(y~pb(x,df=4)+qrt,data=aids,family=PO)
#' m.aids1<-gamlssReport(aids1)
#' score.gamlssReport(m.aids1,centile=c(0.1,0.5,0.9),newdata=data.frame(x=4,qrt=4))
#' score.gamlssReport(m.aids1,centile=0.1,newdata=data.frame(x=c(4,5),qrt= c(4,2)))
#' score.gamlssReport(m.aids1,centile=c(0.1,0.9),newdata=data.frame(x=c(4,5),qrt= c("reference",2) ))


score.gamlssReport<-function(object,centile,newdata,gamlss.prediction=FALSE){

  pred<-predict.gamlssReport(object,newdata,gamlss.prediction  )$fv
  family<-object$family[1]

  f<-paste("q",family,sep="")
  fun <- get(f)
  if (length(pred[[1]])==1&length(centile)>1) pred<-lapply(pred,function(x,y) rep(x,length(y)),centile)


  pred[[names(formals(fun))[1]]]<-centile

  do.call(fun,pred)

}




#' @title Internal predict function
#' @description Internal function used by \code{\link{predict.gamlssReport}}; not intended to be used by the user.
#'
#' @param object the object of class \code{"gamlssReport"} generated by \code{\link{gamlssReport}}.
#' @param xnew.spline a named list containing a vector of values for the variable that is modeled as P-spline in \code{pb()} in a call to \code{gamlss}; needs to be of the same length as the number of parameters (\code{mu}, \code{sigma}, \code{nu}, \code{tau}), \code{NULL} if the parameter does not have a P-spline term.
#' @param z.new a named list containing fixed effects design matrix for each parameter.
#' z.new[[ii]] also needs to be a named list, with the names being equal to the names of the spline var
#' @param gamlss.prediction logical, if \code{TRUE} it uses the spline interpolation approach implemented in \code{predict.gamlss}, otherwise it uses the interpolation of the B-spline basis and the estimated penalized coefficients. Defaults to \code{FALSE}. Note: if set to \code{TRUE}, the \code{\link{gamlssReport}} had to be used with \code{extract.smooth=TRUE}, oterwise the error is returned.
#' @seealso \code{\link{predict.gamlssReport}}
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @export
#' @examples
#' library(gamlss)
#'
#'
#' data(aids)
#' aids1<-gamlss(y~pb(x,df=4)+qrt,data=aids,family=PO)
#' m.aids1<-gamlssReport(aids1)
#' make_prediction(m.aids1,list("mu"=list("x"=c(1,6,45))),
#' list("mu"=rbind(c(1,1,0,0,0),c(1,6,1,0,0),c(1,45,1,0,0))))




make_prediction<-function(object,xnew.spline,z.new,gamlss.prediction=FALSE){

  if (gamlss.prediction==TRUE&is.null(object$term.smooth)) stop("Smooth terms were not extracted by gamlssReport. Use gamlssReport with extract.smooth=TRUE.")

  funi<-function(i,x,s,xnew.spline){ #this can be used for each param, to get smooth preds for ech nonlinear term
    spline(x=x[,i],y=s[,i],xout=xnew.spline[[i]],method="natural")$y
  }
  params<-object$params


  fv<-lp<-vector("list",length=length(params))
  names(fv)<-names(lp)<-params

  for (i in params){

    cfi<-object$coef.beta[[i]]

    if (!is.null(object$terms$splinevar[[i]])){

      if (gamlss.prediction==FALSE){
  cfs<-D<-vector("list",length=length(object$terms$splinevar[[i]]))
	names(cfs)<-names(D)<-object$terms$splinevar[[i]]
	for (ii in object$terms$splinevar[[i]] ){
      cfs[[ii]]<-object$coef.spline[[i]][[ii]]
      D[[ii]]<-make_spline(xnew.spline[[i]][[ii]],object$knots.spline[[i]][[ii]],object$degree.spline[[i]][[ii]],object$range.x[[i]][[ii]])
	}#end for
	} #end if
}
    if (!is.null(object$terms$splinevar[[i]])){
	lp[[i]]<-z.new[[i]]%*%matrix(cfi,ncol=1)
	for (ii in object$terms$splinevar[[i]]){
	  if (gamlss.prediction==FALSE){lp[[i]]<-lp[[i]]+D[[ii]]%*%cfs[[ii]]} else {
	    whichs<-which(object$terms$splinevar[[i]]==ii)
	    s<-spline(x=object$term.smooth$x[[i]][,whichs],y=object$term.smooth$y[[i]][,whichs],xout=xnew.spline[[i]][[ii]],method="natural")$y #!!!!
	    lp[[i]]<-lp[[i]]+s
	  }
	}
    } else {
      lp[[i]]<-z.new[[i]]%*%matrix(cfi,ncol=1)
    }
    fv[[i]]<-make_inverse(object$link[[i]],lp[[i]])

  }

  list(lp=lp,fv=fv)

}



#' @title Aux extract function
#' @description aux function used by \code{\link{gamlssReport}}; not intendent to be used by the user.
#'
#'
#' @param object An object of class \code{"gamlss"}.
#' @seealso \code{\link{gamlssReport}}, \code{\link{print.gamlssReport}}
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @export
#' @examples
#' library(gamlss)
#'
#' data(aids)
#' aids1<-gamlss(y~pb(x,df=4)+qrt,data=aids,family=PO)
#' extract_terms(aids1)


extract_terms<-function(object){
  params<-object$parameters

  splinevar<-fixformula<-vector("list",length=length(params))
  names(splinevar)<-names(fixformula)<-params

  for (ii in params){

    nms<-names(unlist(attributes(object[[paste(ii,".terms",sep="")]])$dataClasses))[-1] #ignore the outcome var

    if (length(nms)==0){
      fixformula[[ii]]<-as.formula(~1)

    } else {

      att<-attributes(object[[paste(ii,".terms",sep="")]])
      logicf<-att$dataClasses[-1]=="factor"
      if (sum(logicf)>0){
        for (jj in which(logicf==TRUE)){
          nmsj<-nms[jj]

          x<-colnames(object[[paste(ii,".x",sep="")]])
          rf<-"\"reference\""
          lvl<-c(rf,unlist(lapply(strsplit(x[grepl(nmsj,x)],nmsj),function(x) x[2]))  )
          levelsj<-paste("c(",paste(lvl,collapse=",") ,")",sep="")
          newnamej<-paste("factor(",nmsj,",levels=",levelsj,")",sep="")
          nms[jj]<-newnamej
        }

      }

      logic<-grepl("pb",nms)

      if (sum(logic)!=0) {
		pbvar<-rep(NA,sum(logic))
		vars<-nms
		zz=0
		for (jj in which(logic==TRUE)){
			zz=zz+1
        		spl<-strsplit(nms[jj],"pb\\(")[[1]][2]

        		if (length(strsplit(spl,",")[[1]])==1) pbvar[zz]<-strsplit(spl,")")[[1]] else  pbvar[zz]<-strsplit( spl  ,",")[[1]][1]
        		if (zz==sum(logic)) splinevar[[ii]]<-pbvar
        		vars[jj]<-pbvar[zz]
			}


        #fixf<-as.formula(paste("~",paste(vars,collapse="+"),sep=""))
        fixf<-paste("~",paste(vars,collapse="+"),sep="")

      } else {
        #fixf<-as.formula(paste("~",paste(nms,collapse="+"),sep=""))
        fixf<-paste("~",paste(nms,collapse="+"),sep="")
      }

      #fixformula[[ii]]<-format(fixf) #it works if we change predict function also
      fixformula[[ii]]<-fixf
    }
  }

  list(fixformula=fixformula,splinevar=splinevar)
}





#' @title Aux function to recreate the B-spline basis
#' @description aux function used by \code{\link{predict.gamlssReport}}; not intendent to be used by the user.
#'
#'
#' @param newx vector containing values of x for which the B-spline bases are to be calculated.
#' @param knots the vector of knots; currently only the \code{length(knots)} matters (ie we assumed that knots in call to \code{gamlss} were uniformly set)
#' @param degree the degree of the spline used in \code{pb()}.
#' @param rng vector with 2 elements, min(x) and max(x)
#'
#' @seealso \code{\link{predict.gamlssReport}}
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @export
#' @examples
#'
#'
#' make_spline(10,seq(from=10,to=20,by=1),3,c(9,21))



make_spline<-function(newx,knots,degree,rng){ #this will only work if in pb they used uniformly set knots!
##taken from pb gamlss
  bbase <- function(x, xl, xr, ndx, deg) #note we dont allow quantile=TRUE option!
  {
    tpower <- function(x, t, p)
      # Truncated p-th power function
      (x - t) ^ p * (x > t)
    # DS xl= min, xr=max, ndx= number of points within
    # Construct B-spline basis
    # if quantiles=TRUE use different bases
    dx <- (xr - xl) / ndx # DS increment

      knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
      P <- outer(x, knots, tpower, deg)# calculate the power in the knots
      n <- dim(P)[2]
      D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg) #
      B <- (-1) ^ (deg + 1) * P %*% t(D)
      attr(B, "knots") <- knots[-c(1:(deg-1), (n-(deg-2)):n)]
      B

  }
##main f starts here

  xr<-rng[2]
  xl<-rng[1]
  ndx<-length(knots)-degree #number of knots-degree!
  xmax<-xr + 0.01 * (xr - xl)
  xmin<-xl - 0.01 * (xr - xl)
  #dt<-(xmax - xmin) / ndx #use if splineDesign
  #knots2<-seq(xmin - degree * dt, xmax + degree * dt, by = dt)

  #B<-splineDesign(knots = knots2, x = newx, ord = degree + 1, derivs = 0,outer.ok = TRUE)
  B<-bbase(x=newx, xl=xmin, xr=xmax, ndx, deg=degree)

  B

}





#' @title Inverse of the link function
#'
#' @param f the name of the link function (a character).
#' @param x the values at which to evaluate the inverse of the link function.
#'
#' @seealso \code{\link{predict.gamlssReport}}
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @export
#' @examples
#' make_inverse("log",1)


make_inverse<-function(f,x){ #add other options!
  if (f=="log") y<-exp(x)
  if (f=="exp") y<-log(x) #check to see if the names are ok!
  if (f=="logit") y<-1/(1+exp(-x))
  if (f=="inverse") y<-1/x
  if (f=="identity") y<-x
  y
}




