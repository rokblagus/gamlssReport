




#' @title Extract information from the GAMLSS fitted objects
#' @description The function extracts the paramaters of the fitted GAMLSS object.
#'
#'
#' @param object an object of class \code{"gamlss"}.
#' Note: only one smooth term per parameter is assumed (this smooth term had to be obtained using pb()).
#'
#'
#' @seealso \code{\link{plot.gamlssReport}}, \code{\link{print.gamlssReport}}, \code{\link{predict.gamlssReport}}, \code{\link{centile.gamlssReport}}, \code{\link{score.gamlssReport}}
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @export
#' @examples
#' library(gamlss)
#' library(splines)
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



gamlssReport<-function(object){

  .family<-object$family

  params<-object$parameters

  #.terms<-
  .link<-coef.beta<-coef.spline<-knots.spline<-range.x<-degree.spline<-vector("list",length=length(params))
  #names(.terms)<-
  names(.link)<-names(coef.beta)<-names(coef.spline)<-names(knots.spline)<-names(range.x)<-names(degree.spline)<-params


  for (i in params){ #currently assuming only one pb term!, assuming there is no offset!
    #.terms[[i]]<-object[[paste(i,".terms",sep="")]]
    .link[[i]]<-object[[paste(i,".link",sep="")]]
    coef.beta[[i]]<-object[[paste(i,".coefficients",sep="")]]
    coef.spline[[i]]<-object[[paste(i,".coefSmo",sep="")]][[1]]$coef
    knots.spline[[i]]<-object[[paste(i,".coefSmo",sep="")]][[1]]$knots
    if (!is.null(object[[paste(i,".coefSmo",sep="")]][[1]]$coef)) degree.spline[[i]]<-3 #not general!

    if (!is.null(object[[paste(i,".coefSmo",sep="")]][[1]]$coef)) {
      x<-object$mu.x[, attributes(object$mu.terms)$term.labels[grepl(object$mu.coefSmo[[1]]$name,attributes(object$mu.terms)$term.labels)] ]#this could be an issue!!
      range.x[[i]]<- c(min(x),max(x))
    }

  }

  .terms<-extract_terms(object)

  res<-list(family=.family,params=params,link=.link,
       coef.beta=coef.beta,
       coef.spline=coef.spline,knots.spline=knots.spline,range.x=range.x,degree.spline=degree.spline,terms=.terms)
class(res)<-"gamlssReport"
res
}





#' @title Print function
#' @description Prints the extracted information that can be used in the publications
#'
#' @param x an object of class \code{"gamlssReport"} obtained as a call to \code{\link{gamlssReport}}.
#' @param digits_range a numeric value determining how many digits to display for range of x; defaults to Inf (no rounding).
#' @param digits_coefs a numeric value determining how many digits to display for coefficients (fixed and penalized); defaults to Inf (no rounding).
#' @param digits_knots a numeric value determining how many digits to display for the knots of the B-spline basis; defaults to Inf (no rounding).
#'
#'
#' @seealso \code{\link{gamlssReport}}, \code{\link{plot.gamlssReport}}
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @export
#' @examples
#' library(gamlss)
#' library(splines)
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

print.gamlssReport<-function(x,digits_range=Inf,digits_coefs=Inf,digits_knots=Inf){
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

    if (!is.null(x$degree.spline[[i]])){
      cat("\n",i," P-spline details: range of x ",paste("from ", round(x$range.x[[i]][1],digits_range)," to ",round(x$range.x[[i]][2],digits_range) ,sep=""),"; P-spline degree ",x$degree.spline[[i]],
          "\n",i," P-spline knots and coefficients\n",sep="")



      resi<-rbind(round(x$knots.spline[[i]],digits_knots),t(round(x$coef.spline[[i]],digits_coefs)))
      rownames(resi)<-c("knots","coefs")
      colnames(resi)<-paste("P-spline ",1:length(x$knots.spline[[i]]),sep="")

      print(t(resi))


    }
    cat("\n",txt2,sep="")

  }



}







#' @title Plot centiles
#'
#' @description Function \code{"plot.gamlssReport"} plots the centile curves for the object created by \code{"gamlssReport"}.
#'
#'
#' @param x An object of class \code{"gamlssReport"} representing the GAMLSS model summarized by the function \code{"gamlssReport"}.
#' @param xname name (a character of length one) of the x variable used on the x-axis of the plot.
#' Note, the name must exactly match the name of the variable that was used when fitting GAMLSS.
#' @param range.x a numeric vector of length 2 specifying the range of x for which to show the centiles.
#' @param centiles a numeric vector with entries in (0,1) containing centiles to be shown. Defaults to \code{"c(0.1,0.5,0.9)"}.
#' @param newdata the data.frame with a single row containig the values of the other variables that were used when fitting the model, set to NULL (default) if the model only has a single x. Needs to be nonnull if there are more xs. The names of the variables must exactly match the names used when fitting GAMLSS.
#' @param return.object logical, if TRUE an object that can be used to make a plot is returned instead of the plot. Defaults to \code{"FALSE"}.
#' @param seq.length length of the sequence for "new x". Defaults to \code{"1e5"}.
#' @param ... other arguments controling the appearance of the plot.
#'
#' @seealso \code{\link{gamlssReport}}, \code{\link{print.gamlssReport}}
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @export
#' @examples
#' library(gamlss)
#' library(splines)
#'
#' data(aids)
#' aids1<-gamlss(y~pb(x,df=4)+qrt,data=aids,family=PO)
#' m.aids1<-gamlssReport(aids1)
#' plot.gamlssReport(m.aids1,"x",c(1,45),seq(from=0.1,to=0.9,by=0.2),data.frame(qrt=c("reference")))
#'
#' data(abdom)
#' mod<-gamlss(y~pb(x),sigma.fo=~pb(x),family=BCT, data=abdom, method=mixed(1,20))
#' m.mod<-gamlssReport(mod)
#' plot.gamlssReport(m.mod,"x",c(13,40),seq(from=0.1,to=0.9,by=0.2))


plot.gamlssReport<-function(x,xname,range.x,centiles=c(0.1,0.5,0.9),newdata=NULL,return.object=FALSE,seq.length=1e5,...){
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
    yf[[ii]]<-score.gamlssReport(object,cet,df)
  }

  if (return.object==FALSE){
    mn<-min(unlist(yf))
    mx<-max(unlist(yf))

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
#' @param x the object of class \code{"gamlssReport"} generated by \code{\link{gamlssReport}}.
#' @param newdata A dataframe containig all the variables (the names must match comepletely) needed to make predictions; see details for factors.
#'
#' @details For factors enter a value; this value must be an existing level of the factor. For the reference value enter \code{"reference"}.
#'
#' @seealso \code{\link{gamlssReport}}
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @export
#' @examples
#' library(gamlss)
#' library(splines)
#'
#' data(aids)
#' aids1<-gamlss(y~pb(x,df=4)+qrt,data=aids,family=PO)
#' m.aids1<-gamlssReport(aids1)
#' new_data<-data.frame(x=c(10,15,2,4),qrt=c("reference",3,4,2) )
#' predict.gamlssReport(m.aids1,newdata=new_data   )


predict.gamlssReport<-function(x,newdata){
object<-x
  splinevar<-object$terms$splinevar
  fixformula<-object$terms$fixformula


  params<-object$params
  list.nl<-list.l<-vector("list",length=length(params))
  names(list.nl)<-names(list.l)<-params

  for (i in params){


    list.l[[i]]<-model.matrix(fixformula[[i]],newdata)

    if (!is.null(object$degree.spline[[i]])){
      list.nl[[i]]<-newdata[,splinevar[[i]]]
    }

  }

  make_prediction(object,list.nl,list.l)


}




#' @title Calculate the centiles
#' @description Using the object generated by \code{\link{gamlssReport}}, the function calculates the centile for a given y and x
#'
#' @param object the object of class \code{"gamlssReport"} generated by \code{\link{gamlssReport}}.
#' @param y the score for which to calculate the centile.
#' @param newdata a dataframe containing the Xs for which to evaluate the centile.
#'
#' @details y can be a vector, in this case newdata need to be a dataframe with single row, or dataframe with the same no of rows as y; if y is a scalar newdata can have as many rows as desired; see examples.
#'
#' @seealso \code{\link{gamlssReport}}, \code{\link{predict.gamlssReport}}
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @export
#' @examples
#' library(gamlss)
#' library(splines)
#'
#' data(aids)
#' aids1<-gamlss(y~pb(x,df=4)+qrt,data=aids,family=PO)
#' m.aids1<-gamlssReport(aids1)
#' centile.gamlssReport(m.aids1,y=c(14,11,20),newdata=data.frame(x=4,qrt=4))
#' centile.gamlssReport(m.aids1,y=14,newdata=data.frame(x=c(4,5),qrt= c(4,2)))
#' centile.gamlssReport(m.aids1,y=c(14,12),newdata=data.frame(x=c(4,5),qrt= c("reference",2) ))


centile.gamlssReport<-function(object,y,newdata){

  pred<-predict.gamlssReport(object,newdata)$fv
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
#'
#' @details centile can be a vector, in this case newdata need to be a dataframe with single row, or dataframe with the same no of rows as centile; if centile is a scalar newdata can have as many rows as desired; see examples.
#'
#' @seealso \code{\link{gamlssReport}}, \code{\link{predict.gamlssReport}}
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @export
#' @examples
#' library(gamlss)
#' library(splines)
#'
#' data(aids)
#' aids1<-gamlss(y~pb(x,df=4)+qrt,data=aids,family=PO)
#' m.aids1<-gamlssReport(aids1)
#' score.gamlssReport(m.aids1,centile=c(0.1,0.5,0.9),newdata=data.frame(x=4,qrt=4))
#' score.gamlssReport(m.aids1,centile=0.1,newdata=data.frame(x=c(4,5),qrt= c(4,2)))
#' score.gamlssReport(m.aids1,centile=c(0.1,0.9),newdata=data.frame(x=c(4,5),qrt= c("reference",2) ))


score.gamlssReport<-function(object,centile,newdata){

  pred<-predict.gamlssReport(object,newdata  )$fv
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
#' @param xnew.spline a named list containing a vector of values for the variable that is modeled as P-spline in pb in a call to gamlss; needs to be of the same length as the number of parameters (mu, sigma, nu, tau), NULL if the parameter does not have a pb term.
#' @param z.new a named list containing fixed effects design matrix for each parameter.
#'
#' @seealso \code{\link{predict.gamlssReport}}
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @export
#' @examples
#' library(gamlss)
#' library(splines)
#'
#' data(aids)
#' aids1<-gamlss(y~pb(x,df=4)+qrt,data=aids,family=PO)
#' m.aids1<-gamlssReport(aids1)
#' make_prediction(m.aids1,list("mu"=c(1,6,45)),
#' list("mu"=rbind(c(1,1,0,0,0),c(1,6,1,0,0),c(1,45,1,0,0))))




make_prediction<-function(object,xnew.spline,z.new){
  params<-object$params


  fv<-lp<-vector("list",length=length(params))
  names(fv)<-names(lp)<-params

  for (i in params){

    cfi<-object$coef.beta[[i]]

    if (!is.null(object$coef.spline[[i]])){
      cfs<-object$coef.spline[[i]]
      D<-make_spline(xnew.spline[[i]],object$knots.spline[[i]],object$degree.spline[[i]],object$range.x[[i]])
    }

    if (!is.null(object$coef.spline[[i]])){
      lp[[i]]<-z.new[[i]]%*%matrix(cfi,ncol=1)+D%*%cfs
    } else {
      lp[[i]]<-z.new[[i]]%*%matrix(cfi,ncol=1)}
    fv[[i]]<-make_inverse(object$link[[i]],lp[[i]])

  }

  list(lp=lp,fv=fv)

}



#' @title Aux extract function
#' @description aux function used by gamlssReport; not intendent to be used by the user.
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

      logic<-grepl("pb",nms) #if TRUE we know that for this mu we have some pb terms, this should go to spline var

      if (sum(logic)!=0) {

        spl<-strsplit(nms[logic],"pb\\(")[[1]][2]

        if (length(strsplit(spl,",")[[1]])==1) pbvar<-strsplit(spl,")")[[1]] else  pbvar<-strsplit( spl  ,",")[[1]][1]
        splinevar[[ii]]<-pbvar
        vars<-nms
        vars[logic]<-pbvar
        fixf<-as.formula(paste("~",paste(vars,collapse="+"),sep=""))

      } else {
        fixf<-as.formula(paste("~",paste(nms,collapse="+"),sep=""))

      }

      fixformula[[ii]]<-fixf

    }
  }

  list(fixformula=fixformula,splinevar=splinevar)
}





#' @title Aux function to recreate the B-spline basis
#' @description aux function used by predict.gamlssReport; not intendent to be used by the user.
#'
#'
#' @param newx vector containing values of x for which the B-spline bases are to be calculated.
#' @param knots the vector of knots; currently only the length(knots) matters (ie we assumed that knots in call to GAMLSS were uniformly set)
#' @param degree the degree of the spline used in pb().
#' @param rng vector with 2 elements, min(x) and max(x)
#'
#' @seealso \code{\link{predict.gamlssReport}}
#' @author Rok Blagus, \email{rok.blagus@@mf.uni-lj.si}
#' @export
#' @examples
#' library(splines)
#'
#' make_spline(10,seq(from=10,to=20,by=1),3,c(9,21))



make_spline<-function(newx,knots,degree,rng){ #this will only work if in pb they used uniformly set knots!

  xr<-rng[2]
  xl<-rng[1]
  ndx<-length(knots)-degree #number of knots-degree!
  xmax<-xr + 0.01 * (xr - xl)
  xmin<-xl - 0.01 * (xr - xl)
  dt<-(xmax - xmin) / ndx
  knots2<-seq(xmin - degree * dt, xmax + degree * dt, by = dt)

  B<-splineDesign(knots = knots2, x = newx, ord = degree + 1, derivs = 0,outer.ok = TRUE)

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




