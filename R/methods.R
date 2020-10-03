#' @rdname lmFilter
#' @importFrom stats filter

summary.spfilter <- function(obj,EV=FALSE){
  #####
  # Print Output
  #####
  # head
  cat("\n\t- Spatial Filtering with Eigenvectors"
      ,paste0("("
              ,paste0(toupper(substring(obj$other$model,1,1))
                      ,substring(obj$other$model,2)
              )
              ," Model)")," -\n\n")

  # estimates & model fit
  signif <- star(p=obj$estimates[,"p-value"])
  estimates <- data.frame(obj$estimates,signif)
  colnames(estimates) <- c(colnames(obj$estimates),"")
  cat("Coefficients",paste0(ifelse(obj$other$model=="linear"
                                   & "condnum" %in% names(obj$other)
                                   ,"(OLS)","(ML)"),":\n"))
  print(estimates)
  if(obj$other$model=="linear"){
    cat("\nAdjusted R-squared:\n")
    print(obj$fit)
  } else {
    cat("\nModel Fit:\n")
    print(obj$fit)
  }

  # additional information on stepwise regression
  cat(paste("\nFiltered for", obj$other$dependence, "spatial autocorrelation"))
  cat(paste(obj$other$nev,"out of",obj$other$ncandidates, "candidate eigenvectors selected\n"))
  if(obj$other$model!="linear" & obj$other$nev>0){
    cat(paste0("Condition Number (Multicollinearity): ",obj$other$condnum,"\n"))
  }
  cat(paste0("Objective Function: \"" ,obj$other$objfn,"\""))
  if(obj$other$objfn=="p"){
    if(obj$other$bonferroni){
      cat(paste0("\ (significance level=",round(obj$other$siglevel*obj$other$ncandidates,5),")\n"))
      cat(paste0("Bonferroni correction: " ,obj$other$bonferroni,""))
      cat(paste0("\ (adjusted significance level=",round(obj$other$siglevel,5),")\n"))
    } else {
      cat(paste0("\ (significance level=",round(obj$other$siglevel,5),")\n"))
      cat(paste0("Bonferroni correction: ",obj$other$bonferroni,"\n"))
    }
  } else cat("\n")

  # optional: information on eigenvectors
  if(EV){
    if(obj$other$nev==0){
      cat("\nNo eigenvectors selected\n")
    } else {
      sigev <- star(p=obj$EV[,"p-value"])
      EV <- data.frame(obj$EV,sigev)
      colnames(EV) <- c(colnames(obj$EV),"")
      cat("\nSummary of selected eigenvectors:\n")
      print(EV)
    }
  }

  # Moran's I
  m_signif <- star(p=obj$moran[,"p-value"])
  moran <- data.frame(obj$moran,m_signif)
  colnames(moran) <- c(colnames(obj$moran),"")
  cat(paste0("\n","Moran's I (",ifelse(obj$other$model!="linear"
                                ,paste0(toupper(substr(obj$other$resid.type,1,1))
                                        ,substr(obj$other$resid.type,2)
                                        ,""),"")
             ," Residuals):\n"))
  print(moran)
}


# print function
print.spfilter <- function(obj){
  cat(paste(obj$other$nev,"out of",obj$other$ncandidates, "candidate eigenvectors selected"))
}


# coef function
coef.spfilter <- function(obj){
  obj$estimates[,"Estimate"]
}

# vcov function
vcov.spfilter <- function(obj){
  obj$varcovar
}

#' @rdname lmFilter
#' @importFrom graphics plot legend polygon abline points
#' @importFrom grDevices rgb

plot.spfilter <- function(obj){
  plot(0,ylim=c(min(obj$evMI),max(obj$evMI)),xlim=c(1,length(obj$evMI))
       ,main="Moran Coefficients for\n all Eigenvectors"
       ,ylab="Moran Coefficient",xlab="Eigenvector",type="n",las=1)
  # area of candidate set
  xstart <- ifelse(obj$other$dependence=="positive",-100,length(obj$evMI)-obj$other$ncandidates)
  xend <- ifelse(obj$other$dependence=="positive",obj$other$ncandidates,length(obj$evMI)*2)
  polygon(x=c(xstart,xend,xend,xstart)
          ,y=c(min(obj$evMI)-1,min(obj$evMI)-1
               ,max(obj$evMI)+1,max(obj$evMI)+1)
          ,col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1)
          ,border=F)
  # not selected EVs
  points(y=obj$evMI[which(!(1:length(obj$evMI) %in% obj$other$sel_id))]
         ,x=which(!(1:length(obj$evMI) %in% obj$other$sel_id)),pch=16,cex=.4,col="gray")
  # selected EVs
  points(y=obj$evMI[obj$other$sel_id],x=obj$other$sel_id,pch=16,cex=.7)
  # legend
  legend("topright",legend=c("selected","other"),pch=16,col=c("black","gray"),cex=.8)
  abline(h=0,lty=2,cex=.5)
}
