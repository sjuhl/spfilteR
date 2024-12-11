#' @rdname lmFilter
#' @importFrom stats filter
#' @export

summary.spfilter <- function(object, EV = FALSE, ...) {
  #####
  # Print Output
  #####
  # head
  mod <- ifelse(object$other$model == "nb", "NegBin", object$other$model)
  cat("\n\t- Spatial Filtering with Eigenvectors"
      ,paste0("(", paste0(toupper(substring(mod, 1, 1)),
                          substring(mod, 2)), " Model)")," -\n\n")

  # estimates & model fit
  signif <- star(p = object$estimates[, "p-value"])
  estimates <- data.frame(object$estimates, signif)
  colnames(estimates) <- c(colnames(object$estimates), "")
  cat("Coefficients",paste0(ifelse(object$other$model == "linear"
                                   & !("condnum" %in% names(object$other)),
                                   "(OLS)", "(MLE)"), ":\n"))
  print(estimates)
  if (object$other$model == "linear") {
    cat("\nAdjusted R-squared:\n")
    print(object$fit)
  } else {
    cat("\nModel Fit:\n")
    print(object$fit)
  }

  # additional information on stepwise regression
  cat(paste("\nFiltered for", object$other$dependence, "spatial autocorrelation\n"))
  cat(paste(object$other$nev, "out of", object$other$ncandidates,
            "candidate eigenvectors selected\n"))
  if (object$other$model != "linear" & object$other$nev > 0) {
    cat(paste0("Condition Number (Multicollinearity): ", object$other$condnum, "\n"))
  }
  cat(paste0("Objective Function: \"", object$other$objfn, "\""))
  if (object$other$objfn == "p") {
    if (object$other$bonferroni) {
      cat(paste0("\ (significance level = ", round(object$other$siglevel * object$other$ncandidates, 5),
                 ")\n"))
      cat(paste0("Bonferroni correction: ", object$other$bonferroni, ""))
      cat(paste0("\ (adjusted significance level = ", round(object$other$siglevel, 5), ")\n"))
    } else {
      cat(paste0("\ (significance level = ", round(object$other$siglevel, 5), ")\n"))
      cat(paste0("Bonferroni correction: ", object$other$bonferroni, "\n"))
    }
  } else {
    cat("\n")
  }

  # optional: information on eigenvectors
  if (EV) {
    if (object$other$nev == 0) {
      cat("\nNo eigenvectors selected\n")
    } else {
      sigev <- star(p = object$EV[, "p-value"])
      EV <- data.frame(object$EV, sigev)
      colnames(EV) <- c(colnames(object$EV), "")
      cat("\nSummary of selected eigenvectors:\n")
      print(EV)
    }
  }

  # Moran's I
  m_signif <- star(p = object$moran[, "p-value"])
  moran <- data.frame(object$moran, m_signif)
  colnames(moran) <- c(colnames(object$moran), "")
  cat(paste0("\n","Moran's I ", ifelse(object$other$model != "linear",
                                paste0("(", toupper(substring(object$other$resid.type, 1, 1)),
                                       substring(object$other$resid.type, 2), ""), "("),
             " Residuals):\n"))
  print(moran)
}


#' @export
print.spfilter <- function(x, ...) {
  cat(paste(x$other$nev, "out of", x$other$ncandidates, "candidate eigenvectors selected\n"))
}


#' @export
coef.spfilter <- function(object, ...) {
  object$estimates[, "Estimate"]
}

#' @export
vcov.spfilter <- function(object, ...) {
  object$varcovar
}


#' @importFrom graphics plot legend polygon abline points
#' @importFrom grDevices rgb
#' @export

plot.spfilter <- function(x, ...) {
  plot(0, ylim = c(min(x$evMI), max(x$evMI)), xlim = c(1, length(x$evMI)),
       main = "Moran Coefficients for\n all Eigenvectors",
       ylab = "Moran Coefficient", xlab = "Eigenvector", type = "n",las = 1, ...)
  # area of candidate set
  xstart <- ifelse(x$other$dependence == "positive", -100, length(x$evMI) - x$other$ncandidates)
  xend <- ifelse(x$other$dependence == "positive", x$other$ncandidates, length(x$evMI) * 2)
  polygon(x = c(xstart, xend, xend, xstart),
          y = c(min(x$evMI) - 1, min(x$evMI) - 1,
                max(x$evMI) + 1, max(x$evMI) + 1),
          col = rgb(red = 0, green = 0, blue = 0, alpha = 0.1),
          border = FALSE)
  # not selected EVs
  points(y = x$evMI[which(!(seq_along(x$evMI) %in% x$other$sel_id))],
         x = which(!(seq_along(x$evMI) %in% x$other$sel_id)), pch = 16, cex = .4, col = "gray")
  # selected EVs
  points(y = x$evMI[x$other$sel_id], x = x$other$sel_id, pch = 16, cex = .7)
  # legend
  legend("topright", legend = c("selected", "other"), pch = 16, col = c("black", "gray"), cex = .8)
  abline(h = 0, lty = 2, cex = .5)
}
