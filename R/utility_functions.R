#' Reorder a variance-covariance matrix
#'
#' @param old The original variance-covariance matrix
#' @param new_order The desired order in terms of the original columns.
#'   Default is the original order: new_order=1:nrow(old).
#'
#'
#' @return A reordered variance-covariance matrix
#' @export
#'
#' @keywords internal
reorder_vcov <- function(old, new_order=1:nrow(old)) {

  new <- diag(nrow(old))

  for (i in 1:nrow(old)) {
    for (j in 1:nrow(old)) {
      new[i, j] = old[new_order[i], new_order[j]]
    }
  }

  new

}

#' Create pickoff matrix
#'
#' @param contrast Contrast for specifying mean of interest
#' @param n_int Number of intercept parameters
#'
#' @return A matrix that can be used to calculate cumulative probabilities
#' @export
#'
#' @keywords internal
create_pickoff_matrix <- function(contrast, n_int){

  # Contrast matrix
  c_mat <- -contrast%*%(t(rep(1,n_int)))

  # matrix for picking off parameters for contrast 1
  A <- cbind(diag(n_int), t(c_mat))

  A
}

#' Create B matrix
#'
#' @param n_int Number of intercept parameters
#'
#' @return A matrix for turning cumulative probabilities into probabilities
#' @export
#'
#' @keywords internal
create_B_matrix <- function(n_int) {

  B=matrix(0, nrow=n_int+1, ncol=n_int)

  for (i in 1:(n_int+1)) {
    for (j in 1:n_int) {
      if (i == j) {
        B[i, j] = 1
      }
      else if (i == j+1) {
        B[i, j] = -1
      }
    }
  }
  B
}

#' Turn cumulative probabilities to probabilities
#'
#' @param cumprobs Vector of cumulative probabilities
#' @param n_int Number of intercept parameters
#' @param B B matrix
#'
#' @return A vector of probabilities
#' @export
#'
#' @keywords internal
cumprobs_to_probs <- function(cumprobs, n_int, B) {

  # calculate prob of last category
  # as 1 - its cumulative prob
  b <- c(rep(0, n_int), 1)

  probs <- B%*%cumprobs + b

  probs

}

#' Calculate z-statistics
#'
#' @param mean Mean value based on multiplying yvalues by their associated probabilities
#' @param se Standard Error. Either based on delta method or bootstrapping
#'
#' @return Vector of mean, SE, z statistic, and p-value
#' @export
#'
#' @keywords internal
calculate_z_statistics <- function(mean, se) {

  z <- mean/se

  pval=2*stats::pnorm(-abs(z))

  cmat <- cbind(mean, se, z, pval)
  colnames(cmat) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")

  cmat
}

#' Calculate a normal-based confidence interval
#'
#' @param output An object of class "ordinalTools"
#' @param level Confidence interval level
#'
#' @return Upper and lower bounds of a confidence interval
#' @export
#'
#' @examples
#' @keywords internal
calculate_normal_ci <- function(output, level){

a <- (1 - level)/2
a <- c(a, 1 - a)
pct <- format.perc(a, 3)
fac <- stats::qnorm(a)

cint <- array(NA, dim = c(nrow(output), 2L),
              dimnames = list(row.names(output), pct))

cint[] <- output[,1] + output[,2] %o% fac

cint
}


#' Print Coefficient Matrices
#'
#' @param object an object of class "ordinalTools"
#'
#'
#' @return a matrix with columns for the estimated coefficient,
#' its standard error, z-statistic and corresponding (two-sided) p-value.
#' @export
#'
#' @examples
print.ordinalTools <- function(object, ...)
{
  printCoefmat(object$output)
}

#' Summarize results of an object of class "ordinalTools"
#'
#' @param object An object of class "ordinalTools" for which a summary is desired
#' @param digits minimum number of significant digits to be used for most numbers
#' @param signif.stars ogical; if TRUE, P-values are additionally encoded
#' visually as ‘significance stars’ in order to help scanning of long coefficient tables.
#' @param dig.tst minimum number of significant digits for the test statistics
#' @param ...
#'
#' @return Parameter estimates, their SEs, test statistic, and p-value
#' @export
#'
#' @examples
summary.ordinalTools <-
  function(object, digits=max(3, getOption("digits") - 2), signif.stars=TRUE, dig.tst = max(1, min(5, digits - 1)), ...) {
    cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    printCoefmat(object$output, digits=digits, signif.stars=signif.stars, dig.tst=dig.tst, cs.ind=1:2, tst.ind=3, Pvalues=TRUE, has.Pvalue=TRUE)
  }

#' Extract model coefficients
#'
#' @param object an object of class "ordinalTools"
#'
#' @return Coefficients extracted from the model object
#' @export
#'
#' @examples
coef.ordinalTools <-
  function(object,...) {
    object$output[,"Estimate"]
  }

#' Confidence Intervals for Means
#'
#' @param object a object of type ordinalTools
#'
#' @return A matrix (or vector) with columns giving
#' lower and upper confidence limits for each parameter.
#' These will be labelled as (1-level)/2 and 1 - (1-level)/2 in %
#' @export
#'
#' @examples
confint.ordinalTools <- function(object) {

  ci <- object$conf.int

  ci

}

#' Format a probability as a percentage
#'
#' @param probs Probabilty value
#' @param digits The number of digits
#'
#' @return Formatted percentages
#' @export
#'
#' @keywords internal
format.perc <- function (probs, digits) {

    paste(format(100 * probs, trim = TRUE,
                 scientific = FALSE, digits = digits), "%")
  }




