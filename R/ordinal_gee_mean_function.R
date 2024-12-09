

#' Calculate means from an ordinal GEE model
#'
#' @param object An object of class "LORgee"
#' @param data A data frame.
#' @param conf.level Confidence interval level
#' @param contrast1 A vector the same length as the number of coefficients in
#   the model. Allows the user to calculate means based on covariate
#   values.
#' @param contrast2 A vector the same length as the number of coefficients in
#   the model. Allows the user to calculate means based on covariate
#   values.
#' @param se.type Method for calculating standard errors. Options are
#' the delta method (se.type="delta") or bootstrapping (se.type="bootstrap")
#' @param R The number of bootstrap samples
#'
#' @return An object of class "ordinalTools. This has components
#'
#' * call the matched call
#' * output Parameter estimates, their standard errors, test statistics, and p-values
#' * conf.int Confidence interval
#' * se.type Method used to calculate standard errors
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- ordinal_means(multgee.fit)
#' summary(fit)
#'
#' # Two coefficient constrasts and their difference
#' fit <- ordinal_means(multgee.fit, contrast1=c(0,0,0), contrast2=c(0,0,0))
#' summary(fit)
#'
#' # Two coefficient constrasts and their difference
#' # Bootstrap standard errors
#' fit <- ordinal_means(multgee.fit, data=data, contrast1=c(0,0,0),
#'        contrast2=c(0,0,0), type="bootstrap")
#' summary(fit)
#' }
#' @details This function uses the parameter estimates from an proportional odds
#' GEE logistic regression model that was fit using the multgee::ordLORgee function to obtain
#' mean values for a given set of coefficient contrasts. Up to two constrasts
#' can be specified. Standard errors are calculated either using the delta
#' method or bootstrapping.
ordinal_gee_means <- function(object, data=data, conf.level = 0.95,
                          contrast1=rep(0, length(coef(object))),
                          contrast2=NULL, se.type="delta", R=500) {

  cl <- match.call()

  # Error handling
  if (class(object) != "LORgee") {
    stop("To use the ordinal_gee_means function, the object must be
         of class LORgee from the multgee package.")
  }

  # Error handling
  if (se.type=="bootstrap") {
    stop("Bootstrap standard errors have not yet been implemented. Please use
          delta method standard errors.")
  }

  # identify the outcome variable
  formula <- formula(object$terms)
  outcome <- all.vars(formula)[1]

  # Number of covariates
  n_cov <- length(all.vars(formula)) - 1

  # Error handling
  if (length(contrast1) != n_cov){
    stop("Length of contrast 1 vector should be equal to ",
         n_cov, " which is the length of the coefficient vector.")
  }

  if (!is.null(contrast2) & length(contrast2) != n_cov){
    stop("Length of contrast 1 vector should be equal to ",
         n_cov, " which is the length of the coefficient vector.")
  }

  # Unique values of the outcome
  # CHANGE SO NAs not counted
  yvals <- sort(as.numeric(unique(data[, outcome])))


  # Number of intercepts
  n_int <- object$categories - 1

  # Full parameter vector
  theta <- object$coefficients


  # B matrix for turning cumulative
  # probabilities to probabilities
  B <- create_B_matrix(n_int)

  if (se.type=="delta" & is.null(contrast2)) {

    c1 <- get_mean_and_se(object, contrast1, n_int=n_int, theta=theta, B=B, yvals=yvals)
    output <- c1
    row.names(output) <- c("Contrast 1")

    cint <- calculate_normal_ci(output, level=conf.level)


  }

  else if (se.type=="bootstrap" & is.null(contrast2)) {
    boot.means <- boot::boot(data=data, statistic=bootstrap_mean,
                             R=R, contrast1=contrast1, object=object)

    mean <- boot.means$t0
    se <- stats::sd(boot.means$t)
    c1 <- calculate_z_statistics(mean, se)
    output <- c1
    row.names(output) <- c("Contrast 1")

    ci1 <- boot::boot.ci(boot.means, conf=conf.level, index=1, type="perc")

    cint <- cbind(ci1$percent[4], ci1$percent[5])
    row.names(cint) <- c("Contrast 1")

    a <- (1 - conf.level)/2
    a <- c(a, 1 - a)
    pct <- format_perc(a, 3)

    dimnames(cint) = list(row.names(output), pct)
  }

  else if (se.type=="delta" & !is.null(contrast2)){
    c1 <- get_mean_and_se(object, contrast1, n_int=n_int, theta=theta, B=B, yvals=yvals)
    c2  <- get_mean_and_se(object, contrast2, n_int=n_int, theta=theta, B=B, yvals=yvals)
    diff <- calculate_difference_and_se(object, contrast1, contrast2, n_int=n_int, theta=theta, B, yvals=yvals)
    output <- rbind(c1, c2, diff)
    row.names(output) <- c("Contrast 1", "Contrast 2", "Difference")
    cint <- calculate_normal_ci(output, level=conf.level)

  }

  else if (se.type=="bootstrap" & !is.null(contrast2)) {
    boot.means <- boot::boot(data=data, statistic=bootstrap_mean,
                             R=R, contrast=contrast1, contrast2=contrast2, object=object)
    mean1 <- boot.means$t0[1]
    mean2 <- boot.means$t0[2]
    diff  <- boot.means$t0[3]

    se1   <-  stats::sd(boot.means$t[,1])
    se2   <-  stats::sd(boot.means$t[,2])
    sediff <- stats::sd(boot.means$t[,3])

    c1 <- calculate_z_statistics(mean1, se1)
    c2 <- calculate_z_statistics(mean2, se2)
    diff <- calculate_z_statistics(diff, sediff)
    output <- rbind(c1, c2, diff)
    row.names(output) <- c("Contrast 1", "Contrast 2", "Difference")

    ci1 <- boot::boot.ci(boot.means, conf=conf.level, index=1, type="perc")
    ci2 <- boot::boot.ci(boot.means, conf=conf.level, index=2, type="perc")
    ci3 <- boot::boot.ci(boot.means, conf=conf.level, index=3, type="perc")

    cint1 <- cbind(ci1$percent[4], ci1$percent[5])
    cint2 <- cbind(ci2$percent[4], ci2$percent[5])
    cint3 <- cbind(ci3$percent[4], ci3$percent[5])

    cint <- rbind(cint1, cint2, cint3)
    row.names(cint) <- c("Contrast 1", "Contrast 2", "Difference")

    a <- (1 - conf.level)/2
    a <- c(a, 1 - a)
    pct <- format_perc(a, 3)

    dimnames(cint) = list(row.names(output), pct)
  }

  coefficients <- object$output[,"Estimate"]

  # Unique values of the outcome
  yvals <- as.numeric(object$lev)

  # Cumulative probabilities
  cumprobs <- calculate_cumprobs(contrast1, n_int, theta)

  if (!is.null(contrast2)) {

    cumprobs2 <- calculate_cumprobs(contrast2, n_int, theta)

    cumprobs <- cbind(cumprobs, cumprobs2)

  }

  object <- list(call = cl, output = output, se.type = se.type,
                 conf.int = cint, coefficients=coefficients,
                 yvals = yvals, cumprobs = cumprobs)

  class(object) <- "ordinalTools"

  object

}




