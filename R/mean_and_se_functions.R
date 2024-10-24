
#' Calculate mean
#'
#' @param object An object of class "polr"
#' @param contrast Contrast vector for calculating a mean based on
#' covariate values
#' @param n_int Number of intercept parameters
#' @param theta Vector of intercept parameters and regression coefficients
#' @param B Matrix for converting cumulative probabilities to probabilities
#'
#' @return The mean associated with the specified contrast
#' @export
#'
#' @keywords internal
calculate_mean <- function(object, contrast, n_int, theta, B) {

  A <- create_pickoff_matrix(contrast, n_int)

  cumprobs <- 1/(1 + exp(-A%*%theta))

  probs <- cumprobs_to_probs(cumprobs, n_int, B)

  # Unique values of the outcome
  yvals <- as.numeric(object$lev)

  # Return mean
  mean <- yvals%*%probs

  mean
}

#' Bootstrap mean and difference in means
#'
#' This function is called when bootstrapping means and differences in means
#'
#' @param data the original data frame
#' @param i indices
#' @param object An object of class "polr"
#' @param contrast1 First contrast vector for calculating a mean based on
#' covariate values
#' @param contrast2 Second contrast vector for calculating a mean based on
#' covariate values
#'
#' @return A mean associated with contrast 1 or two means and their difference
#' based on constrasts 1 and 2
#' @export
#'
#' @keywords internal
bootstrap_mean <- function(data=data, i, object, contrast1, contrast2=NULL) {

  d2 <- data[i,]

  # Fit ordinal logistic regression model
  bootfit <- MASS::polr(formula=formula(object$terms), data=d2)

  # number of intercepts will change in each bootstrap sample

  # Extract betas only
  betas <- as.vector(coef(bootfit))

  # Number of covariates
  n_cov <- length(bootfit)

  # Extract intercepts
  zetas <- as.vector(bootfit$zeta)

  # Number of intercepts
  n_int <- length(zetas)

  # Full parameter vector
  theta <- c(zetas, betas)

  # Number of parameters
  n_params <- n_cov + n_int

  # B matrix for turning cumulative
  # probabilities to probabilities
  B <- create_B_matrix(n_int)

  mean1 <- calculate_mean(bootfit, contrast1, n_int, theta, B)

  if (!is.null(contrast2)) {
    mean2 <- calculate_mean(bootfit, contrast2, n_int, theta, B)
    diff <- mean1 - mean2
    output <- c(mean1, mean2, diff)
  }
  else{
    output <- mean1
  }

  return(output)

}

#' Calculate Delta Method Standard Error
#'
#' @param object An object of class "polr"
#' @param contrast Contrast vector for calculating a mean based on
#' covariate values
#' @param n_int Number of intercept parameters
#' @param theta Vector of intercepts and regression coefficients
#' @param B Matrix for converting cumulative probabilities to probabilities
#'
#' @return The standard error of mean associated with the specified contrast
#' vector
#' @export
#'
#' @keywords internal
calculate_delta_method_se <- function(object, contrast, n_int, theta, B) {

  A <- create_pickoff_matrix(contrast, n_int)

  cumprobs <- 1/(1 + exp(-A%*%theta))

  # Unique values of the outcome
  yvals <- as.numeric(object$lev)

  ## Standard error of mean via delta method ##

  # Calculate derivative

  middle <- diag(length(cumprobs))

  diag(middle) <- diag(cumprobs%*%t((1-cumprobs)))

  deriv <- yvals%*%B%*%middle%*%A

  vcov <- vcov(object)

  # variance covariance matrix has coefficients first
  # reorder so intercepts first
  n_cov <- length(theta) - n_int

  order <- c((n_cov+1):nrow(vcov),1:n_cov)

  new_vcov <- reorder_vcov(vcov, new_order=order)

  se <- sqrt(deriv%*%new_vcov%*%t(deriv))

  se
}

#' Calculate the difference of two means and their standard error
#'
#' @param object An object of class "polr"
#' @param contrast1 First contrast vector for calculating a mean based on
#' covariate values
#' @param contrast2 Second contrast vector for calculating a mean based on
#' covariate values
#' @param n_int Number of intercept parameters
#' @param theta A vector of intercept parameters and regression coefficients
#' @param B Matrix for converting cumulative probabilities to probabilities
#'
#' @return A vector consisting of the difference in means, its SE, Z-statistic,
#' and p-value
#' @export
#'
#' @keywords internal
calculate_difference_and_se <- function(object,
                                        contrast1,
                                        contrast2,
                                        n_int,
                                        theta,
                                        B) {
  # Return mean
  mean1 <- calculate_mean(object, contrast1, n_int, theta, B)
  mean2 <- calculate_mean(object, contrast2, n_int, theta, B)
  diff <- mean1 - mean2

  ## Standard error of mean via delta method ##

  # Calculate derivative

  # matrix for picking off parameters for contrast 1
  A1 <- create_pickoff_matrix(contrast1, n_int)

  # Contrast 1 cumulative probabilities
  cumprobs1 <- 1/(1 + exp(-A1%*%theta))

  # matrix for picking off parameters for contrast 2
  A2 <- create_pickoff_matrix(contrast2, n_int)

  # Contrast 2 cumulative probabilities
  cumprobs2 <- 1/(1 + exp(-A2%*%theta))

  # Unique values of the outcome
  yvals <- as.numeric(object$lev)

  middle1 <- diag(length(cumprobs1))

  diag(middle1) <- diag(cumprobs1%*%t((1-cumprobs1)))

  middle2 <- diag(length(cumprobs2))

  diag(middle2) <- diag(cumprobs2%*%t((1-cumprobs2)))

  middle_diff <- middle1%*%A1 - middle2%*%A2

  deriv <- yvals%*%B%*%middle_diff

  vcov <- vcov(object)

  # variance covariance matrix has coefficients first
  # reorder so intercepts first
  n_cov <- length(theta) - n_int

  order <- c((n_cov+1):nrow(vcov),1:n_cov)

  new_vcov <- reorder_vcov(vcov, new_order=order)

  se <- sqrt(deriv%*%new_vcov%*%t(deriv))

  cmat <- calculate_z_statistics(diff, se)

  cmat

}

#' Get mean and standard error from an ordinal model via
#' the delta method
#'
#' This function is a wrapper function that uses an object of class "polr" and a specified contrast
#' vector to calculate a mean, its SE, and Z-statistics
#'
#' @param object An object of class ``polr"
#' @param contrast Contrast vector for calculating a mean based on
#' covariate values
#' @param n_int Number of intercept parameters
#' @param theta Vector of intercept parameters and regression coefficients
#'
#' @return A data frame consisting of a mean, its SE, z-statistic, and p-value
#' @export
#'
#' @keywords internal
get_mean_and_se <- function(object, contrast, n_int, theta, B) {

  # Calculate mean
  mean <- calculate_mean(object, contrast, n_int, theta, B)

  ## Standard error of mean via delta method ##
  se <- calculate_delta_method_se(object, contrast, n_int, theta, B)

  cmat <- calculate_z_statistics(mean, se)

  cmat


}
