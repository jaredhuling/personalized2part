
#' Main fitting function for group lasso and cooperative lasso penalized two part models
#'
#' @description This function fits penalized two part models with a logistic regression model
#' for the zero part and a gamma regression model for the positive part. Each covariate's effect
#' has either a group lasso or cooperative lasso penalty for its effects for the two consituent
#' models
#'
#' @param x an n x p matrix of covariates for the zero part data, where each row is an observation
#' and each column is a predictor
#' @param z a length n vector of responses taking values 1 and 0, where 1 indicates the response is positive
#' and zero indicates the response has value 0.
#' @param x_s an n_s x p matrix of covariates (which is a submatrix of x) for the positive part data, where
#' each row is an observation and each column is a predictor
#' @param s a length n_s vector of responses taking strictly positive values
#' @param weights a length n vector of observation weights for the zero part data
#' @param weights_s a length n_s vector of observation weights for the positive part data
#' @param offset a length n vector of offset terms for the zero part data
#' @param offset_s a length n_s vector of offset terms for the positive part data
#' @param penalty either \code{"grp.lasso"} for the group lasso penalty or \code{"coop.lasso"} for the
#' cooperative lasso penalty
#' @param penalty_factor a length p vector of penalty adjustment factors corresponding to each covariate.
#' A value of 0 in the jth location indicates no penalization on the jth variable, and any positive value will
#' indicate a multiplicative factor on top of the common penalization amount. The default value is 1 for
#' all variables
#' @param algorithm either \code{"irls"} for an iteratively reweighted least squares type algorithm or
#' \code{"mm"} for a majorization minimization type algorithm. \code{"irls"} is recommended and is the default
#' @param nlambda the number of lambda values. The default is 100.
#' @param lambda_min_ratio Smallest value for \code{lambda}, as a fraction of lambda.max, the data-derived largest lambda value
#' The default depends on the sample size relative to the number of variables.
#' @param lambda a user supplied sequence of penalization tuning parameters. By default, the program automatically
#' chooses a sequence of lambda values based on \code{nlambda} and \code{lambda_min_ratio}
#' @param intercept whether or not to include an intercept in the model. Default is \code{TRUE}.
#' @param maxit_irls maximum number of IRLS iterations
#' @param tol_irls convergence tolerance for IRLS iterations
#' @param maxit_mm maximum number of MM iterations. Note that for \code{algorithm = "irls"}, MM is used within
#' each IRLS iteration, so \code{maxit_mm} applies to the convergence of the inner iterations in this case.
#' @param tol_mm convergence tolerance for MM iterations. Note that for \code{algorithm = "irls"}, MM is used within
#' each IRLS iteration, so \code{tol_mm} applies to the convergence of the inner iterations in this case.
#' @export
#'
#' @examples
#'
#' library(personalized2part)
#'
hd2part <- function(x, z,
                    x_s, s,
                    weights          = rep(1, NROW(x)),
                    weights_s        = rep(1, NROW(x_s)),
                    offset           = NULL,
                    offset_s         = NULL,
                    penalty          = c("grp.lasso", "coop.lasso"),
                    penalty_factor   = NULL,
                    algorithm        = c("irls", "mm"),
                    nlambda          = 100L,
                    lambda_min_ratio = ifelse(n_s < p, 0.05, 0.005),
                    lambda           = NULL,
                    intercept        = TRUE,
                    maxit_irls       = 50,
                    tol_irls         = 1e-5,
                    maxit_mm         = 500,
                    tol_mm           = 1e-5)
{
    p <- NCOL(x)
    n <- NROW(x)
    n_s <- NROW(x_s)

    if (p != ncol(x_s))
    {
        stop("'x' and 'x_s' must have the same number of columns")
    }

    if (length(weights) != n)
    {
        stop("'weights' must be same length as number of observations in 'x'")
    }

    if (length(weights_s) != n_s)
    {
        stop("'weights_s' must be same length as number of observations in 'x_s'")
    }

    is.offset   <- !is.null(offset)
    is.offset.s <- !is.null(offset_s)

    if (is.offset)
    {
        offset <- drop(as.double(offset))
    } else
    {
        offset <- rep(0, n)
    }

    if (is.offset.s)
    {
        offset_s <- drop(as.double(offset_s))
    } else
    {
        offset_s <- rep(0, n_s)
    }

    if (length(offset) != n)
    {
        stop("'offset' must be same length as number of observations in 'x'")
    }

    if (length(offset_s) != n_s)
    {
        stop("'offset_s' must be same length as number of observations in 'x_s'")
    }



    algorithm <- match.arg(algorithm)
    penalty   <- match.arg(penalty)

    if (is.null(penalty_factor))
    {
        penalty_factor <- numeric(0)
    }
    if (is.null(lambda))
    {
        lambda <- numeric(0)
    }

    if (length(penalty_factor) > 0)
    {
        if (length(penalty_factor) != p)
        {
            stop("'penalty_factor' must be same length as the number of observations")
        }
    }

    ## run checks on outcomes
    z <- setup_y(z, "binomial")
    s <- setup_y(s, "gamma")

    groups           <- as.integer(rep(1:NCOL(x), 2))
    unique_groups    <- unique(groups)

    penalty_factor   <- as.double(penalty_factor)
    weights          <- as.double(weights)
    weights_s        <- as.double(weights_s)
    lambda           <- as.double(lambda)
    nlambda          <- as.integer(nlambda[1])
    maxit_mm         <- as.integer(maxit_mm[1])
    tol_mm           <- as.double(tol_mm[1])
    maxit_irls       <- as.integer(maxit_irls[1])
    tol_irls         <- as.double(tol_irls[1])
    lambda_min_ratio <- as.double(lambda_min_ratio[1])
    intercept        <- as.logical(intercept[1])

    if (nlambda <= 0)     stop("'nlambda' must be a positive integer")
    if (maxit_mm <= 0)    stop("'maxit_mm' must be a positive integer")
    if (maxit_irls <= 0)  stop("'maxit_irls' must be a positive integer")
    if (tol_irls <= 0)    stop("'tol_irls' must be a strictly positive number")
    if (tol_mm <= 0)      stop("'tol_mm' must be a strictly positive number")

    if (length(lambda) > 0 && any(lambda) <= 0) stop("every value in 'lambda' must be a strictly positive number")

    if (lambda_min_ratio <= 0 || lambda_min_ratio >= 1) stop("'lambda_min_ratio' must be strictly positive and less than 1")

    lambda <- rev(sort(lambda))

    if (algorithm == "mm")
    {
        res <- mmbcd_twopart_cpp(X = x, Z = z,
                                 Xs = x_s, S = s,
                                 groups = groups,
                                 unique_groups = unique_groups,
                                 group_weights = penalty_factor,
                                 weights = weights,
                                 weights_s = weights_s,
                                 lambda = lambda, nlambda = nlambda,
                                 lambda_min_ratio = lambda_min_ratio,
                                 maxit = maxit_mm, tol = tol_mm,
                                 intercept = intercept,
                                 penalty = penalty)
    } else
    {
        res <- irls_mmbcd_twopart_cpp(X = x, Z = z,
                                      Xs = x_s, S = s,
                                      groups = groups,
                                      unique_groups = unique_groups,
                                      group_weights = penalty_factor,
                                      weights = weights,
                                      weights_s = weights_s,
                                      lambda = lambda, nlambda = nlambda,
                                      lambda_min_ratio = lambda_min_ratio,
                                      maxit = maxit_mm, tol = tol_mm,
                                      maxit_irls = maxit_irls,
                                      tol_irls = tol_irls,
                                      intercept = intercept,
                                      penalty = penalty)
    }

    res$offset   <- is.offset
    res$offset_s <- is.offset.s

    class(res)   <- "hd2part"
    res
}
