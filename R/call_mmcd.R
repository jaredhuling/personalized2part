
#' @export
mmbcd <- function(x, y,
                  groups = NULL,
                  group_weights = NULL,
                  weights = rep(1, NROW(x)),
                  penalty = c("grp.lasso", "coop.lasso"),
                  family = c("gaussian", "binomial", "gamma"),
                  lambda = NULL,
                  nlambda = 100L,
                  lambda_min_ratio = 1e-3,
                  maxit = 500,
                  tol = 1e-5,
                  intercept = TRUE)
{
    p <- NCOL(x)
    n <- NROW(x)

    penalty <- match.arg(penalty)
    family  <- match.arg(family)

    if (is.null(groups))
    {
        groups <- 1:p
    }

    if (is.null(group_weights))
    {
        group_weights <- numeric(0)
    }
    if (is.null(lambda))
    {
        lambda <- numeric(0)
    }

    y <- setup_y(y, family)

    groups <- as.integer(groups)

    unique_groups <- unique(groups)

    group_weights <- as.double(group_weights)
    weights <- as.double(weights)
    lambda <- as.double(lambda)
    nlambda <- as.integer(nlambda)
    maxit <- as.integer(maxit)
    tol <- as.double(tol)
    lambda_min_ratio <- as.double(lambda_min_ratio)
    intercept <- as.logical(intercept)

    res <- mmbcd_cpp(X = x, Y = y, groups = groups,
                     unique_groups = unique_groups,
                     group_weights = group_weights,
                     weights = weights, lambda = lambda, nlambda = nlambda,
                     lambda_min_ratio = lambda_min_ratio,
                     maxit = maxit, tol = tol, intercept = intercept,
                     penalty = penalty, family = family)

    res
}


