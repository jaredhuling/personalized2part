

fit_subgroup_2part <- function(x,
                               y,
                               trt,
                               propensity_func = NULL,
                               match_id = NULL,
                               penalty          = c("grp.lasso", "coop.lasso"),
                               augment_func_zero = NULL,
                               augment_func_positive = NULL,
                               cutpoint              = 1,
                               larger_outcome_better = TRUE,
                               y_eps = 1e-6,
                               ...)
{
    penalty <- match.arg(penalty)
    dims    <- dim(x)
    if (is.null(dims)) stop("x must be a matrix object.")

    y       <- drop(y)
    vnames  <- colnames(x)

    p   <- NCOL(x)
    n   <- NROW(x)

    if (is.null(vnames))
    {
        vnames <- paste0("V", 1:p)
    }

    # check to make sure arguments of augment_func are correct
    if (!is.null(augment_func_zero))
    {
        augmentfunc.names <- sort(names(formals(augment_func_zero)))
        if (length(augmentfunc.names) == 3)
        {
            if (any(augmentfunc.names != c("trt", "x", "y")))
            {
                stop("arguments of augment.func() should be 'trt', 'x', and 'y'")
            }
        } else if (length(augmentfunc.names) == 2)
        {
            if (any(augmentfunc.names != c("x", "y")))
            {
                stop("arguments of augment.func() should be 'x' and 'y'")
            }
            augment.func2 <- augment_func_zero
            augment_func_zero  <- function(trt, x, y) augment.func2(x = x, y = y)
        } else
        {
            stop("augment_func_zero() should only have either two arguments: 'x' and 'y', or three arguments:
                 'trt', 'x', and 'y'")
        }
    }
    if (!is.null(augment_func_positive))
    {
        augmentfunc.names <- sort(names(formals(augment_func_positive)))
        if (length(augmentfunc.names) == 3)
        {
            if (any(augmentfunc.names != c("trt", "x", "y")))
            {
                stop("arguments of augment.func() should be 'trt', 'x', and 'y'")
            }
        } else if (length(augmentfunc.names) == 2)
        {
            if (any(augmentfunc.names != c("x", "y")))
            {
                stop("arguments of augment.func() should be 'x' and 'y'")
            }
            augment.func2 <- augment_func_positive
            augment_func_positive  <- function(trt, x, y) augment.func2(x = x, y = y)
        } else
        {
            stop("augment_func_positive() should only have either two arguments: 'x' and 'y', or three arguments:
                 'trt', 'x', and 'y'")
        }
    }

    if (is.factor(trt))
    {
        # drop any unused levels of trt
        trt         <- droplevels(trt)
        unique.trts <- levels(trt)
        n.trts      <- length(unique.trts)
        trt         <- ifelse(trt == unique.trts[1], 0, 1)
    } else
    {
        unique.trts <- sort(unique(trt))
        n.trts      <- length(unique.trts)
        trt         <- ifelse(trt == unique.trts[1], 0, 1)
    }

    if (n.trts > 2)
    {
        stop("multiple treatments (>2) not currently available.")
    }

    if (n.trts < 2)           stop("trt must have at least 2 distinct levels")
    if (n.trts > dims[1] / 3) stop("trt must have no more than n.obs / 3 distinct levels")


    # defaults to constant propensity score within trt levels
    # the user will almost certainly want to change this
    if (is.null(propensity_func))
    {
        if (is.null(match_id))
        { # No propensity score supplied and no match.id supplied
            if (n.trts == 2)
            {
                mean.trt <- mean(trt == 1)
                propensity_func <- function(trt, x) rep(mean.trt, length(trt))
            } else
            {
                mean.trt <- numeric(n.trts)
                for (t in 1:n.trts)
                {
                    mean.trt[t] <- mean(trt == unique.trts[t])
                }
                propensity_func <- function(trt, x)
                {
                    pi.x <- numeric(length(trt))
                    for (t in 1:n.trts)
                    {
                        which.t       <- trt == unique.trts[t]
                        pi.x[which.t] <- mean(which.t)
                    }
                    pi.x
                }
            }
        } else
        { # No propensity score supplied but match.id supplied
            if (n.trts == 2)
            {
                # default to pct in treatment group
                mean.trt <- mean(trt == 1)
                pf <- function(trt, x, match.id) rep(mean.trt, NROW(trt))
                propensity_func <- pf
            } else
            {
                mean.trt <- numeric(n.trts)
                for (t in 1:n.trts)
                {
                    mean.trt[t] <- mean(trt == unique.trts[t])
                }
                pf <- function(trt, x, match.id)
                {
                    pi.x <- numeric(length(trt))
                    for (t in 1:n.trts)
                    {
                        which.t       <- trt == unique.trts[t]
                        pi.x[which.t] <- mean(which.t)
                    }
                    pi.x
                }
                propensity_func <- pf
            }
        }
    }



    larger_outcome_better <- as.logical(larger_outcome_better[1])

    this.call     <- mget(names(formals()), sys.frame(sys.nframe()))
    this.call$... <- NULL
    this.call     <- c(this.call, list(...))


    # check to make sure arguments of propensity_func are correct
    propfunc.names <- sort(names(formals(propensity_func)))
    if (length(propfunc.names) == 3)
    {
        if (any(propfunc.names != c("match.id", "trt", "x")))
        {
            stop("arguments of propensity_func() should be 'trt','x', and (optionally) 'match.id'")
        }
    } else if (length(propfunc.names) == 2)
    {
        if (any(propfunc.names != c("trt", "x")))
        {
            stop("arguments of propensity_func() should be 'trt','x', and (optionally) 'match.id'")
        }
    } else
    {
        stop("propensity_func() should only have two or three arguments: 'trt' and 'x', or: 'trt', 'x', and 'match.id'")
    }

    # compute propensity scores
    if (is.null(match_id) | length(propfunc.names) == 2)
    {
        pi.x <- drop(propensity_func(x = x, trt = trt))
    } else
    {
        pi.x <- drop(propensity_func(x = x, trt = trt, match.id = match.id))
    }

    # make sure the resulting propensity scores are in the
    # acceptable range (ie 0-1)
    rng.pi <- range(pi.x)

    if (rng.pi[1] <= 0 | rng.pi[2] >= 1) stop("propensity_func() should return values between 0 and 1")

    ## if returned propensity score
    ## is a matrix, then pick out the
    ## right column for each row so we
    ## always get Pr(T = T_i | X = x)

    dim.pi.x <- dim(pi.x)
    if (!is.null(dim.pi.x))
    {
        if (length(dim.pi.x) == 1)
        {
            pi.x <- as.vector(pi.x)
        } else if (length(dim.pi.x) == 2)
        {
            if (ncol(pi.x) != n.trts)
            {
                stop("Number of columns in the matrix returned by propensity_func() is not the same
                     as the number of levels of 'trt'.")
            }
            if (is.factor(trt))
            {
                values <- levels(trt)[trt]
            } else
            {
                values <- trt
            }

            levels.pi.mat <- colnames(pi.x)
            if (is.null(levels.pi.mat))
            {
                levels.pi.mat <- unique.trts
            }

            # return the probability corresponding to the
            # treatment that was observed
            pi.x <- pi.x[cbind(1:nrow(pi.x), match(values, levels.pi.mat))]
        } else
        {
            stop("propensity_func() returns a multidimensional array; it can only return a vector or matrix.")
        }
    }


    ## indicator of zero
    z      <- as.integer(1*(y  > y_eps))
    is_nz  <- y > y_eps
    s      <- y[is_nz]
    x_s    <- x[is_nz,]
    trt_s  <- trt[is_nz]
    pi.x_s <- pi.x[is_nz]

    # construct design matrix to be passed to fitting function
    x.tilde.s <- create.design.matrix(x             = x_s,
                                      pi.x          = pi.x_s,
                                      trt           = trt_s,
                                      method        = "weighting",
                                      reference.trt = 0)

    ## need to re-order data
    reorder <- c(which(is_nz), which(!is_nz))
    x_z <- cbind(1, x[reorder,])
    z   <- z[reorder]
    trt <- trt[reorder]

    pi.x <- pi.x[reorder]


    # construct observation weight vector
    wts     <- create.weights(pi.x   = pi.x,
                              trt    = trt,
                              method = "weighting")

    # construct observation weight vector
    wts_s     <- create.weights(pi.x   = pi.x_s,
                                trt    = trt_s,
                                method = "weighting")


    extra.args <- NULL
    # check to make sure arguments of augment.func are correct
    if (!is.null(augment_func_zero))
    {
        B.x   <- unname(drop(augment_func_zero(trt = trt, x = x_z, y = z)))

        if (NROW(B.x) != NROW(y))
        {
            stop("augment_func_zero() should return the same number of predictions as observations in y")
        }

        extra.args.zero <- list(offset = B.x)

    } else
    {
        extra.args.zero <- list(offset = rep(0, NROW(z)))
    }
    if (!is.null(augment_func_positive))
    {
        B.x   <- unname(drop(augment_func_positive(trt = trt_s, x = x_s, y = s)))

        if (NROW(B.x) != NROW(y))
        {
            stop("augment_func_positive() should return the same number of predictions as observations in y")
        }

        extra.args.pos <- list(offset = B.x)

    } else
    {
        extra.args.pos <- list(offset = rep(0, NROW(s)))
    }


    comparison.trts <- 1

    if (comparison.trts[1] != 1 & comparison.trts[1] != "1")
    {
        trt.name.cur <- comparison.trts[1]
    } else
    {
        trt.name.cur <- "Trt1"
    }
    all.cnames <- c( trt.name.cur[1],
                     vnames )

    colnames(x.tilde.s) <- all.cnames
    colnames(x_z)       <- all.cnames

    fitted.model <- list()

    resid_outcome <- z
    trt_aug <- ifelse(resid_outcome >= 0, trt, 1 - trt)

    fitted.model$model <- cv.hd2part(x_z, trt_aug,  ## zero part data
                                     x.tilde.s, s,  ## positive part data
                                     weights   = wts * (abs(resid_outcome)), ## observation weights for zero part
                                     weights_s = wts_s, ## observation weights for positive part
                                     penalty   = penalty,
                                     algorithm = "irls",
                                     intercept = FALSE, ...)

    fitted.model$call                  <- this.call
    fitted.model$propensity.func       <- propensity_func
    fitted.model$augment_func_zero     <- augment_func_zero
    fitted.model$augment_func_positive <- augment_func_positive
    fitted.model$larger.outcome.better <- larger_outcome_better
    fitted.model$cutpoint              <- cutpoint
    fitted.model$var.names             <- vnames
    fitted.model$n.trts                <- n.trts
    fitted.model$comparison.trts       <- comparison.trts
    fitted.model$reference.trt         <- 0
    fitted.model$trts                  <- unique.trts
    fitted.model$trt.received          <- trt
    fitted.model$pi.x                  <- pi.x
    fitted.model$pi.x_s                <- pi.x_s
    fitted.model$y                     <- y
    fitted.model$z                     <- z
    fitted.model$s                     <- s


    pred_func <- function(x)
    {
        xbeta_zero <- predict(fitted.model$model, newx = cbind(1, x), model = "zero")
        xbeta_pos  <- predict(fitted.model$model, newx = cbind(1, x), model = "positive" )

        risk_ratio_estimate <- exp(2 * xbeta_pos + 1 * xbeta_zero)

        risk_ratio_estimate
    }

    fitted.model$predict <- pred_func

    bene.scores <- fitted.model$predict(x)

    attr(bene.scores, "comparison.trts") <- comparison.trts
    attr(bene.scores, "reference.trt")   <- 0
    attr(bene.scores, "trts")            <- unique.trts

    fitted.model$benefit.scores          <- bene.scores

    if (NROW(fitted.model$benefit.scores) != NROW(y))
    {
        warning("predict function returned a vector of predictions not equal to the number of observations
                when applied to the whole sample. Please check predict function.")
    }

    fitted.model$loss   <- "2part"
    fitted.model$family <- "2part"
    fitted.model$recommended.trts      <- personalized:::predict.subgroup_fitted(fitted.model,
                                                                                 newx = x,
                                                                                 type = "trt.group",
                                                                                 cutpoint = cutpoint)

    # calculate sizes of subgroups and the
    # subgroup treatment effects based on the
    # benefit scores and specified benefit score cutpoint
    fitted.model$subgroup.trt.effects <- subgroup.effects(fitted.model$benefit.scores,
                                                          y, trt,
                                                          pi.x,
                                                          cutpoint,
                                                          larger_outcome_better,
                                                          reference.trt = 0)

    class(fitted.model) <- "subgroup_fitted"

    fitted.model

}
