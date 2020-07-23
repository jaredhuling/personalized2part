
#' Fitting subgroup identification models for semicontinuous positive outcomes
#'
#' @description Fits subgroup identification models
#'
#' @param x The design matrix (not including intercept term)
#' @param y The nonnegative response vector
#' @param trt treatment vector with each element equal to a 0 or a 1, with 1 indicating
#'            treatment status is active.
#' @param propensity.func function that inputs the design matrix x and the treatment vector trt and outputs
#' the propensity score, ie Pr(trt = 1 | X = x). Function should take two arguments 1) x and 2) trt. See example below.
#' For a randomized controlled trial this can simply be a function that returns a constant equal to the proportion
#' of patients assigned to the treatment group, i.e.:
#' \code{propensity.func = function(x, trt) 0.5}.
#' @param propensity.func.positive function that inputs the design matrix x and the treatment vector trt and outputs
#' the propensity score for units with positive outcome values, ie Pr(trt = 1 | X = x, Z = 1). Function should take
#' two arguments 1) x and 2) trt. See example below.
#' For a randomized controlled trial this can simply be a function that returns a constant equal to the proportion
#' of patients assigned to the treatment group, i.e.:
#' \code{propensity.func = function(x, trt) 0.5}.
#' @param match.id a (character, factor, or integer) vector with length equal to the number of observations in \code{x}
#' indicating using integers or levels of a factor vector which patients are
#' in which matched groups. Defaults to \code{NULL} and assumes the samples are not from a matched cohort. Matched
#' case-control groups can be created using any method (propensity score matching, optimal matching, etc). If each case
#' is matched with a control or multiple controls, this would indicate which case-control pairs or groups go together.
#' If \code{match.id} is supplied, then it is unecessary to specify a function via the \code{propensity.func} argument.
#' A quick usage example: if the first patient is a case and the second and third are controls matched to it, and the
#' fouth patient is a case and the fifth through seventh patients are matched with it, then the user should specify
#' \code{match.id = c(1,1,1,2,2,2,2)} or \code{match.id = c(rep("Grp1", 3),rep("Grp2", 4)) }
#' the covariates \code{x}, and \code{trt} and outputs predicted values (on the probability scale) for the response using a model
#' constructed with \code{x}. \code{augment.func.zero()} can also be simply
#' a function of \code{x} and \code{y}. This function is used for efficiency augmentation.
#' When the form of the augmentation function is correct, it can provide efficient estimation of the subgroups. Some examples of possible
#' augmentation functions are:
#'
#' Example 1: \code{augment.func <- function(x, y) {lmod <- glm(y ~ x, family = binomial()); return(fitted(lmod))}}
#'
#' Example 2:
#' \preformatted{
#' augment.func <- function(x, y, trt) {
#'     data <- data.frame(x, y, trt)
#'     lmod <- glm(y ~ x * trt, family = binomial())
#'     ## get predictions when trt = 1
#'     data$trt <- 1
#'     preds_1  <- predict(lmod, data, type = "response")
#'
#'     ## get predictions when trt = -1
#'     data$trt <- -1
#'     preds_n1 <- predict(lmod, data, type = "response")
#'
#'     ## return predictions averaged over trt
#'     return(0.5 * (preds_1 + preds_n1))
#' }
#' }
#'
#' @param augment.func.positive (similar to augment.func.zero) function which inputs the positive part response
#'  (ie all observations in \code{y} which are strictly positive),
#'  the covariates \code{x}, and \code{trt} and outputs predicted values (on the link scale) for the response using a model
#'  constructed with \code{x}. \code{augment.func.positive()} can also be simply
#'  a function of \code{x} and \code{y}. This function is used for efficiency augmentation.
#'  @param augment.func.zero (similar to augment.func.positive) function which inputs the
#'  indicators oof whether each response is positive (\code{1*(y > 0)}),
#'  the covariates \code{x}, and \code{trt} for all samples and outputs predicted values (on the link scale) for the response using a model
#'  constructed with \code{x}. \code{augment.func.positive()} can also be simply
#'  a function of \code{x} and \code{y}. This function is used for efficiency augmentation.
#' @param cutpoint numeric value for patients with benefit scores above which
#'  (or below which if \code{larger.outcome.better = FALSE})
#'  will be recommended to be in the treatment group. Defaults to 1, since the benefit score is a risk ratio
#' @param larger.outcome.better boolean value of whether a larger outcome is better/preferable. Set to \code{TRUE}
#'  if a larger outcome is better/preferable and set to \code{FALSE} if a smaller outcome is better/preferable. Defaults to \code{TRUE}.
#' @param y_eps positive value above which observations in \code{y} will be considered positive
#' @param ... options to be passed to \code{\link[personalized2part]{cv.hd2part}}
#' @export
#'
#' @examples
#'
#' set.seed(42)
#'
#' dat <- sim_semicontinuous_data(250, n.vars = 15)
#' x <- dat$x
#' y <- dat$y
#' trt <- dat$trt
#'
#' prop_func <- function(x, trt)
#' {
#'     propensmod <- glm(trt ~ x, family = binomial())
#'
#'     propens <- unname(fitted(propensmod))
#'     propens
#' }
#'
#' fitted_model <- fit_subgroup_2part(x, y, trt, prop_func, prop_func)
#'
#' fitted_model
#'
#' ## correlation of estimated covariate-conditional risk ratio and truth
#' cor(fitted_model$benefit.scores, dat$treatment_risk_ratio, method = "spearman")
#'
#'
fit_subgroup_2part <- function(x,
                               y,
                               trt,
                               propensity.func = NULL,
                               propensity.func.positive = NULL,
                               match.id = NULL,
                               augment.func.zero = NULL,
                               augment.func.positive = NULL,
                               cutpoint              = 1,
                               larger.outcome.better = TRUE,
                               y_eps = 1e-6,
                               ...)
{
    dims    <- dim(x)
    if (is.null(dims)) stop("x must be a matrix object.")

    y       <- drop(y)
    vnames  <- colnames(x)

    p   <- NCOL(x)
    n   <- NROW(x)

    if (any(y < 0))
    {
        stop("y must be nonnegative")
    }

    if (is.null(vnames))
    {
        vnames <- paste0("V", 1:p)
    }

    # check to make sure arguments of augment_func are correct
    if (!is.null(augment.func.zero))
    {
        augmentfunc.names <- sort(names(formals(augment.func.zero)))
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
            augment.func2 <- augment.func.zero
            augment.func.zero  <- function(trt, x, y) augment.func2(x = x, y = y)
        } else
        {
            stop("augment.func.zero() should only have either two arguments: 'x' and 'y', or three arguments:
                 'trt', 'x', and 'y'")
        }
    }
    if (!is.null(augment.func.positive))
    {
        augmentfunc.names <- sort(names(formals(augment.func.positive)))
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
            augment.func2 <- augment.func.positive
            augment.func.positive  <- function(trt, x, y) augment.func2(x = x, y = y)
        } else
        {
            stop("augment.func.positive() should only have either two arguments: 'x' and 'y', or three arguments:
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
    if (is.null(propensity.func))
    {
        if (is.null(match.id))
        { # No propensity score supplied and no match.id supplied
            if (n.trts == 2)
            {
                mean.trt <- mean(trt == 1)
                propensity.func <- function(trt, x) rep(mean.trt, length(trt))
            } else
            {
                mean.trt <- numeric(n.trts)
                for (t in 1:n.trts)
                {
                    mean.trt[t] <- mean(trt == unique.trts[t])
                }
                propensity.func <- function(trt, x)
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
                propensity.func <- pf
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
                propensity.func <- pf
            }
        }
    }

    if (is.null(propensity.func.positive))
    {
        if (n.trts == 2)
        {
            mean.trt <- mean(trt[which_pos] == 1)
            propensity.func.positive <- function(trt, x) rep(mean.trt, length(trt))
        } else
        {
            mean.trt <- numeric(n.trts)
            for (t in 1:n.trts)
            {
                mean.trt[t] <- mean(trt[which_pos] == unique.trts[t])
            }
            propensity.func <- function(trt, x)
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
    }


    which_pos <- y > y_eps


    if (sum(which_pos) == 0)
    {
        stop("no observations are larger than zero")
    }
    if (sum(!which_pos) == 0)
    {
        stop("no observations are zero")
    }

    ## reorder observations so positive responses are first
    reorder.idx <- c(which(which_pos),
                     which(!which_pos))

    x <- x[reorder.idx,]
    y <- y[reorder.idx]
    trt <- trt[reorder.idx]

    if (!is.null(match.id))
    {
        match.id <- match.id[reorder.idx]
    }



    larger.outcome.better <- as.logical(larger.outcome.better[1])

    this.call     <- mget(names(formals()), sys.frame(sys.nframe()))
    this.call$... <- NULL
    this.call     <- c(this.call, list(...))


    pi.x <- check_compute_propensity_func(x, trt, match.id, propensity.func,
                                          "propensity.func()", n.trts, unique.trts)


    ## indicator of zero
    z      <- as.integer(1*(y  > y_eps))

    ## flip meaning of z
    if (!larger.outcome.better)
    {
        z <- 1 - z
    }

    is_nz  <- y > y_eps
    s      <- y[is_nz]
    x_s    <- x[is_nz,]
    trt_s  <- trt[is_nz]
    #pi.x_s <- pi.x[is_nz]

    pi.x_s <- check_compute_propensity_func(x_s, trt_s, NULL, propensity.func.positive,
                                            "propensity.func.positive()", n.trts, unique.trts)

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
    if (!is.null(augment.func.zero))
    {
        B.x   <- unname(drop(augment.func.zero(trt = trt, x = x_z, y = z)))

        if (NROW(B.x) != NROW(y))
        {
            stop("augment.func.zero() should return the same number of predictions as observations in y")
        }

        extra.args.zero <- list(offset = B.x)

    } else
    {
        extra.args.zero <- list(offset = rep(0, NROW(z)))
    }
    if (!is.null(augment.func.positive))
    {
        B.x   <- unname(drop(augment.func.positive(trt = trt_s, x = x_s, y = s)))

        if (NROW(B.x) != NROW(s))
        {
            stop("augment.func.positive() should return the same number of predictions as observations in y's with positive values")
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

    resid_outcome <- z - extra.args.zero$offset
    trt_aug <- ifelse(resid_outcome >= 0, trt, 1 - trt)

    fitted.model$model <- cv.hd2part(x_z, trt_aug,  ## zero part data
                                     x.tilde.s, s,  ## positive part data
                                     weights   = wts * (abs(resid_outcome)), ## observation weights for zero part
                                     weights_s = wts_s, ## observation weights for positive part
                                     opposite_signs = !larger.outcome.better,
                                     intercept_z = FALSE,
                                     intercept_s = FALSE, ...)

    if (!larger.outcome.better)
    {
        fitted.model$model$beta_z <- -fitted.model$model$beta_z
    }


    fitted.model$call                  <- this.call
    fitted.model$propensity.func       <- propensity.func
    fitted.model$augment.func.zero     <- augment.func.zero
    fitted.model$augment.func.positive <- augment.func.positive
    fitted.model$larger.outcome.better <- larger.outcome.better
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
                                                          larger.outcome.better,
                                                          reference.trt = 0)

    treat.effects.2part <- function(benefit.scores, loss = "2part", method = "weighting", pi.x = NULL)
    {
        loss   <- match.arg(loss)
        method <- match.arg(method)

        benefit.scores <- drop(benefit.scores)
        if (!is.null(pi.x))
        {
            pi.x <- drop(pi.x)
        }

        trt_eff_delta <- trt_eff_gamma <- NA

        trt_eff_gamma <- benefit.scores

        effects <- list(delta = trt_eff_delta,
                        gamma = trt_eff_gamma)
        class(effects) <- c("individual_treatment_effects", class(effects) )

        effects
    }

    fitted.model$individual.trt.effects <- treat.effects.2part(fitted.model$benefit.scores,
                                                               "2part",
                                                               "weighting",
                                                               pi.x)

    class(fitted.model) <- "subgroup_fitted"

    fitted.model

}
