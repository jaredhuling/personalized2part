
setup_y <- function(y, family)
{
    if (family == "binomial")
    {
        yvals <- sort(unique(y))
        if (all(yvals == c(0, 1)))
        {
            y <- 2 * y - 1
        } else
        {
            if (!all(yvals == c(-1, 1)))
            {
                stop("y has invalid values, can only take values 0 and 1")
            }
        }
    } else if (family == "gamma")
    {
        if (any(y <= 0))
        {
            stop("y must be strictly positive")
        }
    }

    y
}



check_compute_propensity_func <- function(x, trt, match.id, propensity.func, propfunc_title, n.trts, unique.trts)
{
    # check to make sure arguments of propensity.func are correct
    propfunc.names <- sort(names(formals(propensity.func)))
    if (length(propfunc.names) == 3)
    {
        if (any(propfunc.names != c("match.id", "trt", "x")))
        {
            stop(paste0("arguments of ", propfunc_title, " should be 'trt','x', and (optionally) 'match.id'"))
        }
    } else if (length(propfunc.names) == 2)
    {
        if (any(propfunc.names != c("trt", "x")))
        {
            stop(paste0("arguments of ", propfunc_title, " should be 'trt','x', and (optionally) 'match.id'"))
        }
    } else
    {
        stop(paste0(propfunc_title, " should only have two or three arguments: 'trt' and 'x', or: 'trt', 'x', and 'match.id'"))
    }

    # compute propensity scores
    if (is.null(match.id) | length(propfunc.names) == 2)
    {
        pi.x <- drop(propensity.func(x = x, trt = trt))
    } else
    {
        pi.x <- drop(propensity.func(x = x, trt = trt, match.id = match.id))
    }

    # make sure the resulting propensity scores are in the
    # acceptable range (ie 0-1)
    rng.pi <- range(pi.x)

    if (rng.pi[1] <= 0 | rng.pi[2] >= 1) stop(paste0(propfunc_title, " should return values between 0 and 1"))

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
                stop(paste0("Number of columns in the matrix returned by ", propfunc_title, " is not the same
                     as the number of levels of 'trt'."))
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
                stop(paste0(propfunc_title, " returns a multidimensional array; it can only return a vector or matrix."))
            }
    }

    pi.x
}
