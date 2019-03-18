

create.design.matrix.binary.trt <- function(x, pi.x, trt, method, reference.trt = NULL)
{
    # trt must be supplied as integer vector
    # where 1 = treatment, 0 = control


    if (is.factor(trt))
    {
        # drop any unused levels of trt
        trt         <- droplevels(trt)
        unique.trts <- levels(trt)
        n.trts      <- length(unique.trts)
    } else
    {
        unique.trts <- sort(unique(trt))
        n.trts      <- length(unique.trts)
    }

    # if not specified, set reference treatment
    # to be the last one
    if (is.null(reference.trt))
    {
        reference.trt <- unique.trts[1L]
    }

    which.reference <- which(unique.trts == reference.trt)

    if (n.trts != 2) stop("two trtment levels only for binary trt design matrix function")

    # construct modified design matrices
    # depending on what method is used
    if (method == "weighting")
    {
        # create 1 and -1 version of treatment vector
        trt2 <- 2 * (trt != reference.trt) - 1

        x.tilde <- trt2 * cbind(1, x)
    } else
    {   # A-learning method
        x.tilde <- (1 * (trt != reference.trt) - pi.x) * cbind(1, x)
    }
    x.tilde
}

create.weights.binary.trt <- function(pi.x, trt, method)
{
    # trt must be supplied as integer vector
    # where 1 = treatment, 0 = control

    if (is.factor(trt))
    {
        # drop any unused levels of trt
        trt         <- droplevels(trt)
        unique.trts <- levels(trt)
        n.trts      <- length(unique.trts)
    } else
    {
        unique.trts <- sort(unique(trt))
        n.trts      <- length(unique.trts)
    }

    if (n.trts != 2) stop("two trtment levels only for binary trt weighting function")

    # construct weights
    # depending on what method is used
    if (method == "weighting")
    {
        wts     <- 1 / (pi.x * (trt == unique.trts[2L]) + (1 - pi.x) * (trt == unique.trts[1L]))
    } else
    {   # A-learning method
        wts     <- rep(1, length(pi.x))
    }
    wts
}






create.weights.mult.trt <- function(pi.x, trt, method)
{
    # trt must be supplied as factor (actually not anymore!)

    # construct weights
    # depending on what method is used
    if (method == "weighting")
    {
        wts     <- 1 / (pi.x)
    } else
    {   # A-learning method
        wts     <- rep(1, length(pi.x))
    }
    wts
}


create.design.matrix <- function(x, pi.x, trt, y, method, reference.trt = NULL)
{
    # check if multiple treatments or not
    if (is.factor(trt))
    {
        n.trts <- length(levels(trt))
    } else
    {
        n.trts <- length(unique(trt))
    }

    is.mult.trt <- n.trts > 2

    if (is.mult.trt)
    {
        # set to factor for multiple trtment trt vector if it isn't already
        if (!is.factor(trt)) trt <- as.factor(trt)

        return( create.design.matrix.mult.trt(x             = cbind(1, x),
                                              pi.x          = pi.x,
                                              trt           = trt,
                                              #y             = y,
                                              method        = method,
                                              reference.trt = reference.trt) )
    } else
    {
        return( create.design.matrix.binary.trt(x      = x,
                                                pi.x   = pi.x,
                                                trt    = trt,
                                                method = method,
                                                reference.trt = reference.trt) )
    }
}

create.weights <- function(pi.x, trt, method)
{
    # check if multiple treatments or not
    if (is.factor(trt))
    {
        n.trts <- length(levels(droplevels(trt)))
    } else
    {
        n.trts <- length(unique(trt))
    }

    is.mult.trt <- n.trts > 2

    if (is.mult.trt)
    {
        # set to factor for multiple trtment trt vector if it isn't already
        if (!is.factor(trt)) trt <- as.factor(trt)

        return( create.weights.mult.trt(pi.x   = pi.x,
                                        trt    = trt,
                                        method = method) )
    } else
    {
        return( create.weights.binary.trt(pi.x   = pi.x,
                                          trt    = trt,
                                          method = method) )
    }
}






# this creates a block matrix of contrasts where the reference
# treatment is treated as a "control".
# For example, if matrix x is supplied and treatment with three
# levels where the third level is the reference treatment,
# create.design.matrix.mult.trt() will return:
# |x_1   0  |
# |0     x_2|
# |-x_3 -x_3|
#
# where x = (x_1', x_2', x_3')' and x_i is the submatrix
# of x of observations with treatment status == i
#
# with binary treatments, this simplifies to the normal
# | x_1|
# |-x_2|
create.design.matrix.mult.trt <- function(x, pi.x, trt, y, reference.trt = NULL,
                                          method = c("weighting", "a_learning"))
{
    trt.levels      <- levels(trt)
    n.trts          <- length(trt.levels)
    trt.idx         <- vector(mode = "list", length = n.trts)
    sample.sizes    <- numeric(n.trts)
    method          <- match.arg(method)

    # if not specified, set reference treatment
    # to be the last one
    if (is.null(reference.trt))
    {
        reference.trt <- trt.levels[1L]
    }
    which.reference <- which(trt.levels == reference.trt)

    # re-order treatments so that
    # the specified reference will be treated as the reference
    trt.idx.vec <- 1:n.trts
    trt.idx.vec <- c(trt.idx.vec[-which.reference], trt.idx.vec[which.reference])
    trt.levels  <- trt.levels[trt.idx.vec]

    for (t in 1:n.trts)
    {
        trt.idx[[t]]    <- which(levels(trt)[trt] == trt.levels[t])
        sample.sizes[t] <- length(trt.idx[[t]])
    }

    n.obs  <- NROW(x)
    n.vars <- NCOL(x)

    # if (sd(x[,1]) != 0) stop("x must have intercept as first column")
    stopifnot(n.obs == length(trt))
    stopifnot(n.obs == sum(sample.sizes))

    x.return <- array(0, dim = c(n.obs, n.vars * (n.trts - 1)))

    var.idx.list        <- vector(mode = "list", length = n.trts - 1)
    names(var.idx.list) <- trt.levels[-n.trts] # remove reference treatment
    n.vars.cumsum       <- c(0, cumsum(rep(n.vars, n.trts)))
    n.obs.cumsum        <- c(0, cumsum(sample.sizes))

    y.return <- numeric(n.obs)
    for (t in 1:(n.trts - 1))
    {
        idx.obs.cur  <- (n.obs.cumsum[t] + 1):n.obs.cumsum[t + 1]
        idx.vars.cur <- (n.vars.cumsum[t] + 1):n.vars.cumsum[t + 1]

        # construct modified design matrices
        # depending on what method is used
        if (method == "weighting")
        {
            x.return[trt.idx[[t]], idx.vars.cur] <- x[trt.idx[[t]],]
        } else
        {
            stop("A-learning not available for multiple treatments")
            # A-learning method

            #x.return[trt.idx[[t]], idx.vars.cur] <- (1 - pi.x[trt.idx[[t]]]) * x[trt.idx[[t]],]
        }


        #y.return[idx.obs.cur] <- y[trt.idx[[t]]] # now we don't want to re-order observations
    }
    t <- n.trts

    for (r in 1:(n.trts - 1))
    {
        # keep observation index static
        idx.obs.cur  <- (n.obs.cumsum[t] + 1):n.obs.cumsum[t + 1]
        # replicate columns
        idx.vars.cur <- (n.vars.cumsum[r] + 1):n.vars.cumsum[r + 1]

        # construct modified design matrices
        # depending on what method is used
        if (method == "weighting")
        {
            x.return[trt.idx[[t]], idx.vars.cur] <- -x[trt.idx[[t]],]
        }# else
        #{
        #stop("A-learning not available for multiple treatments")
        # A-learning method
        #x.return[trt.idx[[t]], idx.vars.cur] <- -(1 - pi.x[trt.idx[[t]]]) * x[trt.idx[[t]],]
        #}

        #if (r == 1)
        #{
        #    y.return[idx.obs.cur] <- y[trt.idx[[t]]]
        #}
    }


    # 'trt.levels' is the treatment levels (re-ordered based off of reference treatment)
    attr(x.return, "trt.levels") <- trt.levels

    x.return
}
