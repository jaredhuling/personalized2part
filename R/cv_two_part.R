

cv.hd2part <- function (x, z,
                        x_s, s,
                        weights = rep(1, NROW(x)),
                        weights_s = rep(1, NROW(x_s)),
                        offset = NULL,
                        lambda = NULL,
                        type.measure = c("mse", "auc", "mae"),
                        nfolds = 10, foldid = NULL,
                        grouped = TRUE, keep = FALSE, parallel = FALSE, ...)
    {
        type.measure <- match.arg(type.measure)
        if (!is.null(lambda) && length(lambda) < 2)
        {
            stop("Need more than one value of lambda for cv.2part")
        }
        n <- nrow(x)
        n_s <- nrow(x_s)
        if (is.null(weights))
        {
            weights <- rep(1, n)
        }
        if (is.null(weights_s))
        {
            weights_s <- rep(1, n_s)
        }

        weights   <- as.double(weights)
        weights_s <- as.double(weights_s)
        y <- drop(y)
        hd2part.call <- match.call(expand.dots = TRUE)
        which <- match(c("type.measure", "nfolds", "foldid", "grouped", "keep"), names(hd2part.call), FALSE)
        if (any(which))
        {
            hd2part.call = hd2part.call[-which]
        }
        hd2part.call[[1]] = as.name("hd2part")
        hd2part.object = hd2part(x = x, z = z,
                                 x = x_s, s = s,
                                 weights = weights, weights_s = weights_s, #offset = offset,
                               lambda = lambda, ...)
        hd2part.object$call = hd2part.call
        subclass <- class(hd2part.object)[[1]]
        type.measure=cvtype(type.measure,subclass)
        is.offset   <- hd2part.object$offset
        is.offset.s <- hd2part.object$offset_s
        ###Next line is commented out so each call generates its own lambda sequence
        # lambda=glmnet.object$lambda
        nz <- sapply(predict(hd2part.object, type = "nonzero", model = "zero"), length)
        nz_s <- sapply(predict(hd2part.object, type = "nonzero", model = "positive"), length)
        if (is.null(foldid))
        {
            foldid <- sample(rep(seq(nfolds), length = n))
        } else
        {
            nfolds <- max(foldid)
        }
        if (nfolds < 2)
        {
            stop("nfolds must be bigger than 2; nfolds=10 recommended")
        }
        outlist <- as.list(seq(nfolds))
        if (parallel)
        {
            outlist = foreach(i = seq(nfolds), .packages = c("personalized2part")) %dopar%
            {
                which   <- foldid == i
                which_s <- which[1:n_s]

                if (length(dim(z))>1)
                {
                    z_sub <- z[!which, ]
                } else
                {
                    z_sub <- z[!which]
                }

                if (length(dim(s))>1)
                {
                    s_sub <- s[!which_s, ]
                } else
                {
                    s_sub <- s[!which_s]
                }

                if (is.offset)
                {
                    offset_sub <- as.matrix(offset)[!which, ]
                } else
                {
                    offset_sub <- NULL
                }

                if (is.offset.s)
                {
                    offset_sub_s <- as.matrix(offset_s)[!which_s, ]
                } else
                {
                    offset_sub_s <- NULL
                }

                hd2part(x[!which, , drop = FALSE], z_sub,
                        x_s[!which_s, , drop = FALSE], s_sub,
                        lambda    = lambda,
                        offset    = offset_sub,
                        offset_s  = offset_sub_s,
                        weights   = weights[!which],
                        weights_s = weights_s[!which_s],
                       ...)
            }
        } else
        {
            for (i in seq(nfolds))
            {
                which   <- foldid == i
                which_s <- which[1:n_s]

                if (length(dim(z))>1)
                {
                    z_sub <- z[!which, ]
                } else
                {
                    z_sub <- z[!which]
                }

                if (length(dim(s))>1)
                {
                    s_sub <- s[!which_s, ]
                } else
                {
                    s_sub <- s[!which_s]
                }

                if (is.offset)
                {
                    offset_sub <- as.matrix(offset)[!which, ]
                } else
                {
                    offset_sub <- NULL
                }

                if (is.offset.s)
                {
                    offset_sub_s <- as.matrix(offset_s)[!which_s, ]
                } else
                {
                    offset_sub_s <- NULL
                }

                outlist[[i]] <- hd2part(x[!which, , drop = FALSE], z_sub,
                                        x_s[!which_s, , drop = FALSE], s_sub,
                                        lambda    = lambda,
                                        offset    = offset_sub,
                                        offset_s  = offset_sub_s,
                                        weights   = weights[!which],
                                        weights_s = weights_s[!which_s],
                                        ...)
            }
        }

        fun <- paste("cv", subclass, sep = ".")
        lambda <- hd2part.object$lambda
        cvstuff <- do.call(fun, list(outlist, lambda, x, y, weights,
                                     offset, foldid, type.measure, grouped, keep))
        cvm = cvstuff$cvm
        cvsd = cvstuff$cvsd
        nas=is.na(cvsd)
        if(any(nas)){
            lambda=lambda[!nas]
            cvm=cvm[!nas]
            cvsd=cvsd[!nas]
            nz=nz[!nas]
        }
        cvname = names(cvstuff$type.measure)
        names(cvname)=cvstuff$type.measure# to be compatible with earlier version; silly, I know
        out = list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm +
                       cvsd, cvlo = cvm - cvsd, nzero = nz, name = cvname, glmnet.fit = glmnet.object)
        if (keep)
            out = c(out, list(fit.preval = cvstuff$fit.preval, foldid = foldid))
        lamin=if(cvname=="AUC")getmin(lambda,-cvm,cvsd)
        else getmin(lambda, cvm, cvsd)
        obj = c(out, as.list(lamin))
        class(obj) = "cv.glmnet"
        obj
    }
