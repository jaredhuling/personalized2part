

cv.hd2part <- function (x,   z,
                        x_s, s,
                        weights      = rep(1, NROW(x)),
                        weights_s    = rep(1, NROW(x_s)),
                        offset       = NULL,
                        offset_s     = NULL,
                        lambda       = NULL,
                        type.measure = c("mae", "mse", "sep-auc-mse", "sep-auc-mae"),
                        nfolds       = 10,
                        foldid       = NULL,
                        grouped      = TRUE,
                        keep         = FALSE,
                        parallel     = FALSE,
                        ...)
{
    type.measure <- match.arg(type.measure)
    if (!is.null(lambda) && length(lambda) < 2)
    {
        stop("Need more than one value of lambda for cv.2part")
    }
    n   <- nrow(x)
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
    s <- drop(s)
    z <- drop(z)
    hd2part.call <- match.call(expand.dots = TRUE)
    which <- match(c("type.measure", "nfolds", "foldid", "grouped", "keep"), names(hd2part.call), FALSE)
    if (any(which))
    {
        hd2part.call = hd2part.call[-which]
    }
    hd2part.call[[1]] = as.name("hd2part")
    hd2part.object = hd2part(x = x, z = z,
                             x_s = x_s, s = s,
                             weights = weights, weights_s = weights_s, #offset = offset,
                           lambda = lambda, ...)
    hd2part.object$call <- hd2part.call
    subclass <- class(hd2part.object)[[1]]
    #type.measure=cvtype(type.measure,subclass)
    is.offset   <- hd2part.object$offset
    is.offset.s <- hd2part.object$offset_s
    ###Next line is commented out so each call generates its own lambda sequence
    # lambda=glmnet.object$lambda
    nz   <- sapply(predict(hd2part.object, type = "nonzero", model = "zero"),     length)
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

    #fun <- paste("cv", subclass, sep = ".")
    lambda <- hd2part.object$lambda

    cv_res_zero <- cv_function(outlist, lambda, x, z,
                               weights, offset, foldid, "auc",
                               "zero",
                               grouped, keep)

    tm_pos <- type.measure
    if (type.measure == "sep-auc-mse")
    {
        tm_pos <- "mse"
    } else if (type.measure == "sep-auc-mae")
    {
        tm_pos <- "mae"
    }

    cv_res_pos <- cv_function(outlist, lambda, x_s, s,
                              weights_s, offset_s, foldid[1:NROW(x_s)], tm_pos,
                              "positive",
                              grouped, keep)



    cvm_z  <- cv_res_zero$cvm
    cvsd_z <- cv_res_zero$cvsd
    nas    <- is.na(cvsd_z)
    if(any(nas))
    {
        lambda_z <- lambda[!nas]
        cvm_z    <- cvm_z[!nas]
        cvsd_z   <- cvsd_z[!nas]
        nz_z     <- nz[!nas]
    }

    cvm_s  <- cv_res_pos$cvm
    cvsd_s <- cv_res_pos$cvsd
    nas    <- is.na(cvsd_s)
    if(any(nas))
    {
        lambda_s <- lambda[!nas]
        cvm_s    <- cvm_s[!nas]
        cvsd_s   <- cvsd_s[!nas]
        nz_s     <- nz_s[!nas]
    }

    lamin_z <- getmin(lambda,-cvm_z, cvsd_z)  ## negative due to AUC
    lamin_s <- getmin(lambda, cvm_s, cvsd_s)

    names(lamin_z) <- c("lambda.min.z", "lambda.1se.z")
    names(lamin_s) <- c("lambda.min.s", "lambda.1se.s")

    if (type.measure %in% c("mae", "mse"))
    {
        cv_res <- cv_function_all(outlist, lambda,
                                  x,   z,
                                  x_s, s,
                                  weights, weights_s,
                                  offset, offset_s,
                                  foldid, type.measure,
                                  grouped, keep = keep)

        cv_res_sep <- cv_function_all(outlist, lambda,
                                      x,   z,
                                      x_s, s,
                                      weights, weights_s,
                                      offset, offset_s,
                                      foldid, type.measure,
                                      grouped, keep = keep,
                                      lambda_s = lamin_s$lambda.min.s)
    } else
    {
        cv_res <- cv_function_all(outlist, lambda,
                                  x,   z,
                                  x_s, s,
                                  weights, weights_s,
                                  offset, offset_s,
                                  foldid, "mae",
                                  grouped, keep = keep)

        cv_res_sep <- cv_function_all(outlist, lambda,
                                      x,   z,
                                      x_s, s,
                                      weights, weights_s,
                                      offset, offset_s,
                                      foldid, "mae",
                                      grouped, keep = keep,
                                      lambda_s = lamin_s$lambda.min.s)
    }


    cvm  <- cv_res$cvm
    cvsd <- cv_res$cvsd
    nas  <- is.na(cvsd)
    if(any(nas))
    {
        lambda <- lambda[!nas]
        cvm    <- cvm[!nas]
        cvsd   <- cvsd[!nas]
        nz     <- nz[!nas]
    }



    cvm_sep  <- cv_res_sep$cvm
    cvsd_sep <- cv_res_sep$cvsd
    nas  <- is.na(cvsd_sep)
    if(any(nas))
    {
        lambda_sep <- lambda[!nas]
        cvm_sep    <- cvm_sep[!nas]
        cvsd_sep   <- cvsd_sep[!nas]
    }

    lamin_sep <- getmin(lambda_sep, cvm_sep, cvsd_sep)

    names(lamin_sep) <- c("lambda.min.z.sep", "lambda.1se.z")

    lamin_sep <- c(lamin_sep, list(lambda.min.s.sep = lamin_s$lambda.min.s))


    cvname <- type.measure

    out = list(lambda   = lambda,
               cvm      = cvm,
               cvsd     = cvsd,
               cvup     = cvm + cvsd,
               cvlo     = cvm - cvsd,
               cvm_sep  = cvm_sep,
               cvsd_sep = cvsd_sep,
               cvup_sep = cvm_sep + cvsd_sep,
               cvlo_sep = cvm_sep - cvsd_sep,
               cvm_z    = cvm_z,
               cvsd_z   = cvsd_z,
               cvup_z   = cvm_z + cvsd_z,
               cvlo_z   = cvm_z - cvsd_z,
               cvm_s    = cvm_s,
               cvsd_s   = cvsd_s,
               cvup_s   = cvm_s + cvsd_s,
               cvlo_s   = cvm_s - cvsd_s,
               nzero_z  = nz_z,
               nzero_s  = nz_s,
               name     = cvname,
               hd2part.fit = hd2part.object)
    if (keep)
        out = c(out, list(fit.preval = cv_res$fit.preval, foldid = foldid))



    lamin   <- getmin(lambda, cvm, cvsd)

    obj = c(out, as.list(lamin),
                 as.list(lamin_sep),
                 as.list(lamin_z),
                 as.list(lamin_s))
    class(obj) = "cv.hd2part"
    obj
}


## from glmnet
cv_function <- function (outlist, lambda, x, y,
                         weights, offset, foldid, type.measure,
                         model = c("zero", "positive"),
                         grouped, keep = FALSE)
{
    prob_min = 1e-05
    prob_max = 1 - prob_min
    model <- match.arg(model)
    nc = dim(y)
    if (is.null(nc))
    {
        if (model == "zero")
        {
            y = as.factor(y)
            ntab = table(y)
            nc = as.integer(length(ntab))
            y = diag(nc)[as.numeric(y), ]
        }
    }
    N = NROW(y)
    nfolds = max(foldid)
    if ((N/nfolds < 10) && type.measure == "auc")
    {
        warning("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; changed to type.measure='deviance'. Alternatively, use smaller value for nfolds",
                call. = FALSE)
        type.measure = "deviance"
    }
    if ((N/nfolds < 3) && grouped) {
        warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold",
                call. = FALSE)
        grouped = FALSE
    }
    if (!is.null(offset))
    {
        is.offset = TRUE
        offset = drop(offset)
    }
    else is.offset = FALSE
    mlami=max(sapply(outlist,function(obj)min(obj$lambda)))
    which_lam=lambda >= mlami

    predmat = matrix(NA, NROW(y), length(lambda))
    nlams = double(nfolds)
    for (i in seq(nfolds))
    {
        which = foldid == i
        fitobj = outlist[[i]]
        if (is.offset)
        {
            off_sub = offset[which]
        } else
        {
            off_sub <- NULL
        }
        preds = predict(fitobj,x[which, , drop = FALSE],
                        s=lambda[which_lam],
                        newoffset = off_sub,
                        model = model,
                        type = "model_response")
        nlami = sum(which_lam)
        predmat[which, seq(nlami)] = preds
        nlams[i] = nlami
    }
    if (type.measure == "auc")
    {
        cvraw = matrix(NA, nfolds, length(lambda))
        good = matrix(0, nfolds, length(lambda))
        for (i in seq(nfolds))
        {
            good[i, seq(nlams[i])] = 1
            which = foldid == i
            for (j in seq(nlams[i]))
            {
                cvraw[i, j] = auc.mat(y[which, ], predmat[which, j], weights[which])
            }
        }
        N = apply(good, 2, sum)
        weights = tapply(weights, foldid, sum)
    } else
    {
        if (model == "zero")
        {
            ywt = apply(y, 1, sum)
            y = y/ywt
            weights = weights * ywt
        }

        N = NROW(y) - apply(is.na(predmat), 2, sum)
        if (model == "zero")
        {
            cvraw = switch(type.measure,
                           mse = (y[, 1] - (1 - predmat))^2 +
                               (y[, 2] - predmat)^2,
                           mae = abs(y[, 1] - (1 - predmat)) +
                               abs(y[, 2] - predmat),
                           deviance = {
                               predmat = pmin(pmax(predmat, prob_min), prob_max)
                               lp = y[, 1] * log(1 - predmat) + y[, 2] * log(predmat)
                               ly = log(y)
                               ly[y == 0] = 0
                               ly = drop((y * ly) %*% c(1, 1))
                               2 * (ly - lp)
                           },
                           class = y[, 1] * (predmat > 0.5) + y[, 2] * (predmat <=0.5)
            )
        } else
        {
            cvraw = switch(type.measure,
                           mse = (y - predmat)^2,
                           deviance = (y - predmat)^2,
                           mae = abs(y - predmat))
        }
        if (grouped)
        {
            cvob = cvcompute(cvraw, weights, foldid, nlams)
            cvraw = cvob$cvraw
            weights = cvob$weights
            N = cvob$N
        }
    }
    cvm  = apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)
    cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean,
                      w = weights, na.rm = TRUE)/(N - 1))
    out  = list(cvm = cvm,
                cvsd = cvsd,
                type.measure = type.measure)
    if (keep)
    {
        out$fit.preval = predmat
    }
    out
}


cv_function_all <- function (outlist, lambda,
                             x,   z,
                             x_s, s,
                             weights, weights_s,
                             offset, offset_s,
                             foldid, type.measure,
                             grouped, keep = FALSE,
                             lambda_s = NULL)
{
    prob_min = 1e-05
    prob_max = 1 - prob_min
    nc = dim(z)
    if (is.null(nc))
    {
        z = as.factor(z)
        ntab = table(z)
        nc = as.integer(length(ntab))
        z = diag(nc)[as.numeric(z), ]
    }
    N   <- NROW(z)
    N_s <- NROW(s)
    nfolds = max(foldid)
    if ((N_s/nfolds < 3) && grouped)
    {
        warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold",
                call. = FALSE)
        grouped = FALSE
    }
    if (!is.null(offset))
    {
        is.offset = TRUE
        offset = drop(offset)
    } else
    {
        is.offset = FALSE
    }
    if (!is.null(offset_s))
    {
        is.offset.s = TRUE
        offset_s = drop(offset_s)
    } else
    {
        is.offset.s = FALSE
    }
    mlami     = max(sapply(outlist,function(obj)min(obj$lambda)))
    which_lam = lambda >= mlami

    predmat = matrix(NA, NROW(z), length(lambda))
    nlams   = double(nfolds)
    for (i in seq(nfolds))
    {
        which   <- foldid == i
        which_s <- which[1:N_s]
        fitobj  <- outlist[[i]]
        if (is.offset)
        {
            off_sub <- offset[which]
        } else
        {
            off_sub <- NULL
        }
        if (is.offset.s)
        {
            off_sub_s <- offset_s[which_s]
        } else
        {
            off_sub_s <- NULL
        }
        preds_z = predict(fitobj,
                          x[which, , drop = FALSE],
                          s=lambda[which_lam],
                          newoffset = off_sub,
                          model = "zero",
                          type = "model_response")

        if (is.null(lambda_s))
        {
            preds_s = predict(fitobj,
                              x[which, , drop = FALSE],
                              s=lambda[which_lam],
                              #newoffset = off_sub,
                              model = "positive",
                              type = "model_response")

            preds <- preds_z * preds_s
        } else
        {
            #lambda_s[which_lam_s]
            preds_s = predict(fitobj,
                              x[which, , drop = FALSE],
                              s=lambda_s[1],
                              #newoffset = off_sub,
                              model = "positive",
                              type = "model_response")
            preds <- preds_z * array(drop(preds_s), dim = dim(preds_z))
        }

        nlami <- sum(which_lam)
        predmat[which, seq(nlami)] = preds
        nlams[i] = nlami
    }




    zwt = apply(z, 1, sum)
    z   = z/zwt
    #weights = weights * ywt


    y <- ifelse(z[,2] == 0, 0, s)

    N = NROW(y) - apply(is.na(predmat), 2, sum)
    cvraw = switch(type.measure,
                   mse = (y - predmat)^2,
                   deviance = (y - predmat)^2,
                   mae = abs(y - predmat))
    if (grouped)
    {
        cvob = cvcompute(cvraw, weights, foldid, nlams)
        cvraw = cvob$cvraw
        weights = cvob$weights
        N = cvob$N
    }
    cvm  = apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)
    cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean,
                      w = weights, na.rm = TRUE)/(N - 1))
    out  = list(cvm = cvm,
                cvsd = cvsd,
                type.measure = type.measure)
    if (keep)
    {
        out$fit.preval = predmat
    }
    out
}



## from glmnet
cvcompute=function(mat, weights, foldid, nlams)
{
    ###Computes the weighted mean and SD within folds, and hence the se of the mean
    wisum=tapply(weights,foldid,sum)
    nfolds=max(foldid)
    outmat=matrix(NA,nfolds,ncol(mat))
    good=matrix(0,nfolds,ncol(mat))
    mat[is.infinite(mat)]=NA#just in case some infinities crept in
    for(i in seq(nfolds)){
        mati=mat[foldid==i,,drop=FALSE]
        wi=weights[foldid==i]
        outmat[i,]=apply(mati,2,weighted.mean,w=wi,na.rm=TRUE)
        good[i,seq(nlams[i])]=1
    }
    N=apply(good,2,sum)
    list(cvraw=outmat,weights=wisum,N=N)
}

## from glmnet
getmin=function(lambda,cvm,cvsd){
    cvmin=min(cvm,na.rm=TRUE)
    idmin=cvm<=cvmin
    lambda.min=max(lambda[idmin],na.rm=TRUE)
    idmin=match(lambda.min,lambda)
    semin=(cvm+cvsd)[idmin]
    idmin=cvm<=semin
    lambda.1se=max(lambda[idmin],na.rm=TRUE)
    list(lambda.min=lambda.min,lambda.1se=lambda.1se)
}
