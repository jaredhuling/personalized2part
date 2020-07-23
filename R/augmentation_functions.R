

#' Fit a penalized gamma augmentation model via cross fitting
#'
#' @description Fits a penalized gamma augmentation model via cross fitting and
#' returns vector of length n of out of sample predictions on the link scale from cross fitting
#'
#' @param x an n x p matrix of covariates for the zero part data, where each row is an observation
#' and each column is a predictor. MUST be ordered such that the first n_s rows align with the
#' observations in x_s and s
#' @param y a length n vector of responses taking positive values
#' @param trt a length n vector of treatment variables with 1 indicating treatment and -1 indicating control
#' @param wts a length n vector of sample weights
#' @param K number of folds for cross fitting
#' @param p tweedie mixing parameter. See \code{\link[HDtweedie]{HDtweedie}} for details
#' @param interactions boolean variable of whether or not to fit model with interactions. For predictions,
#' interactions will be integrated out
#' @export
#' @importFrom HDtweedie cv.HDtweedie
HDtweedie_kfold_aug <- function(x, y, trt, wts = NULL, K = 10, p = 1.5, interactions = FALSE)
{

    if (is.null(wts))
    {
        wts <- rep(1, NROW(x))
    }

    trt.levels <- sort(unique(trt))

    stopifnot("trt must take only values -1 and 1" = all(trt.levels == c(-1, 1)))


    if (interactions)
    {

        ## full model for nonzeroness
        df_all <- data.frame(x, trt = trt)
        df_1   <- data.frame(x, trt = 1)
        df_0   <- data.frame(x, trt = -1)

        mm_all <- model.matrix(~x*trt-1, data = df_all)
        mm_1   <- model.matrix(~x*trt-1, data = df_1)
        mm_0   <- model.matrix(~x*trt-1, data = df_0)

        colnames(mm_1) <- colnames(mm_0) <- colnames(mm_all)
    } else
    {
        mm_all <- x
    }



    n <- NROW(mm_all)

    foldid = sample(rep(seq(K), length = n))

    predvec <- numeric(n)

    for (i in seq(K))
    {
        which <- foldid == i

        ## full model for nonzeroness
        gamfit_pos_main  <- HDtweedie::cv.HDtweedie(y = y[!which], x = mm_all[!which,,drop = FALSE],
                                                    weights = wts[!which],
                                                    p = p, eps = 1e-5,
                                                    pred.loss = "mae")

        if (interactions)
        {
            ## get predictions for trt = 1 & -1
            pred1 <- unname(drop(predict(gamfit_pos_main, newx = mm_1[which,,drop=FALSE], s = "lambda.min", type = "link")))
            pred0 <- unname(drop(predict(gamfit_pos_main, newx = mm_0[which,,drop=FALSE], s = "lambda.min", type = "link")))

            predvec[which] <- 0.5 * (pred1 + pred0)
        } else
        {
            ## get predictions
            pred <- unname(drop(predict(gamfit_pos_main, newx = mm_all[which,,drop=FALSE], s = "lambda.min", type = "link")))

            predvec[which] <- pred
        }

    }

    predvec
}
