#' Plot method for hd2part fitted objects
#'
#' @param x fitted "hd2part" model object
#' @param model either \code{"zero"} for the zero part model or \code{"positive"} for the positive part model
#' @param xvar What is on the X-axis. \code{"norm"} plots against the L1-norm of the coefficients, \code{"lambda"} against the log-lambda sequence, and \code{"dev"}
#' against the percent deviance explained.
#' @param labsize size of labels for variable names. If labsize = 0, then no variable names will be plotted
#' @param xlab label for x-axis
#' @param ylab label for y-axis
#' @param main main title for plot
#' @param ... other graphical parameters for the plot
#' @rdname plot
#' @export
#' @examples
#' set.seed(123)
#'
plot.hd2part <- function(x,
                         model = c("zero", "positive"),
                         xvar = c("loglambda", "norm", "lambda"),
                         labsize = 0.6,
                         xlab = iname, ylab = NULL,
                         main = paste(model, "model"),
                         ...)
{
    model <- match.arg(model)
    xvar <- match.arg(xvar)

    if (model == "zero")
    {
        nbeta <- as.matrix(x$beta_z[-1,]) ## remove intercept
        nzero <- x$nonzero_z
    } else
    {
        nbeta <- as.matrix(x$beta_s[-1,]) ## remove intercept
        nzero <- x$nonzero_s
    }


    remove <- apply(nbeta, 1, function(betas) all(betas == 0) )
    switch(xvar,
           "norm" = {
               index    <- apply(abs(nbeta), 2, sum)
               iname    <- expression(L[1] * " Norm")
               xlim     <- range(index)
               approx.f <- 1
           },
           "lambda" = {
               index    <- x$lambda
               iname    <- expression(lambda)
               xlim     <- rev(range(index))
               approx.f <- 0
           },
           "loglambda" = {
               index    <- log(x$lambda)
               iname    <- expression(log(lambda))
               xlim     <- rev(range(index))
               approx.f <- 1
           },
           "dev" = {
               index    <- x$sumSquare
               iname    <- "Sum of Squares"
               xlim     <- range(index)
               approx.f <- 1
           }
    )
    if (all(remove)) stop("All beta estimates are zero for all values of lambda. No plot returned.")


    cols <- rainbow(sum(!remove))

    ## create sequence that grabs one of ROYGBIV and repeats with
    ## an increment up the rainbow spectrum with each step from 1:7 on ROYGBIV
    n.cols <- 7L
    scramble.seq <- rep(((1:n.cols) - 1) * (length(cols) %/% (n.cols)) + 1, length(cols) %/% n.cols)[1:length(cols)] +
        (((0:(length(cols)-1)) %/% n.cols))

    scramble.seq[is.na(scramble.seq)] <- which(!(1:length(cols) %in% scramble.seq))
    colseq <- cols[scramble.seq]


    matplot(index, t(nbeta[!remove,,drop=FALSE]),
            lty = 1,
            xlab = xlab,
            ylab = "",
            col = colseq,
            xlim = xlim,
            type = 'l', ...)

    if (is.null(ylab))
    {
        mtext(expression(hat(beta)), side = 2, cex = par("cex"), line = 3, las = 1)
    } else
    {
        mtext(ylab, side = 2, cex = par("cex"), line = 3)
        ylab = ""
    }

    atdf <- pretty(index, n = 10L)
    plotnz <- approx(x = index, y = nzero, xout = atdf, rule = 2, method = "constant", f = approx.f)$y
    axis(side=3, at = atdf, labels = plotnz, tick=FALSE, line=0, ...)

    title(main, line = 2.5, ...)



    # Adjust the margins to make sure the labels fit
    labwidth <- ifelse(labsize > 0, max(strwidth(rownames(nbeta[!remove,]), "inches", labsize)), 0)

    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))

    margins <- par("mai")
    par("mai" = c(margins[1:3], max(margins[4], labwidth*1.4)))
    if ( labsize > 0 && !is.null(rownames(nbeta)) )
    {
        take <- which(!remove)
        for (i in 1:sum(!remove)) {
            j <- take[i]
            axis(4, at = nbeta[j, ncol(nbeta)], labels = rownames(nbeta)[j],
                 las=1, cex.axis=labsize, col.axis = colseq[i],
                 lty = (i - 1) %% 5 + 1, col = colseq[i], ...)
        }
    }
    par("mai"=margins)
}




#' @param sign.lambda Either plot against log(lambda) (default) or its negative if \code{sign.lambda = -1}.
#' @rdname plot
#' @method plot cv.hd2part
#' @export
#' @examples
#' set.seed(123)
#'
plot.cv.hd2part <- function(x, sign.lambda = 1, ...)
{
    # modified from glmnet
    object = x
    num.models <- length(object$cvm)

    main.txt <- "Overall Outcomes"

    xlab <- expression(log(lambda))
    if(sign.lambda<0)xlab=paste("-",xlab,sep="")
    plot.args=list(x    = sign.lambda * log(object$lambda),
                   y    = object$cvm,
                   ylim = range(object$cvup, object$cvlo),
                   xlab = xlab,
                   ylab = object$name,
                   type = "n")
    new.args=list(...)
    if(length(new.args))plot.args[names(new.args)]=new.args
    do.call("plot", plot.args)
    error.bars(sign.lambda * log(object$lambda),
               object$cvup,
               object$cvlo, width = 0.005)

    nzero.txt <- paste(object$nonzero_z + object$nonzero_s)
    points(sign.lambda*log(object$lambda), object$cvm, pch=20, col="dodgerblue")
    axis(side=3,at=sign.lambda*log(object$lambda), labels = nzero.txt, tick=FALSE, line=0, ...)
    abline(v = sign.lambda * log(object$lambda.min), lty=2, lwd = 2, col = "firebrick1")
    abline(v = sign.lambda * log(object$lambda.1se), lty=2, lwd = 2, col = "firebrick1")
    title(main.txt, line = 2.5, ...)
    invisible()
}

# taken from glmnet
error.bars <- function(x, upper, lower, width = 0.02, ...)
{
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, col = 8, lty = 5, lwd = 0.5, ...)
    segments(x - barw, upper, x + barw, upper, col = "grey50", lwd = 1, ...)
    segments(x - barw, lower, x + barw, lower, col = "grey50", lwd = 1, ...)
    range(upper, lower)
}
