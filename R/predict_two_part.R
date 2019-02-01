## the code here is largely based on the code
## from the glmnet package

#' Prediction method for two part fitted objects
#'
#' @param object fitted "2part" model object
#' @param newx Matrix of new values for \code{x} at which predictions are to be made. Must be a matrix
#' This argument is not used for \code{type=c("coefficients","nonzero")}
#' @param s Value(s) of the penalty parameter lambda for the zero part at which predictions are required. Default is the entire sequence used to create
#' the model.
#' @param model either \code{"zero"} for the zero part model or \code{"positive"} for the positive part model
#' @param type Type of prediction required. \code{type = "link"} gives the linear predictors;
#' \code{type = "model_response"} gives the fitted probabilities for the zero part and fitted expected values for the positive part.
#' \code{type = "response"} gives the combined response prediction across the two models using the full unconditional expected
#' value of the response. When \code{type = "response"}, argument \code{"model"} is unused.
#' \code{type = "coefficients"} computes the coefficients at the requested values for \code{s}.
#' @param ... not used
#' @return An object depending on the type argument
#' @export
#' @examples
predict.2part <- function(object, newx, s = NULL,
                          model = c("zero", "positive"),
                          type = c("link",
                                   "model_response",
                                   "response",
                                   "coefficients",
                                   "nonzero"), ...)
{
    type  <- match.arg(type)
    model <- match.arg(model)

    if(missing(newx))
    {
        if(!match(type, c("coefficients", "nonzero"), FALSE)) stop("A value for 'newx' must be supplied")
    }

    if (type != "response")
    {
        if (model == "zero")
        {
            nbeta <- object$beta_z
        } else
        {
            nbeta <- object$beta_s
        }

        if(!is.null(s))
        {
            #vnames=dimnames(nbeta)[[1]]
            lambda  <- object$lambda
            lamlist <- lambda.interp(object$lambda, s)
            nbeta   <- nbeta[,lamlist$left,drop=FALSE]*lamlist$frac + nbeta[,lamlist$right,drop=FALSE]*(1-lamlist$frac)
            #dimnames(nbeta)=list(vnames,paste(seq(along=s)))
        }
        if (type == "coefficients") return(nbeta)
        if (type == "nonzero")
        {
            nbeta[1,] <- 0 ## rem intercept
            newbeta <- abs(as.matrix(nbeta)) > 0
            index <- 1:(dim(newbeta)[1])
            nzel <- function(x, index) if(any(x)) index[x] else NULL
            betaList <- apply(newbeta, 2, nzel, index)
            return(betaList)
        }

        newx <- as.matrix(newx)
        # add constant column if needed
        if (ncol(newx) < nrow(nbeta))
        {
            newx <- cbind(rep(1, nrow(newx)), newx)
        }

        eta <- as.matrix(newx %*% nbeta)
        if (type == "link")
        {
            return(eta)
        } else if (type == "model_response")
        {
            if (model == "zero")
            {
                return(1 / (1 + exp(-eta)))
            } else
            {
                return(exp(eta))
            }
        }
    } else
    {
        nbeta_z <- object$beta_z
        nbeta_s <- object$beta_s

        if(!is.null(s))
        {
            #vnames=dimnames(nbeta)[[1]]
            lambda  <- object$lambda
            lamlist <- lambda.interp(object$lambda, s)
            nbeta_z <- nbeta_z[,lamlist$left,drop=FALSE] * lamlist$frac + nbeta_z[,lamlist$right,drop=FALSE] * (1 - lamlist$frac)
            nbeta_s <- nbeta_s[,lamlist$left,drop=FALSE] * lamlist$frac + nbeta_s[,lamlist$right,drop=FALSE] * (1 - lamlist$frac)
            #dimnames(nbeta)=list(vnames,paste(seq(along=s)))
        }

        newx <- as.matrix(newx)
        # add constant column if needed
        if (ncol(newx) < nrow(nbeta_z))
        {
            newx <- cbind(rep(1, nrow(newx)), newx)
        }

        eta_z <- as.matrix(newx %*% nbeta_z)
        eta_s <- as.matrix(newx %*% nbeta_s)

        return(exp(eta_s) * as.vector(  (1 / (1 + exp(-eta_z)))  )  )
    }
}

