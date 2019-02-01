## the code here is largely based on the code
## from the glmnet package

#' Prediction method for two part fitted objects
#'
#' @param object fitted "hd2part" model object
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
#' @param newoffset f an offset is used in the fit, then one must be supplied for making predictions
#' @param ... not used
#' @return An object depending on the type argument
#' @export
#' @examples
predict.hd2part <- function(object, newx, s = NULL,
                            model = c("zero", "positive"),
                            type = c("link",
                                     "model_response",
                                     "response",
                                     "coefficients",
                                     "nonzero"),
                            newoffset = NULL,
                            ...)
{
    type  <- match.arg(type)
    model <- match.arg(model)

    if(missing(newx))
    {
        if(!match(type, c("coefficients", "nonzero"), FALSE)) stop("A value for 'newx' must be supplied")
    }

    if (!match(type, c("coefficients", "nonzero"), FALSE))
    {

        if (is.null(newoffset))
        {
            newoffset <- rep(0, NROW(newx))
        }

        if (NROW(newx) != NROW(newoffset))
        {
            stop("newoffset must be same length as the number of observations in newx")
        }
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
        eta <- eta + array(newoffset, dim = dim(eta))
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

        eta_z <- eta_z + array(newoffset, dim = dim(eta_z))
        eta_s <- eta_s + array(newoffset, dim = dim(eta_s))

        return(exp(eta_s) * as.vector(  (1 / (1 + exp(-eta_z)))  )  )
    }
}



## taken from glmnet

lambda.interp=function(lambda,s){
    ### lambda is the index sequence that is produced by the model
    ### s is the new vector at which evaluations are required.
    ### the value is a vector of left and right indices, and a vector of fractions.
    ### the new values are interpolated bewteen the two using the fraction
    ### Note: lambda decreases. you take:
    ### sfrac*left+(1-sfrac*right)

    if(length(lambda)==1){# degenerate case of only one lambda
        nums=length(s)
        left=rep(1,nums)
        right=left
        sfrac=rep(1,nums)
    }
    else{
        s[s > max(lambda)] = max(lambda)
        s[s < min(lambda)] = min(lambda)
        k=length(lambda)
        sfrac <- (lambda[1]-s)/(lambda[1] - lambda[k])
        lambda <- (lambda[1] - lambda)/(lambda[1] - lambda[k])
        coord <- approx(lambda, seq(lambda), sfrac)$y
        left <- floor(coord)
        right <- ceiling(coord)
        sfrac=(sfrac-lambda[right])/(lambda[left] - lambda[right])
        sfrac[left==right]=1
        sfrac[abs(lambda[left]-lambda[right])<.Machine$double.eps]=1

    }
    list(left=left,right=right,frac=sfrac)
}

## taken from glmnet

nonzeroCoef = function (beta, bystep = FALSE)
{
    ### bystep = FALSE means which variables were ever nonzero
    ### bystep = TRUE means which variables are nonzero for each step
    nr=nrow(beta)
    if (nr == 1) {#degenerate case
        if (bystep)
            apply(beta, 2, function(x) if (abs(x) > 0)
                1
                else NULL)
        else {
            if (any(abs(beta) > 0))
                1
            else NULL
        }
    }
    else {
        beta=abs(beta)>0 # this is sparse
        which=seq(nr)
        ones=rep(1,ncol(beta))
        nz=as.vector((beta%*%ones)>0)
        which=which[nz]
        if (bystep) {
            if(length(which)>0){
                beta=as.matrix(beta[which,,drop=FALSE])
                nzel = function(x, which) if (any(x))
                    which[x]
                else NULL
                which=apply(beta, 2, nzel, which)
                if(!is.list(which))which=data.frame(which)# apply can return a matrix!!
                which
            }
            else{
                dn=dimnames(beta)[[2]]
                which=vector("list",length(dn))
                names(which)=dn
                which
            }

        }
        else which
    }
}


## taken from glmnet

auc <- function(y,prob,w)
{
    if(missing(w)){
        rprob=rank(prob)
        n1=sum(y);n0=length(y)-n1
        u=sum(rprob[y==1])-n1*(n1+1)/2
        exp(log(u) - log(n1) - log(n0))
    }
    else{
        rprob=runif(length(prob))
        op=order(prob,rprob)#randomize ties
        y=y[op]
        w=w[op]
        cw=cumsum(w)
        w1=w[y==1]
        cw1=cumsum(w1)
        wauc = log(sum(w1*(cw[y==1]-cw1)))
        sumw1 = cw1[length(cw1)]
        sumw2  = cw[length(cw)] - sumw1
        exp(wauc - log(sumw1) - log(sumw2))
    }
}

## taken from glmnet

auc.mat <- function(y,prob,weights=rep(1,nrow(y)))
{
    Weights=as.vector(weights*y)
    ny=nrow(y)
    Y=rep(c(0,1),c(ny,ny))
    Prob=c(prob,prob)
    auc(Y,Prob,Weights)
}

deviance.2part <- function(object,...)
{
    dev     <- object$deviance
    nulldev <- object$deviance[1]
    (1 - dev) * nulldev
}
