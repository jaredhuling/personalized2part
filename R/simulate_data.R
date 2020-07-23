

#' Generates data from a two part distribution with a point mass at zero and heterogeneous treatment effects
#'
#' @description Generates semicontinuous data with heterogeneity of treatment effect
#'
#' @param n.obs number of observations
#' @param n.vars number of variables. Must be at least 10
#' @return returns list with values \code{y} for outcome, \code{x} for design matrix, \code{trt} for
#' treatment assignments, \code{betanonzero} for true coefficients for treatment-covariate interactions for model for
#' whether or not a response is nonzero, \code{betapos} for true coefficients for treatment-covariate interactions
#' for positive model, \code{treatment_risk_ratio} for the true covariate-conditional treatment effect risk ratio for
#' each observation, \code{pi.x} for the true underlying propensity score
#' @export
sim_semicontinuous_data <- function(n.obs = 1000, n.vars = 25)
{

    if (n.vars < 10) stop("n.vars must be at least 10")

    x <- matrix(rnorm(n.obs * n.vars, sd = 1), n.obs, n.vars)

    trtbeta <- c(0.75, 0.5, -0.5, -0.25, -0.75, 0, 0, 0, 0, 0, rep(0, n.vars - 10))

    trtprob <- 1 / (1 + exp(-x %*% trtbeta))

    trt <- 2 * rbinom(n.obs, 1, trtprob) - 1

    betanonzero <- 0.5 * c(0.5, 0.35, 0.85, -0.5, -1.45, 0.75, 0, 0, 0, 0, 0, rep(0, n.vars - 10))

    beta     <- 0.5 * c(0.35, -0.5, 0.65, -0.95, -1.25, 1.05, 0.5, 0, 0, 0, 0, rep(0, n.vars - 10))

    ######## simulate response

    ## the treatment-covariate interactions for the zero part and the positive part
    delta     <- drop(cbind(1, x) %*% beta)
    deltazero <- drop(cbind(1, x) %*% betanonzero)

    ## main effects
    xbeta <- x[,1] + 0.5 * x[,10] + 0.5 * x[,9] + 0.25 * x[,8] + 0.5 * x[,7] + 0.15 * x[,7] ^ 2 + 0.5 * x[,2] - 0.5 * x[,5]

    ## add main effects to interactions
    xbeta <- 1 * xbeta + delta * trt

    ## generate the linear predictor portions of the means of the potential outcomes
    xbeta_1  <- 1 * xbeta + delta * 1
    xbeta_n1 <- 1 * xbeta - delta * 1

    ## now repeat the same for the positive part

    ## main effects
    xbetazero <- 0.15 * x[,1] - 0.5 * x[,10] - 0.5 * x[,9] - 0.25 * x[,7] + 0.25 * x[,6] ^ 2 + 0.2 * x[,2] + 0.2 * x[,5]

    zero_int <- -1
    xbetazero <- zero_int + 0.5 * xbetazero + deltazero * trt


    xbetazero_1  <- zero_int + 0.5 * xbetazero + deltazero * 1
    xbetazero_n1 <- zero_int + 0.5 * xbetazero - deltazero * 1



    nonzeroprob <- 1 / (1 + exp(-xbetazero))

    nonzeroprob_t1  <- 1 / (1 + exp(-xbetazero_1))
    nonzeroprob_tn1 <- 1 / (1 + exp(-xbetazero_n1))



    treatment_risk_ratio <- (exp((xbeta_1 - xbeta_n1) )) * (nonzeroprob_t1 / nonzeroprob_tn1)

    ## generate positive indicator
    ynzero <- rbinom(n.obs, 1, prob = nonzeroprob)

    #yzero <- 1 - 1 *(xbetazero + rnorm(NROW(xbetazero)) > 0)


    zeroprob_t1  <- 1 - 1 / (1 + exp(-xbetazero_1))
    zeroprob_tn1 <- 1 - 1 / (1 + exp(-xbetazero_n1))

    yzero_t1  <- 1 - rbinom(n.obs, 1, prob = zeroprob_t1)
    yzero_tn1 <- 1 - rbinom(n.obs, 1, prob = zeroprob_tn1)


    nu <- 25
    gam_mean <- 1
    yg_t1  <- rgamma(n.obs, shape = nu, rate = nu / exp(gam_mean + xbetazero_1) )
    yg_tn1 <- rgamma(n.obs, shape = nu, rate = nu / exp(gam_mean + xbetazero_n1) )

    y_t1  <- ifelse(ynzero == 0, 0, yg_t1)
    y_tn1 <- ifelse(ynzero == 0, 0, yg_tn1)



    ## generate a gamma-distributed outcome for the positive part

    yg <- rgamma(n.obs, shape = nu, rate = nu / exp(gam_mean + xbeta) )
    range(yg)
    #hist(yg, breaks = 50)

    ## create two part outcome
    y <- ifelse(ynzero == 0, 0, yg)

    list(y = y, x = x,
         trt = ifelse(trt == 1, 1, 0),
         betanonzero = betanonzero,
         betapos     = beta,
         treatment_risk_ratio = treatment_risk_ratio,          ### Gamma(X)
         pi.x = trtprob                                        ### Pr(T=1|X)
         )
}
