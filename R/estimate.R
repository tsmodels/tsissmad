arma_conditions <- function(pars, parnames) {
    if (any(grepl("theta",parnames))) {
        ar <- pars[grepl("theta",parnames)]
        if (all(ar == 0)) {
            
        } else {
            if (min(Mod(polyroot(c(1, -ar)))) < 1.01) {
                return(TRUE)
            }
        }
    }
    if (any(grepl("psi",parnames))) {
        ma <- pars[grepl("psi",parnames)]
        if (all(ma == 0)) {
            
        } else {
            if (min(Mod(polyroot(c(1, ma)))) < 1.01) {
                return(TRUE)
            }
        }
    }
    return(FALSE)
}


prepare_inputs_issm <- function(spec)
{
    estimate <- NULL
    parameters <- NULL
    S <- spec$S
    S <- S[matrix %in% c("F0","F1","F2","g","w")]
    S <- S[matrix != "xreg"]
    V <- S$values
    V[which(is.na(V))] <- 0
    findex <- which(!is.na(S$pars)) - 1
    p <- spec$parmatrix
    p <- p[!grepl("kappa",parameters)]
    p <- p[!grepl("lambda",parameters)]
    allpars <- p$initial
    parnames_all <- p$parameters
    pars <- p[estimate == 1]
    parnames_estimate <- pars$parameters
    tmb_names <- rep("pars", length(parnames_estimate))
    lower <- pars$lower
    upper <- pars$upper
    ppindex <- match(pars$parameters, p$parameters) - 1
    pars <- pars$initial
    kappa <- spec$parmatrix[grepl("kappa",parameters)]$initial
    lambda <- spec$parmatrix[grepl("lambda",parameters)]$initial
    fpindex <- match(na.omit(S$pars), p$parameters) - 1
    f0_index <- range(which(S$matrix == "F0"))
    f0_index <- c(f0_index[1] - 1, f0_index[2] - f0_index[1] + 1)
    f1_index <- range(which(S$matrix == "F1"))
    f1_index <- c(f1_index[1] - 1, f1_index[2] - f1_index[1] + 1)
    f2_index <- range(which(S$matrix == "F2"))
    f2_index <- c(f2_index[1] - 1, f2_index[2] - f2_index[1] + 1)
    modeli <- spec$dims[1]
    windex <-  range(which(S$matrix == "w"))
    windex <- c(windex[1] - 1, windex[2] - windex[1] + 1)
    gindex <-  range(which(S$matrix == "g"))
    gindex <- c(gindex[1] - 1, gindex[2] - gindex[1] + 1)
    fshape <- c(f0_index, f1_index, f2_index, windex, gindex)
    y <- spec$target$y
    X <- spec$xreg$xreg
    map <- list()
    if (!spec$transform$include_lambda) {
        map$lambda <- factor(NA)
    } else {
        lower <- c(lower, spec$parmatrix[parameters == "lambda"]$lower)
        upper <- c(upper, spec$parmatrix[parameters == "lambda"]$upper)
        parnames_estimate <- c(parnames_estimate, "lambda")
        tmb_names <- c(tmb_names, "lambda")
    }
    lambda <- spec$transform$lambda
    if (!spec$xreg$include_xreg) {
        map$kappa <- factor(NA)
        kappa <- 0
    } else {
        # check for fixed pars
        if (any(spec$parmatrix[grepl("kappa",parameters)]$estimate == 0)) {
            xn <- NROW(spec$parmatrix[grepl("kappa",parameters)])
            ix <- which(spec$parmatrix[grepl("kappa",parameters)]$estimate == 0)
            kappa_fixed <- rep(1:xn)
            kappa_fixed[ix] <- factor(NA)
            kappa_fixed <- as.factor(kappa_fixed)
            map$kappa <- kappa_fixed
        }
        if (any(spec$parmatrix[grepl("kappa",parameters)]$estimate == 1)) {
            lower <- c(lower, spec$parmatrix[grepl("kappa",parameters) & estimate == 1]$lower)
            upper <- c(upper, spec$parmatrix[grepl("kappa",parameters) & estimate == 1]$upper)
            parnames_estimate <- c(parnames_estimate, spec$parmatrix[grepl("kappa",parameters) & estimate == 1]$parameters)
            tmb_names <- c(tmb_names, rep("kappa",nrow(spec$parmatrix[grepl("kappa",parameters) & estimate == 1])))
        }
        kappa <- spec$parmatrix[grepl("kappa",parameters)]$initial
    }
    # check for fixed lambda
    # check for NA values
    good <- rep(1, NROW(y))
    if (any(is.na(y))) {
        good[which(is.na(y))] <- 0
        y <- na.fill(y, fill = 0)
    }
    modeli <- c(modeli, spec$dims[5])
    # create function for ARMA and non ARMA models
    if (sum(spec$arma$order) > 0) {
        llh_fun <- function(pars, fun, issmenv) {
            if (arma_conditions(pars, issmenv$parnames_estimate)) {
                lik <- issmenv$lik + 0.25 * abs(issmenv$lik)
                issmenv$lik <- lik
            } else {
                names(pars) <- issmenv$tmb_names
                lik <- fun$fn(pars)
                D <- abs(Re(eigen(fun$report()$D, only.values = TRUE)$values))
                if (is.na(lik) | any(D > 1.01) | !is.finite(lik)) {
                    lik <- issmenv$lik + 0.25 * abs(issmenv$lik)
                    issmenv$lik <- lik
                } else {
                    issmenv$lik <- lik
                }
            }
            return(lik)
        }
    } else {
        llh_fun <- function(pars, fun, issmenv) {
            names(pars) <- issmenv$tmb_names
            lik <- fun$fn(pars)
            D <- abs(Re(eigen(fun$report()$D, only.values = TRUE)$values))
            if (is.na(lik) | any(D > 1.01) | !is.finite(lik)) {
                lik <- issmenv$lik + 0.25 * abs(issmenv$lik)
                issmenv$lik <- lik
            } else {
                issmenv$lik <- lik
            }
            return(lik)
        }
    }
    grad_fun <- function(pars, fun, issmenv)
    {
        names(pars) <- issmenv$tmb_names
        fun$gr(pars)
    }
    
    hess_fun <- function(pars, fun, issmenv)
    {
        names(pars) <- issmenv$tmb_names
        fun$he(pars)
    }
    
    data <- list(V = V, X = X, good = good, y = y, allpars = allpars, 
                 findex = findex, fpindex = fpindex, ppindex = ppindex, fshape = fshape, 
                 modeli = modeli)
    par_list <- list(pars = pars, kappa = kappa, lambda = lambda)
    L <- list(data = data, par_list = par_list, map = map, lower = lower, upper = upper, 
              llh_fun = llh_fun, grad_fun = grad_fun, hess_fun = hess_fun, 
              parnames_estimate = parnames_estimate, tmb_names = tmb_names, parnames_all = parnames_all)
    return(L)
}

#' Estimates an ISSM model given a specification object using maximum likelihood and autodiff
#'
#' @param object An object of class tsissm.spec.
#' @param solver Only \dQuote{nlminb} currently supported.
#' @param use_hessian Whether to include the hessian in the calculation. This is currently
#' problematic and we suggest not overriding the default (FALSE) until further investigation.
#' @param control Solver control parameters.
#' @param ... additional parameters passed to the estimation function
#' @return An list of coefficients and other information.
#' @details This function is not expected to be used by itself but rather as a plugin
#' to be called from the estimate method of the tsissm package.
#' @export estimate_ad.tsissm.spec
#' @aliases estimate_ad
#' @export
#' 
estimate_ad.tsissm.spec <- function(object, solver = "nlminb", control = list(trace = 0, eval.max = 300, iter.max = 500), use_hessian = FALSE, ...)
{
    parameters <- NULL
    
    spec_list <- prepare_inputs_issm(object)
    other_opts <- list(...)
    if (!is.null(other_opts$silent)) {
        silent <- other_opts$silent
    } else {
        silent <- TRUE
    }
    if (object$transform$name == "logit") {
        spec_list$data$model <- "issmb"
    } else if (object$transform$name == "box-cox") {
        spec_list$data$model <- "issma"
    } else {
        stop("\nunknown transformation used")
    }
    fun <- try(MakeADFun(data = spec_list$data, hessian = use_hessian, parameters = spec_list$par_list, DLL = "tsissmad_TMBExports", 
                         map = spec_list$map, trace = FALSE, silent = silent, checkParameterOrder = F), silent = FALSE)
    fun$env$tracemgc <- FALSE
    
    if (inherits(fun, 'try-error')) {
        stop("\nestimate_ad found an error. Please use non ad version of estimator and contact developer with reproducible example.")
    }
    issmenv <- new.env()
    issmenv$lik <- 1
    issmenv$grad <- NULL
    issmenv$parameter_names <- spec_list$parnames_all
    issmenv$parnames_estimate <- spec_list$parnames_estimate
    issmenv$tmb_names <- spec_list$tmb_names
    issmenv$parmatrix <- object$parmatrix
    ## if (solver != "nlminb") warning("\nonly nlminb solver currently supported for issm with autodiff. Using nlminb.")
    if (use_hessian) hessian <- spec_list$hess_fun else hessian <- NULL
    # scf <- scale_kappa(object)
    # if (scf$use_scaling) {
    #     tmp <- object$parmatrix
    #     tmp[grepl("kappa",parameters), scale := as.numeric(scf$scale_factors)]
    #     scf <- tmp[estimate == 1]$scale
    # } else {
    #     scf <- rep(1, length(fun$par))
    # }
    if (solver == "nlminb") {
        sol <- nlminb(start = fun$par, objective = spec_list$llh_fun, 
                      gradient = spec_list$grad_fun, hessian = hessian,  
                      lower = spec_list$lower, upper = spec_list$upper,  control = control, 
                      fun = fun, issmenv = issmenv)
    } else {
        sol <- optim(par = fun$par, fn = spec_list$llh_fun, gr = spec_list$grad_fun, 
                      lower = spec_list$lower, upper = spec_list$upper,  control = control, 
                      method = "L-BFGS-B", fun = fun, issmenv = issmenv)
    }
    pars <- sol$par
    names(pars) <- issmenv$tmb_names
    llh <- spec_list$llh_fun(pars, fun, issmenv)
    gradient <- spec_list$grad_fun(pars, fun, issmenv)
    hessian <- spec_list$hess_fun(pars, fun, issmenv)
    names(pars) <- issmenv$parnames_estimate
    parmatrix <- object$parmatrix
    parmatrix[parameters %in% issmenv$estimation_names]$initial <- pars
    colnames(gradient) <- issmenv$estimation_names
    colnames(hessian) <- rownames(hessian) <- issmenv$estimation_names
    xseed <- fun$report()$states[1,,drop = FALSE]
    out <- list(pars = pars, llh = llh, gradient = gradient, hessian = hessian, xseed = xseed, solver_out = sol)
    return(out)
}



scale_kappa <- function(object)
{
    if (object$transform$include_lambda) {
        use_scaling <- FALSE
        scale_factors <- 1
    } else {
        if (object$xreg$include_xreg) {
            use_scaling <- TRUE
            yt <- object$transform$transform(object$target$y_orig)
            scale_factors <- 1/(max(yt)/apply(object$xreg$xreg, 2, max))
        } else {
            use_scaling <- FALSE
            scale_factors <- 1
        }
    }
    return(list(use_scaling = use_scaling, scale_factors = scale_factors))
}
