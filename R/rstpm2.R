## TODO bundle vuniroot
## needs Rcpp

## (c) Mark Clements. From the rstpm2 package, used under GPL.
## https://raw.githubusercontent.com/mclements/rstpm2/refs/heads/master/R/vuniroot.R
## 2025-05-02

vuniroot <- 
    function (f, interval, ..., lower=pmin(interval[,1], interval[,2]),
              upper=pmax(interval[,1], interval[,2]),
              f.lower = f(lower, ...), f.upper = f(upper, ...),
              extendInt = c("no", "yes", "downX", "upX"),
              check.conv = FALSE, tol = .Machine$double.eps^0.25,
              maxiter = 1000, trace = 0, n=NULL)
{
    if (!missing(interval) && ncol(interval) != 2L)
        stop("'interval' must be a matrix with two columns")
    if (all(!is.numeric(lower) | !is.numeric(upper) | lower >= upper))
        stop("lower < upper  is not fulfilled")
    if (is.null(n) && length(lower) == 1 && length(upper) == 1 && length(f.lower) > 1)
        n <- length(f.lower)
    ## warning("Bounds have length 1: possibly replicate for the bounds to have the correct length")
    if (!is.null(n) && length(lower) == 1 && length(upper) == 1) {
        lower = rep(lower,n)
        upper = rep(upper,n)
        f.lower = f(lower, ...)
        f.upper = f(upper, ...)
    }
    if (any(is.na(f.lower)))
        stop("f.lower = f(lower) is NA at least once")
    if (any(is.na(f.upper)))
        stop("f.upper = f(upper) is NA at least once")
    Sig <- switch(match.arg(extendInt), yes = NULL, downX = -1,
        no = 0, upX = 1, stop("invalid 'extendInt'; please report"))
    truncate <- function(x) pmax.int(pmin(x, .Machine$double.xmax),
                                     -.Machine$double.xmax)
    NAF <- function(x) ifelse(is.na(x),FALSE,x)
    f.low. <- truncate(f.lower)
    f.upp. <- truncate(f.upper)
    fun <- function(x) f(x, ...)
    doX <- (is.null(Sig) && any(f.low. * f.upp. > 0, na.rm=TRUE) ||
            (is.numeric(Sig) &&
             (any(Sig * f.low. > 0, na.rm=TRUE) || any(Sig * f.upp. < 0, na.rm=TRUE))))
    if (doX) {
        if (trace)
            cat(sprintf("search in [%g,%g]%s", lower, upper,
                if (trace >= 2)
                  "\n"
                else " ... "))
        Delta <- function(l,u) 0.01 * cbind(pmax(1e-04, abs(l)), pmax(1e-04, abs(u)))
        it <- 0L
        delta <- Delta(lower, upper)
        if (is.null(Sig)) {
            while (any(i <- NAF(f.lower * f.upper > 0)) &&
                   (any((iFl <- is.finite(lower)) | (iFu <- is.finite(upper))))) {
                       if ((it <- it + 1L) > maxiter)
                           stop(gettextf("no sign change found in %d iterations",
                                         it - 1), domain = NA)
                       if (any(j <- iFl & i)) {
                           ol <- lower[j]
                           of <- f.lower[j]
                           lower[j] <- lower[j] - delta[j,1]
                           f.lower[j] <- fun(lower)[j]
                           if (any(k <- is.na(f.lower[j]))) {
                               lower[j][k] <- ol[k]
                               f.lower[j][k] <- of[k]
                               delta[j,1][k] <- delta[j,1][k]/4
                           }
                       }
                       if (any(j <- iFu & i)) {
                           ol <- upper[j]
                           of <- f.upper[j]
                           upper[j] <- upper[j] + delta[j,2]
                           f.upper[j] <- fun(upper)[j]
                           if (any(k <- is.na(f.upper[j]))) {
                               upper[j][k] <- ol[k]
                               f.upper[j][k] <- of[k]
                               delta[j,2][k] <- delta[j,2][k]/4
                           }
                       }
                       if (trace >= 2)
                           cat(sprintf(" .. modified lower,upper: (%15g,%15g)\n",
                                       lower, upper))
                       delta[i & (iFl | iFu)] <- 2 * delta[i & (iFl | iFu)]
                   }
        }
        else {
            while (any(i <- NAF(Sig * f.lower > 0))) {
                if ((it <- it + 1L) > maxiter)
                    stop(gettextf("no sign change found in %d iterations",
                                  it - 1), domain = NA)
                lower[i] <- lower[i]-delta[i,1]
                f.lower[i] <- fun(lower)[i]
                if (trace >= 2)
                    cat(sprintf(" .. modified lower: %g\n", lower))
                delta[i] <- 2 * delta[i]
            }
            while (any(i <- NAF(Sig * f.upper < 0))) {
                if ((it <- it + 1L) > maxiter)
                    stop(gettextf("no sign change found in %d iterations",
                                  it - 1), domain = NA)
                upper[i] <- upper[i] + delta[i,2]
                f.upper[i] <- f(upper, ...)[i]
                if (trace >= 2)
                    cat(sprintf(" .. modified upper: %g\n", upper))
                delta[i] <- 2 * delta[i]
            }
        }
        if (trace && trace < 2)
            cat(sprintf("extended to [%g, %g] in %d steps\n",
                lower, upper, it))
    }
    if (any(!NAF(as.vector(f.lower * f.upper <= 0))))
        stop(if (doX)
                 "did not succeed extending the interval endpoints for f(lower) * f(upper) <= 0"
             else "f() values at end points not of opposite sign")
    if (check.conv) {
        val <- tryCatch(vunirootRcpp(fun, lower, upper, f.lower, f.upper, maxiter, tol),
                        warning = function(w) w)
        if (inherits(val, "warning")) 
            stop("convergence problem in zero finding: ", conditionMessage(val))
    }
    else {
        val <- vunirootRcpp(fun, lower, upper, f.lower, f.upper, as.integer(maxiter), tol)
    }
    iter <- as.integer(val[[2L]])
    if (any(iter < 0)) {
        (if (check.conv) 
             stop
         else warning)(sprintf(ngettext(maxiter, "_NOT_ converged in %d iteration",
                                        "_NOT_ converged in %d iterations"), maxiter), domain = NA)
        iter[iter<0] <- maxiter
    }
    if (doX)
        iter <- iter + it
    else it <- NA_integer_
    list(root = val[[1L]], f.root = f(val[[1L]], ...), iter = iter,
         init.it = it, estim.prec = val[[3L]])
}

vunirootRcpp <- function(f, lower, upper, fa, fb, numiter, tol) {
    .Call(`_survextrap_vunirootRcpp`, f, lower, upper, fa, fb, numiter, tol)
}

## then in src/RcppExports.cpp

## Rcpp::List vunirootRcpp(Rcpp::Function f, Rcpp::NumericVector lower, Rcpp::NumericVector upper, Rcpp::NumericVector fa, Rcpp::NumericVector fb, int numiter, double tol);
## RcppExport SEXP _survextrap_vunirootRcpp(SEXP fSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP faSEXP, SEXP fbSEXP, SEXP numiterSEXP, SEXP tolSEXP) {
## BEGIN_RCPP
##     Rcpp::RObject rcpp_result_gen;
##     Rcpp::RNGScope rcpp_rngScope_gen;
##     Rcpp::traits::input_parameter< Rcpp::Function >::type f(fSEXP);
##     Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lower(lowerSEXP);
##     Rcpp::traits::input_parameter< Rcpp::NumericVector >::type upper(upperSEXP);
##     Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fa(faSEXP);
##     Rcpp::traits::input_parameter< Rcpp::NumericVector >::type fb(fbSEXP);
##     Rcpp::traits::input_parameter< int >::type numiter(numiterSEXP);
##     Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
##     rcpp_result_gen = Rcpp::wrap(vunirootRcpp(f, lower, upper, fa, fb, numiter, tol));
##     return rcpp_result_gen;
## END_RCPP
## }
