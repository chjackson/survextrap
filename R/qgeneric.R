## Generic function to find quantiles of a distribution
##
## Generic function to find the quantiles of a distribution, given the
## equivalent probability distribution function.   Copied from the version
## distributed in the flexsurv package. 
##
## This function is used by default for custom distributions for which a
## quantile function is not provided.
##
## It works by finding the root of the equation \eqn{h(q) = pdist(q) -
## p = 0}.  Starting from the interval \eqn{(-1, 1)}, the interval
## width is expanded by 50\% until \eqn{h()} is of opposite sign at
## either end.  The root is then found using numerical methods
## (\code{vuniroot} from the `rstpm2` package by Mark Clements).
##
## This assumes a suitably smooth, continuous distribution.
##
## @param pdist Probability distribution function, for example,
## \code{\link{pnorm}} for the normal distribution, which must be defined in
## the current workspace.  This should accept and return vectorised parameters
## and values.  It should also return the correct values for the entire real
## line, for example a positive distribution should have \code{pdist(x)==0}
## for \eqn{x<0}.
##
## @param p Vector of probabilities to find the quantiles for.
##
## @param matargs Character vector giving the elements of \code{...} which
## represent vector parameters of the distribution.  Empty by default.  When
## vectorised, these will become matrices.  
##
## @param scalarargs Character vector naming arguments of the distribution function that cannot be vectorised.
##
## @param ...  The remaining arguments define parameters of the distribution
## \code{pdist}.  These MUST be named explicitly.
##
## This may also contain the standard arguments \code{log.p}, and
## \code{lower.tail} (as used in, e.g. \code{\link{qnorm}})
##
## If the distribution is bounded above or below, then this should contain
## arguments \code{lbound} and \code{ubound} respectively, and these will be
## returned if \code{p} is 0 or 1 respectively.  Defaults to \code{-Inf} and
## \code{Inf} respectively.
##
## @return Vector of quantiles of the distribution at \code{p}.
##
## @author Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
##
## @keywords internal distribution
qgeneric <- function(pdist, p, matargs=NULL, scalarargs=NULL, ...)
{
    args <- list(...)
    if (is.null(args$log.p)) args$log.p <- FALSE
    if (is.null(args$lower.tail)) args$lower.tail <- TRUE
    if (is.null(args$lbound)) args$lbound <- -Inf
    if (is.null(args$ubound)) args$ubound <- Inf
    if (args$log.p) p <- exp(p)
    if (!args$lower.tail) p <- 1 - p
    ret <- numeric(length(p))
    ret[p == 0] <- args$lbound
    ret[p == 1] <- args$ubound
    ## args containing vector params of the distribution (e.g. gamma and knots in dsurvspline)
    args.mat <- args[matargs]
    ## Arguments that cannot be vectorised
    args.scalar <- args[scalarargs]
    args[c(matargs,scalarargs,"lower.tail","log.p","lbound","ubound")] <- NULL
    ## Other args assumed to contain vectorisable parameters of the distribution.
    ## Replicate all to their maximum length, along with p
    matlen <- if(is.null(matargs)) NULL else sapply(args.mat, function(x){if(is.matrix(x))nrow(x) else 1})
    veclen <- if (length(args) == 0) NULL else sapply(args, length)
    maxlen <- max(c(length(p), veclen, matlen))
    na_inds <- rep(FALSE, length(ret))
    for (i in seq(along=args)){
        args[[i]] <- rep(args[[i]], length.out=maxlen)
        na_inds <- na_inds | is.na(args[[i]])
    }
    for (i in seq(along=args.mat)){
        if (is.matrix(args.mat[[i]])){
            args.mat[[i]] <- matrix(
              apply(args.mat[[i]], 2, function(x)rep(x, length=maxlen)),
              ncol=ncol(args.mat[[i]]),
              byrow=F
            )
        }
        else args.mat[[i]] <- matrix(args.mat[[i]], nrow=maxlen, ncol=length(args.mat[[i]]), byrow=TRUE)
        na_inds <- na_inds | apply(args.mat[[i]], 1, function(x)any(is.na(x)))
    }
    p <- rep(p, length.out=maxlen)
    ret[p < 0 | p > 1] <- NaN
    ret[na_inds] <- NA
    ind <- (p > 0 & p < 1 & !na_inds)
    if (any(ind)) {
        hind <- seq(along=p)[ind]
        n <- length(p[ind])
        ptmp <- numeric(n)
        interval <- matrix(rep(c(-1, 1), n), ncol=2, byrow=TRUE)
        h <- function(y) {
            args <- lapply(args, function(x)x[hind])
            args.mat <- lapply(args.mat, function(x)x[hind,])
            p <- p[hind]
            args$q <- y
            args <- c(args, args.mat, args.scalar)
            (do.call(pdist, args) - p)
        }
        ptmp <- vuniroot(h, interval, tol=.Machine$double.eps, extendInt="yes", maxiter=10000)$root
        ret[ind] <- ptmp
    }
    if (any(is.nan(ret))) warning("NaNs produced")
    ret
}
