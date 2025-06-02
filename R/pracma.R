## copied from pracma package version 2.4.4 
## by Hans W. Borchers
## used under GPL>=3
gaussLegendre <- function(n, a, b) 
{
    stopifnot(is.numeric(a), length(a) == 1, is.numeric(b), length(b) == 
        1, is.numeric(n), length(n) == 1, n >= 2)
    i <- seq(1, n - 1, by = 1)
    d <- i/sqrt(4 * i^2 - 1)
    E <- eigen(Diag(d, 1) + Diag(d, -1), symmetric = TRUE)
    L <- E$values
    V <- E$vectors
    inds <- order(L)
    x <- L[inds]
    V <- t(V[, inds])
    w <- 2 * V[, 1]^2
    x <- 0.5 * ((b - a) * x + a + b)
    w <- -0.5 * (a - b) * w
    return(list(x = x, w = w))
}

Diag <- function(x, k = 0) 
{
    if (!is.numeric(x) && !is.complex(x)) 
        stop("Argument 'x' must be a real or complex vector or matrix.")
    if (!is.numeric(k) || k != round(k)) 
        stop("Argument 'k' must be an integer.")
    if (is.matrix(x)) {
        n <- nrow(x)
        m <- ncol(x)
        if (k >= m || -k >= n) {
            y <- matrix(0, nrow = 0, ncol = 0)
        }
        else {
            y <- x[col(x) == row(x) + k]
        }
    }
    else {
        if (is.vector(x)) {
            n <- length(x)
            m <- n + abs(k)
            y <- matrix(0, nrow = m, ncol = m)
            y[col(y) == row(y) + k] <- x
        }
        else {
            stop("Argument 'x' must be a real or complex vector or matrix.")
        }
    }
    return(y)
}
