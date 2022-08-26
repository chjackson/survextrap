
# Stop without printing call
stop2    <- function(...) stop(..., call. = FALSE)

# Check whether a vector/matrix/array contains an "(Intercept)"
check_for_intercept <- function(x, logical = FALSE) {
  nms <- if (is.matrix(x)) colnames(x) else names(x)
  sel <- which("(Intercept)" %in% nms)
  if (logical) as.logical(length(sel)) else sel
}

# Drop intercept from a vector/matrix/array of named coefficients
drop_intercept <- function(x) { 
  sel <- check_for_intercept(x)
  if (length(sel) && is.matrix(x)) {
    x[, -sel, drop = FALSE]
  } else if (length(sel)) {
    x[-sel]
  } else {
    x
  }
}

warning2 <- function(...) warning(..., immediate. = TRUE, call. = FALSE)

