
# Return the cutpoints for a specified number of quantiles of 'x'
#
# @param x A numeric vector.
# @param nq Integer specifying the number of quantiles.
# @return A vector of percentiles corresponding to percentages 100*k/m for 
#   k=1,2,...,nq-1.
qtile <- function(x, nq = 2) {
  if (nq > 1) {
    probs <- seq(1, nq - 1) / nq
    return(quantile(x, probs = probs))
  } else if (nq == 1) {
    return(NULL)
  } else {
    stop("'nq' must be >= 1.")
  }
}

# Throw error if parameter isn't a positive scalar
#
# @param x The object to test.
validate_positive_scalar <- function(x, not_greater_than = NULL) {
  nm <- deparse(substitute(x))
  if (is.null(x))
    stop(nm, " cannot be NULL", call. = FALSE)
  if (!is.numeric(x))
    stop(nm, " should be numeric", call. = FALSE)
  if (any(x <= 0)) 
    stop(nm, " should be postive", call. = FALSE)
  if (!is.null(not_greater_than)) {
    if (!is.numeric(not_greater_than) || (not_greater_than <= 0))
      stop("'not_greater_than' should be numeric and postive")
    if (!all(x <= not_greater_than))
      stop(nm, " should less than or equal to ", not_greater_than, call. = FALSE)
  }
}



# Methods for creating linear predictor
#
# Make linear predictor vector from x and point estimates for beta, or linear 
# predictor matrix from x and full posterior sample of beta.
#
# @param beta A vector or matrix or parameter estimates.
# @param x Predictor matrix.
# @param offset Optional offset vector.
# @return A vector or matrix.
linear_predictor <- function(beta, x, offset = NULL) {
  UseMethod("linear_predictor")
}
linear_predictor.default <- function(beta, x, offset = NULL) {
  eta <- as.vector(if (NCOL(x) == 1L) x * beta else x %*% beta)
  if (length(offset))
    eta <- eta + offset
  
  return(eta)
}
linear_predictor.matrix <- function(beta, x, offset = NULL) {
  if (NCOL(beta) == 1L) 
    beta <- as.matrix(beta)
  eta <- beta %*% t(x)
  if (length(offset)) 
    eta <- sweep(eta, 2L, offset, `+`)

  return(eta)
}

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


# Reformulate as LHS of a formula
#
# @param x A character string or expression.
# @return A formula.
reformulate_lhs <- function(x) {
  x <- formula(substitute(LHS ~ 1, list(LHS = x)))
  x
}

# Reformulate as RHS of a formula
#
# @param x A character string or expression.
# @return A formula.
reformulate_rhs <- function(x) {
  x <- formula(substitute(~ RHS, list(RHS = x)))
  x
}


# Return a data frame with NAs excluded
#
# @param formula The parsed model formula.
# @param data The (user-specified) data frame.
# @return A data frame, with only complete cases for the variables that
#   appear in the model formula.
make_model_data <- function(formula, data) {
  mf <- model.frame(formula, data, na.action = na.pass)
  include <- apply(mf, 1L, function(row) !any(is.na(row)))
  data[include, , drop = FALSE]
}

# Paste character vector collapsing with a comma
comma <- function(x) {
  paste(x, collapse = ", ")
}
