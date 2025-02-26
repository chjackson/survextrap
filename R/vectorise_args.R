## Utilities for creating vectorised functions 

#' @param args a named list
#'
#' @param args_mat character vector giving components of args that are matrices
#'
#' @param args_unvectorised character vector giving components of args that we don't want to vectorise
#'
#' @return a list with components
#'
#' args   vectors
#' 
#' args_mat   matrices
#'
#' with (length) (number of rows) set to the maximum over all, and contents replicated
#'
#' @noRd
vectorise_args <- function(args, mat_args, unvectorised_args){
  margs <- args[mat_args]
  uargs <- args[unvectorised_args]
  args[c(mat_args,unvectorised_args)] <- NULL
  vargs <- args
  matlens <- if(is.null(mat_args)) NULL else sapply(margs, function(x){if(is.matrix(x))nrow(x) else 1})
  veclens <- if (length(vargs) == 0) NULL else sapply(vargs, length)
  maxlen <- max(c(veclens, matlens)) 
  na_inds <- rep(FALSE, maxlen)
  for (i in seq(along=args)){
    vargs[[i]] <- rep(vargs[[i]], length.out=maxlen)
    na_inds <- na_inds | is.na(vargs[[i]])
  }
  margs <- rep_matrix_list(margs, maxlen)
  for (i in seq(along=margs))
    na_inds <- na_inds | apply(margs[[i]], 1, function(x)any(is.na(x)))
  list(vargs=vargs, margs=margs, uargs=uargs, na_inds=na_inds)
}

##' For each of a list of matrices `margs`, replicate contents
##' vertically until all elements have `maxlen` rows
##'
##' If any elements are vectors, they are stacked sideways to make
##' a matrix of dimensions  maxlen x length(vector) 
##' 
##' @noRd
rep_matrix_list <- function(margs, maxlen){
  for (i in seq(along=margs)){
    ## replicate matrix contents vertically until length maxlen.  TODO bring this out as a function 
    if (is.matrix(margs[[i]])){
      margs[[i]] <- margs[[i]][rep(seq_len(nrow(margs[[i]])),
                                   length.out=maxlen),,drop=FALSE]
    }
    else margs[[i]] <- matrix(margs[[i]], nrow=maxlen,
                              ncol=length(margs[[i]]), byrow=TRUE)
  }
  margs
}

##' Prepend a new dimension to each of a list of objects, and stretch
##' the contents to a given size.  Used for vectorising functions
##'
##' @param vargs list of vectors
##' 
##' @param margs list of matrices
##' 
##' @param nt size of new dimension
##' 
##' @noRd
stretch_dim_to <- function(vargs, margs, nt){
  for (i in seq_along(vargs))
    vargs[[i]] <- matrix(rep(vargs[[i]], each=nt), nrow=nt)
  for (i in seq_along(margs))
    margs[[i]] <- array(rep(margs[[i]], each=nt),
                           dim=c(nt, dim(margs[[i]])))  
  list(vargs=vargs, margs=margs)
}
