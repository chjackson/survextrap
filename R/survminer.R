## Code from the survminer package, used under the GPL-2
##   Kassambara A, Kosinski M, Biecek P (2021). _survminer: Drawing Survival
##  Curves using 'ggplot2'_. R package version 0.4.9,
##  <https://CRAN.R-project.org/package=survminer>.

## Copied, rather than imported, to reduce dependencies
## tidyr::gather call replaced with stats::reshape

surv_summary <- function (x, data = NULL) 
{
    res <- as.data.frame(.compact(unclass(x)[c("time", "n.risk", 
        "n.event", "n.censor")]))
    if (inherits(x, "survfitms")) {
        surv <- 1 - x$prev
        upper <- 1 - x$upper
        lower <- 1 - x$lower
        res <- cbind(res, surv = c(surv), std.err = c(x$std.err), 
            upper = c(upper), lower = c(lower))
        res$state <- rep(x$states, each = nrow(surv))
    }
    else {
        if (is.matrix(x$surv)) {
            ncurve <- ncol(x$surv)
            res <- data.frame(time = rep(x$time, ncurve), n.risk = rep(x$n.risk, 
                ncurve), n.event = rep(x$n.event, ncurve), n.censor = rep(x$n.censor, 
                ncurve))
            res <- cbind(res, surv = .flat(x$surv), std.err = .flat(x$std.err), 
                upper = .flat(x$upper), lower = .flat(x$lower))
            res$strata <- as.factor(rep(colnames(x$surv), each = nrow(x$surv)))
        }
        else res <- cbind(res, surv = x$surv, std.err = x$std.err, 
            upper = x$upper, lower = x$lower)
    }
    if (!is.null(x$strata)) {
        data <- .get_data(x, data = data)
        res$strata <- rep(names(x$strata), x$strata)
        res$strata <- .clean_strata(res$strata, x)
        variables <- .get_variables(res$strata, x, data)
        for (variable in variables) res[[variable]] <- .get_variable_value(variable, 
            res$strata, x, data)
    }
    structure(res, class = c("data.frame", "surv_summary"))
    attr(res, "table") <- as.data.frame(summary(x)$table)
    res
}

.compact <- function (x) 
{
    Filter(Negate(is.null), x)
}

.flat <- function (x) 
{
    if (is.null(x)) 
        return(NA)
    x <- as.data.frame(x)
#    x <- tidyr::gather_(x, key_col = "key", value_col = "value", 
#        gather_cols = colnames(x))
    x <- reshape(x, direction="long",
                 varying=colnames(x), timevar="key", v.names="value")
    x$value
}

.get_data <- function (fit, data = NULL, complain = TRUE) 
{
    if (is.null(data)) {
        if (complain) 
            warning("The `data` argument is not provided. Data will be extracted from model fit.")
        data <- eval(fit$call$data)
        if (is.null(data)) 
            stop("The `data` argument should be provided either to ggsurvfit or survfit.")
    }
    data
}

.get_variables <- function (strata, fit, data = NULL) 
{
    variables <- sapply(as.vector(strata), function(x) {
        x <- unlist(strsplit(x, "=|,\\s+", perl = TRUE))
        x[seq(1, length(x), 2)]
    })
    variables <- unique(as.vector(variables))
    variables <- intersect(variables, colnames(.get_data(fit, 
        data)))
    variables
}

.clean_strata <- function (strata, fit) 
{
    is_dollar_sign <- grepl("$", as.character(strata)[1], fixed = TRUE)
    if (is_dollar_sign) {
        strata <- as.character(strata)
        data_name <- unlist(strsplit(strata[1], "$", fixed = TRUE))[1]
        strata <- gsub(paste0(data_name, "$"), "", strata, fixed = TRUE)
        strata <- as.factor(strata)
    }
    else if (!missing(fit)) 
        strata <- factor(strata, levels = names(fit$strata))
    return(strata)
}

.get_variable_value <- function (variable, strata, fit, data = NULL) 
{
    res <- sapply(as.vector(strata), function(x) {
        x <- unlist(strsplit(x, "=|(\\s+)?,\\s+", perl = TRUE))
        index <- grep(paste0("^", variable, "$"), x)[1]
        .trim(x[index + 1])
    })
    res <- as.vector(res)
    var_levels <- levels(.get_data(fit, data)[, variable])
    if (!is.null(var_levels)) 
        res <- factor(res, levels = var_levels)
    else res <- as.factor(res)
    res
}

.trim <- function (x) 
{
    gsub("^\\s+|\\s+$", "", x)
}
