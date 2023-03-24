## internal function in stats
## Used to mimic the names returned by quantile()
## NO DO NT NEED THIS AN Y MORE

format_perc <- function (x, digits = max(2L, getOption("digits")), probability = TRUE, 
                         use.fC = length(x) < 100, ...) 
{
    if (length(x)) {
        if (probability) 
            x <- 100 * x
        ans <- paste0(if (use.fC) 
            formatC(x, format = "fg", width = 1, digits = digits)
        else format(x, trim = TRUE, digits = digits, ...), "%")
        ans[is.na(x)] <- ""
        ans
    }
    else character(0)
}
