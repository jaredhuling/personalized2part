
setup_y <- function(y, family)
{
    if (family == "binomial")
    {
        yvals <- sort(unique(y))
        if (all(yvals == c(0, 1)))
        {
            y <- 2 * y - 1
        } else
        {
            if (!all(yvals == c(-1, 1)))
            {
                stop("y has invalid values, can only take values 0 and 1")
            }
        }
    } else if (family == "gamma")
    {
        if (any(y <= 0))
        {
            stop("y must be strictly positive")
        }
    }

    y
}
