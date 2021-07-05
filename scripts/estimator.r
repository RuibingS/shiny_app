# ################################################## estimator.r

err.abs <- function(e) {
    if(!is.finite(e)) {
      # NaN/Inf -> Inf for undefined outcomes
      return(Inf)
    }
    return(abs(e))
}

err.lsq <- function(e) {
    if(!is.finite(e)) {
      # NaN/Inf -> Inf for undefined outcomes
      return(Inf)
    }
    return(e^2)
}

err.huber <- function(e, k = 1.345, sigma = 1) {
    if(!is.finite(e) || !is.finite(k) || !is.finite(sigma)) {
      # NaN/Inf -> Inf for undefined outcomes
      return(Inf)
    }
    if(is.infinite(sigma)) {
      return(Inf)
    }
    k <- k * sigma
    if (abs(e) <= k) {
        return((1/2) * (e^2))
    } else {
        return(k * abs(e) - (1/2) * (k^2))
    }
}

err.tukey <- function(e, k = 4.685, sigma = 1) {
    if(is.infinite(e) || is.infinite(k) || is.infinite(sigma)) {
      # NaN -> Inf for undefined outcomes
      return(Inf)
    }
    k <- k * sigma
    c <- (k^2)/6
    if (abs(e) <= k) {
        return(c * (1 - (1 - (e/k)^2)^3))
    } else {
        return(c)
    }
}

est.sigma <- function(xy, model, args) {
    x <- xy[, 1]
    y <- xy[, 2]
    s <- xy[, 1]
    for (i in 1:length(s)) {
        ei <- model.wrapper(f = "f.", model = model, x = x[i], args = args) - y[i]
        if(is.nan(ei)) {
          # NaN -> Inf for undefined outcomes
          return(Inf)
        }
        s[i] <- abs(ei)
    }
    sigma <- (median(s)/0.6745)
    return(sigma)
}

# ################################################## End of Document
