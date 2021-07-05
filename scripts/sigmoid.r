# ################################################## sigmoid.r

# ################################################## Model: Probit

f.probit <- function(x, beta1 = 0, beta2 = 1, minv = 0, maxv = 1) {
    return(pnorm(beta1 + beta2 * x) * (maxv - minv) + minv)
}

dif1.probit <- function(x, beta1 = 0, beta2 = 1, minv = 0, maxv = 1) {
    return((beta2*exp(-(1/2)*(beta1+beta2*x)^2)*(maxv-minv))/(sqrt(2*pi)))
    #term <- beta1 + beta2 * x
    #return(dnorm(term) * beta2 * (maxv - minv) + minv)
}

dif2.probit <- function(x, beta1 = 0, beta2 = 1, minv = 0, maxv = 1) {
    return(((beta2^2)*exp(-(1/2)*(beta1+beta2*x)^2)*(maxv-minv)*(beta1+beta2*x))/(sqrt(2*pi)))
    #term <- beta1 + beta2 * x
    #return(-(beta2^2) * (maxv - minv) * term * dnorm(term))
}

inv.probit <- function(x, beta1 = 0, beta2 = 1, minv = 0, maxv = 1, absval = F) {
    if (absval) {
        return((qnorm((x - minv)/(maxv - minv)) - beta1)/beta2)
    } else {
        return((qnorm(x) - beta1)/beta2)
    }
}

init.probit <- function(xy, c = 2) {
    x <- xy[, 1]
    y <- xy[, 2]
    a <- min(x)
    b <- max(x)
    if (cor(x, y) < 0) {
        aux <- a
        a <- b
        b <- aux
    }
    beta1 <- c * (a + b)/(a - b)
    beta2 <- 2 * c/(b - a)
    
    #' minv <- min(y)
    #' maxv <- max(y)
    minv <- unname(quantile(y, .10))
    maxv <- unname(quantile(y, .90))
    
    return(c(beta1, beta2, minv, maxv))
}

obj.probit <- function(pars, xy, err = "lsq", minva = NULL, maxva = NULL) {
    err.f <- get(paste("err.", err, sep = ""))
    
    x <- xy[, 1]
    y <- xy[, 2]
    
    beta1 <- pars[1]
    beta2 <- pars[2]
    minv <- minva
    maxv <- maxva
    if (is.null(minva)) {
        minv <- pars[3]
        if (is.null(maxva)) {
            maxv <- pars[4]
        }
    } else {
        if (is.null(maxva)) {
            maxv <- pars[3]
        }
    }
    s <- 0
    sigma <- 1
    if (err == "huber" || err == "tukey") {
        pvecc = c(beta1, beta2, minv, maxv)
        sigma <- est.sigma(xy = xy, model = "probit", args = pvecc)
    }
    for (i in 1:length(x)) {
        ei <- f.probit(x[i], beta1 = beta1, beta2 = beta2, minv = minv, maxv = maxv) - y[i]
        if (err == "huber" || err == "tukey") {
            rho.ei <- err.f(e = ei, sigma = sigma)
        } else {
            rho.ei <- err.f(e = ei)
        }
        s <- s + rho.ei
    }
    return(s)
}

conf.probit <- function(mat, xmin, xmax, no.intervals = 100, conf.level = 0.95, lwd = 2, lty = 1) {
    x <- ((0:no.intervals) * (xmax - xmin)/no.intervals) + xmin
    yl <- rep(NA, no.intervals + 1)
    yu <- rep(NA, no.intervals + 1)
    y <- rep(NA, nrow(mat))
    for (i in 0:no.intervals) {
        for (j in 1:nrow(mat)) {
            y[j] <- f.probit(x[i + 1], beta1 = mat[j, 1], beta2 = mat[j, 2], minv = mat[j, 3], maxv = mat[j, 4])
        }
        yl[i + 1] <- quantile(y, (1 - conf.level)/2)
        yu[i + 1] <- quantile(y, (1 + conf.level)/2)
    }
    return(cbind(x, yl, yu))
}

# ################################################## Model: Logit

f.logit <- function(x, beta1 = 0, beta2 = 1, minv = 0, maxv = 1) {
    return((1/(1 + exp(-(beta1 + beta2 * x)))) * (maxv - minv) + minv)
}

dif1.logit <- function(x, beta1 = 0, beta2 = 1, minv = 0, maxv = 1) {
    term <- exp(beta1 + beta2 * x)
    return((beta2 * (maxv - minv) * term)/((term + 1)^2))
}

dif2.logit <- function(x, beta1 = 0, beta2 = 1, minv = 0, maxv = 1) {
    term <- exp(beta1 + beta2 * x)
    return(((beta2^2) * (maxv - minv) * term * (term - 1))/((term + 1)^3))
}

inv.logit <- function(x, beta1 = 0, beta2 = 1, minv = 0, maxv = 1, absval = F) {
    y <- x
    if (absval) {
        y <- (x - minv)/(maxv - minv)
    }
    k <- y/(1 - y)
    return((log(k) - beta1)/beta2)
}

init.logit <- function(xy, c = 2) {
    x <- xy[, 1]
    y <- xy[, 2]
    a <- min(x)
    b <- max(x)
    if (cor(x, y) < 0) {
        aux <- a
        a <- b
        b <- aux
    }
    beta1 <- c * (a + b)/(a - b)
    beta2 <- 2 * c/(b - a)
    
    #' minv <- min(y)
    #' maxv <- max(y)
    minv <- unname(quantile(y, .10))
    maxv <- unname(quantile(y, .90))
    
    return(c(beta1, beta2, minv, maxv))
}

obj.logit <- function(pars, xy, err = "lsq", minva = NULL, maxva = NULL) {
    err.f <- get(paste("err.", err, sep = ""))
    
    x <- xy[, 1]
    y <- xy[, 2]
    beta1 <- pars[1]
    beta2 <- pars[2]
    minv <- minva
    maxv <- maxva
    if (is.null(minva)) {
        minv <- pars[3]
        if (is.null(maxva)) {
            maxv <- pars[4]
        }
    } else {
        if (is.null(maxva)) {
            maxv <- pars[3]
        }
    }
    s <- 0
    sigma <- 1
    if (err == "huber" || err == "tukey") {
        pvecc = c(beta1, beta2, minv, maxv)
        sigma <- est.sigma(xy = xy, model = "logit", args = pvecc)
    }
    for (i in 1:length(x)) {
        ei <- f.logit(x[i], beta1 = beta1, beta2 = beta2, minv = minv, maxv = maxv) - y[i]
        if (err == "huber" || err == "tukey") {
            rho.ei <- err.f(e = ei, sigma = sigma)
        } else {
            rho.ei <- err.f(e = ei)
        }
        s <- s + rho.ei
    }
    return(s)
}

conf.logit <- function(mat, xmin, xmax, no.intervals = 100, conf.level = 0.95, lwd = 2, lty = 1) {
    x <- ((0:no.intervals) * (xmax - xmin)/no.intervals) + xmin
    yl <- rep(NA, no.intervals + 1)
    yu <- rep(NA, no.intervals + 1)
    y <- rep(NA, nrow(mat))
    for (i in 0:no.intervals) {
        for (j in 1:nrow(mat)) {
            y[j] <- f.logit(x[i + 1], beta1 = mat[j, 1], beta2 = mat[j, 2], minv = mat[j, 3], maxv = mat[j, 4])
        }
        yl[i + 1] <- quantile(y, (1 - conf.level)/2)
        yu[i + 1] <- quantile(y, (1 + conf.level)/2)
    }
    return(cbind(x, yl, yu))
}

# ################################################## Model: Weibull

f.weibull <- function(x, beta1 = 0, beta2 = 1, minv = 0, maxv = 1) {
    return((1 - exp(-exp(beta1 + beta2 * x))) * (maxv - minv) + minv)
}

dif1.weibull <- function(x, beta1 = 0, beta2 = 1, minv = 0, maxv = 1) {
    term <- beta1 + beta2 * x
    return(beta2 * (maxv - minv) * exp(-(exp(term)) + term))
}

dif2.weibull <- function(x, beta1 = 0, beta2 = 1, minv = 0, maxv = 1) {
    term <- beta1 + beta2 * x
    return((beta2 * (maxv - minv) * exp(-(exp(term)) + term)) * (beta2 - beta2 * exp(term)))
}

inv.weibull <- function(x, beta1 = 0, beta2 = 1, minv = 0, maxv = 1, absval = F) {
    y <- x
    if (absval) {
        y <- (x - minv)/(maxv - minv)
    }
    k <- -log(1 - y)
    return((log(k) - beta1)/beta2)
}

init.weibull <- function(xy, c = 2) {
    x <- xy[, 1]
    y <- xy[, 2]
    a <- min(x)
    b <- max(x)
    if (cor(x, y) < 0) {
        aux <- a
        a <- b
        b <- aux
    }
    beta1 <- c * (a + b)/(a - b)
    beta2 <- 2 * c/(b - a)
    
    #' minv <- min(y)
    #' maxv <- max(y)
    minv <- unname(quantile(y, .10))
    maxv <- unname(quantile(y, .90))
    
    return(c(beta1, beta2, minv, maxv))
}

obj.weibull <- function(pars, xy, err = "lsq", minva = NULL, maxva = NULL) {
    err.f <- get(paste("err.", err, sep = ""))
    
    x <- xy[, 1]
    y <- xy[, 2]
    beta1 <- pars[1]
    beta2 <- pars[2]
    minv <- minva
    maxv <- maxva
    if (is.null(minva)) {
        minv <- pars[3]
        if (is.null(maxva)) {
            maxv <- pars[4]
        }
    } else {
        if (is.null(maxva)) {
            maxv <- pars[3]
        }
    }
    s <- 0
    sigma <- 1
    if (err == "huber" || err == "tukey") {
        pvecc = c(beta1, beta2, minv, maxv)
        sigma <- est.sigma(xy = xy, model = "weibull", args = pvecc)
    }
    for (i in 1:length(x)) {
        ei <- f.weibull(x[i], beta1 = beta1, beta2 = beta2, minv = minv, maxv = maxv) - y[i]
        if (err == "huber" || err == "tukey") {
            rho.ei <- err.f(e = ei, sigma = sigma)
        } else {
            rho.ei <- err.f(e = ei)
        }
        s <- s + rho.ei
    }
    return(s)
}

conf.weibull <- function(mat, xmin, xmax, no.intervals = 100, conf.level = 0.95, lwd = 2, lty = 1) {
    x <- ((0:no.intervals) * (xmax - xmin)/no.intervals) + xmin
    yl <- rep(NA, no.intervals + 1)
    yu <- rep(NA, no.intervals + 1)
    y <- rep(NA, nrow(mat))
    for (i in 0:no.intervals) {
        for (j in 1:nrow(mat)) {
            y[j] <- f.weibull(x[i + 1], beta1 = mat[j, 1], beta2 = mat[j, 2], minv = mat[j, 3], maxv = mat[j, 4])
        }
        yl[i + 1] <- quantile(y, (1 - conf.level)/2)
        yu[i + 1] <- quantile(y, (1 + conf.level)/2)
    }
    return(cbind(x, yl, yu))
}

# ################################################## Model: Generalised Logit I

f.glogitI <- function(x, beta1 = 0, beta2 = 1, beta3 = 1, minv = 0, maxv = 1) {
    return((1/(1 + exp(-(beta1 + beta2 * x)))^beta3) * (maxv - minv) + minv)
}

dif1.glogitI <- function(x, beta1 = 0, beta2 = 1, beta3 = 1, minv = 0, maxv = 1) {
    term <- exp(beta1 + beta2 * x)
    return((beta2 * beta3 * (maxv - minv) * ((exp(-beta1 - beta2 * x) + 1)^-beta3))/(term + 1))
}

dif2.glogitI <- function(x, beta1 = 0, beta2 = 1, beta3 = 1, minv = 0, maxv = 1) {
    term <- exp(beta1 + beta2 * x)
    return((beta2^2 * beta3 * ((exp(-beta1 - beta2 * x) + 1)^-beta3) * (term - beta3))/((term + 1)^2))
}

inv.glogitI <- function(x, beta1 = 0, beta2 = 1, beta3 = 1, minv = 0, maxv = 1, absval = F) {
    y <- x
    if (absval) {
        y <- (x - minv)/(maxv - minv)
    }
    k <- (1/y)^(1/beta3) - 1
    return((-log(k) - beta1)/beta2)
}

init.glogitI <- function(xy, c = 2) {
    x <- xy[, 1]
    y <- xy[, 2]
    a <- min(x)
    b <- max(x)
    if (cor(x, y) < 0) {
        aux <- a
        a <- b
        b <- aux
    }
    beta1 <- c * (a + b)/(a - b)
    beta2 <- 2 * c/(b - a)
    beta3 <- 1
    
    #' minv <- min(y)
    #' maxv <- max(y)
    minv <- unname(quantile(y, .10))
    maxv <- unname(quantile(y, .90))
    
    return(c(beta1, beta2, beta3, minv, maxv))
}

obj.glogitI <- function(pars, xy, err = "lsq", minva = NULL, maxva = NULL) {
    err.f <- get(paste("err.", err, sep = ""))
    
    x <- xy[, 1]
    y <- xy[, 2]
    beta1 <- pars[1]
    beta2 <- pars[2]
    beta3 <- pars[3]
    minv <- minva
    maxv <- maxva
    if (is.null(minva)) {
        minv <- pars[4]
        if (is.null(maxva)) {
            maxv <- pars[5]
        }
    } else {
        if (is.null(maxva)) {
            maxv <- pars[4]
        }
    }
    s <- 0
    sigma <- 1
    if (err == "huber" || err == "tukey") {
        pvecc = c(beta1, beta2, beta3, minv, maxv)
        sigma <- est.sigma(xy = xy, model = "glogitI", args = pvecc)
    }
    for (i in 1:length(x)) {
        ei <- f.glogitI(x[i], beta1 = beta1, beta2 = beta2, beta3 = beta3, minv = minv, maxv = maxv) - y[i]
        if (err == "huber" || err == "tukey") {
            rho.ei <- err.f(e = ei, sigma = sigma)
        } else {
            rho.ei <- err.f(e = ei)
        }
        s <- s + rho.ei
    }
    return(s)
}

conf.glogitI <- function(mat, xmin, xmax, no.intervals = 100, conf.level = 0.95, lwd = 2, lty = 1) {
    x <- ((0:no.intervals) * (xmax - xmin)/no.intervals) + xmin
    yl <- rep(NA, no.intervals + 1)
    yu <- rep(NA, no.intervals + 1)
    y <- rep(NA, nrow(mat))
    for (i in 0:no.intervals) {
        for (j in 1:nrow(mat)) {
            y[j] <- f.glogitI(x[i + 1], beta1 = mat[j, 1], beta2 = mat[j, 2], beta3 = mat[j, 3], minv = mat[j, 4], maxv = mat[j, 5])
        }
        yl[i + 1] <- quantile(y, (1 - conf.level)/2)
        yu[i + 1] <- quantile(y, (1 + conf.level)/2)
    }
    return(cbind(x, yl, yu))
}

# ################################################## Model: Generalised Logit II

f.glogitII <- function(x, beta1 = 0, beta2 = 1, beta3 = 1, minv = 0, maxv = 1) {
    return((1 - (1/(1 + exp(beta1 + beta2 * x))^beta3)) * (maxv - minv) + minv)
}

dif1.glogitII <- function(x, beta1 = 0, beta2 = 1, beta3 = 1, minv = 0, maxv = 1) {
    term <- exp(beta1 + beta2 * x)
    return(beta2 * beta3 * (maxv - minv) * term * ((term + 1)^(-beta3 - 1)))
}

dif2.glogitII <- function(x, beta1 = 0, beta2 = 1, beta3 = 1, minv = 0, maxv = 1) {
    term <- exp(beta1 + beta2 * x)
    return(-(beta2^2) * beta3 * (maxv - minv) * term * ((term + 1)^(-beta3 - 2)) * (beta3 * term - 1))
}

inv.glogitII <- function(x, beta1 = 0, beta2 = 1, beta3 = 1, minv = 0, maxv = 1, absval = F) {
    y <- x
    if (absval) {
        y <- (x - minv)/(maxv - minv)
    }
    k <- (1/(1 - y))^(1/beta3) - 1
    return((log(k) - beta1)/beta2)
}

init.glogitII <- function(xy, c = 2) {
    x <- xy[, 1]
    y <- xy[, 2]
    a <- min(x)
    b <- max(x)
    if (cor(x, y) < 0) {
        aux <- a
        a <- b
        b <- aux
    }
    beta1 <- c * (a + b)/(a - b)
    beta2 <- 2 * c/(b - a)
    beta3 <- 1
    
    #' minv <- min(y)
    #' maxv <- max(y)
    minv <- unname(quantile(y, .10))
    maxv <- unname(quantile(y, .90))
    
    return(c(beta1, beta2, beta3, minv, maxv))
}

obj.glogitII <- function(pars, xy, err = "lsq", minva = NULL, maxva = NULL) {
    err.f <- get(paste("err.", err, sep = ""))
    
    x <- xy[, 1]
    y <- xy[, 2]
    beta1 <- pars[1]
    beta2 <- pars[2]
    beta3 <- pars[3]
    minv <- minva
    maxv <- maxva
    if (is.null(minva)) {
        minv <- pars[4]
        if (is.null(maxva)) {
            maxv <- pars[5]
        }
    } else {
        if (is.null(maxva)) {
            maxv <- pars[4]
        }
    }
    s <- 0
    sigma <- 1
    if (err == "huber" || err == "tukey") {
        pvecc = c(beta1, beta2, beta3, minv, maxv)
        sigma <- est.sigma(xy = xy, model = "glogitII", args = pvecc)
    }
    for (i in 1:length(x)) {
        ei <- f.glogitII(x[i], beta1 = beta1, beta2 = beta2, beta3 = beta3, minv = minv, maxv = maxv) - y[i]
        if (err == "huber" || err == "tukey") {
            rho.ei <- err.f(e = ei, sigma = sigma)
        } else {
            rho.ei <- err.f(e = ei)
        }
        s <- s + rho.ei
    }
    return(s)
}

conf.glogitII <- function(mat, xmin, xmax, no.intervals = 100, conf.level = 0.95, lwd = 2, lty = 1) {
    x <- ((0:no.intervals) * (xmax - xmin)/no.intervals) + xmin
    yl <- rep(NA, no.intervals + 1)
    yu <- rep(NA, no.intervals + 1)
    y <- rep(NA, nrow(mat))
    for (i in 0:no.intervals) {
        for (j in 1:nrow(mat)) {
            y[j] <- f.glogitII(x[i + 1], beta1 = mat[j, 1], beta2 = mat[j, 2], beta3 = mat[j, 3], minv = mat[j, 4], maxv = mat[j, 5])
        }
        yl[i + 1] <- quantile(y, (1 - conf.level)/2)
        yu[i + 1] <- quantile(y, (1 + conf.level)/2)
    }
    return(cbind(x, yl, yu))
}

# ################################################## Model: Aranda-Ordaz

f.ao <- function(x, beta1 = 0, beta2 = 1, beta3 = 1, minv = 0, maxv = 1) {
    return((1 - (1/(1 + exp(beta1 + beta2 * x)/beta3)^beta3)) * (maxv - minv) + minv)
}

dif1.ao <- function(x, beta1 = 0, beta2 = 1, beta3 = 1, minv = 0, maxv = 1) {
    term <- exp(beta1 + beta2 * x)
    return(beta2 * (maxv - minv) * term * (((term + beta3)/(beta3))^(-beta3 - 1)))
}

dif2.ao <- function(x, beta1 = 0, beta2 = 1, beta3 = 1, minv = 0, maxv = 1) {
    term <- exp(beta1 + beta2 * x)
    return(((beta2^2) * (beta3^2) * (maxv - minv) * term * (term - 1) * (((term + beta3)/(beta3))^-beta3))/((term + beta3)^2))
}

inv.ao <- function(x, beta1 = 0, beta2 = 1, beta3 = 1, minv = 0, maxv = 1, absval = F) {
    y <- x
    if (absval) {
        y <- (x - minv)/(maxv - minv)
    }
    k <- beta3 * ((1/(1 - y))^(1/beta3) - 1)
    return((log(k) - beta1)/beta2)
}

init.ao <- function(xy, c = 2) {
    x <- xy[, 1]
    y <- xy[, 2]
    a <- min(x)
    b <- max(x)
    if (cor(x, y) < 0) {
        aux <- a
        a <- b
        b <- aux
    }
    beta1 <- c * (a + b)/(a - b)
    beta2 <- 2 * c/(b - a)
    beta3 <- 1
    
    #' minv <- min(y)
    #' maxv <- max(y)
    minv <- unname(quantile(y, .10))
    maxv <- unname(quantile(y, .90))
    
    return(c(beta1, beta2, beta3, minv, maxv))
}

obj.ao <- function(pars, xy, err = "lsq", minva = NULL, maxva = NULL) {
    err.f <- get(paste("err.", err, sep = ""))
    
    x <- xy[, 1]
    y <- xy[, 2]
    beta1 <- pars[1]
    beta2 <- pars[2]
    beta3 <- pars[3]
    minv <- minva
    maxv <- maxva
    if (is.null(minva)) {
        minv <- pars[4]
        if (is.null(maxva)) {
            maxv <- pars[5]
        }
    } else {
        if (is.null(maxva)) {
            maxv <- pars[4]
        }
    }
    s <- 0
    sigma <- 1
    if (err == "huber" || err == "tukey") {
        pvecc = c(beta1, beta2, beta3, minv, maxv)
        sigma <- est.sigma(xy = xy, model = "ao", args = pvecc)
    }
    for (i in 1:length(x)) {
        ei <- f.ao(x[i], beta1 = beta1, beta2 = beta2, beta3 = beta3, minv = minv, maxv = maxv) - y[i]
        if (err == "huber" || err == "tukey") {
            rho.ei <- err.f(e = ei, sigma = sigma)
        } else {
            rho.ei <- err.f(e = ei)
        }
        s <- s + rho.ei
    }
    return(s)
}

conf.ao <- function(mat, xmin, xmax, no.intervals = 100, conf.level = 0.95, lwd = 2, lty = 1) {
    x <- ((0:no.intervals) * (xmax - xmin)/no.intervals) + xmin
    yl <- rep(NA, no.intervals + 1)
    yu <- rep(NA, no.intervals + 1)
    y <- rep(NA, nrow(mat))
    for (i in 0:no.intervals) {
        for (j in 1:nrow(mat)) {
            y[j] <- f.ao(x[i + 1], beta1 = mat[j, 1], beta2 = mat[j, 2], beta3 = mat[j, 3], minv = mat[j, 4], maxv = mat[j, 5])
        }
        yl[i + 1] <- quantile(y, (1 - conf.level)/2)
        yu[i + 1] <- quantile(y, (1 + conf.level)/2)
    }
    return(cbind(x, yl, yu))
}

# ################################################## Model: Gompertz

f.gompertz <- function(x, beta1 = 0, beta2 = 1, minv = 0, maxv = 1) {
    return(exp(-exp(beta2 * (x - beta1))) * (maxv - minv) + minv)
}

dif1.gompertz <- function(x, beta1 = 0, beta2 = 1, minv = 0, maxv = 1) {
    term <- beta2 * (x - beta1)
    return(beta2 * (minv - maxv) * exp(term - exp(term)))
}

dif2.gompertz <- function(x, beta1 = 0, beta2 = 1, minv = 0, maxv = 1) {
    term <- beta2 * (x - beta1)
    return(beta2^2 * (maxv - minv) * exp(term * -(exp(term))) * (exp(term) - 1))
}

inv.gompertz <- function(x, beta1 = 0, beta2 = 1, minv = 0, maxv = 1, absval = F) {
    if (absval) {
        return((beta2 * beta1 - log(abs(1/(log((minv - x)/(minv - maxv))))))/beta2)
    } else {
        return((beta2 * beta1 - log(abs(1/(log(x)))))/beta2)
    }
}

init.gompertz <- function(xy, c = 2) {
  x <- xy[, 1]
  y <- xy[, 2]
  a <- min(x)
  b <- max(x)
  if (cor(x, y) < 0) {
    aux <- a
    a <- b
    b <- aux
  }
  beta1 <- c * (a + b)/(a - b)
  beta2 <- -2 * c/(b - a)
  
  #' minv <- min(y)
  #' maxv <- max(y)
  minv <- unname(quantile(y, .10))
  maxv <- unname(quantile(y, .90))
  
  return(c(beta1, beta2, minv, maxv))
}

obj.gompertz <- function(pars, xy, err = "lsq", minva = NULL, maxva = NULL) {
    err.f <- get(paste("err.", err, sep = ""))
    
    x <- xy[, 1]
    y <- xy[, 2]
    beta1 <- pars[1]
    beta2 <- pars[2]
    minv <- minva
    maxv <- maxva
    if (is.null(minva)) {
        minv <- pars[3]
        if (is.null(maxva)) {
            maxv <- pars[4]
        }
    } else {
        if (is.null(maxva)) {
            maxv <- pars[3]
        }
    }
    s <- 0
    sigma <- 1
    if (err == "huber" || err == "tukey") {
        pvecc = c(beta1, beta2, minv, maxv)
        sigma <- est.sigma(xy = xy, model = "gompertz", args = pvecc)
    }
    for (i in 1:length(x)) {
        ei <- f.gompertz(x[i], beta1 = beta1, beta2 = beta2, minv = minv, maxv = maxv) - y[i]
        if (err == "huber" || err == "tukey") {
            rho.ei <- err.f(e = ei, sigma = sigma)
        } else {
            rho.ei <- err.f(e = ei)
        }
        s <- s + rho.ei
    }
    return(s)
}

conf.gompertz <- function(mat, xmin, xmax, no.intervals = 100, conf.level = 0.95, lwd = 2, lty = 1) {
    x <- ((0:no.intervals) * (xmax - xmin)/no.intervals) + xmin
    yl <- rep(NA, no.intervals + 1)
    yu <- rep(NA, no.intervals + 1)
    y <- rep(NA, nrow(mat))
    for (i in 0:no.intervals) {
        for (j in 1:nrow(mat)) {
            y[j] <- f.gompertz(x[i + 1], beta1 = mat[j, 1], beta2 = mat[j, 2], minv = mat[j, 3], maxv = mat[j, 4])
        }
        yl[i + 1] <- quantile(y, (1 - conf.level)/2)
        yu[i + 1] <- quantile(y, (1 + conf.level)/2)
    }
    return(cbind(x, yl, yu))
}

# ################################################## End of Document
