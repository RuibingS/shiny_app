# ################################################## model.r

model.fit <- function(xy, model = "probit", minva = NULL, maxva = NULL, err = "lsq", c = 2) {
    init.f <- get(paste("init.", model, sep = ""))
    obj.f <- get(paste("obj.", model, sep = ""))
    par <- c()
    if (err == "abs" || err == "lsq") {
      np <- length(init.f(xy, c = c))
      parinds <- 1:(np - 2)
      if (is.null(minva)) {
        parinds <- c(parinds, np - 1)
      }
      if (is.null(maxva)) {
        parinds <- c(parinds, np)
      }
      par <- init.f(xy, c = c)[parinds]
    } else if (err == "huber" || err == "tukey") {
      fit <- model.fit(xy, model = model, minva = minva, maxva = maxva, err = "abs", c = c)
      par <- fit$par
    }
    return(optim(par = par, fn = obj.f, xy = xy, minva = minva, maxva = maxva, err = err))
}

model.bootstrap <- function(xy, model = "probit", err = "lsq", k = 1000, minva = NULL, maxva = NULL, c = 2) {
    x <- xy[, 1]
    y <- xy[, 2]
    
    init.f <- get(paste("init.", model, sep = ""))
    obj.f <- get(paste("obj.", model, sep = ""))
    
    np <- length(init.f(xy, c = c))
    parinds <- 1:(np - 2)
    if (is.null(minva)) {
        parinds <- c(parinds, np - 1)
    }
    if (is.null(maxva)) {
        parinds <- c(parinds, np)
    }
    n <- length(x)
    mat <- matrix(NA, ncol = length(parinds), nrow = k)
    for (i in 1:k) {
        inds <- sample(1:n, replace = T)
        xx <- x[inds]
        yy <- y[inds]
        optr <- optim(par = init.f(cbind(xx, yy), c = c)[parinds], fn = obj.f, xy = cbind(xx, yy), minva = minva, maxva = maxva, err = err)
        mat[i, ] <- t(optr$par)
    }
    return(mat)
}

model.loo <- function(xy, model = "probit", minva = NULL, maxva = NULL, err = "lsq", c = 2) {
    f.f <- get(paste("f.", model, sep = ""))
    err.f <- get(paste("err.", err, sep = ""))
    
    x <- xy[, 1]
    y <- xy[, 2]
    
    sumerr <- 0
    for (i in 1:length(x)) {
        xx <- x[-i]
        yy <- y[-i]
        optr <- model.fit(cbind(xx, yy), model = model, minva = minva, maxva = maxva, err = err, c = c)
        
        nop <- model.params(model)
        
        pvec <- optr$par[1:nop]
        indmax <- nop + 2
        
        minv <- minva
        maxv <- maxva
        if (is.null(minva)) {
            minv <- optr$par[nop + 1]
            if (is.null(maxva)) {
                maxv <- optr$par[nop + 2]
            }
        } else {
            if (is.null(maxva)) {
                maxv <- optr$par[nop + 1]
            }
        }
        
        pvecc <- c(pvec, minv, maxv)
        
        sigma <- 1
        if (err == "huber" || err == "tukey") {
            sigma <- est.sigma(xy = cbind(xx, yy), model = model, args = pvecc)
        }
        
        #erri <- 0
        ei <- model.wrapper(f = "f.", model = model, x = x[i], args = pvecc) - y[i]
        
        if (err == "huber" || err == "tukey") {
            rho.ei <- err.f(e = ei, sigma = sigma)
        } else {
            rho.ei <- err.f(e = ei)
        }
        
        sumerr <- sumerr + rho.ei
    }
    
    return(sumerr)
}

model.params <- function(model) {
    if (model == "probit") {
        return(2)
    }
    if (model == "logit") {
        return(2)
    }
    if (model == "weibull") {
        return(2)
    }
    if (model == "glogitI") {
        return(3)
    }
    if (model == "glogitII") {
        return(3)
    }
    if (model == "ao") {
        return(3)
    }
    if (model == "gompertz") {
        return(2)
    }
}

model.select <- function(xy, models = c("probit", "logit", "weibull", "glogitI", "glogitII", "ao", "gompertz"), minva = NULL, maxva = NULL, 
    err = "lsq", criterion = "model.crit.aicc", c = 2) {
    modevals <- rep(NA, length(models))
    for (i in 1:length(models)) {
        erri <- model.loo(xy, model = models[i], minva = minva, maxva = maxva, err = err, c = c)
        nop <- model.params(models[i])
        if (!is.null(minva)) {
            nop <- nop + 1
        }
        if (!is.null(maxva)) {
            nop <- nop + 1
        }
        modevals[i] <- model.crit.aicc(nop, nrow(xy), erri/nrow(xy))
    }
    return(models[order(modevals)])
}

model.crit.aicc <- function(k, n, s2) {
    return(n * log2(s2) + 2 * k + (2 * k * (k + 1))/(n - k - 1))
}

model.crit.aic <- function(k, n, s2) {
  return(n * log2(s2) + 2 * k)
}

model.crit.bic <- function(k, n, s2) {
  return(n * log2(s2) + k * log2(n))
}

model.wrapper <- function(f, model, x, args) {
    wrapper <- get(paste(f, model, sep = ""))
    nop <- model.params(model)
    if (nop == 2) {
        return(wrapper(x = x, beta1 = args[1], beta2 = args[2], minv = args[3], maxv = args[4]))
    }
    if (nop == 3) {
        return(wrapper(x = x, beta1 = args[1], beta2 = args[2], beta3 = args[3], minv = args[4], maxv = args[5]))
    }
}

# ################################################## End of Document
