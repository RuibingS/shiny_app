# ################################################## ic.r

curve.mic <- function(model, args, ic_left = -1, ic_right = 1) {
    minv <- args[model.params(model) + 1]
    
    t.n <- curve.inflection(model = model, args = args, a = ic_left, b = ic_right)
    t.m <- model.wrapper(f = "dif1.", model = model, x = t.n, args = args)
    t.y <- model.wrapper(f = "f.", model = model, x = t.n, args = args)
    
    t.b <- t.y - (t.m * t.n)
    return((minv - t.b)/t.m)
}

curve.nic <- function(model, args, ic_left = -1, ic_right = 1) {
    maxv <- args[model.params(model) + 2]
    
    t.n <- curve.inflection(model = model, args = args, a = ic_left, b = ic_right)
    t.m <- model.wrapper(f = "dif1.", model = model, x = t.n, args = args)
    t.y <- model.wrapper(f = "f.", model = model, x = t.n, args = args)
    
    t.b <- t.y - (t.m * t.n)
    return((maxv - t.b)/t.m)
}

curve.inflection <- function(model, args, a, b, prec.x = 1e-06, prec.y = 1e-06) {
    aa <- a
    bb <- b
    f.a <- model.wrapper(f = "dif2.", model = model, x = a, args = args)
    f.b <- model.wrapper(f = "dif2.", model = model, x = b, args = args)
    s.a <- sign(f.a)
    s.b <- sign(f.b)
    
    root <- (bb + aa)/2
    while (bb - aa > prec.x & abs(f.a - f.b) > prec.y) {
        f.root <- model.wrapper(f = "dif2.", model = model, x = root, args = args)
        s.root <- sign(f.root)
        if (s.root == s.a) {
            aa <- root
            f.a <- f.root
        } else {
            bb <- root
            f.b <- f.root
        }
        root <- (bb + aa)/2
    }
    
    return(root)
}

# ################################################## End of Document
