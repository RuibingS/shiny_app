# ################################################## cplot.r

#' Dependencies: > library('ggplot2') must be imported!
#' Dependencies: > library('scales') must be imported!

#' windows() unter MS Windows
#' quartz() unter osX
#' x11() unter Linux

cplot.base <- function(p.data) {
    #' p.y.breaks <- seq(0, 1.5, 0.25)
    #' p.x.breaks <- seq(-10, 10, 1)
    #' p.x.labels <- parse(text = paste('1E', seq(-10, 10, 1), sep = ''))
    breaks <- length(unique(p.data$concentration))
    #' breaks_x <- round(seq(min(breaks_x), max(breaks_x), by = 0.5),1)
   base <- ggplot(p.data) + theme_minimal() + scale_x_continuous(breaks = pretty_breaks(n = breaks)) + scale_y_continuous(breaks = pretty_breaks(n = breaks))
    #' base <- ggplot() + theme_bw()
    #' cplot <- cplot + geom_hline(yintercept = c(0.5, 0.9), linetype = 'dotted')
    #' base <- base + scale_x_continuous(breaks = p.x.breaks, labels = p.x.labels)
    #' base <- base + scale_y_continuous(breaks = p.y.breaks, labels = percent)
    base <- base + labs(title = "Title", x = "concentration", y = "value")
   return(base)
}





cplot.layer.grid <- function(p.data) {
    layer <- scale_x_continuous(breaks = pretty_breaks(n = 10))
    layer <- layer + scale_y_continuous(breaks = pretty_breaks(n = 10))
    return(layer)
}

cplot.layer.data <- function(p.data) {
    layer <- geom_point(data = p.data, aes(x = concentration, y = value), stat = "identity", colour = "darkgray")
    return(layer)
}

cplot.layer.point <- function(p.data, symbol = NULL) {
    layer <- geom_point(data = p.data, aes(x = concentration, y = value), stat = "identity", colour = "black", fill = "black", shape = symbol, size = 2)
    return(layer)
}

cplot.layer.sd <- function(p.data) {
  layer <- geom_errorbar(data = p.data, aes(x = concentration, y = value, ymin = value - sd, ymax = value + sd), colour = "red", width = 0.25)
  return(layer)
}

cplot.layer.sdp <- function(p.data) {
  layer <- geom_point(data = p.data, aes(x = concentration, y = value), stat = "identity", colour = "red")
  return(layer)
}

cplot.layer.model <- function(p.data, fun = NULL, args = NULL) {
    layer <- geom_path(data = p.data, aes(x = concentration, y = value), colour = "blue", stat = "function", fun = fun, args = args)
    return(layer)
}

cplot.layer.conf <- function(p.data) {
    layer <- geom_ribbon(data = p.data, aes(x = concentration, ymin = value_l, ymax = value_u), colour = "darkgray", alpha = 0.25)
    return(layer)
}

# ################################################## End of Document
