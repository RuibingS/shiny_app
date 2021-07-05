library(shiny)

shinyServer(function(input, output, session) {
  
  df <- reactive({
    req(input$file1)
    read.table(input$file1$datapath,sep = ",", header = TRUE)
   
    })
  #' load and process data
  #' load and process data
  get.dat <- reactive({
    dataset <- df()
    compound <- input$compound
    sub_blank <- input$sub_blank
    
      dat <- df()
      blank <- subset(dat,is.na(concentration))
      blank <- blank$value
      if(sub_blank && length(blank) > 0){
        blank <- mean(blank)
        dat$value <- dat$value - blank
      }
      dat <- dat[dat$compound %in% compound,]
      dat <- subset(dat, is.finite(concentration))
      return(dat)
    
  })
  #' set or calculate model
  get.mod <- reactive({
    dat <- get.dat()
    model <- input$model
    method <- input$method
    if(is.null(dat)) {
      return(NULL)
    } else if(model == "none") {
      return(NULL)
    }
    else if(model == "auto") {
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Selecting best model.", detail = "Please wait...", value = .5)
      xy <- cbind(x = dat$concentration, y =dat$value)
      model <- model.select(xy=xy, err = method)[1]
      progress$set(message = "Selecting best model.", detail = "Please wait...", value = 1)
      return(model)
    } else {
      return(model)
    }
  })
  
  #' fit data to a model
  get.fit <- reactive({
    dat <- get.dat()
    mod <- get.mod()
    if(is.null(mod) || is.null(dat)) {
      return(NULL)
    } else if(mod!="none" ) {
      fit <- model.fit(xy = dat, model = mod, minva = NULL, maxva = NULL, err = input$method)
      return(fit)
    } else {
      return(NULL)
    }
  })
  
  #' calculate confidence band
  get.band <- reactive({
    mod <- get.mod()
    dat <- get.dat()
    show <- input$show_band
    if(is.null(mod) || is.null(dat)) {
      return(NULL)
    } else if(mod!="none" && show) {
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Rendering conf-band.", detail = "Please wait...", value = .5)
      f.conf <- get(paste("conf.", mod, sep=""))
      mat <- model.bootstrap(xy = dat, model = mod, err = input$method, k = 10)
      conf <- f.conf(mat, min(dat[, 1]), max(dat[, 1]), no.intervals = 10, conf.level = 0.95, lwd = 2, lty = 2)
      band <- data.frame(concentration = conf[,1], value_l = conf[,2], value_u = conf[,3])
      progress$set(message = "Rendering conf-band.", detail = "Please wait...", value = 1)
      return(band)
    } else {
      return(NULL)
    }
  })
  
  #' calculate ic values
  get.ic <- reactive({
    mod <- get.mod()
    dat <- get.dat()
    fit <- get.fit()
    if(is.null(mod) || is.null(dat) || is.null(fit)) {
      return(NULL)
    }
    #' calculate ic50 and ic90
    minv <- fit$par[model.params(mod) + 1]
    maxv <- fit$par[model.params(mod) + 2]
    ic01.y <- ((maxv - minv) * 0.01 + minv)
    ic50.y <- ((maxv - minv) * 0.50 + minv)
    ic90.y <- ((maxv - minv) * 0.90 + minv)
    ic99.y <- ((maxv - minv) * 0.99 + minv)
    ic01.x <- model.wrapper(f = "inv.", model = mod, x = 0.01, args = fit$par)
    ic50.x <- model.wrapper(f = "inv.", model = mod, x = 0.50, args = fit$par)
    ic90.x <- model.wrapper(f = "inv.", model = mod, x = 0.90, args = fit$par)
    ic99.x <- model.wrapper(f = "inv.", model = mod, x = 0.99, args = fit$par)
    #' estimate interval that includes inflection point
    a <- -1
    b <- 1
    if (ic01.x < ic99.x) {
      a <- ic01.x
      b <- ic99.x
    } else {
      a <- ic99.x
      b <- ic01.x
    }
    #' calculate mic and nic
    curve.mic.x <- curve.mic(model = mod, args = fit$par, ic_left = a, ic_right = b)
    curve.nic.x <- curve.nic(model = mod, args = fit$par, ic_left = a, ic_right = b)
    curve.mic.y <- model.wrapper(f = "f.", model = mod, x = curve.mic.x, args = fit$par)
    curve.nic.y <- model.wrapper(f = "f.", model = mod, x = curve.nic.x, args = fit$par)
    #' bind concentrations
    ic01 <- cbind(concentration = ic01.x, value = ic01.y)
    ic50 <- cbind(concentration = ic50.x, value = ic50.y)
    ic90 <- cbind(concentration = ic90.x, value = ic90.y)
    ic99 <- cbind(concentration = ic99.x, value = ic99.y)
    mic <- cbind(concentration = curve.mic.x, value = curve.mic.y)
    nic <- cbind(concentration = curve.nic.x, value = curve.nic.y)
    #' return concentrations
    return(list(mic=mic,nic=nic,ic50=ic50,ic90=ic90))
  })
  
  render.plot <- reactive({
    mod <- get.mod()
    dat <- get.dat()
    fit <- get.fit()
    ic <- get.ic()
    band <- get.band()
    if(is.null(mod) || is.null(dat) || is.null(fit) || is.null(ic)) {
      return(ggplot() + ggtitle("No data selected")) ## <-- ###### hier
    }
    
    cplot <- cplot.base(p.data = dat)
    if(input$show_band) {
      if(!is.null(band)){
        cplot <- cplot + cplot.layer.conf(band)
      }
    }
    if(input$show_sd) {
      sd_dat <- data.frame(concentration = unique(dat$concentration), value = unique(dat$concentration), sd = unique(dat$concentration))
      for(i in 1:nrow(sd_dat)){
        val <- subset(dat, concentration == sd_dat$concentration[i])$value
        sd_dat$value[i] <- mean(val)
        sd_dat$sd[i] <- sd(val)
      }
      cplot <- cplot + cplot.layer.sd(p.data = sd_dat)
      cplot <- cplot + cplot.layer.sdp(p.data = sd_dat)
    }
    if(input$show_data) {
      cplot <- cplot + cplot.layer.data(p.data = dat)
    }
    if(input$show_curve) {
      cplot <- cplot + cplot.layer.model(p.data = dat, fun = get(paste("f.",mod,sep="")), args = fit$par)
    }
    if(input$show_mic) {
      cplot <- cplot + cplot.layer.point(p.data = data.frame(ic$mic), symbol = 24)
    }
    if(input$show_nic) {
      cplot <- cplot + cplot.layer.point(p.data = data.frame(ic$nic), symbol = 25)
    }
    if(input$show_ic50) {
      cplot <- cplot + cplot.layer.point(p.data = data.frame(ic$ic50), symbol = 23)
    }
    if(input$show_ic90) {
      cplot <- cplot + cplot.layer.point(p.data = data.frame(ic$ic90), symbol = 26)
    }
    ######################## HIER
    plot.title <- ggtitle(dat$compound)
    center <- theme(plot.title = element_text(hjust = 0.5, size=22, color = "blue"))
    
    
    #grid.x <- scale_y_continuous(name = "value", labels = math_format(10^.x)) 
    # y labels and name changed by Raimo
    grid.x <- scale_y_continuous(name = input$y_axis)#"growth OD600") 
    grid.y <- scale_x_continuous(name = input$x_axis, labels = math_format(10^.x))  #expression(paste("conc [", mu, "M]"))
    
    return(cplot + plot.title + center + annotation_logticks() + grid.x + grid.y)
    ##########################
  })
  render.summary <- reactive({
    mod <- get.mod()
    dat <- get.dat()
    fit <- get.fit()
    ic <- get.ic()
    if(is.null(mod) || is.null(dat) || is.null(fit) || is.null(ic)) {
      return(NULL)
    }
    
    summary <- list(data=c(),model=c(),parameters=c(),concentration=c(),settings=c())
    
    summary$data <- c(df(), input$compound)
    names(summary$data) <- c("dataset","compound")
    
    summary$settings <- c(input$model,input$method,input$sub_blank)
    names(summary$settings) <- c("model","method","sub_blank")
    
    summary$model <- c(paste("f.",mod,sep=""),fit$value)
    names(summary$model) <- c("name","rss")
    
    if(model.params(model = mod) == 2) {
      summary$parameters <- c(fit$par[1], fit$par[2], NA, fit$par[3], fit$par[4])
    } else {
      summary$parameters <- c(fit$par[1], fit$par[2], fit$par[3], fit$par[4], fit$par[5])
    }
    names(summary$parameters) <- c("beta1","beta2","beta3","min","max")
    
    # summary$concentration <- matrix(c(ic$mic[1], ic$nic[1], ic$ic50[1], ic$ic90[1], ic$mic[2], ic$nic[2], ic$ic50[2], ic$ic90[2]), nrow=2, ncol=4, byrow = TRUE)
    # dimnames(summary$concentration) <- list(c("concentration [10^]","value [10^]"),c("MIC","NIC","IC50","IC90"))
    # #changes made by Raimo with 10^
    summary$concentration <- matrix(c(10^ic$mic[1], 10^ic$nic[1], 10^ic$ic50[1], 10^ic$ic90[1], ic$mic[2], ic$nic[2], ic$ic50[2], ic$ic90[2]), nrow=2, ncol=4, byrow = TRUE)
    dimnames(summary$concentration) <- list(c("concentration","value"),c("MIC","NIC","IC50","IC90"))
    summary$show <- c(input$show_data,input$show_band,input$show_curve,input$show_sd,input$show_mic,input$show_nic,input$show_ic50,input$show_ic90)
    names(summary$show) <- c("show_data","show_band","show_curve","show_sd","show_mic","show_nic","show_ic50","show_ic90")
    
    return(summary)
  })
  
  
  observe({
    comp <- as.vector(unique(df()$"compound"))
    updateSelectInput(session, "compound",
                      choices = comp[!is.na(comp)])
  })
  observe({
    if(input$Button>0){
      output$plot <-renderPlot(isolate(render.plot()))
      output$summary <-renderPrint(isolate(render.summary()))
      }
    })
  
  observe({
  output$show <-renderDataTable(isolate(data.load("../db", input$inspect))) #, options = list(paging = TRUE, searching = FALSE)
  })
  
  observe({
    if (is.null(input$file) || is.null(input$filename)) {
      return()
    } else{
      data.save("../db",
                paste(input$filename, "dat", sep = "."),
                data.extract(input$file$datapath))
      updateSelectInput(session, "dataset",
                        choices = data.list("../db"))
      updateSelectInput(session, "inspect",
                        choices = data.list("../db"))
    }
  })
  
  observe({

    output$download <- downloadHandler(
      filename = function() {
        paste("plot", ".png", sep = "")
      },
      content = function(file) {
        ggsave(file,  render.plot())


      })
 })
    
  
  
})