# ################################################## data.r

data.extract <- function(fpath) {
   # wb <- #XLConnect::loadWorkbook(fpath)
    
    value.xls <- read_excel(path = fpath, sheet = 1, skip = 0,col_names=F,.name_repair = "minimal")
    concentration.xls <- read_excel(path = fpath, sheet = 2, skip = 0,col_names=F,.name_repair = "minimal")
    compound.xls <- read_excel(path = fpath, sheet = 3, skip = 0,col_names=F,.name_repair = "minimal")
    names(value.xls) <- NULL
    names(concentration.xls) <- NULL
    names(compound.xls) <- NULL
    value<-unlist(value.xls)
    concentration<-unlist(concentration.xls)
    compound<-unlist( compound.xls)
    
    
    concentration <- as.numeric(concentration)
    concentration[is.na(concentration)] <- 0
    concentration <- sapply(concentration, log10)

    df <- data.frame(concentration, value, compound)
    df <- df[complete.cases(df), ]
    df$concentration[which(df$compound == "BLANK")] <- NA
    df$compound[which(df$compound == "BLANK")] <- NA
    #' df$value[is.na(df$compound)] <- df$value[is.na(df$compound)] * -1
    return(df)
}

data.save <- function(fdir, fname, fdata) {
    write.table(x = fdata, file = file.path(fdir, fname), sep = ",", row.names = FALSE, col.names = TRUE)
}

data.load <- function(fdir, fname) {
    data = read.table(file = file.path(fdir, fname), sep = ",", header = TRUE)
    return(data)
}

data.list <- function(fdir) {
    files <- list.files(path = fdir)
    return(files)
}


# ################################################## End of Document
