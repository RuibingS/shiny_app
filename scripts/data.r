# ################################################## data.r

data.extract <- function(fpath) {
    wb <- XLConnect::loadWorkbook(fpath)
    
    value <- XLConnect::readWorksheet(object = wb, sheet = 1, header = FALSE, simplify = TRUE)
    concentration <- XLConnect::readWorksheet(object = wb, sheet = 2, header = FALSE, simplify = TRUE)
    compound <- XLConnect::readWorksheet(object = wb, sheet = 3, header = FALSE, simplify = TRUE)
    
    #' Convert concentration to numeric.
    #' We need to replace na with 0, so the
    #' BLANK values dont get deleted in the last step.
    concentration <- as.numeric(concentration)
    concentration[is.na(concentration)] <- 0
    concentration <- sapply(concentration, log10)
    #' concentration[is.infinite(concentration)] <- 0
    #' value <- sapply(value, data.divh)
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
