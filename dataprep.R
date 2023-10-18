## 10.17.2023
## get predictor variable names for sample location
getptdata <- function(df1, latname, lonname) {
    require(maps)
    require(mapproj)
    require(mgcv)

    ## load PRISM data
    load("grid.data.rda")

    ## function to interpolate values from gridded data
    ## from PRISM
    interpairtemp <- function(xlat, xlon, pdat.all) {
        lat <- pdat.all$lat
        lon <- pdat.all$lon
        prismdat <- pdat.all$prismdat
        airtemp <- rep(NA, times = length(xlat))
        for (i in 1:length(xlat)) {
            ilat <- which.min(abs(xlat[i]-lat))
            ilon <- which.min(abs(xlon[i]-lon))
            airtemp[i] <- prismdat[ilat, ilon]
        }
        return(airtemp)
    }
    ## calculate mean air temperature and precip at new locations
    df1$airtemp  <- interpairtemp(df1[, latname], df1[, lonname], pdat.all)
    df1$ppt <- interpairtemp(df1[, latname], df1[, lonname], ppt)

    ## project and map new locations
    dev.new()
    par(mar = c(1,1,1,1))
    map("state", proj = "albers", par=c(30,40))
    pout <- mapproject(df1[,lonname], df1[,latname], proj = "")
    points(pout$x, pout$y)

    ## calculate BFI at new locations
    load("mod.bfi.rda")
    new.data <- data.frame(xloc = pout$x, yloc = pout$y,
                           area = log(df1$wsareasqkm))
    df1$bfi.gam  <- predict(mod.bfi, new.data, type = "response")

    return(df1)
}

## download data from StreamCat
getscdata <- function(df1, comidname) {
    require(StreamCatTools)

    ## load names of predictor variables
    load("predvar.sc.rda")

    ## predictor variables
    predvar2 <- unique(unlist(predvar.sc))

    ## correct field names from 2018-19 data to
    ## field names in streamcat
    orig.name <- c("msst.2019", "inorgnwetdep.2008ws", "elevation")
    rev.name <- c("msst_2014", "inorgnwetdep_2008ws", "elevcat")
    for (i in 1:length(orig.name)) {
        ip <- which(predvar2 == orig.name[i])
        predvar2[ip] <- rev.name[i]
    }

    ## drop variables missing from streamcat
    droplist <- c("airtemp", "ppt","bfi.gam", "wsareasqkm")

    for (i in droplist) {
        ip <- which(i == predvar2)
        if (length(ip) == 1) {
            predvar2 <- predvar2[-ip]
        }
    }

    ## find stem of variable name for querying streamcat
    ## this drops 'ws' from the names of certain variables
    predvar.q <- predvar2
    stem.ws <- regexpr("ws", predvar.q) != -1
    predvar.q[stem.ws] <- substring(predvar.q[stem.ws], 1,
                                    nchar(predvar.q[stem.ws])-2)
    stem.cat <- regexpr("cat", predvar.q) != -1
    predvar.q[stem.cat] <- substring(predvar.q[stem.cat], 1,
                                     nchar(predvar.q[stem.cat])-3)
    names(stem.ws) <- predvar.q
    names(stem.cat) <- predvar.q

    ## download predictors from streamcat
    for (i in predvar.q) {
        df <- sc_get_data(metric=i, comid=df1[, comidname])
        if (!is.null(df)) {
            names(df) <- tolower(names(df))
            if (stem.ws[i]) {
                ip <- which(paste(i, "ws", sep = "")== names(df))
                ## wsarea is provided for all metrics so just
                ## save it once for pctcrop
                if (i == "pctcrop2019") {
                    ipw <- which("wsareasqkm" == names(df))
                    df1 <- merge(df1, df[, c(1,ip, ipw)], by = "comid",
                                 all.x = T)
                }
                else
                    df1 <- merge(df1, df[, c(1,ip)], by = "comid",
                                 all.x = T)
            }
            else {
                if (stem.cat[i]) {
                    ip <- which(paste(i, "cat", sep = "")== names(df))
                    df1 <- merge(df1, df[, c(1,ip)], by = "comid",
                                 all.x = T)
                }
                else
                    df1 <- merge(df1, df[, c("comid", i)], by = "comid", all.x = T)
            }
        }
        else {
            print("not found")
        }
    }
    ## check for missing predictors in df1
    checkmiss <- FALSE
    if (checkmiss) {
        for (i in predvar.q) {
            ip <- which(regexpr(i, names(df1)) != -1)
            ismiss <- is.na(df1[, ip])
            if (sum(ismiss) > 0) {
                cat("Number of missing samples for",i, ":", sum(ismiss), "\n")
            }
        }
    }

    ## rename variables to be consistent with 18-19
    for (i in 1:length(rev.name)) {
        ip <- which(names(df1) == rev.name[i])
        names(df1)[ip] <- orig.name[i]
    }
    return(df1)
}

dat <- read.delim("example.dat.txt")
dfenv <- getscdata(dat, comidname = "comid")
dfenv <- getptdata(dfenv, latname = "lat.dd83", lonname = "lon.dd83")

