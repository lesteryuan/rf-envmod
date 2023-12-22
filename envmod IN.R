## Lester Yuan
## Model local environmental conditions with national Random Forest model

## 12.21.2023
## rerun models, calibrating with IN observations

envmod_IN <- function(envdata, pred, var = "cond") {
    require(ranger)
    varlist <- names(pred)[5:length(names(pred))]

    df1 <- envdata

    ## drop forest...correlated with human activity
#    varlist <- varlist[-which(varlist == "pctdecid2019ws")]
#    varlist <- varlist[-which(varlist == "pctmxfst2019ws")]
#    varlist <- varlist[-which(varlist=="pctsilicicws")]
    ## drop kffact for now because of incomplete coverage
    varlist <- varlist[-which(varlist=="kffactws")]


    if (var == "cond") {
        df1$cond <- log(df1$Conductance..uS.cm.)
        ## drop 6 outliers
        incvec <- df1$cond > 3.5
        cond.mn <- tapply(df1$cond, df1$STATION_NAME, mean)
        cond.sd <- tapply(df1$cond, df1$STATION_NAME, sd)

        nsamp <- tapply(df1$cond, df1$STATION_NAME, function(x) sum(!is.na(x)))
        incvec <- nsamp > 0
        cond.mn <- cond.mn[incvec]

        df2 <- data.frame(STATION_NAME = names(cond.mn),
                          cond = as.vector(cond.mn))
    }
    else {
        if (var =="tss") {
            ##set everything with a "<" to NA
            incvec <- regexpr("<", df1$TSS..mg.l.) != -1
            print(sum(incvec))
            df1$TSS..mg.l.[incvec] <- NA
            ## set every measurement less than 10 to NA
            df1$TSS..mg.l. <- as.numeric(df1$TSS..mg.l.)
            print(summary(df1$TSS..mg.l.))
            df1 <- df1[!is.na(df1$TSS..mg.l.),]
            df1$tss <- log(df1$TSS..mg.l.)

            require(lme4)
            mod <- lmer(tss ~ 1 + (1|STATION_NAME), data = df1)
            print(summary(mod))

            tss.mn <- tapply(df1$tss, df1$STATION_NAME, mean)
            df2 <- data.frame(STATION_NAME = names(tss.mn),
                              tss = as.vector(tss.mn))
#            print(sum(!is.na(df1$TSS..mg.l)))
#            hist(log(df1$TSS..mg.l))
        }
    }
    pred <- unique.data.frame(pred[, c("LSite", varlist)])

    docluster <- F
    if (docluster) {
        ## cluster predictor variables
        dftemp <- na.omit(pred[, varlist])
        cor0 <- 1-abs(cor(dftemp[, varlist], method = "spearman"))

        require(cluster)
        clust0 <- agnes(cor0, diss = T, method = "average")

        lab0 <- varlist
        dev.new()
        par(mar = c(2,1,1,4), mgp = c(2.3,1,0))
        plot(clust0,labels = lab0,  which.plots = 2,
         main = "", ylab = "" , xlab = "", sub = "", axes
             = F)
        axis(4, at = seq(0, 1, by = 0.2), lab = seq(1, 0, by = -0.2))
        mtext(expression(abs(italic(r[s]))),side = 4, line = 2.3)
    }

    print(varlist)
    print(nrow(df2))
    df1 <- merge(df2, pred, by.y = "LSite", by.x = "STATION_NAME")
    print(nrow(df1))

    ## drop missing kfact
#    incvec <- !is.na(df1$kffactws)
#    df1 <- df1[incvec,]
#    print(summary(df1))

    ## latitude stratification
    plot(df1$LATITUDE_MEASURE, df1$rckdepws, xlab = "Latitude",
         ylab = "Bedrock depth", pch = 21, col = "grey39", bg = "white")
    abline(v = 39.6)
    incvec <- df1$LATITUDE_MEASURE < 39.6

    dflist <- split(df1, incvec)
    yall <- c(dflist[[1]][, var], dflist[[2]][, var])
    predlist <- as.list(rep(NA, times = length(dflist)))
    modlist <- as.list(rep(NA, times = length(dflist)))

    for (i in 1:length(dflist)) {
        modlist[[i]] <- ranger(data = dflist[[i]][, c(var, varlist)],
                      dependent.variable.name = var,
                      num.trees = 5000, importance = "permutation")

        predlist[[i]] <- modlist[[i]]$predictions
    }
    predall <- c(predlist[[1]], predlist[[2]])

    dev.new()
    par(mar = c(4,4,1,1), mgp = c(2.3,1,0), bty = "l")
    plot(predall, yall, pch = 21, col = "grey", bg = "white", axes = F,
         xlab = "Predicted", ylab ="Observed")
    logtick.exp(0.001, 10, c(1,2), c(T,T))
    abline(0,1)
    print(mean((predall - yall)^2))

#    vimport <- mod$variable.importance
#    vimport <- sort(vimport)
#    print(vimport)
    ## data frame with human activities set to zero
    varadj <- c("pctimp2019ws", "pcturbop2019ws", "pctcrop2019ws")
    predref <- as.list(rep(NA,times=2))
    for (i in 1:length(dflist)) {
        new.data2 <- dflist[[i]]
        for (j in varadj) new.data2[,j] <- 0
        ## plot "ref" vs. observed difference
        predref[[i]] <- predict(modlist[[i]], new.data2)$predictions
    }
    predref.all <- c(predref[[1]], predref[[2]])

    require(maps)
    require(mapproj)
    dev.new()
    par(mar = c(1,1,1,1))
    map("state", region = "indiana", proj = "albers", par = c(30,40))
    delt <- exp(predall- predref.all) - 1
    cex0 <- delt
    cex0[cex0 < 0] <- 0
    cex0 <- cex0*4 + 0.1

    pout <- mapproject(c(dflist[[1]]$LONGITUDE_MEASURE,
                         dflist[[2]]$LONGITUDE_MEASURE),
                       c(dflist[[1]]$LATITUDE_MEASURE,
                         dflist[[2]]$LATITUDE_MEASURE), proj = "")
    points(pout$x, pout$y, pch = 21, col = "grey39", bg = "white",
           cex = cex0)
    stop()

    ## plot partial dependence relationship
    ## set up dataframe for calculating partial dependence relationship
    new.data <- data.frame(matrix(NA, ncol = length(varlist), nrow = 40))
    names(new.data) <- varlist
    for (i in varlist) {
        new.data[, i] <- median(df1[, i])
    }
    varpick <- "pctcrop2019ws"
    new.data[, varpick]<- seq(min(df1[,varpick]),
                                  max(df1[,varpick]), length = 40)

    dev.new()
    predout <- predict(mod, new.data)
    plot(new.data[, varpick], predout$predictions, xlab = varpick)



}

envmod_IN(condpred, pred, var = "cond")
#envmod_IN(tsspred, pred, var = "tss")

envmod <- function(dfenv, idname) {
    require(ranger)
    load("predvar.sc.rda")
    load("envdata1819.rda")
    envdata <- envdata1819

    ## select wadeable sites from calibration data
    envdata <- envdata[envdata$protocol == "WADEABLE",]

    ## logtransforms of modeled envioronmental variables
    varlist <- c("doc", "ptl.diss",  "tss.result", "cond")
    for (i in varlist) {
        incvec <- envdata[,i] <= 0
        incvec[is.na(incvec)] <- F
        minval <- min(envdata[!incvec,i], na.rm = T)

        envdata[incvec,i] <- minval*0.5
        envdata[, i] <- log(envdata[,i])
    }

    ## log transforms of some of the predictors
    ## to make the extrapolation test more robust
    predlog <- c("caows", "elevation", "hydrlcondws", "nws", "omws",
                 "permws", "ppt", "runoffws", "wsareasqkm")
    for (i in predlog) {
        if (sum(envdata[,i] <=0, na.rm = T) > 0) {
            incvec <- envdata[,i] <= 0
            incvec[is.na(incvec)] <- F
            minval <- min(envdata[!incvec,i], na.rm = T)
            envdata[incvec,i] <- 0.5*minval
        }
        envdata[,i] <- log(envdata[,i])
        if (sum(dfenv[,i] <=0, na.rm = T) > 0) {
            incvec <- dfenv[,i] <= 0
            incvec[is.na(incvec)] <- F
            minval <- min(dfenv[!incvec,i], na.rm = T)
            dfenv[incvec,i] <- 0.5*minval
        }
        dfenv[,i] <- log(dfenv[,i])
    }

    ## scale variables that are used for prediction for assessing extrapolation
    name0 <- unique(unlist(predvar.sc))
    for (i in name0) {
        if (substring(i, 1,3) == "pct") {
            envdata[,i] <- envdata[,i]/100
            dfenv[,i] <- dfenv[,i]/100
        }
        else {
            min0 <- min(envdata[,i], na.rm = T)
            r0 <- range(envdata[,i], na.rm = T)
            envdata[,i] <- (envdata[,i] - min0)/diff(r0)
            dfenv[,i] <- (dfenv[,i]-min0)/diff(r0)
        }
    }

    svar <- c("pct.safn", "cond", "tss.result", "ptl.diss")
    svarlab <- c("% sand/fines", "Conductivity", "TSS",  "Pdiss")

    ## random forest model for environment
    envmod <- function(df1,y,dfout.ref, dfout.test) {
        mod <- ranger(data = df1,
                      dependent.variable.name = y,
                      num.trees = 5000, importance = "permutation",
                      mtry = floor(ncol(df1)/3))
        mse <- mod$prediction.error
#        print(sqrt(mse)/diff(range(df1[,y])))
        predout <- mod$predictions

        pred.test <- predict(mod, dfout.test)$predictions
        pred.ref <- predict(mod, dfout.ref)$predictions
        imp <- mod$variable.importance
        imp <- imp/max(imp)

        ## scale predictors by importance
        for (i in names(imp)) {
            df1[,i] <- df1[,i]*imp[i]
            dfout.ref[,i] <- dfout.ref[,i]*imp[i]
            dfout.test[,i] <- dfout.test[,i]*imp[i]
        }
        ## get distances
        n1 <- nrow(df1)
        n2 <- nrow(dfout.ref)
        n3 <- nrow(dfout.test)
        d <- as.matrix(dist(rbind(df1[, 2:ncol(df1)], dfout.ref, dfout.test)))
        mat1 <- d[1:n1, 1:n1]
        meand <- mean(mat1[lower.tri(mat1)])
        mat2 <- d[1:n1, (n1+1):(n1+n2)]
        dref <- apply(mat2, 2, min)/meand
        mat3 <- d[1:n1, (n1+n2+1):(n1+n2+n3)]
        dtest <- apply(mat2, 2, min)/meand

        dfout <- data.frame(pred.ref = pred.ref,
                            pred.test = pred.test,
                            dref = dref,
                            dtest = dtest)
        return(dfout)
    }

    ## generate environmental predictions for new sites
    set.seed(10)
    envtest <- dfenv
    envref <- dfenv
    ## set human land use variables to a low percentile
    ## of all sites in the new data set
    selref <- envref$rt.nrsa == "R"
    humpct <- c("pctcrop2019ws", "pctimp2019ws", "pcthay2019ws", "npdesdensws",
                "pcturbop2019ws", "pctag2006slp10ws",
                "pctnonagintrodmanagvegws")

    pvals <- 0.25
    md0 <- apply(envref[, humpct], 2,quantile, prob = pvals,  na.rm = T)
    cat("Reference levels for human disturbance variables:\n")
    print(md0)
    for (j in 1:length(md0)) {
        envref[, names(md0)[j]] <- md0[j]
    }

    for (i in 1:length(svar)) {
        varlist <- predvar.sc[[svar[i]]]

        ## drop missing variables for response
        envdata0 <- na.omit(envdata[, c("uid", varlist, svar[i])])

        envtest0 <- na.omit(envtest[, c(idname, varlist)])
        envref0 <- na.omit(envref[, c(idname, varlist)])

        xout <- envmod(envdata0[, c(svar[i], varlist)],svar[i],
                       envref0[, varlist], envtest0[, varlist])

        xout$id <- envtest0[, idname]
        ## put meaningful names on xout
        names(xout) <- c(paste(svar[i], "ref", sep = "_"),
                         paste(svar[i], "test", sep = "_"),
                         paste("d",svar[i], "ref", sep = "_"),
                         paste("d", svar[i], "test", sep = "_"), idname)
        if (i==1) {
            dfout <- xout
        }
        else {
            dfout <- merge(dfout, xout, by = idname, all.x = T, all.y = T)
        }
    }
    return(dfout)

}

#envout <-  envmod(dfenv, idname ="uid")

