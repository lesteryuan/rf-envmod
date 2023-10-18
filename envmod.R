## Lester Yuan
## Model local environmental conditions with national Random Forest model

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

envout <-  envmod(dfenv, idname ="uid")

