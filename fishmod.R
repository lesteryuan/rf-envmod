## explore biological metrics values
## questions for IN:
## How is biological data used currently in water programs?

## model relationship between richness and environmental
## variables

richmod <- function(envdat.tss, envdat.cond, pred, ibi) {


    varpred <- names(pred)[-c(1:4)]

#    mname <- "CATCHPERUNITEFFORT"
#    mname2 <-"CATCHPERUNITEFFORT_METRIC"
#    mname <- "SPECIESCOUNT"
#    mname2 <- "SPECIES_COUNT_METRIC"
#    mname <- "TOLERANT_PERCENT"
#    mname2 <- "TOLERANT_PERCENT_METRIC"
    mname <- "SIMPLELITHOPHIL_PERCENT"
    mname2 <- "SIMPLELITHOPHIL_PERCENT_METRIC"
    dftemp <- subset(ibi, metric.name == mname)
    names(dftemp)[4] <- mname
#    dftemp <- subset(ibi, metric.name == "TOLERANT_PERCENT")
    print(nrow(dftemp))
    print(length(unique(dftemp$event.num)))
    dftemp <- merge(dftemp, unique.data.frame(fishcnt[,c("event.num",
                                                   "station.name")]),
                 by = "event.num")
    print(nrow(dftemp))
    dftemp2 <- subset(ibi, metric.name == mname2)
    names(dftemp2)[4] <- mname2
    dftemp <- merge(dftemp, dftemp2, by = "event.num")
    print(summary(dftemp))

#    dftemp[, mname] <- log(dftemp[, mname])

    ## merge in envdata
    names(envdat.tss) <- paste(names(envdat.tss), "tss", sep = ".")
    names(envdat.cond) <- paste(names(envdat.cond), "cond", sep = ".")
    df1 <- merge(dftemp,envdat.tss,
                 by.y = "LSite.tss", by.x = "station.name", all.x = F,
                 all.y = F)
    df1 <- merge(df1, envdat.cond, by.x = "station.name", by.y = "LSite.cond")
    print(nrow(df1))
    df1 <- merge(df1, unique.data.frame(pred[, c("LSite", varpred)]),
                 by.y = "LSite", by.x = "station.name")
    print(nrow(df1))

    ## drop forest variables and kf factor
    varpred <- varpred[-which(varpred == "kffactws")]
    varpred <- varpred[-which(varpred == "pctdecid2019ws")]
    varpred <- varpred[-which(varpred == "pctconif2019ws")]
    varpred <- varpred[-which(varpred == "pctmxfst2019ws")]
    ## add local env variables
    varpred <- c(varpred, "predicted.tss", "predicted.cond")

    print(varpred)

    ## names of anthropogenic variables to plot for each subset
#    varplot <- list(c("inorgnwetdep.2008ws", "pctcrop2019ws",
#                      "pctnonagintrodmanagvegws"),
#                    c("tss.pred", "pctconif2019ws"))
#    lab0 <- list(c("Wet N deposition", "Percent crop",
#                   "Percent introduced managed vegetation"),
#                 c("Predicted TSS", "Percent coniferous forest"))
    varplot <- c("predicted.cond","predicted.tss")
    lab0 <- c("Conductivity", "TSS")

    require(pdp)
    dev.new()
#    par(mar = c(4,4,2,1), mfrow = c(1,2), mgp = c(2.3,1,0),
#        bty = "l")
    par(mar = c(4,4,1,1), mfrow = c(1,3), mgp = c(2.3,1,0),
        bty = "l")

    ## set up bootstrap
    isamp <- sample(nrow(df1))
    print(nrow(df1))
    nval <- nrow(df1)/10
    ip <- round(seq(1, nrow(df1), length = 11))
    ip2 <- ip-1
    ip2 <- ip2[-1]
    ip2[length(ip2)] <- nrow(df1)
    print(ip)
    print(ip2)

    ## set up strata
    cutp <- quantile(df1$wsareasqkm,prob = seq(0, 1, length = 6))
    cutf <- cut(df1$wsareasqkm, cutp, include.lowest = T)
    dflist <- split(df1, cutf)

    ## roll back all human variables to zero
    new.data <- df1
    varset <- c("pctimp2019ws", "pcthay2019ws","pctcrop2019ws",
                "pctbl2019ws", "pcturbop2019ws", "inorgnwetdep.2008ws",
                "pctnonagintrodmanagvegws", "pctag2006slp10ws",
                "npdesdensws")
    for (i in varset) {
        new.data[,i] <- 0
    }
    new.data$predicted.tss <- new.data$predref.tss
    new.data$predicted.cond <- new.data$predref.cond

    pred <- rep(NA, times = nrow(df1))
    predref <- rep(NA, times = nrow(df1))
    set.seed(10)
    require(ranger)
    require(pdp)
    for (i in 1:10) {
        samp0 <- isamp[ip[i]:ip2[i]]
        mod <- ranger(data = df1[-samp0, c(varpred, mname)],
                      dependent.variable.name = mname,
                      num.trees = 5000, importance = "permutation")
        pred[samp0] <- predict(mod, df1[samp0,])$predictions
        predref[samp0] <- predict(mod, new.data[samp0,])$predictions
    }
    print(mod)
    print(rev(sort(mod$variable.importance))[1:20])
    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(2,2))

    plot(log(df1$wsareasqkm), df1[,mname], pch = 21, col = "grey39",
         bg = "white")
    levs <- c(1,5)
    col0 <- c("red", "green")
    for (i in 1:2) {
        incvec <- df1[, mname2] == levs[i]
        points(log(df1$wsareasqkm)[incvec], df1[incvec, mname],
               pch = 21, col ="grey39", bg = col0[i])
    }

    plot(pred, df1[, mname])
    abline(0,1)
    rms0<- sqrt(sum((log(pred) - log(df1[,mname]))^2)/length(pred))

    plot(predref, pred)
    abline(0, 1, lty = "dashed")

    x <- df1$predicted.tss - new.data$predicted.tss
#    x <- log(df1$npdesdensws)
#    x <- df1$predicted.cond - new.data$predicted.cond
    plot(x, pred-predref)
    stop()


    y <- pred - predref
    print(summary(x))
    print(summary(y))

    plot(x,y)
    require(mgcv)
    mod <- gam(y ~ s(x, k = 3))
    predout <- predict(mod, type = "response", se.fit = T)
    iord <- order(x)
    up <- predout$fit + 2*predout$se.fit
    dn <- predout$fit - 2*predout$se.fit

    lines(x[iord], predout$fit[iord])
    lines(x[iord], up[iord], lty ="dashed")
    lines(x[iord], dn[iord], lty ="dashed")
    abline(h = 0, lty = "dashed")
    incvec <- df1[, mname2] == 1
    points(x[incvec], y[incvec], pch = 16, col = "red")

    plot(df1$predicted.tss, df1[, mname], pch = 21,
         col ="grey38", bg = "white")

#    points(df1$predicted.tss[incvec], df[incvec, mname],
#           pch = 16, col= "red")
    stop()


    mod <- ranger(data = df1[, c(varpred, "metric.value")],
                  dependent.variable.name = "metric.value",
                  num.trees = 5000, importance = "permutation")
    print(mod)
    err0 <- mod$prediction.error
    print(rev(sort(mod$variable.importance))[1:10])
    varsel <- names(rev(sort(mod$variable.importance)))
    pout <- partial(mod, "predicted.tss", plot = F)
    plot(pout[,1], pout[,2], type = "l")
    stop()


    qval <- c(0.1, 0.3, 0.5, 0.7, 0.9)
    predout <- matrix(NA, ncol = length(qval), nrow = 40)
    xnew <- matrix(NA, ncol = length(qval), nrow = 40)
    for (k in 1:length(qval)) {
        xnew[,k] <- seq(quantile(dflist[[k]]$predicted.tss, prob = 0.05),
                    quantile(dflist[[k]]$predicted.tss, prob = 0.95),
                    length = 40)
        dftemp <- data.frame(predicted.tss = xnew[,k])
        varrest <- varpred
        varrest <- varrest[-which("predicted.tss" == varrest)]
        matnew <- matrix(NA, ncol = length(varrest), nrow = nrow(xnew))
        for (i in 1:length(varrest)) {
            if (varrest[i] == "wsareasqkm") {
                matnew[,i] <- rep(quantile(df1[, varrest[i]], prob = qval[k]),
                                  times = nrow(xnew))
            }
            else {
                matnew[,i] <- rep(median(df1[, varrest[i]]),
                                  times = nrow(xnew))
            }
        }
        dimnames(matnew)[[2]] <- varrest
        dftemp <- cbind(dftemp, matnew)
        predout[,k] <- predict(mod, dftemp)$predictions
        predout[,k] <- predout[,k] - mean(predout[,k])
    }
    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(1,1))
    xlim <- range(xnew)
    ylim <- range(predout)
    plot(xnew[,1], predout[,1], xlim = xlim, ylim = ylim, type = "l")
    for (i in 2:5) {
        lines(xnew[,i], predout[,i])
    }

    stop()


    stop()

#    png(width = 4,height = 4, pointsize = 10, units = "in", res = 600,
#        file = "plot1.png")
    dev.new()
    par(mar = c(4,4,1,1), mgp = c(2.3,1,0), bty = "l")
    plot(predref, pred, xlab = "EPT reference prediction",
         ylab ="EPT current", pch = 21, col = "grey39", bg = "white")
    abline(0,1, lty = "dashed")
    incvec <- (predref-pred) > 1
    points(predref[incvec], pred[incvec], pch = 21, col = "grey39", bg = "blue")
#    hist(pred/predref, breaks = 10)
    dev.off()


    ## set plotpred to TRUE to plot predicted vs. observed
    ## otherwise, will plot partial dependence relationships
    ## for anthropogenic variables
    plotpred <- F
    if (plotpred) {
        png(width = 4, height = 4, pointsize = 10, units = "in",
            res = 600, file = "plot.png")
        par(mar = c(4,4,1,1), mgp = c(2.3,1,0), bty = "l")
        plot(pred, df1$ept,
             xlab = "EPT(predicted)", ylab = "EPT(observed)",
             pch = 21, col = "grey", bg = "white")
        abline(0,1, lty = "dashed")
        dev.off()
    }
    else {
        png(width = 6, height = 2.5, pointsize = 10, units = "in",
            res = 600, file = "plot.png")
        par(mar = c(4,4,1,1), mgp = c(2.3,1,0), bty = "l", mfrow = c(1,2))
        for (k in 1:length(varplot)) {
            pout <- partial(mod, pred.var = varplot[k], plot = F,
                            approx = T)
            print(names(pout))
            plot(pout[,1], pout[,2] - mean(pout[,2]), type = "l",
                 xlab = lab0[k], ylab = "Change in EPT", axes = F)
            logtick.exp(0.001, 10, c(1), c(T,F))
            axis(2)
        }
        dev.off()
    }



}

richmod(envdat.tss = dftss, envdat.cond = dfcond, pred = pred,
        ibi = ibi1)
