## explore biological metrics values
## questions for IN:
## How is biological data used currently in water programs?
exploremet <- function() {
    ## figure out unique sample id
    richness.in$idnew <- paste(richness.in$lsite, richness.in$date, sep = "-")

    ## need to ask mitchell about repeat samples on the same day

    ## pick first sample
    richness.in <- richness.in[!duplicated(richness.in$idnew),]

    metrics$idnew <- paste(metrics$lsite, metrics$date, sep = "-")
    nsamp <- table(metrics$idnew)
    idp <- names(nsamp)[nsamp > 1]
    incvec <- rep(F, times = nrow(metrics))

    for (i in idp) incvec <- incvec|metrics$idnew == i

    df1 <- metrics[incvec,]
    print(summary(tapply(df1$individuals.count, df1$idnew, function(x) diff(range(x)))))
    stop()
    incvec <- metrics$idnew == idp[1]
    print(metrics[incvec,])

    plot(log(metrics$individuals.count), log(metrics$taxa.count))
    stop()

}
#exploremet()

## model relationship between richness and environmental
## variables

richmod <- function(envdat.tss, envdat.cond, pred, richdat) {


    varpred <- names(pred)[-c(1:4)]

    ## calculate mean values of data from the same site and merge data
    trich <- tapply(richdat$TAXA_COUNT, richdat$LSite, mean)
    pintol <- tapply(metrics$pct.intolerant, metrics$lsite, mean)
    pintol.sd <- tapply(metrics$pct.intolerant, metrics$lsite, sd)
    ept <- tapply(metrics$ept.index, metrics$lsite, mean)


    tol <- tapply(metrics$pct.tolerant, metrics$lsite, mean)

    dftemp <- data.frame(LSite = names(trich), trich=as.vector(trich))
    print(nrow(dftemp))
    dftemp <- merge(dftemp, data.frame(LSite = names(pintol),
                                       pintol = as.vector(pintol),
                                       pintol.sd = as.vector(pintol.sd),
                                       ept = as.vector(ept)),
                    by = "LSite")
    dftemp <- merge(dftemp, data.frame(LSite = names(tol),
                                       tol = as.vector(tol)),
                    by = "LSite")

    ## merge in envdata
    names(envdat.tss) <- paste(names(envdat.tss), "tss", sep = ".")
    names(envdat.cond) <- paste(names(envdat.cond), "cond", sep = ".")
    df1 <- merge(dftemp,envdat.tss,
                 by.y = "LSite.tss", by.x = "LSite")
    print(nrow(df1))
    df1 <- merge(df1, envdat.cond, by.x = "LSite", by.y = "LSite.cond")
    print(nrow(df1))
    df1 <- merge(df1, unique.data.frame(pred[, c("LSite", varpred)]),
                 by = "LSite")
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

    ## roll back all human variables to zero
    new.data <- df1
    varset <- c("pctimp2019ws", "pcthay2019ws",#"pctcrop2019ws",
                "pctbl2019ws", "pcturbop2019ws",# "inorgnwetdep.2008ws",
                "pctnonagintrodmanagvegws", #"pctag2006slp10ws",
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
    for (i in 1:10) {
        samp0 <- isamp[ip[i]:ip2[i]]
        mod <- ranger(data = df1[-samp0, c(varpred, "ept")],
                      dependent.variable.name = "ept",
                      num.trees = 5000, importance = "permutation")
#        print(mod)
#        err0 <- mod$prediction.error
#        varsel <- names(rev(sort(mod$variable.importance)))
#    print(varsel[1:20])

        pred[samp0] <- predict(mod, df1[samp0,])$predictions
        predref[samp0] <- predict(mod, new.data[samp0,])$predictions
    }


    png(width = 4,height = 4, pointsize = 10, units = "in", res = 600,
        file = "plot1.png")
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
        richdat = richness.in)
