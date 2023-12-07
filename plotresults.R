## plot environmental model results

plotresults <- function() {

#    mcodat <- read.csv("MCO 2022 sites model test.csv")
    mcodat <- read.csv("MCO ref22 model test.csv")
    print(names(mcodat))

    dev.new()
    par(mar = c(4,4,1,1), mgp = c(2.3,1,0), bty = "l")
    plot( log(envdata1819$tss.result), log(envdata1819$ptl),
         pch = 21, col = "grey", bg = "white", axes = F,
         xlab = "TSS (mg/L)", ylab = expression(Total~P~(mu*g/L)))
    points( log(allchem$tss),log(allchem$tp*1000), col = "red")
    logtick.exp(0.001, 10, c(1,2), c(F,F))
    stop()
          
#    png(width = 6, height = 3, pointsize = 10, units = "in",
#        res = 600, file = "pred.obs.png")
    dev.new()
    par(mar = c(4,4,1,1), mfrow= c(1,2), mgp = c(2.3,1,0))
    plot(mcodat$tss.result_test, log(mcodat$Avg.TSS), axes = F,
         pch = 21, col = "grey39", bg = "white",
         xlab ="TSS (predicted)",ylab = "TSS (observed)")
    logtick.exp(0.001, 10, c(1,2), c(T,T))
    abline(0,1)
    plot(mcodat$cond_test, log(mcodat$Avg.SC), axes = F,
         pch = 21, col = "grey39", bg = "white",
         xlab = "Cond (predicted)", ylab = "Cond (observed)")
    logtick.exp(0.001, 10, c(1,2), c(T,T))
    abline(0,1)
    stop()
    dev.off()
    
}

plotresults()
