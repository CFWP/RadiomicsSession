###############################################################################
###############################################################################
## 
## Source file containing some convenience functions to FMradio 1.1
## Statistics for Omics course: Radiomics
##
## Code:   Carel F.W. Peeters
##         Department of Epidemiology & Biostatistics
##         VU University medical center Amsterdam
##         Amsterdam, the Netherlands
##         cf.peeters@vumc.nl
## Date:   13/06/2019, Amsterdam, VUmc
##
###############################################################################
###############################################################################


#'
#'**---------------------------------**\
#'**Convenience functions*\
#'**---------------------------------**\

Or2 <- function(this){
  ##############################################################################
  # Convenience function for calculating (overall) explained residual variation 
  # this > crps object from pec
  #
  # Notes:
  # - Function tailored towards the application/comparison at hand.
  #   Can be further generalized.
  ##############################################################################
  
  # Dependencies:
  # require("base")
  # require("pec")
  
  Ar2 <- 1 - this[,1]/this[1,1]
  Cr2 <- 1 - this[,2]/this[1,2]
  Tab <- cbind(Ar2,Cr2)
  colnames(Tab)    <- c("AppR2", "crossvalR2")
  rownames(Tab)[2] <- "FMradio"
  return(Tab)
}


plotR2Table <- function(R2Table, type = "CV", Ylab = NULL, Xlab = "Time"){
  ##############################################################################
  # Convenience function for plotting explained residual variation 
  # R2Table > R2 object from pec
  # type    > character indicating if the cross-validated or apparent plot
  #           must be returned. Must be one of "CV" or "AE"
  # Ylab    > control over the y-axis label
  #
  # Notes:
  # - Function tailored towards the application/comparison at hand.
  #   Can be further generalized.
  ##############################################################################
  
  # Dependencies:
  # require("base")
  # require("pec")
  
  if (type == "CV"){
    if (!is.null(Ylab)){Ylab = Ylab}
    if (is.null(Ylab)){Ylab = expression(Averaged ~ cross-validated ~ R^2)}
    MIN <- min(R2Table$crossvalErr[,c(-1)])
    MAX <- max(R2Table$crossvalErr[,c(-1)])
    COL <- 2:3
    plot(R2Table$crossvalErr$time, 
         R2Table$crossvalErr$RR.MetaCox, 
         type = "l",
         axes = FALSE,
         xlab = Xlab,
         ylab = Ylab,
         col = COL[1],
         lwd = 1.5,
         ylim = c(MIN,MAX)#,
         #main = "Explained residual variation under cross-validated error")
    )
    axis(2, col = "black", lwd = 1)
    axis(1, col = "black", lwd = 1)
    abline(h = 0)
    lines(R2Table$crossvalErr$time,
          R2Table$crossvalErr$RR.RforestCox,
          type = "l",
          col = COL[2],
          lwd = 1.5)
    legend("bottomright", inset = .03,
           legend = c("FMradio",
                      "Random survival forest"),
           col = COL, 
           lty = 1,
           cex = .9,
           box.lty = 0,
           lwd = 2,
           y.intersp = 1.7)
  }
  if (type == "AE"){
    if (!is.null(Ylab)){Ylab = Ylab}
    if (is.null(Ylab)){Ylab = expression(Apparent ~ R^2)}
    MIN <- min(R2Table$AppErr[,c(-1)])
    MAX <- max(R2Table$AppErr[,c(-1)])
    COL <- 2:5
    plot(R2Table$AppErr$time, 
         R2Table$AppErr$RR.MetaCox, 
         type = "l",
         axes = FALSE,
         xlab = Xlab,
         ylab = Ylab,
         col = COL[1],
         lwd = 1.5,
         ylim = c(MIN,MAX)#,
         #main = "Explained residual variation under apparent error")
    )
    axis(2, col = "black", lwd = 1)
    axis(1, col = "black", lwd = 1)
    abline(h = 0)
    lines(R2Table$AppErr$time,
          R2Table$AppErr$RR.RforestCox,
          type = "l",
          col = COL[2],
          lwd = 1.5)
    legend("bottomright", inset = .03,
           legend = c("FMradio",
                      "Random survival forest"),
           col = COL, 
           lty = 1,
           cex = .9,
           box.lty = 0,
           lwd = 2,
           y.intersp = 1.7)
  }
}
