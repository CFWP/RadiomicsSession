###############################################################################
###############################################################################
## 
## Radiomics prediction Exercises
## On the basis of FMradio 1.1
## Course: Statistics for Omics: Radiomics
## Place:  Amsterdam, 28/06/2019
##
## Code:   Carel F.W. Peeters
##         Department of Epidemiology & Biostatistics
##         VU University medical center Amsterdam
##         Amsterdam, the Netherlands
##         cf.peeters@vumc.nl
## Date:   22/06/2019, Amsterdam, VUmc
##
###############################################################################
###############################################################################



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Preliminaries**
#' **------------------------------------------------------------------------**

## Set working directory
setwd("C:/Users/cf.peeters/Desktop/Teaching/Courses_Advanced/2019_Statistics_for_Omics/2. Radiomics/02. Exercises")

## Needed library
library(FMradio)

## Other requirements
require("Biobase")
require("rags2ridges")
require("DandEFA")
require("randomForestSRC")
require("pec")
require("survival")

## Convenience functions
source("Convenience.R")




#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 1: Get acquainted with the data**
#' **------------------------------------------------------------------------**

## Data packaged as expressionset
## Will invoke basic functions for looking at data

## Load data, get to know objects
load("AnonymousRadio.Rdata")

## Look at radiomic measurements
exprs(AnonymousRadio)[10:30,1:5]

## Look at Clinical (sample) information
head(pData(AnonymousRadio))




###############################################################################
###############################################################################
## 
## Section 1: Factor Analytic Projection
##
###############################################################################
###############################################################################

#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 2: Assess redundancy in the correlation matrix**
#' **------------------------------------------------------------------------**

## Scale data and get raw correlation matrix
DATAscaled <- scale(t(exprs(AnonymousRadio)))
R          <- cor(DATAscaled)

## Redundancy visualization
radioHeat(R, diag = FALSE, 
          threshold = TRUE, 
          threshvalue = .95,
          labelsize = .01)



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 3: Redundancy-filter the correlation matrix**
#' **------------------------------------------------------------------------**

## Redundancy filtering
## And subsetting data
## 124 features remain
RFs         <- RF(R, t = .95)
DATAscaledS <- subSet(DATAscaled, RFs)



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 4: Find an optimal regularized correlation matrix**
#' **------------------------------------------------------------------------**

## Optimal penalty
set.seed(303)
OPT <- regcor(DATAscaledS, fold = 5)

## Look at optimal penalty-value
## Obtain regularized correlation matrix
## Conditioning can again be assessed with, e.g., CNplot from rag2ridges
OPT$optPen
Re = OPT$optCor



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 5: Perform factor analysis on the regularized correlation matrix**
#' **------------------------------------------------------------------------**

## Assess dimensionality factor solution
## 13 considered upper bound
## Variance explained would suggest 8 factors
dimGB(Re)
dimVAR(Re, 15, graph = TRUE)

## Assessing solutions around 8
## 9th factor seems weak
## Will keep solution at 8
## ML factor analysis with Varimax rotation
fito <- mlFA(R = Re, m = 8)
print(fito$Loadings, digits = 2, cutoff = .3, sort = TRUE)

## Visualizing solution using Dandelion plot
dandpal <- rev(rainbow(100, start = 0.4, end = 0.6))
dandelion(fito$Loadings, bound = .3, mcex = c(1,1), palet = dandpal)

## Export to pdf for inspection
pdf(file = "Dandelion.pdf", width = 11, height = 11)
dandelion(fito$Loadings, bound = .3, mcex = c(1,1), palet = dandpal)
dev.off()



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 6: Obtain factor scores**
#' **------------------------------------------------------------------------**

## Factor scores
Lambda <- fito$Loadings
Psi    <- fito$Uniqueness
Scores <- facScore(DATAscaledS, Lambda, Psi)



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 7: Assess factor scores**
#' **------------------------------------------------------------------------**

## Determinacy factor scores
## Highly determinate
DF <- facSMC(Re, Lambda); DF




###############################################################################
###############################################################################
## 
## Section 2: Prediction
##
###############################################################################
###############################################################################

#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 8: Concatenate original and projected data**
#' **------------------------------------------------------------------------**

## Combine original (scaled) data with projected meta-features
DAT <- cbind(DATAscaled, Scores)

## Include the survival information
Status  <- as.numeric(AnonymousRadio$Death_yesno) - 1
time    <- AnonymousRadio$Death_followuptime_months
DAT     <- cbind(time, Status, DAT)



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 9: Set up model comparisons**
#' **------------------------------------------------------------------------**

## Formulating the model formula's
FitRSF     <- as.formula(paste("Surv(time, Status)~", 
                               paste(colnames(DAT)[3:434], collapse="+")))
FitMetaCox <- as.formula(paste("Surv(time, Status) ~", 
                               paste(colnames(DAT[,c(435:442)]), collapse="+")))

models <- list("MetaCox" = coxph(FitMetaCox, data = DAT, x = TRUE, y = TRUE),
               "RforestCox" = rfsrc(FitRSF, data = DAT))



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 10: Compare models w.r.t. prediction error**
#' **------------------------------------------------------------------------**

## Assessing prediction error
## Median follow-up time = 25.7
## (Averaged) repeated 5-fold cross-validation
set.seed(446464)
PredError <- pec(object = models,
                 formula = Surv(time, Status) ~ 1,
                 data = DAT,
                 exact = TRUE,
                 maxtime = median(time),
                 cens.model = "marginal",
                 splitMethod = "cv5",
                 B = 50,
                 verbose = TRUE)

## Summary results:
## Apparent and cross-validated prediction error and R2
crps(PredError)
Or2(crps(PredError))



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 11: Visualize the results**
#' **------------------------------------------------------------------------**

## Visualize apparent prediction error
plot(PredError, what = "AppErr",
     xlab = "Time (months)",
     ylab = "Apparent prediction error",
     legend.cex = .9,
     legend.lty = 1,
     legend.lwd = 2,
     legend.legend = c("Reference model", 
                       "FMradio",
                       "Random survival forest"),
     add.refline = TRUE,
     lwd = 1.5,
     legend.y.intersp = 1.7)

## Visualize cross-validated prediction error
plot(PredError, what = "crossvalErr",
     xlab = "Time (months)",
     ylab = "Averaged cross-validated prediction error",
     legend.cex = .9,
     legend.lty = 1,
     legend.lwd = 2,
     legend.legend = c("Reference model", 
                       "FMradio",
                       "Random survival forest"),
     add.refline = TRUE,
     lwd = 1.5,
     legend.y.intersp = 1.7)

## Visualize apparent residual explained variation 
R2Table <- R2(PredError, times = seq(0,median(time),.01), reference = 1)
plotR2Table(R2Table, "AE", Xlab = "Time (months)")

## Visualize cross-validated residual explained variation 
plotR2Table(R2Table, "CV", Xlab = "Time (months)")


## Exporting
pdf("Internal.pdf", width = 15, height = 15)
par(mfrow=c(2,2))
## Visualize apparent prediction error
plot(PredError, what = "AppErr",
     xlab = "Time (months)",
     ylab = "Apparent prediction error",
     legend.cex = .9,
     legend.lty = 1,
     legend.lwd = 2,
     legend.legend = c("Reference model", 
                       "FMradio",
                       "Random survival forest"),
     add.refline = TRUE,
     lwd = 1.5,
     legend.y.intersp = 1.7)

## Visualize cross-validated prediction error
plot(PredError, what = "crossvalErr",
     xlab = "Time (months)",
     ylab = "Averaged cross-validated prediction error",
     legend.cex = .9,
     legend.lty = 1,
     legend.lwd = 2,
     legend.legend = c("Reference model", 
                       "FMradio",
                       "Random survival forest"),
     add.refline = TRUE,
     lwd = 1.5,
     legend.y.intersp = 1.7)

## Visualize apparent residual explained variation 
R2Table <- R2(PredError, times = seq(0,median(time),.01), reference = 1)
plotR2Table(R2Table, "AE", Xlab = "Time (months)")

## Visualize cross-validated residual explained variation 
plotR2Table(R2Table, "CV", Xlab = "Time (months)")
dev.off()




###############################################################################
###############################################################################
## 
## Section 3: Hidden Gems
##
###############################################################################
###############################################################################

## Really?
FMradio:::.Airwolf()

## You betcha!
FMradio:::.Airwolf2()


