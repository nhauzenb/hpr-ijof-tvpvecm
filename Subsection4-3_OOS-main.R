###--------------------------------------------------------------------------###
###-- Sparse TVP VECMs with an application to modeling electricity prices ---###
###------------------ Hauzenberger, Pfarrhofer, & Rossini -------------------###
###----------------- International Journal of Forecasting -------------------###
###--------------------------------------------------------------------------###
###------------------ Reproducing out-of-sample results ---------------------###
###-------------------------------- Tab. 1 ----------------------------------###
###--------------------------------------------------------------------------###
rm(list=ls())

###--------------------------------------------------------------------------###
###------------------- Load packages and set directories --------------------###
###--------------------------------------------------------------------------###
# Packages
library(Hmisc)

# Directories
wdir <- "OOS-data/"   # Working directory  
pdir <- "OOS-res/"   # Plotting directory
dir.create(pdir)

###--------------------------------------------------------------------------###
###------- Compute average forecast performance based on raw scores ---------###
###--------------------------------------------------------------------------###
# Load raw scores 
load(file = paste0(wdir, "Tab1_OOS-RawScores.rda"))

TT <- nrow(MSE.avg.scores) # 540 hold out observations
benchmark <- "ARX-2_ARp"   # Benchmark

# Marginal scores
RMSE.bench <- sqrt(apply(MSE.scores, c(2,3), cumsum)[,benchmark,]/(1:TT)) # RMSE
CRPS.bench <- apply(CRPS.scores, c(2,3), cumsum)[,benchmark,]/(1:TT)      # CRPS
    
RMSE.scores <- sqrt(apply(MSE.scores, c(2,3), cumsum)/(1:TT)) # RMSE
CRPS.scores <- apply(CRPS.scores, c(2,3), cumsum)/(1:TT)      # CRPS
    
# Total scores
RMSE.bench.avg <- sqrt(apply(MSE.avg.scores, 2, cumsum)[,benchmark]/(1:TT)) # RMSE
CRPS.bench.avg <- apply(CRPS.avg.scores, 2, cumsum)[,benchmark]/(1:TT)      # CRPS
    
RMSE.avg.scores <- sqrt(apply(MSE.avg.scores, 2, cumsum)/(1:TT)) # RMSE
CRPS.avg.scores <- apply(CRPS.avg.scores, 2, cumsum)/(1:TT)      # CRPS
    
# Marginal scores
for(kk in 1:dim(RMSE.scores)[[3]]){
RMSE.scores[,,kk] <- RMSE.scores[,,kk]/RMSE.scores[,benchmark,kk] # RMSE
CRPS.scores[,,kk] <- CRPS.scores[,,kk]/CRPS.scores[,benchmark,kk] # CRPS
}
  
# Total scores
RMSE.avg.scores <- RMSE.avg.scores/RMSE.avg.scores[,benchmark]    # RMSE
CRPS.avg.scores <- CRPS.avg.scores/CRPS.avg.scores[,benchmark]    # CRPS
  
# Focus on average at the final OOS period
RMSE.tab <- RMSE.scores[TT,,]
CRPS.tab <- CRPS.scores[TT,,]
RMSE.tab.avg <- RMSE.avg.scores[TT,]
CRPS.tab.avg <- CRPS.avg.scores[TT,]

save(file = paste0(wdir, "Tab1_OOS-AvgScores.rda"), list = c("RMSE.tab", "CRPS.tab", "RMSE.tab.avg", "CRPS.tab.avg"))

###--------------------------------------------------------------------------###
###------------------------ Prepare final table -----------------------------###
###--------------------------------------------------------------------------###
# Column-bind total scores and marginals 
RMSE.final.tab <- as.data.frame(cbind(RMSE.tab.avg, RMSE.tab))
CRPS.final.tab <- as.data.frame(cbind(CRPS.tab.avg, CRPS.tab))
colnames(RMSE.final.tab) <- c("Tot.", colnames(RMSE.tab))
colnames(CRPS.final.tab) <- c("Tot.", colnames(CRPS.tab))

RMSE.final.tab$Class <- CRPS.final.tab$Class <-  c(rep("VECM (sparse)", 6), rep("VECM", 6), "TVP-VAR", "TIV-VAR", 
                                                   "AR(7)-X(all)", "AR(2)-X(all)", "AR(7)-X(RES)", "AR(2)-X(RES)", "AR(7)", "AR(2)")
RMSE.final.tab$TVP <- CRPS.final.tab$TVP     <- c(rep(c("TVP", "TVP", "TVP", "TIV", "TIV", "TIV"), 2), "TVP", rep("TIV", 7))
RMSE.final.tab$SV <- CRPS.final.tab$SV       <- c("t", "n", "hom", "t", "n", "hom", "t", "n", "hom", "t", "n", "hom", rep("n", 8))

n.eval <- 13

for(ii in 1:n.eval)
{
  ind.1 <- which(as.numeric(RMSE.final.tab[,ii]) == min(as.numeric(RMSE.final.tab[,ii]), na.rm = TRUE)) #MSE
  ind.2 <- which(as.numeric(CRPS.final.tab[,ii]) == min(as.numeric(CRPS.final.tab[,ii]), na.rm = TRUE)) #MAE

  RMSE.final.tab[,ii]  <- format(round(RMSE.final.tab[,ii], 3), nsmall = 3)
  CRPS.final.tab[,ii]  <- format(round(CRPS.final.tab[,ii] , 3), nsmall = 3)
  
  RMSE.final.tab[ind.1,ii] <- paste0("\\textbf{", gsub(RMSE.final.tab[ind.1,ii], pattern = " ", replacement = ""), "}")
  CRPS.final.tab[ind.2,ii] <- paste0("\\textbf{", gsub(CRPS.final.tab[ind.2,ii], pattern = " ", replacement = ""), "}")
}

# Load absolute, re-scaled benchmark scores 
# For estimation, we standardized the data and then rescale the forecasts 
# and the scores to the original scale
load(file = paste0(wdir, "Tab1_OOS-AR2resc.rda"))
RMSE.final.tab[benchmark,1:13] <- format(round(c(RMSEraw.bench.avg, RMSEraw.bench), 2), nsmall = 2)
CRPS.final.tab[benchmark,1:13] <- format(round(c(CRPSraw.bench.avg, CRPSraw.bench), 2), nsmall = 2)

full.tab <- data.frame(matrix("", nrow(RMSE.final.tab)*2, ncol(RMSE.final.tab)))
rownames(full.tab)[(0:(nrow(RMSE.final.tab)-1))*2+1] <- rownames(RMSE.final.tab)
colnames(full.tab) <- colnames(RMSE.final.tab)

full.tab[(0:(nrow(RMSE.final.tab)-1))*2+1,] <- RMSE.final.tab
full.tab[(1:(nrow(RMSE.final.tab)))*2,1:n.eval] <- apply(CRPS.final.tab[,1:n.eval], c(2), function(x){paste0("(", gsub(x, pattern = " ", replacement = ""), ")")})

full.tab <- full.tab[,c(14:16,1:n.eval)]

rgroup <- c("", "", "", "", "")
cgroup <-  c("Specification", "1-day-ahead")
colheads <- colnames(full.tab)
n.rgroup <- c(12,12,4,8,4) # Row groups
n.cgroup <- c(3,13)      # Column groups

latex(full.tab, file = paste0(pdir, "Tab1_OOS-overview.tex"), title = "", ctable = FALSE,  rgroupTexCmd = "scshape", n.cgroup = n.cgroup, cgroup = cgroup, colheads = colheads, rgroup = rgroup, n.rgroup = n.rgroup,  
      caption =  "Forecast performance for point and density forecasts", label = "", size = "tiny", booktabs = TRUE, rowname = "",  numeric.dollar = FALSE, caption.loc = "bottom", landscape = TRUE)
