###--------------------------------------------------------------------------###
###-- Sparse TVP VECMs with an application to modeling electricity prices ---###
###------------------ Hauzenberger, Pfarrhofer, & Rossini -------------------###
###----------------- International Journal of Forecasting -------------------###
###--------------------------------------------------------------------------###
###-------------------- Reproducing out-of-sample MSC -----------------------###
###--------------------- Tabs. 2, B.2, B.3, and B.4 -------------------------###
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
###------------------------- Set up MCS tables ------------------------------###
###--------------------------------------------------------------------------###

# Select table that should be produced
tab.slct <- "Tab-2"
#tab.slct <- "Tab-B2"
#tab.slct <- "Tab-B3"
#tab.slct <- "Tab-B4"

# Full model set
model.set <- c("VECM-TVP-iSV-t_2_sps", "VECM-TVP-iSV-n_2_sps", "VECM-TVP-iHOMO_2_sps", 
               "VECM-TIV-iSV-t_2_sps", "VECM-TIV-iSV-n_2_sps", "VECM-TIV-iHOMO_2_sps",
               "VECM-TVP-iSV-t_2"    ,     "VECM-TVP-iSV-n_2", "VECM-TVP-iHOMO_2", 
               "VECM-TIV-iSV-t_2"    ,     "VECM-TIV-iSV-n_2", "VECM-TIV-iHOMO_2",
               "VARl-TVP-iSV-n_2"    ,    
               "VARl-TIV-iSV-n_2"    ,    
               "ARX-7_all","ARX-2_all", "ARX-7_RES", "ARX-2_RES", 
               "ARX-7_ARp", "ARX-2_ARp")

if(tab.slct == "Tab-2"){
  # Tab. 2: Model confidence set (MCS) at a 10% 
  # significance level using the TR test statistic
  mcs.alph <- "010" # 10% vs. 25% signficance
  mcs.stat <- "TR"  # TR vs. Tmax 
  
}else if(tab.slct == "Tab-B2"){
  # Tab. B.2: Model confidence set (MCS) at a 10% 
  # significance level using the Tmax test statistic
  mcs.alph <- "010" # 10% vs. 25% signficance
  mcs.stat <- "Tmax" # TR vs. Tmax 

}else if(tab.slct == "Tab-B3"){
  # Tab. B.3: Model confidence set (MCS) at a 25% 
  # significance level using the TR test statistic
  mcs.alph <- "025" # 10% vs. 25% signficance
  mcs.stat <- "TR" # TR vs. Tmax 
  
}else if(tab.slct == "Tab-B4"){
  # Tab. B.4: Model confidence set (MCS) at a 25% 
  # significance level using the Tmax test statistic
  mcs.alph <- "025" # 10% vs. 25% signficance
  mcs.stat <- "Tmax" # TR vs. Tmax 
  
}
  
mcs.rank <- "Rank_R"  # M vs. R
load(paste0(wdir, "MCS_", mcs.alph, "_", mcs.stat, ".rda"))

MCS.tab <- as.data.frame(matrix("", length(model.set), 13), stringsAsFactors = FALSE)
rownames(MCS.tab) <- model.set
colnames(MCS.tab) <- c("Total", paste0(8:11, "a.m."), "12 noon", paste0(1:6, "p.m."), "Night")
MCS.tab[rownames(MCS.eval.avg),1] <- as.character(MCS.eval.avg[,mcs.rank])
for(ii in 2:ncol(MCS.tab)) MCS.tab[rownames(MCS.eval[[ii-1]]),ii] <- as.character(MCS.eval[[ii-1]][,mcs.rank])

MCS.tab[] <- apply(MCS.tab, 2, function(x) gsub(gsub(gsub(x, pattern = "^1$", replacement = "\\\\textbf{1}"), 
                                                     pattern = "^2$", replacement = "\\\\textbf{2}"),
                                                pattern = "^3$", replacement = "\\\\textbf{3}"))

MCS.tab$Class <- c(rep("VECM (sparse)", 6), rep("VECM", 6), "TVP-VAR", "TIV-VAR", 
                   "AR(7)-X(all)", "AR(2)-X(all)", "AR(7)-X(RES)", "AR(2)-X(RES)", "AR(7)", "AR(2)")
MCS.tab$TVP   <- c(rep(c("TVP", "TVP", "TVP", "TIV", "TIV", "TIV"), 2), "TVP", rep("TIV", 7))
MCS.tab$SV    <- c("t", "n", "hom", "t", "n", "hom", "t", "n", "hom", "t", "n", "hom", rep("n", 8))


MCS.tab <- MCS.tab[,c(14:16,1:13)]

rgroup <- c("", "", "", "", "")
cgroup <-  c("Specification", "1-day-ahead")
colheads <- colnames(MCS.tab)
n.rgroup <- c(6,6,2,4,2) # Row groups
n.cgroup <- c(3,13)      # Column groups

latex(MCS.tab, file = paste0(pdir, tab.slct, "_MCS_", mcs.stat, "_", mcs.alph, ".tex"), title = "", ctable = FALSE,  rgroupTexCmd = "scshape", n.cgroup = n.cgroup, cgroup = cgroup, colheads = colheads, rgroup = rgroup, n.rgroup = n.rgroup, 
      caption = paste0("Model confidence set (MCS) for density forecasts: ", as.numeric(mcs.alph), " percent significance using ", mcs.stat, " test statistic"), label = "", size = "tiny", booktabs = TRUE, rowname = "", numeric.dollar = FALSE, caption.loc = "bottom", landscape = TRUE)

