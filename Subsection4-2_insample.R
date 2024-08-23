###--------------------------------------------------------------------------###
###-- Sparse TVP VECMs with an application to modeling electricity prices ---###
###------------------ Hauzenberger, Pfarrhofer, & Rossini -------------------###
###----------------- International Journal of Forecasting -------------------###
###--------------------------------------------------------------------------###
###-------------------- Reproducing in-sample results -----------------------###
###-------------------- Figs. 1, 2, and 3 + Tab. B.1 ------------------------###
###--------------------------------------------------------------------------###
rm(list=ls())

###--------------------------------------------------------------------------###
###------------------- Load packages and set directories --------------------###
###--------------------------------------------------------------------------###
# Packages
library(xtable)
require(tidyr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(scales)

# Directories
wdir <- "INSAMPLE-data/"   # Working directory  
pdir <- "INSAMPLE-res/"   # Plotting directory
dir.create(pdir)

###--------------------------------------------------------------------------###
###------------------- Reproduce Figure 1 and Figure 2 ----------------------###
###--------------------------------------------------------------------------###
load(paste0(wdir, "Figures1-2_PIPs-TVP-VECM-SV-t.rda"))

# Fig. 1: Posterior probability of the rank (PPR) over time
rank_tab <- melt(rank_tab)
colnames(rank_tab) <- c("Periods", "Rank", "Prob")
l_col <- "white"
h_col <- "red3"
rank_tab$Periods <- as.numeric(rank_tab$Periods)

rank_plot <- ggplot(rank_tab, aes(x=Periods, y= Rank, fill=Prob)) +
  geom_tile() + geom_raster()  +
  labs(x="",y="") +
  scale_x_continuous(breaks = seq(2017, 2020, 1), label = c("2017", "2018", "2019","2020")) + 
  scale_fill_gradientn(colours=c(l_col, alpha(h_col, seq(0.3,1,0.1))),
                       values=rescale(c(-3, -2, -1, 0, 1, 2, 3)),name= "PPR", na.value = "grey80",limits = c(0,1),
                       guide="colorbar") +  
  theme_bw() + theme( panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y=element_blank(), strip.background = element_rect(colour="white", fill="grey99"), legend.text=element_text(size=14), strip.text.x = element_text(size = 14), legend.title=element_text(size=14),axis.text=element_text(size=14))


pdf(paste0(pdir, "Figure1_rank.pdf"), width = 15, height = 4)
print(rank_plot)
dev.off()

# Fig. 2 (a) : Posterior inclusion probability (PIP) over time of coefficients 
A_tab <- melt(coeff_pip, id = "Periods")
colnames(A_tab) <- c("Periods", "Coeff", "Prob")
vars <- t(as.data.frame(strsplit(as.character(A_tab$Coeff), ":")))
vars[,2] <- paste0(gsub(vars[,2], pattern = "_", replacement = "("), ")")

A_tab$Endog <- vars[,1]
A_tab$Coeff <- vars[,2]

l_col<- "white"
h_col <- "red3"
A_tab$Periods <- as.numeric(A_tab$Periods)
A_tab$Coeff <- as.factor(A_tab$Coeff)

A_pip_plot <-  ggplot(A_tab, aes(x= Periods, y= Coeff, fill=Prob)) +
  geom_tile() + geom_raster()  +
  labs(x="",y="") +
  scale_y_discrete(limits = rev(levels(A_tab$Coeff))) +
  scale_x_continuous(breaks = seq(2017, 2020, 1), label = c("2017", "2018", "2019", "2020")) + 
  scale_fill_gradientn(colours=c(l_col, alpha(h_col, seq(0.3,1,0.1))),
                       values=rescale(c(-3, -2, -1, 0, 1, 2,3)),name= "PIP", na.value = "grey80",limits = c(0,1),
                       guide="colorbar") + 
  facet_wrap(~Endog, ncol = 5) +
  theme_bw() + theme( panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y=element_blank(), strip.background = element_rect(colour="white", fill="grey99"), legend.text=element_text(size=14), strip.text.x = element_text(size = 14), legend.title=element_text(size=14),axis.text=element_text(size=14))

pdf(paste0(pdir, "Figure2-a_pip-coeff.pdf"), width = 22, height = 10)
print(A_pip_plot)
dev.off()


# Fig. 2 (b) : Posterior inclusion probability (PIP) over time of covariances 
sig_tab <- melt(cov_pip, id = "Periods")
colnames(sig_tab) <- c("Periods", "Cov.", "Prob")

vars <- t(as.data.frame(strsplit(as.character(sig_tab$Cov.), "-")))
sig_tab$Endog <- vars[,1]
sig_tab$Cov. <- vars[,2]

l_col <- "white"
h_col <- "red3"
sig_tab$Periods <- as.numeric(sig_tab$Periods)
sig_tab$Cov. <- as.factor(sig_tab$Cov.)

cov_pip_plot <-  ggplot(sig_tab, aes(x= Periods, y= Cov., fill=Prob)) +
  geom_tile() + geom_raster()  +
  labs(x="",y="") +
  scale_y_discrete(limits = rev(levels(sig_tab$Cov.))) +
  scale_x_continuous(breaks = seq(2017, 2020, 1), label = c("2017", "2018", "2019", "2020")) + 
  scale_fill_gradientn(colours=c(l_col, alpha(h_col, seq(0.3,1,0.1))),
                       values=rescale(c(-3, -2, -1, 0, 1, 2,3)),name= "PIP", na.value = "grey",limits = c(0,1),
                       guide="colorbar") +  
  facet_wrap(~Endog, ncol = 4) +
  theme_bw() + theme( panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y=element_blank(), strip.background = element_rect(colour="white", fill="grey99"), legend.text=element_text(size=14), strip.text.x = element_text(size = 14), legend.title=element_text(size=14),axis.text=element_text(size=14))


pdf(paste0(pdir, "Figure2-b_pip-cov.pdf"), width = 18, height = 8)
print(cov_pip_plot)
dev.off()


###--------------------------------------------------------------------------###
###------------- Reproduce panels of Figure 3 and Table B.1 -----------------###
###--------------------------------------------------------------------------###
cols <- c(grey(c(0, 0.5)), "deepskyblue2", "royalblue", "royalblue4")

# Fig. 3 (a) Posterior median of SV with t-distributed errors

# Load specification: t-distributed errors ("SV-t") 
load(paste0(wdir, "Figure3-TableB1_TVP-VECM-SV-t.rda"))

# First set of countries: 
pdf(paste0(pdir, "Figure3-a-1_SV-t.pdf"), width = 9, height = 6)
matplot(x = seq(2017+3/365, 2019+364/365, 1/365), 
        y = cbind(logVolas[,c("BALT", "DK", "FI", "NO", "SE")], logVolas_pca), ylim = c(-6, 2),
        type = "l", col = c(cols, "darkred"), lty = c(1,1,1,1,1,1), lwd = c(rep(1.5,5),2), 
        main ="", xlab = "", ylab ="", xaxt = "n")
axis(side = 1, at= c(2017, 2018, 2019, 2020), labels= c("2017", "2018", "2019", "2020"))
legend("bottom", legend = c("BALT", "DK", "FI", "NO", "SE", "PCA") , pch = NA, col = c(cols, "darkred"), lty = 1, lwd = c(rep(2,9),2), bty = "n", ncol = 6) 
dev.off()

# Second set of countries:
pdf(paste0(pdir, "Figure3_a-2_SV-t.pdf"), width = 9, height = 6)
matplot(x = seq(2017+3/365, 2019+364/365, 1/365), 
        y = cbind(logVolas[,c("CH", "DE",  "FR", "IT")], logVolas_pca), ylim = c(-6, 2),
        type = "l", col = c(cols[c(1,2,4,5)], "darkred"), lty = c(1,1,1,1,1), lwd = c(rep(1.5,4),2), 
        main ="", xlab = "", ylab ="", xaxt = "n")
axis(side = 1, at= c(2017, 2018, 2019, 2020), labels= c("2017", "2018", "2019", "2020")) 
legend("bottom", legend = c("CH", "DE",  "FR", "IT", "PCA"), pch = NA, col = c(cols[c(1,2,4,5)], "darkred"), lty = 1, lwd = c(rep(2,9),2), bty = "n", ncol = 6) 
dev.off()

# Upper panel of Tab. B.1: Posterior estimates for the state equations of the error variances
SVpara <- cbind("Type" = rownames(SVpara), SVpara)
SVpara[c(1,3,5,7), "Type"] <- paste0("$\\", SVpara[c(1,3,5,7), "Type"], "$")


print(xtable(SVpara, caption = "Posterior estimates for SV with t-distributed errros"), 
      include.colnames = TRUE, 
      include.rownames = FALSE,
      sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      file = paste0(pdir, "Tab-B1-a_SV-t.tex"))

# Fig. 3 (b) Posterior median of SV with Gaussian-distributed errors

# Load specification: Gaussian-distributed errors ("SV-n") 
load(paste0(wdir, "Figure3-TableB1_TVP-VECM-SV-n.rda"))

# First set of countries: 
pdf(paste0(pdir, "Figure3-a-1_SV-n.pdf"), width = 9, height = 6)
matplot(x = seq(2017+3/365, 2019+364/365, 1/365), 
        y = cbind(logVolas[,c("BALT", "DK", "FI", "NO", "SE")], logVolas_pca), ylim = c(-6, 2),
        type = "l", col = c(cols, "darkred"), lty = c(1,1,1,1,1,1), lwd = c(rep(1.5,5),2), 
        main ="", xlab = "", ylab ="", xaxt = "n")
axis(side = 1, at= c(2017, 2018, 2019, 2020), labels= c("2017", "2018", "2019", "2020"))
legend("bottom", legend = c("BALT", "DK", "FI", "NO", "SE", "PCA") , pch = NA, col = c(cols, "darkred"), lty = 1, lwd = c(rep(2,9),2), bty = "n", ncol = 6) 
dev.off()

# Second set of countries:
pdf(paste0(pdir, "Figure3_a-2_SV-n.pdf"), width = 9, height = 6)
matplot(x = seq(2017+3/365, 2019+364/365, 1/365), 
        y = cbind(logVolas[,c("CH", "DE",  "FR", "IT")], logVolas_pca), ylim = c(-6, 2),
        type = "l", col = c(cols[c(1,2,4,5)], "darkred"), lty = c(1,1,1,1,1), lwd = c(rep(1.5,4),2), 
        main ="", xlab = "", ylab ="", xaxt = "n")
axis(side = 1, at= c(2017, 2018, 2019, 2020), labels= c("2017", "2018", "2019", "2020")) 
legend("bottom", legend = c("CH", "DE",  "FR", "IT", "PCA"), pch = NA, col = c(cols[c(1,2,4,5)], "darkred"), lty = 1, lwd = c(rep(2,9),2), bty = "n", ncol = 6) 
dev.off()

# Lower panel of Tab. B.1: Posterior estimates for the state equations of the error variances
SVpara <- cbind("Type" = rownames(SVpara), SVpara)
SVpara[c(1,3,5), "Type"] <- paste0("$\\", SVpara[c(1,3,5), "Type"], "$")

print(xtable(SVpara, caption = "Posterior estimates for SV with Gaussian errros"), 
      include.colnames = TRUE, 
      include.rownames = FALSE,
      sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      file = paste0(pdir, "Tab-B1-b_SV-n.tex")
)
