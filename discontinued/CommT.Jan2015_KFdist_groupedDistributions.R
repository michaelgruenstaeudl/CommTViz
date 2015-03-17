#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2014-2015 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl")
#email = "gruenstaeudl.1@osu.edu"
#version = "2015.02.27.1300"


####################
## Load libraries ##
####################
library(ape)
library(ggplot2)
library(gridExtra)
library(phangorn)
library(reshape)


######################
## Global variables ##
######################
setwd("~/Desktop/")
n_genes = 10
n_sims = sprintf("%02d", 1:20)


####################
## Start sim-loop ##
####################
for (s in n_sims) {
cat(paste("\n", "Analyzing simulation", s, "\n"))


####################################################
## Loading data sets of treatment group "swapped" ##
####################################################
BEAST_prefix = paste("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/01_input/02_BEAST_GeneTrees/CommT.Jan2015.fitonly.sim.", s, ".gene", sep="")
starBEAST_prefix = paste("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/01_input/01_starBEAST_geneTrees/CommT.Jan2015.g1swap.sim.", s, ".gene", sep="")


####################################################
## Loading data sets of treatment group "swapped" ##
####################################################
#BEAST_prefix = paste("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/01_input/02_BEAST_GeneTrees/CommT.Jan2015.noswap.sim.", s, ".gene", sep="")
#starBEAST_prefix = paste("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/01_input/01_starBEAST_geneTrees/CommT.Jan2015.noswap.sim.", s, ".gene", sep="")


######################
## Initialize lists ##
######################
BEAST_postGTdistr = list()
starBEAST_postGTdistr = list()


###################################################
## Loading the posterior gene tree distributions ##
###################################################
for (g in 1:n_genes) {
g_lz = sprintf("%03d", g)
BEAST_postGTdistr[[g]] = read.nexus(paste(BEAST_prefix, g_lz, ".trees", sep=""))
starBEAST_postGTdistr[[g]] = read.nexus(paste(starBEAST_prefix, g_lz, ".trees", sep=""))
}


##########################
## Initialize outmartix ##
##########################
out_df = matrix(nrow=length(BEAST_postGTdistr[[1]]), ncol=n_genes, NA)


######################################
## Calculate various tree distances ##
######################################
# Note: [2] = Kuhner-Felsenstein distance
for (g in 1:n_genes) {
for (i in 1:length(BEAST_postGTdistr[[g]])) {
out_df[i,g] = treedist(BEAST_postGTdistr[[g]][[i]], starBEAST_postGTdistr[[g]][[i]], check.labels = TRUE)[2]
}}


######################
## Add column names ##
######################
columnnames = lapply(sprintf("%03d", 1:10), function(x){paste("gene", x, sep="")})
colnames(out_df) = columnnames


####################
## Stack the data ##
####################
results = melt(out_df, id.vars=c(columnnames))


###########################
## Add grouping variable ##
###########################
results[,4] = rep("gene002-gene010", length(results[,3]))
results[which(results[,2]=="gene001"),4] = "gene001"


##################
## Add colnames ##
##################
colnames(results) = c('generation', 'gene', 'KFdist', 'swapped')


######################################
## Save results to file as R-object ##
######################################
outFn=paste("BEAST.to.starBEAST_KFdist_g1swap_sim.", s, "_data", sep="")
assign(outFn, results)
save(list=outFn, file=paste("BEAST.to.starBEAST_KFdist_g1swap_sim.", s, "_data.rda", sep=""))
#outFn=paste("BEAST.to.starBEAST_KFdist_noswap_sim.", s, "_data", sep="")
#assign(outFn, results)
#save(list=outFn, file=paste("BEAST.to.starBEAST_KFdist_noswap_sim.", s, "_data.rda", sep=""))


###################
## Visualization ##
###################
my_plot = ggplot(data=results, aes(x=KFdist, group=gene, color=factor(swapped))) +
geom_density(line=2) +
ggtitle(paste("sim.", s, sep="")) +
theme_bw() +
theme(legend.position = "none",
plot.title = element_text(size = rel(1.5)))


###############
## Save plot ##
###############
assign(paste("sim.", s, sep=""), my_plot)


##################
## End sim-loop ##
##################
}


###############################
## Setup layout for plotting ##
###############################
svg("~/Desktop/BEAST.to.starBEAST_KFdist_g1swap.svg", width=20, height=15)
#svg("~/Desktop/BEAST.to.starBEAST_KFdist_noswap.svg", width=20, height=15)
grid.arrange(
sim.01, sim.02, sim.03, sim.04, sim.05,
sim.06, sim.07, sim.08, sim.09, sim.10,
sim.11, sim.12, sim.13, sim.14, sim.15,
sim.16, sim.17, sim.18, sim.19, sim.20, ncol=5)
dev.off()
