#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2014-2015 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl")
#email = "gruenstaeudl.1@osu.edu"
#version = "2015.02.01.1800"

#############
# Libraries #
#############
require(ggplot2)

####################
# Global Variables #
####################
setwd("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/02_BEAST_GeneTrees/07_parsed_results/")

#setwd("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/01_starBEAST_GeneTrees/07_parsed_results/sorted/")


n_sims = 20
n_genes_fixed = 10
n_sum_stats = 3
nums = sprintf("%02d", c(1:n_sims))
#fn_prefixes = c("CommT.Jan2015.noswap.sim.", "CommT.Jan2015.g1swap.sim.")
fn_prefixes = c("CommT.Jan2015.fitonly.sim.")

## Load .rda-files
for (fns in fn_prefixes) {
    for (num in nums) {
      load(paste(fns, num, ".rda", sep=""))
    }
}

####################
# Helper Functions #
####################

customStack = function(inData, simNum, nLoci, special=FALSE) {
# Custom function which 
# (a) converts results matrices into presence/absence matrices,
# (b) stacks the matrix columns,
# (c) adds identifier information

  handle = inData
  # Order of stats must be alphabetic because input data
  # also sorted alphabetically
  colnames(handle) = c("gtp", "ndc", "ray")
  handle[grepl("\\*", handle)] = 1
  handle[grepl("n.s.", handle)] = 0
  handle[grepl(" 0", handle)] = 0
  handle = stack(data.frame(handle, stringsAsFactors=FALSE))
  colnames(handle)[1] = "value"
  colnames(handle)[2] = "stat"
  if (special) {
    handle[,3] = rep(c("sum", "mean", "median", "mode", "cv"), n_sum_stats)
  }
  else {handle[,3] = rep(c(1:nLoci), n_sum_stats)}
  colnames(handle)[3] = "gene"
  handle[,4] = simNum
  colnames(handle)[4] = "sim"
  return(handle)
}

wrapper = function(inData, nums, nLoci, special=FALSE) {
  out_list = list()
  for (num in nums) {
    out_list[[num]] = customStack(inData[[num]], num, nLoci, special)
  }
  return(out_list)
}

########################################
# STEP1. Load data and save into lists #
######################################## 

# perGene - regular #
dataHandle_reg_a0.01_perGene = list()
for (num in nums) {
  name_handle = paste(fn_prefixes[[1]], num, "$results$alpha0.01$perGene", sep="")
  dataHandle_reg_a0.01_perGene[[num]] = eval(parse(text = name_handle))
}
# acrGenes - regular #
dataHandle_reg_a0.01_acrGene = list()
for (num in nums) {
  name_handle = paste(fn_prefixes[[1]], num, "$results$alpha0.01$acrGene", sep="")
  dataHandle_reg_a0.01_acrGene[[num]] = eval(parse(text = name_handle))
}
## perGene - swapped #
#dataHandle_swap_a0.01_perGene = list()
#for (num in nums) {
#  name_handle = paste(fn_prefixes[[2]], num, "$results$alpha0.01$perGene", sep="")
#  dataHandle_swap_a0.01_perGene[[num]] = eval(parse(text = name_handle))
#}
## acrGenes - swapped #
#dataHandle_swap_a0.01_acrGene = list()
#for (num in nums) {
#  name_handle = paste(fn_prefixes[[2]], num, "$results$alpha0.01$acrGene", sep="")
#  dataHandle_swap_a0.01_acrGene[[num]] = eval(parse(text = name_handle))
#}

######################
# STEP2. Format data #
#######################

# perGene - 5 loci at alpha=0.01 #
dataHandle_reg_pG_0.01 = do.call("rbind", wrapper(dataHandle_reg_a0.01_perGene, nums, n_genes_fixed))
dataHandle_reg_pG_0.01[, ncol(dataHandle_reg_pG_0.01)+1] = "regular"
colnames(dataHandle_reg_pG_0.01)[ncol(dataHandle_reg_pG_0.01)] = "sim_type"

# acrGenes - 5 loci at alpha=0.01 #
dataHandle_reg_aG_0.01 = do.call("rbind", wrapper(dataHandle_reg_a0.01_acrGene, nums, n_genes_fixed, special=T))
dataHandle_reg_aG_0.01[, ncol(dataHandle_reg_aG_0.01)+1] = "regular"
colnames(dataHandle_reg_aG_0.01)[ncol(dataHandle_reg_aG_0.01)] = "sim_type"

## perGene - 10 loci at alpha=0.01 #
#dataHandle_swap_pG_0.01 = do.call("rbind", wrapper(dataHandle_swap_a0.01_perGene, nums, n_genes_fixed))
#dataHandle_swap_pG_0.01[, ncol(dataHandle_swap_pG_0.01)+1] = "swapped"
#colnames(dataHandle_swap_pG_0.01)[ncol(dataHandle_swap_pG_0.01)] = "sim_type"

## acrGenes - 10 loci at alpha=0.01 #
#dataHandle_swap_aG_0.01 = do.call("rbind", wrapper(dataHandle_swap_a0.01_acrGene, nums, n_genes_fixed, special=T))
#dataHandle_swap_aG_0.01[, ncol(dataHandle_swap_aG_0.01)+1] = "swapped"
#colnames(dataHandle_swap_aG_0.01)[ncol(dataHandle_swap_aG_0.01)] = "sim_type"

dataHandle_0.01 = rbind(dataHandle_reg_pG_0.01, dataHandle_reg_aG_0.01)#,
                        #dataHandle_swap_pG_0.01, dataHandle_swap_aG_0.01)

#####################
# STEP4. Make plots #
#####################

plot_0.01 = ggplot(data=dataHandle_0.01, aes(x=sim, y=gene)) +
    geom_point(aes(colour=value), size=3, alpha=1.0) +
    scale_colour_manual(values=c(NA, 'black')) +
    facet_grid(stat ~ sim_type) +
    ggtitle(paste("Regular and Gene1-Swapped Simulations",
                  "alpha = 0.01", "\n",
                  sep="\n")) +
    theme_bw() +
    scale_x_discrete(breaks=c(nums), labels=c(nums)) +
    scale_y_discrete(limits=c("cv", "mode", "median", "mean", "sum", c(n_genes_fixed:1))) +
    theme(axis.text = element_text(size=6),
          strip.background=element_rect(fill="white")) +
    xlab("\nSimulations") + 
    ylab("Genes\n")

svg("/home/michael/Desktop/CommTJan2015_0.01.svg", width=7, height=6)
plot_0.01
dev.off()
