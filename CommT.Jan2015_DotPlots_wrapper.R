#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2014-2015 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl")
#email = "gruenstaeudl.1@osu.edu"
#version = "2015.03.17.1800"

#############
# Libraries #
#############
require(ggplot2)

####################
# Global Variables #
####################
in_folder = "/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/01_starBEAST_GeneTrees/07_parsed_results/"
out_folder = "/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/01_starBEAST_GeneTrees/08_visualization/"

n_sims = 20
n_genes_fixed = 10
n_sum_stats = 3
nums = sprintf("%02d", c(1:n_sims))
#fn_prefixes = c("noswap", "g1swap", "g1a2swap", "g1t4swap")


dataHandle_0.05_noswap = dotPlots(n_genes_fixed, n_sum_stats, nums, "noswap")
dataHandle_0.05_noswap[, ncol(dataHandle_0.05_noswap)+1] = "noswap"
colnames(dataHandle_0.05_noswap)[ncol(dataHandle_0.05_noswap)] = "sim_type"

dataHandle_0.05_g1swap = dotPlots(n_genes_fixed, n_sum_stats, nums, "g1swap")
dataHandle_0.05_g1swap[, ncol(dataHandle_0.05_g1swap)+1] = "g1swap"
colnames(dataHandle_0.05_g1swap)[ncol(dataHandle_0.05_g1swap)] = "sim_type"

dataHandle_0.05_g1a2swap = dotPlots(n_genes_fixed, n_sum_stats, nums, "g1a2swap")
dataHandle_0.05_g1a2swap[, ncol(dataHandle_0.05_g1a2swap)+1] = "g1a2swap"
colnames(dataHandle_0.05_g1a2swap)[ncol(dataHandle_0.05_g1a2swap)] = "sim_type"

dataHandle_0.05_g1t4swap = dotPlots(n_genes_fixed, n_sum_stats, nums, "g1t4swap")
dataHandle_0.05_g1t4swap[, ncol(dataHandle_0.05_g1t4swap)+1] = "g1t4swap"
colnames(dataHandle_0.05_g1t4swap)[ncol(dataHandle_0.05_g1t4swap)] = "sim_type"

dataHandle_0.05 = rbind(dataHandle_0.05_noswap, dataHandle_0.05_g1swap,
                        dataHandle_0.05_g1a2swap, dataHandle_0.05_g1t4swap)

dataHandle_0.05[, 'sim_type'] = factor(dataHandle_0.05[, 'sim_type'], levels = c("noswap", "g1swap", "g1a2swap", "g1t4swap"))

my_plot = ggplot(data=dataHandle_0.05, aes(x=sim, y=gene)) +
    geom_point(aes(colour=value), size=3, alpha=1.0) +
    scale_colour_manual(values=c(NA, 'black')) +
    facet_grid(sim_type ~ .) +
    ggtitle(paste("Combined, alpha = 0.05", "\n", sep="")) +
    theme_bw() +
    scale_x_discrete(breaks=c(nums), labels=c(nums)) +
    #scale_y_discrete(limits=c("cv", "mode", "median", "mean", "sum", c(n_genes_fixed:1))) +
    scale_y_discrete(limits=as.character(c(n_genes_fixed:1))) +
    theme(legend.position = "none",
          axis.text = element_text(size=8),
          strip.background=element_rect(fill="white")) +
    xlab("\nSimulations") + 
    ylab("Genes\n")

ggsave(filename=paste(out_folder, "Combined_0.05.svg", sep=""),
       plot=my_plot, width=5, height=10)
