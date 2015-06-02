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
#in_folder = "/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/02_BEAST_GeneTrees/07_parsed_results/"
#out_folder = "/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/02_BEAST_GeneTrees/08_visualizations/"

#in_folder = "/home/michael/Desktop/CommT_Visualization/02_BEAST_GeneTrees/"
#out_folder = "/home/michael/Desktop/CommT_Visualization/output/"

in_folder = "~/TEMP/input/"
out_folder = "~/TEMP/output/"

selected_metric = "ndc"
alpha_level = "0.05"
n_sims = 20
n_genes_fixed = 10
n_sum_stats = 3
nums = sprintf("%02d", c(1:n_sims))
#fn_prefixes = c("noswap", "g1swap", "g1a2swap", "g1t4swap", "gAllShuf")


dataHandle_noswap = dotPlots(n_genes_fixed, n_sum_stats, nums, "noswap", alpha_level, selected_metric)
dataHandle_noswap[, ncol(dataHandle_noswap)+1] = "noswap"
colnames(dataHandle_noswap)[ncol(dataHandle_noswap)] = "sim_type"

dataHandle_g1swap = dotPlots(n_genes_fixed, n_sum_stats, nums, "g1swap", alpha_level, selected_metric)
dataHandle_g1swap[, ncol(dataHandle_g1swap)+1] = "g1swap"
colnames(dataHandle_g1swap)[ncol(dataHandle_g1swap)] = "sim_type"

dataHandle_g1a2swap = dotPlots(n_genes_fixed, n_sum_stats, nums, "g1a2swap", alpha_level, selected_metric)
dataHandle_g1a2swap[, ncol(dataHandle_g1a2swap)+1] = "g1a2swap"
colnames(dataHandle_g1a2swap)[ncol(dataHandle_g1a2swap)] = "sim_type"

dataHandle_g1t4swap = dotPlots(n_genes_fixed, n_sum_stats, nums, "g1t4swap", alpha_level, selected_metric)
dataHandle_g1t4swap[, ncol(dataHandle_g1t4swap)+1] = "g1t4swap"
colnames(dataHandle_g1t4swap)[ncol(dataHandle_g1t4swap)] = "sim_type"

dataHandle_gAllShuf = dotPlots(n_genes_fixed, n_sum_stats, nums, "gAllShuf", alpha_level, selected_metric)
dataHandle_gAllShuf[, ncol(dataHandle_gAllShuf)+1] = "gAllShuf"
colnames(dataHandle_gAllShuf)[ncol(dataHandle_gAllShuf)] = "sim_type"

dataHandle = rbind(dataHandle_noswap, dataHandle_g1swap,
                   dataHandle_g1a2swap, dataHandle_g1t4swap,
                   dataHandle_gAllShuf)

dataHandle[, 'sim_type'] = factor(dataHandle[, 'sim_type'], levels = c("noswap", "g1swap", "g1a2swap", "g1t4swap", "gAllShuf"))

my_plot = ggplot(data=dataHandle, aes(x=sim, y=gene)) +
    geom_point(aes(colour=value), size=3, alpha=1.0) +
    scale_colour_manual(values=c(NA, 'black')) +
    facet_grid(sim_type ~ .) +
    ggtitle(paste(selected_metric, ", alpha = ", alpha_level, "\n", sep="")) +
    theme_bw() +
    scale_x_discrete(breaks=c(nums), labels=c(nums)) +
    #scale_y_discrete(limits=c("cv", "mode", "median", "mean", "sum", c(n_genes_fixed:1))) +
    scale_y_discrete(limits=as.character(c(n_genes_fixed:1))) +
    theme(legend.position = "none",
          axis.text = element_text(size=8),
          strip.background=element_rect(fill="white")) +
    xlab("\nSimulations") + 
    ylab("Genes\n")

ggsave(filename=paste(out_folder, selected_metric, "_", alpha_level,".svg", sep=""),
       plot=my_plot, width=5, height=10)
