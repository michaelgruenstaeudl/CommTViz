#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2014-2015 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl")
#email = "gruenstaeudl.1@osu.edu"
#version = "2015.03.17.1800"


####################
# Global Variables #
####################
#in_folder = "/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/01_starBEAST_GeneTrees/07_parsed_results/"
#out_folder = "/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/01_starBEAST_GeneTrees/08_visualization/"

#n_sims = 20
#n_genes_fixed = 10
#n_sum_stats = 3
#nums = sprintf("%02d", c(1:n_sims))
#fn_prefixes = c("noswap", "g1swap", "g1a2swap", "g1t4swap")


dotPlots = function (n_genes_fixed, n_sum_stats, nums, fn_prefix) {

    #############
    # Libraries #
    #############
    require(ggplot2)

    ## Load .rda-files
    for (num in nums) {
          load(paste(in_folder, "CommT.Jan2015.", fn_prefix, ".sim.", num, ".rda", sep=""))
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

    # perGene #
    dataHandle_a0.05_perGene = list()
    for (num in nums) {
      name_handle = paste("CommT.Jan2015.", fn_prefix, ".sim.", num, "$results$alpha0.05$perGene", sep="")
      dataHandle_a0.05_perGene[[num]] = eval(parse(text = name_handle))
    }
    # acrGenes #
    dataHandle_a0.05_acrGene = list()
    for (num in nums) {
      name_handle = paste("CommT.Jan2015.", fn_prefix, ".sim.", num, "$results$alpha0.05$acrGene", sep="")
      dataHandle_a0.05_acrGene[[num]] = eval(parse(text = name_handle))
    }

    ######################
    # STEP2. Format data #
    ######################

    # perGene - 5 loci at alpha=0.05 #
    dataHandle_pG_0.05 = do.call("rbind", wrapper(dataHandle_a0.05_perGene, nums, n_genes_fixed))
    # acrGenes - 5 loci at alpha=0.05 #
    dataHandle_aG_0.05 = do.call("rbind", wrapper(dataHandle_a0.05_acrGene, nums, n_genes_fixed, special=T))

    #dataHandle_0.05 = rbind(dataHandle_pG_0.05, dataHandle_aG_0.05)

    ####################################
    # STEP3. Extract relevant sections #
    ####################################

    # Only observe pergene values
    dataHandle_0.05 = dataHandle_pG_0.05

    # Extract ndc values
    out_data = dataHandle_0.05[which(dataHandle_0.05[,2]=="ndc"),]

    # convert to characters for better plotting
    out_data[,'gene'] = as.character(out_data[,'gene'])

    #####################
    # STEP4. Make plots #
    #####################

    my_plot = ggplot(data=out_data, aes(x=sim, y=gene)) +
        geom_point(aes(colour=value), size=3, alpha=1.0) +
        scale_colour_manual(values=c(NA, 'black')) +
        facet_grid(stat ~ .) +
        ggtitle(paste(fn_prefix, ", alpha = 0.05", "\n", sep="")) +
        theme_bw() +
        scale_x_discrete(breaks=c(nums), labels=c(nums)) +
        #scale_y_discrete(limits=c("cv", "mode", "median", "mean", "sum", c(n_genes_fixed:1))) +
        scale_y_discrete(limits=as.character(c(n_genes_fixed:1))) +
        theme(legend.position = "none",
              axis.text = element_text(size=8),
              strip.background=element_rect(fill="white")) +
        xlab("\nSimulations") + 
        ylab("Genes\n")

    ggsave(filename=paste(out_folder, fn_prefix, "_0.05.svg", sep=""),
           plot=my_plot, width=5, height=3)

    ######################
    # STEP5. Return data #
    ######################
    return(out_data)

}
