#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2014-2015 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl")
#email = "gruenstaeudl.1@osu.edu"
#version = "2015.03.02.2000"


#############
# Load data #
#############

library(PMCMR)

# Global variables
setwd("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/02a_output_KF/")
fn_prefixes = c("BEAST.to.starBEAST_KFdist_noswap_sim.")
#fn_prefixes = c("BEAST.to.starBEAST_KFdist_g1swap_sim.")

n_sims = 20
nums = sprintf("%02d", c(1:n_sims))


# Generate infile names
inFns = c()
for (fns in fn_prefixes) {
    for (num in nums) {inFns = append(inFns, paste(fns, num, sep=""))}
}

# Load infiles
for (inFn in inFns) {load(paste(inFn, "_data.rda", sep=""))}

########################
# Perform calculations #
########################

out_handle = list()
for (inFn in inFns) {

    # Extract the difference distributions
    handle = eval(parse(text = paste(inFn, "_data", sep="")))
    
    # Conduct Kruskal-Wallis test
    KW_result = kruskal.test(KFdist ~ gene, data = handle)
    # Conduct posthoc pairwise Kruskal-Wallis tests
    PH_KW_result = posthoc.kruskal.nemenyi.test(x=handle$KFdist, g=handle$gene, method="Chisquare")

    # Conduct one-way ANOVA
    av = aov(KFdist ~ gene, data = handle)
    AV_result = summary(av)
    # Calculate Tukey Honest Significant Differences
    TukeyHSD_result = TukeyHSD(av)
    # Conduct two-way ANOVA
    twav = aov(KFdist ~ (gene+swapped), data = handle)
    TWAV_result = summary(twav)
    
    # Conduct pairwise t-tests with Bonferroni correction
    PTT_BF_result = pairwise.t.test(handle$KFdist, handle$gene, p.adj = "bonf")
    # Conduct pairwise t-tests with Holm correction
    PTT_HLM_result = pairwise.t.test(handle$KFdist, handle$gene, p.adj = "holm")

    ###############
    # Save output #
    ###############
    # Save result to out_handle
    out_handle = list("KruskalWallis" = KW_result,
                      "Posthoc_KruskalWallis" = PH_KW_result,
                      "ANOVA" = AV_result,
                      "Tukey_HonestSignDiff" = TukeyHSD_result,
                      "TwoWayANOVA" = TWAV_result,
                      "PairwiseTtests_BonferroniCorr" = PTT_BF_result,
                      "PairwiseTtests_HolmCorr" = PTT_HLM_result)
    outFn = file(paste(inFn, "_stats.txt", sep=""))
    sink(outFn, append=F)
    print(out_handle)
    cat("\n")
    sink()

#################
# Visualization #
#################
    ggplot(data=handle, aes(x=KFdist, group=gene)) +
    geom_histogram(binwidth=0.0025) +
    facet_grid(gene ~ .) +
    theme_bw()    
    ggsave(paste(inFn, ".svg", sep=""))
    ggsave(paste(inFn, ".eps", sep=""))
    
}
