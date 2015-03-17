#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2014 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl")
#email = "gruenstaeudl.1@osu.edu"
#version = "2015.02.20.1200"

#############
# Load data #
#############

# Global variables
#setwd("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/02_BEAST_GeneTrees/07_parsed_results/")
setwd("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/02_BEAST_GeneTrees/07_parsed_results/no_swap/")

n_sims = 20
nums = sprintf("%02d", c(1:n_sims))
#fn_prefixes = c("CommT.Jan2015.fitonly.sim.")
fn_prefixes = c("CommT.Jan2015.noswap.sim.")

# Generate infile names
inFns = c()
for (fns in fn_prefixes) {
    for (num in nums) {inFns = append(inFns, paste(fns, num, sep=""))}
}

# Load infiles
for (inFn in inFns) {load(paste(inFn, ".rda", sep=""))}

########################
# Perform calculations #
########################

out_handle = list()
for (inFn in inFns) {

    # Extract the difference distributions
    handle = eval(parse(text = paste(inFn, "$rawStats$RAY$dif", sep="")))

    # Change column names to numbers
    colnames(handle) = 1:length(colnames(handle))

    # Stack the different columns
    handle = stack(data.frame(handle, stringsAsFactors=FALSE))

    # Convert the gene names ("X1", "X2", ...) into integers
    handle[,2] = as.integer(handle[,2])

    # Provide appropriate column names
    colnames(handle) = c("values", "gene")

    # Conduct Kruskal-Wallis-Test
    KW_result = kruskal.test(values ~ gene, data = handle)

    # Conduct one-way ANOVA
    av = aov(values ~ gene, data = handle)
    AV_result = summary(av)

    # Save result to out_handle
    out_handle[[inFn]] = list("KruskalWallis" = KW_result, "ANOVA" = AV_result)
}

###############
# Save output #
###############

outFn = file(paste(fn_prefixes, "DistrComparison.results", sep=""), "a")
sink(outFn, append=F)
print(out_handle)
cat("\n")
sink()
