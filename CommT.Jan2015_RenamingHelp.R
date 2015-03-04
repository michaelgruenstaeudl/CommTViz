setwd("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/02_output/")

n_sims = 20
nums = sprintf("%02d", c(1:n_sims))
fn_prefixes = c("BEAST.to.starBEAST_KFdist_")

# Generate infile names
inFns = c()
for (fns in fn_prefixes) {
    for (num in nums) {inFns = append(inFns, paste(fns, "noswap_sim.", num, sep=""))}
    #for (num in nums) {inFns = append(inFns, paste(fns, "g1swap_sim.", num, sep=""))}
}

# Load infiles
for (fns in fn_prefixes) {
    for (num in nums) {
        inVn=paste(fns, "regular_sim.", num, "_data", sep="")
        outVn=paste(fns, "noswap_sim.", num, "_data", sep="")
        load(paste(fns, "noswap_sim.", num, "_data.rda", sep=""))
        assign(outVn, eval(parse(text = inVn)))
        save(list=outVn, file=paste(fns, "noswap_sim.", num, "_data.rda", sep=""))
    }
    #for (num in nums) {
    #    inVn=paste(fns, "swapped_sim.", num, "_data", sep="")
    #    outVn=paste(fns, "g1swap_sim.", num, "_data", sep="")
    #    load(paste(fns, "g1swap_sim.", num, "_data.rda", sep=""))
    #    assign(outVn, eval(parse(text = inVn)))
    #    save(list=outVn, file=paste(fns, "g1swap_sim.", num, "_data.rda", sep=""))
    #}
}
