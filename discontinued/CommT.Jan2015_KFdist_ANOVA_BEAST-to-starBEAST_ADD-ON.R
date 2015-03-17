
setwd("/run/media/michael/6974AEDF6DB4DD0C/05_Project_CommT_Jan2015/03_KFdist_distributions/02a_output_KF/")
#fn_prefixes = c("BEAST.to.starBEAST_KFdist_noswap_sim.")
fn_prefixes = c("BEAST.to.starBEAST_KFdist_g1swap_sim.")
n_sims = 20
nums = sprintf("%02d", c(1:n_sims))
# Generate infile names
inFns = c()
for (fns in fn_prefixes) {
  for (num in nums) {inFns = append(inFns, paste(fns, num, sep=""))}
}
# Load infiles
for (inFn in inFns) {load(paste(inFn, "_data.rda", sep=""))}
out_handle = list()
for (inFn in inFns) {
  # Extract the difference distributions
  handle = eval(parse(text = paste(inFn, "_data", sep="")))
  handle[,4] = 0
  tmp_list = list()
  for (i in sprintf("%02d", 1:10)) {
    target = paste("gene0", i, sep="")
    handle[which(handle[,2]==target),4] = target
    handle[which(handle[,2]!=target),4] = "blocked"
    aov_results = summary(aov(KFdist ~ swapped + Error(gene), data = handle))
    #tmp_list[[i]] = TukeyHSD(aov(KFdist ~ swapped + gene, data = handle))
    tmp_list[[paste("gene0", i, sep="")]] = aov_results$'Error: gene'[[1]][1,5]
  }
  out_handle[[inFn]] = t(as.data.frame(tmp_list))
  out_handle[[inFn]] = round(out_handle[[inFn]], digits=5)
  colnames(out_handle[[inFn]]) = "Pr(>F)"
}
#outFn = file("~/Desktop/blockedDesign_noswap.txt")
outFn = file("~/Desktop/blockedDesign_g1swap.txt")
sink(outFn, append=F)
print(out_handle)
cat("\n")
sink()
