# remove(error.per.sample)
error.per.sample <- ""
for (date in unique(time)) {
  target.names <-  names(  counts.by.samling.date.subset.gene[ grep(date, names(counts.by.samling.date.subset.gene)) ]  )
  temp <- counts.by.samling.date.subset.gene[,c(target.names)]
  temp <- melt(temp)
  temp$group <- "blah"
  head(temp)
  source("SUB-summarySE.R")
  the.stats <- summarySE(temp, measurevar="value", groupvars=c("group"))
  temp$se <- the.stats$se
  unique.error.value <- unique(temp$se) # unique so that we only add one
  error.per.sample <- c(error.per.sample,unique.error.value)
}
# add to time variable
error.per.sample <- error.per.sample[-1] # remove empty
error.per.sample

error.per.sample[is.na(error.per.sample)] <- 0 # replace NA with 0 for the error bar
class(error.per.sample)
error.per.sample <- as.numeric(error.per.sample)



