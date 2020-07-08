idTop <- targetgenes

pathwayset = FUN.create.pathwayset.database(pathwayset.gmt.path, all.query.vector = idAll, min.threshold = 1, max.threshold = 10000)
temp = FUN.hyperGTest(idTop = idTop, idAll = idAll, pathwayset.database.list = pathwayset$pathwayset,pvalue = 0.05)
pathwayset = FUN.create.pathwayset.database(pathwayset.gmt.path, all.query.vector = idAll, min.threshold = 1, max.threshold = 10000)
pathwayset.database.list = pathwayset$pathwayset
#
times = 5000
p.value = 0.05
enrichment = FUN.hyperGTest(idTop, idAll, pathwayset.database.list, pvalue = p.value)$summary


# if only one hit, script will fail, so let us just duplicate
if(length(enrichment[,2]) == 1)
enrichment <- rbind(enrichment,enrichment)  





#
all.name = names(pathwayset.database.list)
num = length(idTop)

result = list()

progbar = utils::txtProgressBar(style = 3)
for (i in 1:times) {
  hold = FUN.hyperGTest(sample(idAll, size = num, replace = F), idAll, pathwayset.database.list, pvalue = 1)$summary
  hold2 = hold$"P Value"
  names(hold2) = hold$"Pathway Name"
  result[[i]] = hold2[all.name]
  utils::setTxtProgressBar(progbar, i/times)
}
close(progbar)

result = do.call(cbind, result)
result = result[as.character(enrichment$"Pathway Name"),]
bootstrap.p = c()
for (i in 1:nrow(enrichment)) {
  bootstrap.p[i] = sum(result[i,] <= enrichment[i, "P Value"]) / times
}
output = cbind.data.frame(enrichment, "Bootstrap P Value" = bootstrap.p)



# output genes!
outputnumber <- length(output[,1])
genes.result = list()
for (i in 1:outputnumber) 
{
  genes.result[[i]] = (   temp$raw.data$id.present[[i]]  ) # log'ed lifspan info
  print(i)  # prints how many genes have been looped through     
}
genes.resultUNLIST <- (unlist(lapply(genes.result, paste, collapse=","))) # separate using comma
# output table
write.table(output, paste("./output/tables/",output.name,"_",the.database,"-table.txt",sep=""), row.names=T, col.names=T, sep="\t", quote=F,append=FALSE)
# list of genes
write.table(genes.resultUNLIST,paste("./output/genes/",output.name,"_",the.database,"-genes.txt",sep=""), row.names=T, col.names=T, sep="\t", quote=F,append=FALSE)

