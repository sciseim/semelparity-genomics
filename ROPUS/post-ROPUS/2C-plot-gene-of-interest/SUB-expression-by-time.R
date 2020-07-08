exprs <- as.matrix(counts.by.samling.date.subset.gene)

exprs <- rbind(exprs,exprs) # bind


dim(exprs) # yup, 8 samples and 9,927 genes for male liver
# and for each time value containing replicates, we calculate the TPM means
# if there are no replicates, we keep the initial TPM
count=1
for ( i in unique(time) ){
  if ( dim(as.data.frame(exprs[,which(time==i)]))[2] == 1 ){
    mean_rpkm=data.frame(exprs[,which(time==i)])
  } else {
    mean_rpkm=data.frame(rowMeans(exprs[,which(time==i)]))
  }
  colnames(mean_rpkm)=i
  if (count == 1){
    mean_rpkm_ok=mean_rpkm
  } else {
    mean_rpkm_ok=merge(mean_rpkm_ok,mean_rpkm,by="row.names")
    rownames(mean_rpkm_ok)=mean_rpkm_ok[,1]
    mean_rpkm_ok=mean_rpkm_ok[,-1]
  }
  count=count+1
}
# here we have an expression matrix containing one column per time value (and not one column per sample)
exprs_with_time=as.matrix(mean_rpkm_ok, header=TRUE, sep="\t",row.names=1,as.is=TRUE)
dim(exprs_with_time)
# 9927     5 # since we have e.g. two replicates 25th Sep and 2nd of Oct
colnames(exprs_with_time)


exprs_with_time <- exprs_with_time[1,]
# row.names(exprs) <- targetgene