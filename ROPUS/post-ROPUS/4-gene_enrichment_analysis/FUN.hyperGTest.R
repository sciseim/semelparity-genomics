
workingdirectory <- getwd()

source(paste(workingdirectory, "/MIT-msigdb/","FUN.create.pathwayset.database.R", sep=""))

FUN.hyperGTest <- function(
  idTop,
  idAll,
  pathwayset.database.list,
  pathwayset.gmt.path=NULL,
  min.threshold=1,
  max.threshold=length(idAll),
  pvalue=0.05,
  SUBGROUP=NULL
) {
  
  require(Biobase)
  
  ##########################################
  ########## create pathway database #######
  ##########################################
  
  ### outputs: $pathwayset, $anno, $original.no, $filter.no  
  
  if (!is.null(pathwayset.gmt.path)) {
    pathwayset.database <- FUN.create.pathwayset.database(
      pathwayset.gmt.path,
      all.query.vector=idAll,
      min.threshold, 
      max.threshold,
      SUBGROUP
    )
    pathwayset.database.list <- pathwayset.database$pathwayset
  }
  
  ##########################################
  ########## hyperGTest ####################
  ##########################################
  
  output <- FUN.hyperGTest.functions(
    idTop,
    idAll,
    pathwayset.database.list
  )
  
  ##########################################
  ########## result output #################
  ##########################################
  
  top <- output$statistics
  
  top <- cbind.data.frame(top, AdjPvalue=p.adjust(top$pval, method="BH"))

  top <- top[which(top$pval <= pvalue),]
  rownames(top) <- NULL
  top <- cbind.data.frame(
    "P Value"=top$pval,
    "Adj Pvalue"=top$AdjPvalue,
    "Odds Ratio"=top$odds,
    "Expected Number"=top$expected,
    "Observed Number"=top$numWDrawn,
    "Pathway Size"=top$numW,
    "Pathway Name"=top$Pathway.ID
  )
  
  return(list(summary=top, raw.data=output))
}

############################################
############################################

FUN.hyperGTest.functions <- function(
  idTop,
  idAll,
  pathwayset.database.list
) {
  require("Biobase")
  
  psl <- pathwayset.database.list
  
  psl.n <- length(psl)
  
  ### remove id not represented by the pathwayset
  psl.all <- unique(unlist(psl))
  id.all <- idAll[idAll %in% psl.all]
  id.top <- idTop[idTop %in% id.all]
  
  ### "ct": count table
  ct <- data.frame(matrix(NA,nrow=psl.n, ncol=7))
  colnames(ct) <- c("numW", "numB", "numDrawn", "numWDrawn", "pval", "odds", "expected")
  
  ### white: metabolites in a pathway
  ### black: metabolites not in a pathway
  
  ct[,"numW"] <- listLen(psl)
  ct[,"numB"] <- length(id.all) - ct[,"numW"]
  
  ### numDrawn: number of top metabolites
  ct[,"numDrawn"] <- length(id.top)
  
  ### numWDrawn: number of pathway metabolite among the top metabolites
  ct[,"numWDrawn"] <- unlist(lapply(psl, function(x) sum(x %in% id.top)))
  
  ### pvalue
  ### based on ".doHyperGInternal" of Category
  ct[,"pval"] <- apply(ct[,1:4], 1, function(x) {
    phyper(x["numWDrawn"] - 1, x["numW"], x["numB"], x["numDrawn"], lower.tail = FALSE)
  })
  
  ### odd ratio
  ### based on ".doHyperGInternal" of Category
  ct[,"odds"] <- apply(ct[,1:4], 1, function(x) {
    n21 <- x["numW"] - x["numWDrawn"]
    n12 <- x["numDrawn"] - x["numWDrawn"]
    n22 <- x["numB"] - n12
    odds.ratio <- (x["numWDrawn"] * n22)/(n12 * n21)
    return(odds.ratio)
  })
  ct[,"expected"] <- apply(ct[,1:4], 1, function(x) {
    n21 <- x["numW"] - x["numWDrawn"]
    n12 <- x["numDrawn"] - x["numWDrawn"]
    n22 <- x["numB"] - n12
    expected <- (x["numWDrawn"] + n12) * (x["numWDrawn"] + n21)
    expected <- expected/(x["numWDrawn"] + n12 + n21 + n22)
    return(expected)
  })
  
  id.present <- lapply(psl, function(x) x[x %in% id.top])
  id.pathway <- psl
  
  result <- list()
  
  result$statistics <- cbind.data.frame("Pathway.ID"=names(psl), ct)
  rownames(result$statistics) <- names(psl)
  result$statistics <- result$statistics[order(result$statistics$pval, decreasing=F), ]
  
  result$id.present <- id.present[rownames(result$statistics)]
  result$id.pathway <- id.pathway[rownames(result$statistics)]
  
  return(result)  
}
