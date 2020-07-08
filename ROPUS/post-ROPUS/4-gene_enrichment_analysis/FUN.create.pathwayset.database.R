# compare to all.gene list and get rid of any genes that are missing = in effect, sizing down the db

# pathwayset.gmt.path = folder with .gmt (db) file
# SUBGROUP ... IF NEED TO FILTER BY CC. MF BP 

# TERM CLASS GENES

FUN.create.pathwayset.database = function(
  pathwayset.gmt.path, 
  all.query.vector,    # ALL GENE NAMES TESTED 
  min.threshold,       # e.g. 10
  max.threshold,
  SUBGROUP=NULL
) {
  
  if (!file.exists(pathwayset.gmt.path)) {
    stop
    print("ERROR! Pathwayset file does not exist!")
  }
  
  if (!grepl(".gmt",pathwayset.gmt.path)) {
    stop
    print("ERROR! Pathwayset file needs to be in .gmt format!")
  }
  
  if (!is.character(all.query.vector)) {
    stop
    print("ERROR! all.query.vector needs to be all characters!")
  }
  
  
  ### "ps: pathwayset"; "op": output; "len": length
  
  ### original pathwayset
  ps0 = readLines(pathwayset.gmt.path)
  ps0 = strsplit(ps0, "\t")
  
  if (!is.null(SUBGROUP)) {
    keep = unlist(lapply(ps0, function(x) {as.character(x[2]) %in% SUBGROUP}))
    ps0 = ps0[keep]
  }
  
  anno = data.frame(matrix(NA, nrow=length(ps0), ncol=3))
  colnames(anno) = c("Original.length", "Filter.length", "Annotation")
  
  names(ps0) = rownames(anno) = unlist(lapply(ps0, function(x) x[1]))
  anno[,"Annotation"] = unlist(lapply(ps0, function(x) x[2]))
  
  ps0 = lapply(ps0, function(x) unique(x[-c(1:2)]))
  anno[,"Original.length"] = listLen(ps0)
  
  
  ### keep only those present in all.query.vector   # ALL GENES TESTED
  ps1 = lapply(ps0, function(x) x[x %in% all.query.vector])
  anno[,"Filter.length"] = listLen(ps1)
  
  
  ### remove according to threshold
  index.keep = which(anno[, "Filter.length"] >= min.threshold & anno[, "Filter.length"] <= max.threshold)
  
  print(paste("original number of pathwaysets:", length(ps0)))
  print(paste("filtered number of pathwaysets:", length(index.keep)))
  

  ### output in list
  output = list(
    pathwayset=ps1[index.keep], 
    anno=anno[index.keep,], 
    original.no=length(ps0), 
    filter.no=length(index.keep)
  )
  
  return(output)
}    
