# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@




# ASSIGN TO GENE CHARACTER VECTOR
# idTop <- targetgenes




# START OF BIG LOOP
for (db in gmt.databases) {
  # the.database <- "CGP"  # testing
  the.database <- db
  
  
  remove(pathwayset.gmt.path)
  
  
  
  
  # CHOOSE A DATABASE
  #########################################################################################################
  #
  # Molecular Signatures Database (MSigDB) 
  #
  #########################################################################################################
  # pathwayset.gmt.path = "custom_pathways.gmt"
  # 
  # http://www.broadinstitute.org/gsea/msigdb/collections.jsp
  #   databases with identifiers corresponding to "gene symbols" and "entrez gene ids" are available
  #   obviously, gene symbols are usually preferred!
  # 
  
  
  
  
  
  
  
  # ############################################################################################## 
  # ############################################################################################## 
  # ############################################################################################## 
  # H: hallmark gene sets
  # Hallmark gene sets summarize and represent specific well-defined biological states or processes and display coherent expression. These gene sets were generated by a computational methodology based on identifying overlaps between gene sets in other MSigDB collections and retaining genes that display coordinate expression. details
  if(the.database == "hallmark")
  {
  
    pathwayset.gmt.path = "./MIT-msigdb/db/h.all.v7.0.symbols.gmt"
    cat("hallmark")
    
    tryCatch(source("SUB-bootstrap.R"))
   
  }
      
  # ############################################################################################## 
  # ############################################################################################## 
  # ############################################################################################## 

  
  
  
  
    
  
  #  C2: curated gene sets
  
  if(the.database == "KEGG")
  {
    pathwayset.gmt.path = "./MIT-msigdb/db/C2: curated gene sets/CP:KEGG: KEGG gene sets/c2.cp.kegg.v7.0.symbols.gmt"
    tryCatch(source("SUB-bootstrap.R"))
  }
  
  if(the.database == "reactome")
  {
    pathwayset.gmt.path = "./MIT-msigdb/db/C2: curated gene sets/CP:REACTOME: Reactome gene sets/c2.cp.reactome.v7.0.symbols.gmt"
    tryCatch(source("SUB-bootstrap.R"))
  }
  
  # combo of KEGG, BIOCARTA, etc
  # if(the.database == "canonical.pathways")
  {
    # pathwayset.gmt.path = "./MIT-msigdb/db/C2: curated gene sets/CP: Canonical pathways/c2.cp.v4.0.symbols.gmt"
   #  tryCatch(source("SUB-bootstrap.R"))
  }
  
  
  
  # CGP: chemical and genetic perturbations
  # Gene sets represent expression signatures of genetic and chemical perturbations. A number of these gene sets come in pairs: xxx_UP (and xxx_DN) gene set representing genes induced (and repressed) by the perturbation
  if(the.database == "CGP")
  {
    pathwayset.gmt.path = "./MIT-msigdb/db/C2: curated gene sets/c2.cgp.v7.0.symbols.gmt" # all curated
     tryCatch(source("SUB-bootstrap.R"))
  }
  
  
  
  
  # C5: GO gene sets
  if(the.database == "GO.BP")
  {
    pathwayset.gmt.path = "./MIT-msigdb/db/C5: GO gene sets//BP: GO biological process/c5.bp.v7.0.symbols.gmt" # GO term BP
     tryCatch(source("SUB-bootstrap.R"))
  }
  
  if(the.database == "GO.CC")
  {
    pathwayset.gmt.path = "./MIT-msigdb/db/C5: GO gene sets//CC: GO cellular component/c5.cc.v7.0.symbols.gmt" # GO term CC
     tryCatch(source("SUB-bootstrap.R"))
  }
  
  # if(the.database == "GO.all")
  {
    # pathwayset.gmt.path = "./MIT-msigdb/db/C5: GO gene sets/c5.all.v7.0.symbols.gmt" # all GO term
   #  tryCatch(source("SUB-bootstrap.R"))
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # DATABASES THAT WILL PROVE MORE USEFUL FOR OTHER DATA SETS!  
  #
  

  
  
  # C3: motif gene sets
  # Gene sets that contain genes that share a cis-regulatory motif that is conserved across the human, mouse, rat, and dog genomes. 
  # ... TFT: transcription factor targets
  if(the.database == "cis.motifs")
  {
    pathwayset.gmt.path = "./MIT-msigdb/db/C3 motif gene sets/c3.tft.v7.0.symbols.gmt" # C3: motif gene sets
     tryCatch(source("SUB-bootstrap.R"))
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # C7: immunologic signatures
  if(the.database == "immuno.sig")
  {
    pathwayset.gmt.path = "./MIT-msigdb/db/C7 immunologic signatures/c7.all.v7.0.symbols.gmt"
     tryCatch(source("SUB-bootstrap.R"))
  }
    # My custom GMT (GenAge, senescenece signatures from Current Biology and NAR papers)
  if(the.database == "custom")
  {
    pathwayset.gmt.path = "./MIT-msigdb/db/MyCustom.gmt"
     tryCatch(source("SUB-bootstrap.R"))
  }
  } # END OF BIG LOOP

