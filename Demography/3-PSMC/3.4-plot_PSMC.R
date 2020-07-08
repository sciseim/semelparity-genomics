rm(list=ls()) # clear

setwd("~/Downloads/PSMC_plotting/no-ChrX/")


# scientific  notation OFF
options(scipen=999)
# scientific  notation ON
options(scipen=0)


# Load libraries needed
library(scales)
library(zoo)


# want to draw PDFs using Arial so that we can add extra text in e.g. PowerPoint without it looking wonky
# see https://blog.revolutionanalytics.com/2012/09/how-to-use-your-favorite-fonts-in-r-charts.html
# extrafont
library(extrafont)
# font_import()
loadfonts()







# Open PDF file to write to
# pdf("psmc.plot.pdf", height = 7, width = 10)
pdf("psmc.plot.pdf", height = 7, width = 22,family="Arial")


# Adjust plot margins
par(mar = c(4.1, 4.1, 2.1, 4.1))

y.limit <- 100 # 100*10e4=1e+07=10M animals
x.limit <- 2*10^6 # in years  # 2M years
# x.limit <- 10^7 # in years
# 1000000 = 1M = 10^6



x.min <- 20000

# Start new empty plot
plot.new()


# log looks better, but warps x-axis a bit... see 
# https://stackoverflow.com/questions/13238273/uneven-spacing-on-axis-with-r
# log(20*100000)=14.5 on log-scale
# log(10*100000)=13.8 on log-scale
 plot.window(xlim = c(x.min, x.limit), ylim = c(0, y.limit), log = "x") # log
# plot.window(xlim = c(x.min, x.limit), ylim = c(0, y.limit))


box()

# axis(1, tck = -0.01, at = c(10000, 100000,1000000,10000000), labels = expression(10^4, 10^5, 10^6, 10^7))
# axis(1, tck = -0.01, at = c(20000, 100000,1000000,), labels = expression(10^4, 10^5, 10^6))
# 10*1000000 # to 10M

# 100k years ago
# 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20
# axis(1, tck = -0.01, at = c(10000, 100000,1000000,10000000), labels = expression(10^4, 10^5, 10^6, 10^7))

# (0*100000)
# 0,

axis(1, tck = -0.01, at = c((0.2*100000),(0.5*100000),(1*100000),(2*100000),(3*100000),(4*100000),(5*100000),(6*100000),(7*100000),(8*100000),(9*100000),(10*100000),(11*100000),(12*100000),(13*100000),(14*100000),(15*100000),(16*100000),(17*100000),(18*100000),(19*100000),(20*100000)), labels = expression(0.2,0.5,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))





# 20*100000 # 2M

# axis(1, tck = -0.01, at = c(20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 2000000,3000000,4000000,5000000,6000000,7000000,8000000,9000000,20000000,30000000,40000000), labels = F)

# Pleistocene AKA ice age : 2,588,000 to 11,700
# Holocene: 11,700 to now ...
# Pliocene 5.333 million to 2.58
# ...Piacenzian (3.600–2.58 Ma) 
# ...Zanclean (5.333–3.600 Ma)
# The Piacenzian is sometimes referred to as the Late Pliocene, whereas the Zanclean is referred to as the Early Pliocene.
# draw shaded lines
# panel.first = rect(c(1,7), -1e6, c(3,10), 1e6, col='green', border=NA)
# The first two arguments c(1,7) are the starting values for the shaded rectangle, and following arguments c(3,10) are where the shading ends. This creates a shaded region from 1-3 and 7-10
# panel.first = rect(c(0,7), -1e6, c(11700), 1e6, col='white', border=NA) # Holocene
# panel.first = rect(c(11700,7), -1e6, c(2588000), 1e6, col='#E9F4FB', border=NA) # Pleistocene
# panel.first = rect(c(2588000,7), -1e6, c(5333000,10), 1e6, col='white', border=NA) # Pliocene
# panel.first = rect(c(5333000,7), -1e6, c(23030000,10), 1e6, col='#E9F4FB', border=NA) # Miocene

# Last glacial period (LGP): 115,000 – c. 11,700
# last ice age 
# AKA 'Tarantian'
# Late Pleistocene: 0.0117	0.126 ... so (0.126*1000000) (0.0117*1000000) 
# Early Pleistocene
# AKA  the Lower Pleistocene
# 2.588 to 0.781 ... The Early Pleistocene consists of the Gelasian and the Calabrian ages ... so separated by 'Chibanian'	0.126	0.781
panel.first = rect(c((0.0117*1000000),7), -1e6, c((0.126*1000000) ,10), 1e6, col='#e1f6fd', border=NA) # 
panel.first = rect(c((0.781*1000000) ,7), -1e6, c((2.588*1000000) ,10), 1e6, col='#e1f6fd', border=NA) 
# mid is 
# panel.first = rect(c((0.126*1000000) ,7), -1e6, c((0.781*1000000) ,10), 1e6, col='#dfdfdf', border=NA) 

# The informal term "Late Quaternary" refers to the past 0.5–1.0 million years
# panel.first = rect(c((0.500*1000000) ,7), -1e6, c((1.000*1000000) ,10), 1e6, col='#dfdfdf', border=NA) 


# Late Pleistocene: 0.0117	0.126 ... so (0.126*1000000) (0.0117*1000000) 
# panel.first = rect(c((0.0117*1000000) ,7), -1e6, c((0.126*1000000) ,10), 1e6, col='#f0f9bf', border=NA) 

# Last glacial period (LGP): 115,000 – c. 11,700
# panel.first = rect(c((0.0117*1000000) ,7), -1e6, c((0.115*1000000) ,10), 1e6, col='#bed0ff', border=NA) 


# fossil of extant 
# no fossils of extant dasyurid genera are known from earlier than the Pliocene (Wroe 2003).
# dasyurid cladogenesis 
# 10-13
#  dasyurid cladogenesis time point: emergence of Murexia and Antechinus 
# ref: PHYLOGENETIC RELATIONSHIPS OF THE DASYURID MARSUPIAL GENUS MUREXIA
# abline(v = 10000000, lty = "dashed", col = "#ff4d4d") # earliest fossil record of 
# axis(1, tck = -0.01, at = c(10000000), labels = F,col="#ff4d4d",lwd=5) # emergence of Murexia and Antechinus 
# abline(v = 48800, lty = "dashed", col = "black") # earliest fossil record of 


# estimated extinction times of 16 megafaunal genera in mainland Australia
# points(x=500000, col="#ff4d4d", pch=19)
# FosSahul is the first database compiling the ages of nonhuman vertebrate fossils from the Middle Pleistocene to the present in the Sahul region. It includes comprehensive metadata with ratings of reliability allocated to each fossil age. Because ecological and evolutionary phenomena are time-dependent, the entire range of archaeological and palaeontological research disciplines benefit from the availability of this data.
# The initial arrival of humans in Australia and New Guinea (then connected by a land bridge and hereafter referred to collectively as ‘Sahul’

# latest extinction time -- from https://www.nature.com/articles/ncomms10511/tables/1
axis(1, tck = -0.01, at = c(63*1000), labels = F,col="red",lwd=5) # Congruus
axis(1, tck = -0.01, at = c(44.2*1000), labels = F,col="red",lwd=5) # Diprotodon
axis(1, tck = -0.01, at = c(35*1000), labels = F,col="red",lwd=5) # Genyornis
axis(1, tck = -0.01, at = c(43.3*1000), labels = F,col="red",lwd=5) # Macropus
axis(1, tck = -0.01, at = c(56*1000), labels = F,col="red",lwd=5) # Megalibgwilia	
axis(1, tck = -0.01, at = c(56*1000), labels = F,col="red",lwd=5) # Metasthenurus
axis(1, tck = -0.01, at = c(56*1000), labels = F,col="red",lwd=5) # Palorchestes
axis(1, tck = -0.01, at = c(46*1000), labels = F,col="red",lwd=5) # Phascolonus
axis(1, tck = -0.01, at = c(40*1000), labels = F,col="red",lwd=5) # Procoptodon
axis(1, tck = -0.01, at = c(40.8*1000), labels = F,col="red",lwd=5) # Protemnodon
axis(1, tck = -0.01, at = c(46*1000), labels = F,col="red",lwd=5) # Sarcophilus
axis(1, tck = -0.01, at = c(44.9*1000), labels = F,col="red",lwd=5) # Simosthenurus
axis(1, tck = -0.01, at = c(52*1000), labels = F,col="red",lwd=5) # Sthenurus
# axis(1, tck = -0.01, at = c(3.5*1000), labels = F,col="red",lwd=5) # Thylacinus
# Tasmanian tiger is too close and beyond the reliable PSMC, so do not draw
axis(1, tck = -0.01, at = c(44.9*1000), labels = F,col="red",lwd=5) # Thylacoleo
axis(1, tck = -0.01, at = c(44.9*1000), labels = F,col="red",lwd=5) # Zygomaturus









# Humans arrive
# abline(v = 48800, lty = "dashed", col = "black")
# abline(v = 10000, lty = "dashed", col = "red") # HOLOCENE
# abline(v = 22000, lty = "dashed", col = "red") # Last glacial maximum
# abline(v = 10000, lty = "dashed", col = "blue") # Last glacial period 10-120ka
# abline(v = 12000, lty = "dashed", col = "red") # Last glacial period 10-120ka
# PLEISTOCENE ~10?
# PLIOCENE
#@ axis(1, tck = -0.01, at = c(48800), labels = F,col="black",lwd=5) # first humans
# abline(v = 48800, lty = "dashed", col = "grey") # 
# ?axis
# panel.first = rect(c((0.0593*1000000) ,0.2), -1e6, c((0.065*1000000) ,0.2), 1e6, col='grey', border=NA) 
abline(v = (0.0593*1000000), lty = "dashed", col = "grey") # 
abline(v = (0.065*1000000), lty = "dashed", col = "grey") # 

# axis(1, tck = -0.01, at = c(6000), labels = F,col="#00cc66",lwd=5) # The Torres Strait Islands were formed when the land separating Australia and New Guinea was flooded by rising sea levels around 6000 BCE.
# cannot reliably see here...

# As the Australia-New Guinea tectonic plate ebbed to the north over millions of years, a group of mammals called marsupials evolved. When this plate collided with the Eurasian plate about 25 million years ago, New Guinea pushed up from under the seas and it was colonized by Australian marsupials. 
# Much later, a rise in sea level separated New Guinea from Australia and allowed mammal migrants from the south to evolve independently. 
# During the Ice Age, sea levels dropped again and a land connection was re-established to Australia. While some mammals crossed this divide, the unfavourable ecological conditions in this intermediate zone prevented most from doing so. 4



y.axis = seq(0, y.limit, by = 1)
axis(2, las = 1)
axis(2, las = 1, tck = -0.01, at = y.axis, labels = F)








library(ggsci)
# Science colours
# https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html#npg
# 2d357f
# e60006
# 107b35
# 4f0268
# 0e706d

col.black <- "#000000"
col.orange <- "#E69F00"
col.sky.blue <-"#56B4E9"
col.blueish.green <- "#009E73"
col.yellow <- "#F0E442"
col.blue <- "#0072B2"
col.vermillion <- "#D55E00"
reddish.purple <- "#CC79A7"


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SAMPLENAME <- "AdamAnt"

remove(bootstrap.list)
bootstrap.list <- list()
# Plot the overall PSMC results as black line
a = read.table(paste(SAMPLENAME,".plot.0.txt",sep=""), header = F)
# Plot each bootstrapped PSMC file in gray color. Remember 1 to 100 are boostraps
no.of.bootstraps <- length(list.files(path = ".",pattern=SAMPLENAME ))-1 # remember 0 is not bt
for (i in 1:no.of.bootstraps){
  b = read.table(paste(SAMPLENAME ,".plot.", i, ".txt", sep = ""), header = F)
  # remove <= 10k and >2M
  #  which(a$V1 <= 10000)
  # class(a) # df, so
  b <- b[which(b$V1 >= x.min),]
  b <- b[which(b$V1 <= x.limit),]
  bootstrap.list[[i]] <- min(b$V1)
  points(b$V1 + 0.0001, b$V2, col = "#bcbfdc", type = "l", lwd = 1)
}
# min(a$V1) # here 19,997 -- why appears cut-off
min.b.x <- min(unlist(bootstrap.list))
# a <- a[which(a$V1 >= min.b.x),]
a <- a[which(a$V1 >= x.min),]
a <- a[which(a$V1 <= x.limit),]
points(a$V1 + 0.0001, a$V2, col = "#2d357f", type = "l", lwd = 3) # draw here so on top
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SAMPLENAME <- "Ant2"

remove(bootstrap.list)
bootstrap.list <- list()
# Plot the overall PSMC results as black line
a = read.table(paste(SAMPLENAME,".plot.0.txt",sep=""), header = F)
# Plot each bootstrapped PSMC file in gray color. Remember 1 to 100 are boostraps
no.of.bootstraps <- length(list.files(path = ".",pattern=SAMPLENAME ))-1 # remember 0 is not bt
for (i in 1:no.of.bootstraps){
  b = read.table(paste(SAMPLENAME ,".plot.", i, ".txt", sep = ""), header = F)
  # remove <= 10k and >2M
  #  which(a$V1 <= 10000)
  # class(a) # df, so
  b <- b[which(b$V1 >= x.min),]
  b <- b[which(b$V1 <= x.limit),]
  bootstrap.list[[i]] <- min(b$V1)
  points(b$V1 + 0.0001, b$V2, col = "#d4ad8e", type = "l", lwd = 1)
}
# min(a$V1) # here 19,997 -- why appears cut-off
min.b.x <- min(unlist(bootstrap.list))
# a <- a[which(a$V1 >= min.b.x),]
a <- a[which(a$V1 >= x.min),]
a <- a[which(a$V1 <= x.limit),]
points(a$V1 + 0.0001, a$V2, col = "#e60006", type = "l", lwd = 3) # draw here so on top
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SAMPLENAME <- "AA100A"

remove(bootstrap.list)
bootstrap.list <- list()
# Plot the overall PSMC results as black line
a = read.table(paste(SAMPLENAME,".plot.0.txt",sep=""), header = F)
# Plot each bootstrapped PSMC file in gray color. Remember 1 to 100 are boostraps
no.of.bootstraps <- length(list.files(path = ".",pattern=SAMPLENAME ))-1 # remember 0 is not bt
for (i in 1:no.of.bootstraps){
  b = read.table(paste(SAMPLENAME ,".plot.", i, ".txt", sep = ""), header = F)
  # remove <= 10k and >2M
  #  which(a$V1 <= 10000)
  # class(a) # df, so
  b <- b[which(b$V1 >= x.min),]
  b <- b[which(b$V1 <= x.limit),]
  bootstrap.list[[i]] <- min(b$V1)
  points(b$V1 + 0.0001, b$V2, col = "#88c79e", type = "l", lwd = 1)
}
# min(a$V1) # here 19,997 -- why appears cut-off
min.b.x <- min(unlist(bootstrap.list))
# a <- a[which(a$V1 >= min.b.x),]
a <- a[which(a$V1 >= x.min),]
a <- a[which(a$V1 <= x.limit),]
points(a$V1 + 0.0001, a$V2, col = "#107b35", type = "l", lwd = 3) # draw here so on top
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SAMPLENAME <- "BD-17-5A"

remove(bootstrap.list)
bootstrap.list <- list()
# Plot the overall PSMC results as black line
a = read.table(paste(SAMPLENAME,".plot.0.txt",sep=""), header = F)
# Plot each bootstrapped PSMC file in gray color. Remember 1 to 100 are boostraps
no.of.bootstraps <- length(list.files(path = ".",pattern=SAMPLENAME ))-1 # remember 0 is not bt
for (i in 1:no.of.bootstraps){
  b = read.table(paste(SAMPLENAME ,".plot.", i, ".txt", sep = ""), header = F)
  # remove <= 10k and >2M
  #  which(a$V1 <= 10000)
  # class(a) # df, so
  b <- b[which(b$V1 >= x.min),]
  b <- b[which(b$V1 <= x.limit),]
  bootstrap.list[[i]] <- min(b$V1)
  points(b$V1 + 0.0001, b$V2, col = "#d698ea", type = "l", lwd = 1)
}
# min(a$V1) # here 19,997 -- why appears cut-off
min.b.x <- min(unlist(bootstrap.list))
# a <- a[which(a$V1 >= min.b.x),]
a <- a[which(a$V1 >= x.min),]
a <- a[which(a$V1 <= x.limit),]
points(a$V1 + 0.0001, a$V2, col = "#4f0268", type = "l", lwd = 3) # draw here so on top
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SAMPLENAME <- "46020A"

remove(bootstrap.list)
bootstrap.list <- list()
# Plot the overall PSMC results as black line
a = read.table(paste(SAMPLENAME,".plot.0.txt",sep=""), header = F)
# Plot each bootstrapped PSMC file in gray color. Remember 1 to 100 are boostraps
no.of.bootstraps <- length(list.files(path = ".",pattern=SAMPLENAME ))-1 # remember 0 is not bt
for (i in 1:no.of.bootstraps){
  b = read.table(paste(SAMPLENAME ,".plot.", i, ".txt", sep = ""), header = F)
  # remove <= 10k and >2M
  #  which(a$V1 <= 10000)
  # class(a) # df, so
  b <- b[which(b$V1 >= x.min),]
  b <- b[which(b$V1 <= x.limit),]
  bootstrap.list[[i]] <- min(b$V1)
  points(b$V1 + 0.0001, b$V2, col = "#a2d9d7", type = "l", lwd = 1)
}
# min(a$V1) # here 19,997 -- why appears cut-off
min.b.x <- min(unlist(bootstrap.list))
# a <- a[which(a$V1 >= min.b.x),]
a <- a[which(a$V1 >= x.min),]
a <- a[which(a$V1 <= x.limit),]
points(a$V1 + 0.0001, a$V2, col = "#0e706d", type = "l", lwd = 3) # draw here so on top
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~












# Add lines for each major geomagnetic reversal
# abline(v = 780000, lty = "dashed", col = "black")
# abline(v = 2580000, lty = "dashed", col = "black")
# abline(v = 3580000, lty = "dashed", col = "black")


mtext(side = 2, expression(paste("Effective Population Size (x", "10"^"4",")")), line = 2.5, cex = 1.25)
mtext(side = 1, "Time (100 ka ago)", line = 2.5, cex = 1.25)






# Write plot to PDF file
dev.off()


# how many Antechinus flavipes?
# 80*(10000) = 800,000 # almost a million at its peak
# 5*(10000) = 50,000 # 50k about 10,000 years ago

