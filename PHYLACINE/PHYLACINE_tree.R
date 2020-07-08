# install.packages("pacman", repos="https://cloud.r-project.org")
pacman::p_load(ggplot2,
               dplyr,
               stringr,
               gridExtra,
               viridisLite,
               raster,
               rasterVis,
               rgdal,
               maptools,
               ape,
               ggtree, update = F)

setwd("~/Downloads/PHYLACINE/PHYLACINE_1.2-master/")
setwd("~/Dropbox/seimlab at NNU/antechinus genome project/--MS--/-Figure workings/PHYLACINE/PHYLACINE_1.2-master/")


forest <- read.nexus("Data/Phylogenies/Complete_phylogeny.nex")
names(forest) <- NULL
set.seed(42)
forest <- forest[sample(1:1000, 30)]

# Load world map and subset to Australia
library(maptools)
data(wrld_simpl)
# australia <- wrld_simpl[wrld_simpl$NAME == "Australia", ]
australia <- wrld_simpl[wrld_simpl$NAME == "Australia" | wrld_simpl$NAME == "Papua New Guinea", ]

# The studied species is native to South America. We can subset the data and collect the occurrence points from South America only. The extent is defined as (1) the westernmost longitude, (2) the easternmost longitude, (3) the southernmost latitude, and (4) the northern most latitude.


# Subset data for South America
southamerica <- wrld_simpl[   which(wrld_simpl$LON <=-50 & wrld_simpl$LON <= -30 & wrld_simpl$LAT -60 & wrld_simpl$LAT < 30) , ]

all.marsupials <- wrld_simpl[   which(wrld_simpl$LON <=-50 & wrld_simpl$LON <= -30 & wrld_simpl$LAT -60 | wrld_simpl$LAT < 30 & wrld_simpl$NAME == "Australia" | wrld_simpl$NAME == "Papua New Guinea" )  , ]



# Load trait data. Remember to always use UTF-8 encoding with PHYLACINE files to avoid weird characters 
mam <- read.csv("Data/Traits/Trait_data.csv", fileEncoding = "UTF-8", stringsAsFactors = F)

# Set factor levels for IUCN status. "EP" is a new status we added to designate species that went extinct in prehistory like Diprotodon  
mam$IUCN.Status.1.2 <- factor(mam$IUCN.Status.1.2, levels=c("EP", "EX", "EW", "CR", "EN", "VU", "NT", "LC", "DD"))

# Subset to species that are in Australian marsupial orders, over 20 kg, and over 90 % herbivorous
marsupial.orders <- c("Dasyuromorphia")
# marsupial.orders <- c("Dasyuromorphia", "Peramelemorphia",
# "Notoryctemorphia", "Diprotodontia")

# South American
marsupial.orders <- c("Didelphimorphia")

# both
marsupial.orders <- c("Dasyuromorphia","Didelphimorphia")


marsupials <- mam[mam$Order.1.2 %in% marsupial.orders, ]
#@ marsupials <- marsupials[marsupials$Mass.g > 20000, ]
#@ marsupials <- marsupials[marsupials$Diet.Plant >= 90, ]

# Current maps show where species live today. Present natural maps represent a counterfactual scenario that shows where species would live without Anthropogenic pressures. These ranges include extinct species.
# Load maps of these marsupial species
maps.current <- paste0("Data/Ranges/Current/", marsupials$Binomial.1.2, ".tif")
r.current <- stack(maps.current)
maps.pres.nat <- paste0("Data/Ranges/Present_natural/", marsupials$Binomial.1.2, ".tif")
r.pres.nat <- stack(maps.pres.nat)

# Project Australia map to the range map projection
australia <- spTransform(australia, crs(r.current))
southamerica <- spTransform(southamerica, crs(r.current))
all.marsupials <- spTransform(all.marsupials, crs(r.current))



# Crop range maps to just the extent of Australia for a cleaner plot
ext <- extent(australia)
ext <- extent(southamerica)
ext <- extent(all.marsupials)
#@ ext[2] <- 15000000 # Reduce eastern extent
#@ ext[3] <- -5200000 # Reduce southern extent
#@ ext[4] <- -1370000 # Reduce northern extent

# want PNG
# 10898232, 15391382, -5986986, -140109.6  (xmin, xmax, ymin, ymax)
#@ ext[2] <- 15000000 # Reduce eastern extent
#@ ext[3] <- -5200000 # Reduce southern extent
#@ ext[4] <- -1370000 # Reduce northern extent



r.current <- crop(r.current, ext)
r.pres.nat <- crop(r.pres.nat, ext)

# Create a blank raster of the region
blank.r <- r.current[[1]]
blank.r[] <- NA
names(blank.r) <- NA

# Load all the current raster data into a matrix for faster handling
m.current <- matrix(NA, nrow=nrow(marsupials), ncol=dim(r.current)[1]*dim(r.current)[2])
rownames(m.current) <- marsupials$Binomial.1.2
for(i in 1:nrow(marsupials)) {
  m.current[i, ] <- getValues(r.current[[i]])
}
# Current species list
sp.current <- marsupials$Binomial.1.2[rowSums(m.current) > 0]

# Load all the present natural raster data into a matrix for faster handling
m.pres.nat <- matrix(NA, nrow=nrow(marsupials), ncol=dim(r.current)[1]*dim(r.current)[2])
rownames(m.pres.nat) <- marsupials$Binomial.1.2
for(i in 1:nrow(marsupials)) {
  m.pres.nat[i, ] <- getValues(r.pres.nat[[i]])
}
# Present natural species list
sp.pres.nat <- marsupials$Binomial.1.2[rowSums(m.pres.nat) > 0]

# Create rasters of taxonomic richness
current.div <- blank.r
current.div[] <- colSums(m.current)
names(current.div) <- "Current diversity"
pres.nat.div <- blank.r
pres.nat.div[] <- colSums(m.pres.nat)
names(pres.nat.div) <- "Present natural diversity"
div <- stack(current.div, pres.nat.div)
# Change 0 to NA, for nicer plotting
div[] <- ifelse(div[] == 0, NA, div[])

# Set up the ranges for plotting.
# Plot the two maps on top of each other, with a nice legend
p.map <- 
  levelplot(div,
            layout = c(1,2),
            colorkey=list(
              space='left',                   
              labels=list(at = 1:max(div[], na.rm = T), font=4),
              axis.line=list(col = 'black'),
              width=0.75
            ),
            par.settings = list(
              strip.border = list(col='transparent'),
              strip.background = list(col='transparent'),
              axis.line = list(col = 'transparent')
            ),
            scales = list(draw = FALSE),            
            col.regions = viridis,                   
            at = 1:max(div[], na.rm = T),
            names.attr = str_replace_all(names(div), "\\.", " "))

# Add Australia polygon outline
p.map <- p.map + layer(sp.polygons(australia, col = "darkgrey"))
p.map <- p.map + layer(sp.polygons(southamerica, col = "darkgrey"))
p.map <- p.map + layer(sp.polygons(all.marsupials, col = "darkgrey"))



# Trim the phylogenies down to just the Australian megafauna we are looking at and prepare them for plotting. We plot a random tree as the focal tree and then 29 other trees behind it to show uncertainty.
# Trim tree down to the species list of Australian marsupials created earlier
pruned.forest <- lapply(forest, 
                        function(x) {
                          drop.tip(x, x$tip.label[!x$tip.label %in% sp.pres.nat])
                        }
)
class(pruned.forest) <- "multiPhylo"

# Pick a single tree for overplotting and labelling
tree <- pruned.forest[[1]]
# Group tree by extant species for marking them later
tree <- groupOTU(tree, tree$tip.label %in% sp.current)
# Turn tree to dataframe
tree <- fortify(tree)
# Add trait info
tree <- left_join(tree, marsupials, by = c("label" = "Binomial.1.2"))
# Reverse age scale
tree$x <- tree$x - max(tree$x)

# Prepare mass based legend
mass.breaks <- marsupials$Mass.g/1000
# mass.breaks <- ceiling(c(min(mass.breaks), median(mass.breaks), max(mass.breaks))) # fix...
# mass.breaks <- (c(0.1, median(mass.breaks), max(mass.breaks))) # try for small marsupials

median(mass.breaks) # 0.03415 = 0.03415*1000  = 34g 
mean(mass.breaks) # 0.7475068 kg

# Prepare IUCN status color legend. We use a modified version of IUCN's official color palette
status.colors <- c("EP" = "#87421F", "EX" = "#8F47B3", "EW" = "#8F47B3", 
                   "CR" = "#D81E05", "EN" = "#FC7F3F", "VU" = "#F9E814", 
                   "NT" = "#CCE226", "LC" = "#60C659", "DD" = "#D1D1C6")

# Convert multiPhylo to data.frame
pruned.forest <- fortify(pruned.forest)
# Reverse ages to show time BP
pruned.forest$x <- pruned.forest$x - max(pruned.forest$x)

# Plot forest of uncertainty (based on only 30 out of 1000 trees for speed)
p.tree <- ggplot(pruned.forest) +
  geom_tree(col = "lightblue", alpha = .3, multiPhylo = TRUE) +
  theme_tree2() +
  scale_x_continuous("Time BP (Ma)",
                     limits = c(min(tree$x), 23), breaks = seq(-50, 0, 10),
                     label = abs(seq(-50, 0, 10)))

# Overlay a plot of a single tree with more info
p.tree <- p.tree + 
  geom_tree(data = tree, aes(x, y, lty = group), size = .8) +
  geom_tiplab(data = tree, aes(label = paste0('italic(', str_replace(label, "_", "~"), ')')),
              offset = .8, parse = T, size= 3.5) +
  scale_linetype_manual(values = c("dashed", "solid"), guide = F) +
  geom_tippoint(data = tree, aes(x, y, color = IUCN.Status.1.2, size = Mass.g/1000)) +
  scale_color_manual("IUCN Status",
                     values = status.colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        panel.background = element_blank(),
        axis.line.x = element_line()) +
  scale_size_continuous("Mass (kg)",
                        breaks = mass.breaks,
                        trans = "log10")

# Finally, put everything together into one big plot and save it.
# Plot it all

# Australian
p <- arrangeGrob(
  p.map, p.tree,
  widths = c(1, 2),
  layout_matrix = rbind(c(1, 2)),
  top = "Diversity of Dasyuromorphia marsupials in Australia and Papua New Guinea\n(30 trees randomly drawn from posterior distribution)"
)


# South American
p <- arrangeGrob(
  p.map, p.tree,
  widths = c(1, 2),
  layout_matrix = rbind(c(1, 2)),
  top = "Diversity of Didelphimorphia in South America\n(30 trees randomly drawn from posterior distribution)"
)

# South American
p <- arrangeGrob(
  p.map, p.tree,
  widths = c(1, 2),
  layout_matrix = rbind(c(1, 2)),
  top = "Diversity of Didelphimorphia and Dasyuromorphia\n(30 trees randomly drawn from posterior distribution)"
)




# do not include any traits here -- just show them



# If you just want to print the plot without saving it, you can use grid.arrange(p)
ggsave("Australian_PNG_megafauna.pdf", plot = p, units = "cm", width = 28, height = 40)

ggsave("SouthAmerica_megafauna.pdf", plot = p, units = "cm", width = 28, height = 40)

ggsave("Marsupial_megafauna.pdf", plot = p, units = "cm", width = 28, height = 40)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Overlay a plot of a single tree with more info
# Plot forest of uncertainty (based on only 30 out of 1000 trees for speed)
p.tree <- ggplot(pruned.forest) +
  geom_tree(col = "lightblue", alpha = .3, multiPhylo = TRUE) +
  theme_tree2() +
  scale_x_continuous("Time BP (Ma)",
#                     limits = c(min(tree$x), 23), breaks = seq(-90, 0, 10),
limits = c(-85, 23), breaks = seq(-90, 0, 10),
label = abs(seq(-90, 0, 10)))

table(tree$group)
status.colors

# set size = 1 below so all are the same
p.tree <- p.tree + 
  geom_tree(data = tree, aes(x, y, lty = group), size = .8) +
  geom_tiplab(data = tree, aes(label = paste0('italic(', str_replace(label, "_", "~"), ')')),
              offset = .8, parse = T, size= 3.5) +
  scale_linetype_manual(values = c("dashed", "solid"), guide = F) +
  geom_tippoint(data = tree, aes(x, y, color = IUCN.Status.1.2)) +
  scale_color_manual("IUCN Status",
                     values = status.colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        panel.background = element_blank(),
        axis.line.x = element_line()) 
ggsave("tree-only.pdf", plot = p.tree, units = "cm", width = 30, height = 55)
