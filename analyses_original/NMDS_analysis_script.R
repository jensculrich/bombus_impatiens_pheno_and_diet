#---------------------------------------------
#  Small set up
#---------------------------------------------

rm(list=ls()) # clears the workspace
setwd("C:/Users/mplat/Documents/thesis")

#---------------------------------------------
# Setting up for NMDS stats
#---------------------------------------------

library(picante)

#read .csv file
comm <- read.csv("bom_plant_interactions_workers only.csv", header = TRUE, row.names = 1)

#check the class of my data
class(comm)

#check the dimensions of my data
dim(comm)

# If data is read in correctly, each row name will be a plant species
rownames(comm)

#If data is read in correctly, each column name will be a bumble bee species
head(colnames(comm))

#just produce a head for comm
head(comm)

#read in metadata so points can be coloured appropriately
metadata <- read.csv("bb_metadata.csv", header = TRUE, row.names = 1)
head(metadata)

ls()

#Check that metadata is in the correct format
all.equal(rownames(comm), rownames(metadata))
metadata <- metadata[rownames(comm), ]
print(metadata)

#--------------------------------------------
# Producing the NMDS graph with plant and bumble bee labels

library(ggplot2)

#create nmds coordinates with the data
comm.nmds <- metaMDS(comm)

# Extract NMDS coordinates
nmds_coords <- as.data.frame(scores(comm.nmds, display = "sites"))
print(nmds_coords)

#plot just bee_site
ordiplot(comm.nmds, display = "species", type = 'text')

#plot nmds data into a graph with labels for both
ordipointlabel(comm.nmds)

ordiplot(comm.nmds, display = "sites", type = 'text')

#create a blank ordination plot
mds.fig <- ordiplot(comm.nmds, type = "none")

# Specify the X and Y axis ranges
ordiplot(comm.nmds, type = "none", xlim = c(-1.5, 1.5), ylim = c(-1.75, 2))


# plot just the samples, colour by bee species, pch=19 means plot a circle
points(mds.fig, "sites", pch = 19, col = "#fb9a99", select = metadata == "bom_imp")
points(mds.fig, "sites", pch = 19, col = "#a6cee3", select = metadata == "bom_fla")
points(mds.fig, "sites", pch = 19, col = "#1f78b4", select = metadata == "bom_mel")
points(mds.fig, "sites", pch = 19, col = "#b2df8a", select = metadata == "bom_mix")
points(mds.fig, "sites", pch = 19, col = "#33a02c", select = metadata == "bom_vos")


# Define colors for each species
species_colors <- c("bom_fla" = "#a6cee3","bom_imp" = "#fb9a99", "bom_mel" = "#1f78b4", "bom_mix" = "#b2df8a", "bom_vos" = "#33a02c")

# plot just the samples, colour by bee species, pch=19 means plot a circle
points(mds.fig, "sites", pch = 19, col = species_colors[metadata])

# Add confidence ellipses around specific bee species
ordiellipse(comm.nmds, groups = metadata, conf = 0.95, label = FALSE, col = species_colors)

# Define the mapping between original species names and desired labels
label_mapping <- c("bom_imp" = "B. impatiens", "bom_fla" = "B. flavifrons", "bom_mel" = "B. melanopygus", "bom_mix" = "B. mixtus", "bom_vos" = "B. vosnesenskii")

# Replace original species names with desired labels
metadata_labels <- label_mapping[metadata]


# Create a legend with the new labels
legend("topright", legend = unique(metadata_labels), fill = unique(species_colors))




#TO DO
perm_anova <- adonis2(as.matrix(comm.nmds) ~ factor(metadata_labels), permutations = 999)
