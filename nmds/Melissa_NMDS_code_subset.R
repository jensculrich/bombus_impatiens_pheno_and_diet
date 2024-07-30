#This code follows Steve Kembel's online guide for biodiversity analyses
#https://kembellab.ca/r-workshop/biodivR/SK_Biodiversity_R.html


#install.packages("vegan")
#install.packages("picante")
#install.packages("dplyr")

library(vegan)
library(picante)
library(dplyr)


comm=read.csv("./nmds/data/Melissa_commun_matrix_no0.csv",header=T, row.names=1)
mini_comm_test= comm[which(rownames(comm)=="th_bom_imp"),]
mini_comm_imp=comm[grep("*bom_imp$",rownames(comm)),]
mini_comm_fla=comm[grep("*bom_fla$",rownames(comm)),]
mini_comm_mix=comm[grep("*bom_mix$",rownames(comm)),]
mini_comm_mel=comm[grep("*bom_mel$",rownames(comm)),]
mini_comm_vos=comm[grep("*bom_vos$",rownames(comm)),]
mini_comm=bind_rows(mini_comm_imp,mini_comm_fla,mini_comm_mel,mini_comm_mix,mini_comm_vos)
rownames(mini_comm)

#mini comm analysis (full below)
dim(mini_comm)
rownames(mini_comm)
head(colnames(mini_comm))
mini_comm[1:5,1:5]
apply(mini_comm,1,sum)
mini_comm=decostand(mini_comm,method="total")
apply(mini_comm,1,sum)

mini_comm[1:15,1:15]

metadata_mini=read.csv("./nmds/data/Community_meta_data_mini.csv", header=T,row.names=1)
###Needed a mini metadata with irrelevant rows removed.. 



mini.comm.bc.mds <- metaMDS(mini_comm, dist = "bray")
mini.comm.bc.mds$stress
#spits out a stress value, lower than 0.1 = great, less than 0.2=good, higher than 0.3
#indicates a lack of 2 dimensional arrangement and should be treated with suspect
stressplot(mini.comm.bc.mds)

##Look at Harmon paper for better ideas - they connect rather than draw elipses, then use MANOVA...

ordiplot(mini.comm.bc.mds,display="sites",type="text", xlim=c(-2,2),ylim=c(-2,2))

# ordination plots are highly customizable set up the plotting area but
# don't plot anything yet
mini.mds.fig <- ordiplot(mini.comm.bc.mds, type = "none",xlim=c(-2,2),ylim=c(-2,2))
# plot just the samples, colour by Species, pch=19 means plot a circle
points(mini.mds.fig, "sites", pch = 19, col = "green", 
       select = metadata_mini$Species == "B_imp")
points(mini.mds.fig, "sites", pch = 19, col = "blue", 
       select = metadata_mini$Species == "B_mel")
points(mini.mds.fig,"sites",pch=19, col="brown", 
       select = metadata_mini$Species == "B_vos")
points(mini.mds.fig,"sites",pch=19, col="pink", 
       select = metadata_mini$Species == "B_mix")
points(mini.mds.fig,"sites",pch=19, col="orange", 
       select = metadata_mini$Species == "B_fla")

## redo with shared colour palette for other figures 
library(RColorBrewer)
my_palette <- hcl.colors(5, "Inferno")

mini.mds.fig <- ordiplot(mini.comm.bc.mds, type = "none",xlim=c(-2,2),ylim=c(-2,2))
# plot just the samples, colour by Species, pch=19 means plot a circle
points(mini.mds.fig, "sites", pch = 19, col = my_palette[1], 
       select = metadata_mini$Species == "B_imp")
points(mini.mds.fig, "sites", pch = 19, col = my_palette[5], 
       select = metadata_mini$Species == "B_mel")
points(mini.mds.fig,"sites",pch=19, col= my_palette[4], 
       select = metadata_mini$Species == "B_vos")
points(mini.mds.fig,"sites",pch=19, col= my_palette[2], 
       select = metadata_mini$Species == "B_mix")
points(mini.mds.fig,"sites",pch=19, col= my_palette[3], 
       select = metadata_mini$Species == "B_fla")

# add confidence ellipses - need to pare down to just the relevant sp.
#use a vector for this?

# this colours all of the ellipses the same colour
# ordiellipse(mini.comm.bc.mds, metadata_mini$Species, conf = 0.95, label = F)

##you can play with frame size above

with(mini_comm, ordiellipse(mini.comm.bc.mds, metadata_mini$Species, conf=0.95,
                           lwd=2, col=my_palette[1], show.groups = "B_imp"))
with(mini_comm, ordiellipse(mini.comm.bc.mds, metadata_mini$Species, conf=0.95,
                            lwd=2, col=my_palette[5], show.groups = "B_mel"))
with(mini_comm, ordiellipse(mini.comm.bc.mds, metadata_mini$Species, conf=0.95,
                            lwd=2, col=my_palette[4], show.groups = "B_vos"))
with(mini_comm, ordiellipse(mini.comm.bc.mds, metadata_mini$Species, conf=0.95,
                            lwd=2, col=my_palette[2], show.groups = "B_mix"))
with(mini_comm, ordiellipse(mini.comm.bc.mds, metadata_mini$Species, conf=0.95,
                            lwd=2, col=my_palette[3], show.groups = "B_fla"))

## ggplot version
library(ggordiplots)

test <- metadata_mini
  
test$Species <- factor(test$Species, 
    levels = c("B_imp", "B_mix", "B_fla", "B_vos", "B_mel"))
labs <- c("B. impatiens", "B. mixtus", "B. flavifrons", "B. vosnesenskii", "B. melanopygus")

p <- gg_ordiplot(mini.comm.bc.mds, 
            groups = test$Species, pt.size = 3, size = 3, 
            kind = "sd", conf = 0.95)
plot <- p$plot
plot + theme_light() +
  scale_y_continuous(limits = c(-2, 2)) +
  scale_color_manual(name="Species", 
                     values=my_palette, 
                     labels = labs) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size = 16),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")
        )


# MANOVA to test for diffs among groups
mini.comm.bc.dist=vegdist(mini_comm, method="bray")
mini.comm.bc.clust <- hclust(mini.comm.bc.dist, method = "average")

# plot cluster diagram
plot(mini.comm.bc.clust, ylab = "Bray-Curtis dissimilarity")

#plot shows which site/bee communities have most similar
#visitation patterns
#Mantel below shows if the similarities are explained by
#the species of bee

adonis2(mini.comm.bc.dist ~ metadata_mini$Species, method = "bray")



#ORIGINAL CODE FOR ALL SPECIES
#this doesn't work because it assumes that is the full row name... need a function like grep (?
#that can search by partial row names...)
dim(comm)
rownames(comm)
head(colnames(comm))
comm[1:5,1:5]
apply(comm,1,sum)
comm=decostand(comm,method="total")
apply(comm,1,sum)

comm[1:15,1:15]

metadata=read.csv("Community_meta_data.csv", header=T,row.names=1)

#note, R spits out error for rows that are all zero but it seems to still be working
#Scratch that, it's not working - delete empty rows? 
#(Jul10, new file w no 0 rows uploaded)


comm.bc.mds <- metaMDS(comm, dist = "bray")
comm.bc.mds$stress
#spits out a stress value, lower than 0.1 = great, less than 0.2=good, higher than 0.3
#indicates a lack of 2 dimensional arrangement and should be treated with suspect
stressplot(comm.bc.mds)

##Look at Harmon paper for better ideas - they connect rather than draw elipses, then use MANOVA...

ordiplot(comm.bc.mds,display="sites",type="text")

# ordination plots are highly customizable set up the plotting area but
# don't plot anything yet
mds.fig <- ordiplot(comm.bc.mds, type = "none")
# plot just the samples, colour by Species, pch=19 means plot a circle
points(mds.fig, "sites", pch = 19, col = "green", select = metadata$Species == 
         "B_imp")
points(mds.fig, "sites", pch = 19, col = "blue", select = metadata$Species == 
         "B_mel")
points(mds.fig,"sites",pch=19, col="brown", select = metadata$Species == "B_vos")
points(mds.fig,"sites",pch=19, col="pink", select = metadata$Species == "B_mix")
points(mds.fig,"sites",pch=19, col="orange", select = metadata$Species == "B_fla")
points(mds.fig,"sites",pch=19, col="yellow", select = metadata$Species == "B_sit")
points(mds.fig,"sites",pch=19, col="grey", select = metadata$Species == "B_cal")
points(mds.fig,"sites",pch=19, col="red", select = metadata$Species == "B_occ")
points(mds.fig,"sites",pch=19, col="purple", select = metadata$Species == "B_ruf")
points(mds.fig,"sites",pch=19, col="black", select = metadata$Species == "B_flv")
#Is this a real species or a mistake?
points(mds.fig,"sites",pch=19, col="darkgreen", select = metadata$Species == "B_hun")
points(mds.fig,"sites",pch=19,col="lightblue", select=metadata$Species == "B_ins")

# add confidence ellipses - need to pare down to just the relevant sp.
#use a vector for this?

ordiellipse(comm.bc.mds, select = metadata$Species =="B_vos", conf = 0.95, label = F)

# overlay the cluster results we calculated earlier
#ordicluster(comm.bc.mds, comm.bc.clust, col = "gray")

# MANOVA to test for diffs among groups
comm.bc.dist=vegdist(comm, method="bray")
comm.bc.clust <- hclust(comm.bc.dist, method = "average")

# plot cluster diagram
plot(comm.bc.clust, ylab = "Bray-Curtis dissimilarity")

#plot shows which site/bee communities have most similar
#visitation patterns
#Mantel below shows if the similarities are explained by
#the species of bee

adonis2(comm.bc.dist ~ metadata$Species, method = "bray")
