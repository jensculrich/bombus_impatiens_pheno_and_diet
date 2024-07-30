#This code follows Steve Kembel's online guide for biodiversity analyses
#https://kembellab.ca/r-workshop/biodivR/SK_Biodiversity_R.html


#install.packages("vegan")
#install.packages("picante")

library(vegan)
library(picante)


comm=read.csv("./nmds/data/Melissa_commun_matrix_no0.csv",header=T, row.names=1)
dim(comm)
rownames(comm)
head(colnames(comm))
comm[1:5,1:5]
apply(comm,1,sum)
comm=decostand(comm,method="total")
apply(comm,1,sum)

comm[1:15,1:15]

metadata=read.csv("./nmds/data/Community_meta_data.csv", header=T,row.names=1)

#note, R spits out error for rows that are all zero but it seems to still be working
#Scratch that, it's not working - delete empty rows? 
#(Jul10, new file w no 0 rows uploaded)
comm.bc.mds <- metaMDS(comm, dist = "bray")
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
points(mds.fig,"sites",pch=19,col="lightblue", select=metadata$Species = "B_ins")

# add confidence ellipses - need to pare down to just the relevant sp.
#use a vector for this?
ordiellipse(comm.bc.mds, metadata$Species, conf = 0.95, label = F)



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

