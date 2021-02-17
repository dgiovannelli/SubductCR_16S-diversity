##########################################################
## Analysis of the metagenome results from the Biology Meet
## Subduction project leg 1 Samples from Norther and Central
## Costa Rica collected in 2017.
## This script and associated data are available from my github repository
## and from the appropriate database (SRA). For questions please contact:
## Donato Giovannelli - donato.giovannelli@unina.it - dgiovannelli.github.io
##########################################################

### Load required libraries
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(dplyr) # data handling
library(network) # networks
library(intergraph)  # networks
library(ggnet)   # network plotting with ggplot
library(igraph)  # networks
library(phyloseq) # ASV ecological analysis package
library(ggplot2) # plotting library
library(gridExtra) # gridding plots
library(ape) # importing and handling phylogenetic trees
library(ggthemes) # additional themes fro ggplot2
library(magrittr) #
library(rioja) # plotting poackages for tabular bubbleplots
library(ggpubr)
library(ggtern) # ternary plots for geochemistry
library(plyr)
library(coda.base)
library(tydiverse)
library(vegan) # Multivariate ecological analysis
library(propr)
library(missForest) # Imputing missing values in dataframes using Random Forests
library(VSURF) # Random Forests approach to variable importance identification
library(car) #for scatterplot

# Importing all the mi-faser csv files as separate objects
temp = list.files(pattern="*.csv")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))

mifaser<-ec_list.csv

test <- merge(x=ec_list.csv,y=BQF_mifaser.csv[,c("ec","readcount")], by="ec",  all.x=TRUE)

temp<-list(mifaser,ARS_mifaser.csv, BQF_mifaser.csv, BQS_mifaser.csv, BQS1_mifaser.csv,BR1F_mifaser.csv,BRF1_mifaser.csv,BRF2_mifaser.csv,BRS1_mifaser.csv,BRS2_mifaser.csv,CYF_mifaser.csv,CYS_mifaser.csv,EPF_mifaser.csv,EPS_mifaser.csv,ESF9_mifaser.csv,ETS_mifaser.csv,FAS_mifaser.csv,mifaser,MTF_mifaser.csv,PBS_mifaser.csv,PFF_mifaser.csv,PFS_mifaser.csv,PGF_mifaser.csv,PGS_mifaser.csv,PLS_mifaser.csv,QH2F_mifaser.csv,QHS1_mifaser.csv,QHS2_mifaser.csv,QNF_mifaser.csv,QNS_mifaser.csv,RSF_mifaser.csv,RSS_mifaser.csv,RVF_mifaser.csv,SIF_mifaser.csv,SIS_mifaser.csv,SLF_mifaser.csv,SLS_mifaser.csv,TCF_mifaser.csv,TCS_mifaser.csv)

dataset <- join_all(temp, "ec")



# Selecting only the readcount data
dataset <- cbind (dataset[,1], dataset[ , grepl( "readcount" , names( dataset ) ) ])

colnames(dataset)<-c('EC', 'ARS', 'BQF', 'BQS', 'BQS1', 'BR1F', 'BRF1', 'BRF2', 'BRS1', 'BRS2', 'CYF', 'CYS', 'EPF', 'EPS', 'ESF9', 'ETS', 'FAS', 'MTF', 'PBS', 'PFF', 'PFS', 'PGF', 'PGS', 'PLS', 'QH2F', 'QHS1', 'QHS2', 'QNF', 'QNS', 'RSF', 'RSS', 'RVF', 'SIF', 'SIS', 'SLF', 'SLS', 'TCF', 'TCS')

colnames(dataset)

dataset_clean<-subset (dataset, select=-c(ARS,PFF,PFS,PGS,PGF,PLS))

dataset_clean

ecdata<-as.matrix(dataset_clean[,-1])

row.names(ecdata)<-dataset_clean[,1]

row.names(ecdata)

ecdata <- phyloseq(otu_table(ecdata, taxa_are_rows = T))

ecdata #inspect the object

# Normalize counts to relative abundance first and then multiply for the median library abundance
ecdata_n <- transform_sample_counts(ecdata, function(x) ((x / sum(x, na.rm=T))*median(colSums(ecdata, na.rm=T))))

median(colSums(ecdata, na.rm=T))

ecdata_n

# extract only the relevanc EC number for carbon metabolism
ec_list<-read.csv("ec_list_carbon.csv", header=F)

ec_list<-as.matrix(ec_list)

rownames(otu_table(ecdata_n))

dataset_carbon <- subset(otu_table(ecdata_n), rownames(otu_table(ecdata_n)) %in% c(ec_list))

dataset_carbon


#Building the co-occurrence network
bac.cor <- cor(t(dataset_carbon),  use="complete.obs", method="spearman")   # use="complete.obs" because we have some missing data
bac.cor[bac.cor < 0.5] = 0
bac.cor.ig <- graph.adjacency(bac.cor, mode='undirected', add.rownames = TRUE, weighted = TRUE)
bac.cor.ig <- igraph::simplify(bac.cor.ig)
plot(bac.cor.ig, layout=layout_with_fr, vertex.size = 8, rescale=T)

# Network statistics
ecount(bac.cor.ig)
vcount(bac.cor.ig)
assortativity_degree(bac.cor.ig, directed=F)

# Identify cliques using Louvain
bac.cor.louvain <- cluster_louvain(bac.cor.ig)
modularity(bac.cor.louvain)
V(bac.cor.ig)$color=bac.cor.louvain$membership

plot(bac.cor.ig, layout=layout_with_fr, col = bac.cor.louvain, vertex.size = 8, rescale=T)

# Identify cliques using Greedy clustering
bac.cor.greedy <- cluster_fast_greedy(bac.cor.ig)
modularity(bac.cor.greedy)
V(bac.cor.ig)$color=bac.cor.greedy$membership
plot(bac.cor.ig, layout=layout_with_fr, col = bac.cor.greedy, vertex.size = 8, rescale=T)

svg("carbon_netwrok_cliques.svg", height = 4, width=4)
plot(bac.cor.ig, layout=layout_with_fr, col = bac.cor.greedy, vertex.size = 8, rescale=T)
dev.off()

# Identify cliques using Walk trap
bac.cor.walk <- cluster_walktrap(bac.cor.ig)
modularity(bac.cor.walk)
V(bac.cor.ig)$color=bac.cor.walk$membership
plot(bac.cor.ig, layout=layout_with_fr, col = bac.cor.walk, vertex.size = 8, rescale=T)

sizes(bac.cor.louvain)
sizes(bac.cor.greedy)
sizes(bac.cor.walk)

### Extracting the list of ASV for the main cliques identified by the Louvain function
V(bac.cor.ig)[bac.cor.louvain$membership == 1] -> bac.cor.group1
V(bac.cor.ig)[bac.cor.louvain$membership == 2] -> bac.cor.group2
V(bac.cor.ig)[bac.cor.louvain$membership == 3] -> bac.cor.group3
V(bac.cor.ig)[bac.cor.louvain$membership == 4] -> bac.cor.group4
V(bac.cor.ig)[bac.cor.louvain$membership == 5] -> bac.cor.group5

bac.cor.group1
bac.cor.group2
bac.cor.group3
bac.cor.group4
bac.cor.group5

#subset cliques
clique1 <- subset(dataset_carbon, rownames(otu_table(dataset_carbon)) %in% c(names(bac.cor.group1)))
clique2 <- subset(dataset_carbon, rownames(otu_table(dataset_carbon)) %in% c(names(bac.cor.group2)))
clique3 <- subset(dataset_carbon, rownames(otu_table(dataset_carbon)) %in% c(names(bac.cor.group3)))
clique4 <- subset(dataset_carbon, rownames(otu_table(dataset_carbon)) %in% c(names(bac.cor.group4)))
clique5 <- subset(dataset_carbon, rownames(otu_table(dataset_carbon)) %in% c(names(bac.cor.group5)))


envdata<-read.csv("subductCR_final_dataset.csv", header=T)

envdata[,1]

colnames(clique1)

colSums(clique4, na.rm=T)

plot(envdata$dic, log(colSums(clique4, na.rm=T)))
cor.test(envdata$dic, colSums(clique4, na.rm=T), method = "spearman")

## Test each clique abundance against environmental predictors using Person, Spearman and scatterplots
#Clique 1
par(mfrow=c(3,3))
for (i in 1:length(envdata)) {
    plot(envdata[,i], log(colSums(otu_table(clique1), na.rm=T)), xlab=colnames(envdata)[i])
}

message("Test with Pearson correlation:")
# Pearson
for (i in 16:length(envdata)) {
    a <- cor.test(envdata[,i], colSums(otu_table(clique1), na.rm=T))
       if (a$p.value<0.05) {
           print(paste(i,colnames(envdata)[i],a$estimate, a$parameter, a$p.value))
       }
}
message("Test with Spearman correlation:")
# Spearman
for (i in 16:length(envdata)) {
    a <- suppressWarnings(cor.test(envdata[,i], colSums(otu_table(clique1), na.rm=T), method = "spearman"))
       if (a$p.value<0.05) {
           print(paste(i,colnames(envdata)[i],a$estimate, a$parameter, a$p.value))
       }
}

#Clique 2
par(mfrow=c(3,3))
for (i in 1:length(envdata)) {
    plot(envdata[,i], log(colSums(otu_table(clique2), na.rm=T)), xlab=colnames(envdata)[i])
}

message("Test with Pearson correlation:")
# Pearson
for (i in 16:length(envdata)) {
    a <- cor.test(envdata[,i], colSums(otu_table(clique2), na.rm=T))
       if (a$p.value<0.05) {
           print(paste(i,colnames(envdata)[i],a$estimate, a$parameter, a$p.value))
       }
}
message("Test with Spearman correlation:")
# Spearman
for (i in 16:length(envdata)) {
    a <- suppressWarnings(cor.test(envdata[,i], colSums(otu_table(clique2), na.rm=T), method = "spearman"))
       if (a$p.value<0.05) {
           print(paste(i,colnames(envdata)[i],a$estimate, a$parameter, a$p.value))
       }
}

#Clique 4
par(mfrow=c(3,3))
for (i in 1:length(envdata)) {
    plot(envdata[,i], log(colSums(otu_table(clique4), na.rm=T)), xlab=colnames(envdata)[i])
}

message("Test with Pearson correlation:")
# Pearson
for (i in 16:length(envdata)) {
    a <- cor.test(envdata[,i], colSums(otu_table(clique4), na.rm=T))
       if (a$p.value<0.05) {
           print(paste(i,colnames(envdata)[i],a$estimate, a$parameter, a$p.value))
       }
}
message("Test with Spearman correlation:")
# Spearman
for (i in 16:length(envdata)) {
    a <- suppressWarnings(cor.test(envdata[,i], colSums(otu_table(clique4), na.rm=T), method = "spearman"))
       if (a$p.value<0.05) {
           print(paste(i,colnames(envdata)[i],a$estimate, a$parameter, a$p.value))
       }
}

#Clique 5
par(mfrow=c(3,3))
for (i in 1:length(envdata)) {
    plot(envdata[,i], log(colSums(otu_table(clique5), na.rm=T)), xlab=colnames(envdata)[i])
}

message("Test with Pearson correlation:")
# Pearson
for (i in 16:length(envdata)) {
    a <- cor.test(envdata[,i], colSums(otu_table(clique5), na.rm=T))
       if (a$p.value<0.05) {
           print(paste(i,colnames(envdata)[i],a$estimate, a$parameter, a$p.value))
       }
}
message("Test with Spearman correlation:")
# Spearman
for (i in 16:length(envdata)) {
    a <- suppressWarnings(cor.test(envdata[,i], colSums(otu_table(clique5), na.rm=T), method = "spearman"))
       if (a$p.value<0.05) {
           print(paste(i,colnames(envdata)[i],a$estimate, a$parameter, a$p.value))
       }
}

sizes(bac.cor.greedy)

### Extracting the list of ASV for the main cliques identified by the Louvain function
V(bac.cor.ig)[bac.cor.greedy$membership == 1] -> bac.cor.gr.group1
V(bac.cor.ig)[bac.cor.greedy$membership == 2] -> bac.cor.gr.group2
V(bac.cor.ig)[bac.cor.greedy$membership == 3] -> bac.cor.gr.group3

bac.cor.gr.group1
bac.cor.gr.group2
bac.cor.gr.group3

#subset cliques
clique.gr1 <- subset(dataset_carbon, rownames(otu_table(dataset_carbon)) %in% c(names(bac.cor.gr.group1)))
clique.gr2 <- subset(dataset_carbon, rownames(otu_table(dataset_carbon)) %in% c(names(bac.cor.gr.group2)))
clique.gr3 <- subset(dataset_carbon, rownames(otu_table(dataset_carbon)) %in% c(names(bac.cor.gr.group3)))

## Test each clique abundance against environmental predictors using Person, Spearman and scatterplots
#Clique 1
par(mfrow=c(3,3))
for (i in 1:length(envdata)) {
    plot(envdata[,i], log(colSums(otu_table(clique.gr1), na.rm=T)), xlab=colnames(envdata)[i])
}

message("Test with Pearson correlation:")
# Pearson
for (i in 16:length(envdata)) {
    a <- cor.test(envdata[,i], colSums(otu_table(clique.gr1), na.rm=T))
       if (a$p.value<0.05) {
           print(paste(i,colnames(envdata)[i],a$estimate, a$parameter, a$p.value))
       }
}
message("Test with Spearman correlation:")
# Spearman
for (i in 16:length(envdata)) {
    a <- suppressWarnings(cor.test(envdata[,i], colSums(otu_table(clique.gr1), na.rm=T), method = "spearman"))
       if (a$p.value<0.05) {
           print(paste(i,colnames(envdata)[i],a$estimate, a$parameter, a$p.value))
       }
}

## Test each clique abundance against environmental predictors using Person, Spearman and scatterplots
#Clique 2
par(mfrow=c(3,3))
for (i in 1:length(envdata)) {
    plot(envdata[,i], log(colSums(otu_table(clique.gr2), na.rm=T)), xlab=colnames(envdata)[i])
}

message("Test with Pearson correlation:")
# Pearson
for (i in 16:length(envdata)) {
    a <- cor.test(envdata[,i], colSums(otu_table(clique.gr2), na.rm=T))
       if (a$p.value<0.05) {
           print(paste(i,colnames(envdata)[i],a$estimate, a$parameter, a$p.value))
       }
}
message("Test with Spearman correlation:")
# Spearman
for (i in 16:length(envdata)) {
    a <- suppressWarnings(cor.test(envdata[,i], colSums(otu_table(clique.gr2), na.rm=T), method = "spearman"))
       if (a$p.value<0.05) {
           print(paste(i,colnames(envdata)[i],a$estimate, a$parameter, a$p.value))
       }
}

## Test each clique abundance against environmental predictors using Person, Spearman and scatterplots
#Clique 1
par(mfrow=c(3,3))
for (i in 1:length(envdata)) {
    plot(envdata[,i], log(colSums(otu_table(clique.gr3), na.rm=T)), xlab=colnames(envdata)[i])
}

message("Test with Pearson correlation:")
# Pearson
for (i in 16:length(envdata)) {
    a <- cor.test(envdata[,i], colSums(otu_table(clique.gr3), na.rm=T))
       if (a$p.value<0.05) {
           print(paste(i,colnames(envdata)[i],a$estimate, a$parameter, a$p.value))
       }
}
message("Test with Spearman correlation:")
# Spearman
for (i in 16:length(envdata)) {
    a <- suppressWarnings(cor.test(envdata[,i], colSums(otu_table(clique.gr3), na.rm=T), method = "spearman"))
       if (a$p.value<0.05) {
           print(paste(i,colnames(envdata)[i],a$estimate, a$parameter, a$p.value))
       }
}


plot_data<-data.frame("clique2"=colSums(otu_table(clique.gr2), na.rm=T), "trench"=envdata$trench, "temp"=envdata$temp, "ph"=envdata$ph, "dic"=envdata$dic, "dic_d13"=envdata$dic_d13, "doc"=envdata$doc, "doc_d13"=envdata$doc_d13)

plot_data

cor.test(plot_data$temp, plot_data$clique2)
cor.test(plot_data$temp, plot_data$clique2, method = "spearman")

cor.test(plot_data$ph, plot_data$clique2)
cor.test(plot_data$ph, plot_data$clique2, method = "spearman")

cor.test(plot_data$dic, plot_data$clique2)
cor.test(plot_data$dic, plot_data$clique2, method = "spearman")

cor.test(plot_data$dic_d13, plot_data$clique2)
cor.test(plot_data$dic_d13, plot_data$clique2, method = "spearman")

cor.test(plot_data$doc, plot_data$clique2)
cor.test(plot_data$doc, plot_data$clique2, method = "spearman")

cor.test(plot_data$doc_d13, plot_data$clique2)
cor.test(plot_data$doc_d13, plot_data$clique2, method = "spearman")

names(envdata)

for (i in 11:length(envdata)) {
    a <- suppressWarnings(cor.test(envdata[,i], colSums(otu_table(clique.gr2), na.rm=T), method = "spearman"))
           print(paste(i,colnames(envdata)[i],a$estimate, a$parameter, a$p.value))
}

ggplot(plot_data, aes(temp, log(clique2))) +
geom_jitter(aes(color=trench), size=3) +
geom_smooth(method = lm, se = F) +
theme_bw()

ggplot(plot_data, aes(ph, log(clique2))) +
geom_jitter(aes(color=trench), size=3) +
geom_smooth(method = lm, se = F) +
theme_bw()

ggplot(plot_data, aes(log(dic), log(clique2))) +
geom_jitter(aes(color=trench), size=3) +
geom_smooth(method = lm, formula= y~x, se = F) +
ylim(NA, 11) +
theme_bw()

ggsave("clique2_dic.svg", height = 4, width=5)

ggplot(plot_data, aes(dic_d13, log(clique2))) +
geom_jitter(aes(color=trench), size=3) +
geom_smooth(method = lm, formula= y~x, se = F) +
ylim(NA, 11) +
theme_bw()

ggsave("clique2_dic_d13.svg", height = 4, width=5)

ggplot(plot_data, aes(log(doc), log(clique2))) +
geom_jitter(aes(color=trench), size=3) +
geom_smooth(method = lm, formula= y~x, se = F) +
ylim(NA, 11) +
theme_bw()

ggsave("clique2_doc.svg", height = 4, width=5)

ggplot(plot_data, aes(doc_d13, log(clique2))) +
geom_jitter(aes(color=trench), size=3) +
geom_smooth(method = lm, formula= y~x, se = F) +
ylim(NA, 11) +
theme_bw()

ggsave("clique2_doc_d13.svg", height = 4, width=5)


# Extract environmental predictor and remove non numeric variables (i.e. factors) and unwanted variables (air temperature)
envdata_rf <- data.frame(envdata[,10:74])
names(envdata_rf)

# Contruct the vetors for each clique abundance response variable and normalize them as Z-scores
cliqueA <- scale(data.frame(colSums(otu_table(clique.gr1), na.rm = T)), center = TRUE, scale = TRUE)
cliqueB <- scale(data.frame(colSums(otu_table(clique.gr2), na.rm = T)), center = TRUE, scale = TRUE)
cliqueC <- scale(data.frame(colSums(otu_table(clique.gr3), na.rm = T)), center = TRUE, scale = TRUE)

#Check if the predictors contain NAs
if (grep("NA", envdata_rf)!=0){
  print("There are NAs! Need to impute variables")
} else {
  print("No NAs in the predictors, good to go directly to the Random Forests")
} # The predictors dataset c ontains NAs

# Impute missing data using a iterative random forest approach with the missForest package
envdata_rf2 <- missForest(envdata_rf, maxiter=10, variablewise = TRUE, verbose = TRUE, ntree=5000)
envdata_rf2$OOBerror
envdata_rf2$ximp # dataframe with imputed missing values to be used in downstrewam analysis

# Run the the Random Forests model for variable selection using VSURF
cliqueA_RF <- VSURF(x = data.frame(envdata_rf2$ximp), y = cliqueA[,1], parallel=T, ntree=2000, mtry=23)
summary(cliqueA_RF)
plot(cliqueA_RF)
colnames(envdata_rf2$ximp)[cliqueA_RF$varselect.thres]
colnames(envdata_rf2$ximp)[cliqueA_RF$varselect.interp]

# Repeat for each cliques
cliqueB_RF <- VSURF(x = data.frame(envdata_rf2$ximp), y = cliqueB[,1], parallel=T, ntree=5000, mtry=50)
plot(cliqueB_RF)
colnames(envdata_rf2$ximp)[cliqueB_RF$varselect.thres]
colnames(envdata_rf2$ximp)[cliqueB_RF$varselect.interp]

cliqueC_RF <- VSURF(x = data.frame(envdata_rf2$ximp), y = cliqueC[,1], parallel=T, ntree=2000, mtry=23)
plot(cliqueC_RF)
colnames(envdata_rf2$ximp)[cliqueC_RF$varselect.thres]
colnames(envdata_rf2$ximp)[cliqueC_RF$varselect.interp]

save.image() # Save all this good work!
