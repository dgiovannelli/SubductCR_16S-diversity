##########################################################
## Analysis of the 16S rRNA diversity results from the Biology Meet Subduction project leg 1
## Samples from Norther and Central Costa Rica collected in 2017
## Donato  Giovannelli - donato.giovannelli@unina.it - dgiovannelli.github.io
## This script and associated data and materials are available from my github repository
## and from the approapriate database (SRA)
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
library(vegan) # Multivariate ecological analysis
library(propr)
library(missForest) # Imputing missing values in dataframes using Random Forests
library(VSURF) # Random Forests approach to variable importance identification


## Few Functions for data conversion to be used later
# convert the sample_data() within a phyloseq object to a vegan compatible data object
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


theme_set(theme_bw()) #Set Global theme for ggplot2




####################################################################################
## Importing the diversity files obtained from the ASVs analysis into a phyloseq object
## Count data have already been normalyzed applying the following function
## transform_sample_counts(bac_data, function(x) ((x / sum(x))*median(readcounts(bac_data))))

# Bacteria
bac_count <- as.matrix(read.csv("../16S/bac_normalized_count.csv", header=T, sep=",", row.names=1))
bac_tax <- as.matrix(read.csv("../16S/bac_taxonomy.csv", header=T, sep=",", row.names=1))
bac_tree <- read_tree("../16S/bac_tree.tre")
bac_sample <- read.csv("../16S/bac_sample_table.csv", header=T, sep=",", row.names=1)
colnames(bac_count) # Check the sample names
length(colnames(bac_count)) # Check the number of samples

bac_data <- phyloseq(otu_table(bac_count, taxa_are_rows = TRUE), phy_tree(bac_tree), tax_table(bac_tax), sample_data(bac_sample))
bac_data # Inspect the object to get stat on number of taxa and samples
readcount(bac_data)





################################################################################
## Data clean up and preprocessing
# Clean up unwanted sequences from mitochrondria e chloroplast
bac_data <- subset_taxa(bac_data, (Order!="Chloroplast") | is.na(Order))
bac_data <- subset_taxa(bac_data, (Family!="Mitochondria") | is.na(Family))
bac_data
readcount(bac_data) # Check how many reads have been lost

# Removing unwanted samples
bac_data <- subset_samples(bac_data, sample_names(bac_data) != "DCO_LLO_Bv4v5..PFF_PF170224" & sample_names(bac_data) !="DCO_LLO_Bv4v5..PFS_PF170222" & sample_names(bac_data) !="DCO_LLO_Bv4v5..PGF_PG170224" & sample_names(bac_data) !="DCO_LLO_Bv4v5..PGS_PG170224")
bac_data = filter_taxa(bac_data, function(x) sum(x) > 0, TRUE) # After removing samples filter the taxa left with zero global abundance
bac_data # get stats on the dataset

# Filter ASV with less than 5 reads
bac_data <- prune_taxa(taxa_sums(bac_data) > 5, bac_data) # remove taxa with less than X reads across samples
bac_data

# Removing the potential human pathogens and contaminants
bac_data <- subset_taxa(bac_data, (Genus != "Acinetobacter") |  is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Pseudomonas") |  is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Abiotrophia") |  is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Achromobacter") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Acinetobacter") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Actinobacillus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Arcanobacterium") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Arcobacter") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Babesia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Bacillus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Bartonella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Bordetella") |  is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Borrelia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Brodetella") |  is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Brucella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Burkholderia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Campylobacter") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Capnocytophaga") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Chlamydia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Clostridium") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Comamonas") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Corynebacterium") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Coxiella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Cronobacter") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Deinococcus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Dermatophilus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Ehrlichia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Enterococcus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Erysipelothrix") |  is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Escherichia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Escherichia/Shigella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Flavobacterium") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Francisella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Gardnerella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Granulicatella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Haemophilus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Hafnia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Halomonas") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Helicobacter") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Klebsiella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Kocuria") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Lawsonia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Legionella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Leptospira") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Listeria") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Merkel_cell") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Micrococcus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Morganella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Mycobacterium") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Mycoplasma") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Neisseria") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Nocardia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Pasteurella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Photobacterium") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Plesiomonas") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Propionibacterium") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Proteus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Providencia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Pseudomonas") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Rhodococcus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Rickettsiae") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Roseomonas") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Rothia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Salmonella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Serratia") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Shewanella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Shigella") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Sphaerophorus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Staphylococcus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Stenotrophomonas") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Streptococcus") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Treponema") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Vibrio") | is.na(Genus))
bac_data <- subset_taxa(bac_data, (Genus != "Yersinia") | is.na(Genus))

bac_data # check number of lost ASVs
readcount(bac_data) # Check number of lost reads for each samples

# Remove potential Eukaryotic domain assignments
bac_data <- subset_taxa(bac_data,  (Domain != "Eukaryota") | is.na(Domain))

bac_data # check number of lost ASVs
readcount(bac_data) # Check number of lost reads for each samples

## Checking for number of ASVs and reads at each steps helps to troubleshoot potential problems

## Make prevalence plots of the phyla
bac.prev0 = apply(X = otu_table(bac_data),
              MARGIN = ifelse(taxa_are_rows(bac_data), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
bac.prevdf = data.frame(Prevalence = bac.prev0,
                    TotalAbundance = taxa_sums(bac_data),
                    tax_table(bac_data))
keepPhyla = table(bac.prevdf$Phyla)[(table(bac.prevdf$Phyla) > 1)]
bac.prevdf1 = subset(bac.prevdf, Phyla %in% names(keepPhyla))

## Define prevalence threshold as 10% of total samples, i.e. taxa present in at least 10% of the samples investigated.
prevalenceThreshold = 0.1 * nsamples(bac_data)
prevalenceThreshold

# Plot prevalence
#png('prevalence_plot.png')
ggplot(bac.prevdf1, aes(TotalAbundance, Prevalence, color = Class)) +
  geom_hline(yintercept = prevalenceThreshold, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.3) +
  scale_y_log10() + scale_x_log10() +
  xlab("Prevalence plots on absolute abindance by Phyla") +
  facet_wrap(~Phyla) +
  theme(
  legend.position = "none"
  )
#dev.off()


##############################################################################
## Inspecting and plotting 16s rRNA diversity

## Subset and aggreate the dataset based on taxonomic level to ease computational work downstream
# Check the number of aggregated taxa
length(get_taxa_unique(bac_data, taxonomic.rank = "Genus"))
length(get_taxa_unique(bac_data, taxonomic.rank = "Family"))
length(get_taxa_unique(bac_data, taxonomic.rank = "Order"))
length(get_taxa_unique(bac_data, taxonomic.rank = "Class"))
length(get_taxa_unique(bac_data, taxonomic.rank = "Phyla"))

# Transform Bacteria abundance to relative abundances for plotting and some stats
bac_ra = transform_sample_counts(bac_data, function(x){x / sum(x)})

# Agglomerate at a specific taxonomic level at the Genus level
bac_ra_genus = tax_glom(bac_ra, "Genus", NArm = TRUE)
bac_ra_family = tax_glom(bac_ra, "Family", NArm = TRUE)
bac_ra_order = tax_glom(bac_ra, "Order", NArm = TRUE)
bac_ra_class = tax_glom(bac_ra, "Class", NArm = TRUE)
bac_ra_phyla = tax_glom(bac_ra, "Phyla", NArm = TRUE)

write_phyloseq(bac_data, type = "all")
write_phyloseq(bac_data, type = "TAXONOMY")

 # Define the taxa level to be plotted
bac_pha <- data.frame(t(otu_table(bac_ra_phyla))) # Phyla
bac_cla <- data.frame(t(otu_table(bac_ra_class))) # Class
bac_fam <- data.frame(t(otu_table(bac_ra_family))) # Family

# Define environmental vectors to reorder plot
ph <- as.vector(sample_data(bac_ra_phyla))$ph
temp <- as.vector(sample_data(bac_ra_phyla))$temp
trench <- as.vector(sample_data(bac_ra_phyla))$trench
fe_ac <- as.vector(sample_data(bac_ra_phyla))$fe_ac
fe_s <- as.vector(sample_data(bac_ra_phyla))$fe_s
cu_ac <- as.vector(sample_data(bac_ra_phyla))$cu_ac
cu_s <- as.vector(sample_data(bac_ra_phyla))$cu_s
h2s <- as.vector(sample_data(bac_ra_phyla))$h2s
so4 <- as.vector(sample_data(bac_ra_phyla))$so4_ac
dic <- as.vector(sample_data(bac_ra_phyla))$dic
dic_d13 <- as.vector(sample_data(bac_ra_phyla))$dic_d13
doc <- as.vector(sample_data(bac_ra_phyla))$doc
na <- as.vector(sample_data(bac_ra_phyla))$na
ca <- as.vector(sample_data(bac_ra_phyla))$ca
si <- as.vector(sample_data(bac_ra_phyla))$si
ca_na <- as.vector(sample_data(bac_ra_phyla))$ca/as.vector(sample_data(bac_ra_phyla))$na

# Plot the bubbleplot for selected taxa and vector
#svg("inkspot_phyla_ph.svg")
inkspot(bac_pha, ph, cex.max = 8,
          cex.axis = 0.3,
          legend.values =NULL,
          use.rank = T,
          x.axis.top = T,
          spec.names = get_taxa_unique(bac_ra_phyla, "Phyla"),
          site.names = data.frame(sample_data(bac_ra_phyla))$station,
          main = "Phyla abundance w/ sites ordered by ph on top axis"
        ) # Phyla with sites ordered by ph on top axis
#dev.off()

#svg("inkspot_phyla_temp.svg")
inkspot(bac_pha, temp, cex.max = 8,
          cex.axis = 0.3,
          legend.values =NULL,
          use.rank = T,
          x.axis.top = T,
          spec.names = get_taxa_unique(bac_ra_phyla, "Phyla"),
          site.names = data.frame(sample_data(bac_ra_phyla))$station,
          main = "Phyla abundance w/ sites ordered by temperature on top axis"
        )
#dev.off()

#svg("inkspot_phyla_na.svg")
inkspot(bac_pha, na, cex.max = 8,
          cex.axis = 0.3,
          legend.values =NULL,
          use.rank = T,
          x.axis.top = T,
          spec.names = get_taxa_unique(bac_ra_phyla, "Phyla"),
          site.names = data.frame(sample_data(bac_ra_phyla))$station,
          main = "Phyla abundance w/ sites ordered by [Na+] on top axis"
        )
#dev.off()

svg("inkspot_legend.svg")
inkspot(bac_pha, trench, cex.max = 8,
          cex.axis = 0.3,
      #    legend.values =NULL,
          use.rank = T,
          x.axis.top = T,
          spec.names = get_taxa_unique(bac_ra_phyla, "Phyla"),
          site.names = data.frame(sample_data(bac_ra_phyla))$station,
          main = "Phyla abundance w/ sites ordered by distance from trench on top axis"
        )
dev.off()

#svg("inkspot_class_na.svg")
inkspot(bac_cla, na, cex.max = 8,
          cex.axis = 0.3,
          legend.values =NULL,
          use.rank = T,
          x.axis.top = T,
          spec.names = get_taxa_unique(bac_ra_class, "Class"),
          site.names = data.frame(sample_data(bac_ra_class))$station,
          main = "Class abundance w/ sites ordered by [Na+] on top axis"
        )
#dev.off()

#svg("inkspot_class_dic.svg")
inkspot(bac_cla, dic, cex.max = 8,
          cex.axis = 0.3,
          legend.values =NULL,
          use.rank = T,
          x.axis.top = T,
          spec.names = get_taxa_unique(bac_ra_class, "Class"),
          site.names = data.frame(sample_data(bac_ra_class))$station,
          main = "Class abundance w/ sites ordered by [DIC] on top axis"
        )
#dev.off()

#svg("inkspot_class_trench.svg")
inkspot(bac_cla, trench, cex.max = 8,
          cex.axis = 0.3,
          legend.values =NULL,
          use.rank = T,
          x.axis.top = T,
          spec.names = get_taxa_unique(bac_ra_class, "Class"),
          site.names = data.frame(sample_data(bac_ra_class))$station,
          main = "Class abundance w/ sites ordered by distance from trench on top axis"
        )
#dev.off()

# Combined abundance plot of potential metabolic genera
# Subset specific taxa involved in selected metabolisms for downstream analysis
bac_iob <- subset_taxa(bac_ra_genus, Genus == "Gallionella" | Genus == "Ferriphaselus" | Genus == "Syderoxidans" | Genus == "Gallionellaceae_unclassified" |
                      Family == "Geobacteraceae" |
                      Genus == "Acidovorax" | Genus == "Leptothrix" | Genus == "Sphaerotilus" | Genus == "Cupriavidus" |   Genus == "Leptospirillum" |
                      Family == "Acidithiobacillaceae" |
                      Genus == "Rhodobacter" |  Genus == "Paracoccus" |  Genus == "Rhodovulum" | Genus == "Rhodoferax" |
                      Family == "Ferrovaceae" |
                      Genus == "Ferritrophicum" |
                      Order == "Mariprofundales" |
                       Genus == "Geothrix"| Genus == "Ferrimicrobium"| Genus == "Rhodopseudomonas"| Genus == "Rhodomicrobium"|
                      Genus == "Thiodictyon"
                    ) # Subsetting only the Iron Oxidizing genera

bac_sob <- subset_taxa(bac_ra_genus, Genus == "Beggiatoa" | Genus == "Acidithiobacillus" | Genus == "Paracoccus"
                       | Genus == "Sulfurovum" |
                      Family == "Sulfurimonas" |
                      Genus == "Caminibacter" | Genus == "Nautilia" | Genus == "Cetia" | Genus == "Thiovolum" |   Genus == "Leptospirillum" |
                      Family == "Thiotrix" |
                      Genus == "Sulfuritalea" |  Genus == "Sulfuricurvum" |  Genus == "Sulfurihydrogenibium" |
                       Genus == "Persephonella" |
                      Family == "Hydrogenobacter" |
                      Genus == "Aquifex" |
                      Order == "Hydrogenobaculum" |
                       Genus == "Chlorobium"| Genus == "Arcobacter"| Genus == "Magnetococcus"| Genus == "Chloroflexi"|
                      Genus == "Sulfurisoma" |  Genus == "Thioalkalivibrio" |  Genus == "Thiorhodovibrio" |
                        Genus == "Thiobacillus"
                    ) # Subsetting only the Sulfur Oxidizing genera

# Iron Oxidizers
plot_bar(bac_iob, fill="Genus", x="province", facet_grid = ~type, title="Abundance of iron oxidizers in fluid and sediments") +
  theme_hc() # Plot bar of Iron Oxidizing Bacteria between fluid and sediments

# Sulfur Oxidizers
plot_bar(bac_sob, fill="Genus", x="province", facet_grid = ~type, title="Abundance of sulfur oxidizers in fluid and sediments") +
  theme_hc() # Plot bar of Iron Oxidizing Bacteria between fluid and sediments

# Iron Oxidizers
plot_bar(bac_iob, fill="Genus", x="province", facet_grid = ~plate, title="Abundance of iron oxidizers between plate types") +
  theme_hc() # Plot bar of Iron Oxidizing Bacteria between the two plates
ggsave("iron_oxidizers.png", width=6, height=6)

plot_bar(bac_sob, fill="Genus", x="province", facet_grid = ~plate, title="Abundance of sulfur oxidizers between plate types") +
  theme_hc() # Plot bar of Iron Oxidizing Bacteria  the two plates
ggsave("sulfur_oxidizers.png", width=6, height=6)


##################################################################################
## Geochemistry plots

## Load the environmental dataset
# Complete observations with respect to the environmental data
# loaded for the coupled microbial-environmental analysis
dataset <- read.csv("../subductCR_final_dataset.csv", header = T, sep=",")
dataset$geol_prov<-as.factor(dataset$geol_prov) # correct import problem

# Concentration in our dataset are in µmol/L and needs to be converted in mg/Kg for the Giggenplot after Giggenbach (1988)
# Conversion can be made on the fly with ggtern
# Giggenbach plot of major water types in geothermal environemnts
ggtern(data=dataset, aes(x=so4*96,y=cl*35.45, z=dic*61)) +
geom_point(aes(color=anions), size=3) +
geom_text(aes(label= code), size=3, hjust = 1, vjust = 0.5, check_overlap = F)
ggsave("ternary_anions_fluids.svg", width=6, height=6)

# Giggenbach plot of major cations in geothermal fluids
ggtern(data=dataset, aes(x=ca*40.1,y=mg*24.3, z=(na*23)+(k*39.1))) +
geom_point(aes(color=cations), size=3) +
geom_text(aes(label= code), size=3, hjust = 1, vjust = -1, check_overlap = F)
ggsave("ternary_cations_fluids.svg", width=6, height=6)

# Equilibrium geothermometer after Giggenbach (1988)
ggtern(data=dataset, aes(x=(k*39.1)/10,y=(na*23)/400, z=(mg*24.3)^0.5)) +
geom_point(aes(color=province), size=3) +
geom_text(aes(label= code), size=3, hjust = -0.5, vjust = 1, check_overlap = F)
ggsave("ternary_geothermometer_fluids.svg", width=6, height=6)




###########################################################################
## Plotting of additional diversity plots for the supplementary materials

# Alpha diversity estimates
plot_richness(bac_data, x="province", measures =c("Simpson"), color="plate") +
   geom_boxplot() +
   geom_jitter() +
   theme_hc()

# Plot Relative Abundance by Phyla for each station
plot_bar(bac_ra_phyla, fill="Phyla", facet_grid = ~type, x="station") +
  theme_hc()
ggsave("diversity_phyla.svg", width=10, height=8)

# Plot Relative Abundance by Phyla for each province
plot_bar(bac_ra_phyla, fill="Phyla", facet_grid = ~type, x="province") +
  theme_hc()
ggsave("diversity_phyla_province.svg", width=10, height=8)

# Plot Relative Abundance by Phyla for each plate type
plot_bar(bac_ra_phyla, fill="Phyla", facet_grid = ~type, x="plate", title ="Relative Abundance by Phyla for each plate type") +
  theme_hc()
ggsave("diversity_phyla_plate.svg", width=10, height=8)

# Plot Relative Abundance by Phyla for each geology
plot_bar(bac_ra_phyla, fill="Phyla", facet_grid = ~type, x="geol_prov", title ="Relative Abundance by Phyla for each geology") +
  theme_hc()
ggsave("diversity_phyla_geology.svg", width=10, height=8)

# Plot Relative Abundance by Phyla for each volcanic province
plot_bar(bac_ra_phyla, fill="Phyla", facet_grid = ~type, x="volcano", title ="Relative Abundance by Phyla for each volcanic province") +
  theme_hc()
ggsave("diversity_phyla_volcano.svg", width=10, height=8)

# Plot Relative Abundance by Phyla for each geochemistry regime obtained from the geochemical analysis detailed below
plot_bar(bac_ra_phyla, fill="Phyla", facet_grid = ~type, x="cations", title ="Relative Abundance by Phyla for each cation geochemistry") +
  theme_hc()
ggsave("diversity_phyla_cations.svg", width=10, height=8)

# Find top 5% abundant Phyla apprearing at least in 40% of the samples (i.e. domiant phylotypes)
bac_ra_top5 <- prune_taxa(genefilter_sample(bac_ra, filterfun_sample(topp(0.05)), A = round(0.4 * nsamples(bac_ra))),bac_ra) #top 5% ASV in 40% samples
get_taxa_unique(bac_ra_top5, "Phyla")

plot_bar(subset_taxa(bac_ra_class, Phyla == "Acidobacteria"), fill="Class", facet_grid = ~type,
         x="station", title = "Acidobacteria class level") +
         theme_hc()
ggsave("diversity_barplot_acidobacteria.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_class, Phyla == "Actinobacteria"), fill="Class", facet_grid = ~type,
         x="station", title = "Actinobacteria class level") +
         theme_hc()
ggsave("diversity_barplot_actinobacteria.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_family, Phyla == "Aquificae"), fill="Family", facet_grid = ~type,
         x="station", title = "Aquificae famiy level") +
         theme_hc()
ggsave("diversity_barplot_aquificae.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_family, Phyla == "Armatimonadetes"), fill="Family", facet_grid = ~type,
         x="station", title = "Armatimonadetes famiy level") +
         theme_hc()
ggsave("diversity_barplot_armatimonadetes.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_family, Phyla == "Bacteria_unclassified"), fill="Family", facet_grid = ~type,
         x="station", title = "Unclassified Bacteria famiy level") +
         theme_hc()
ggsave("diversity_barplot_unclassified.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_class, Phyla == "Bacteroidetes"), fill="Class", facet_grid = ~type,
         x="station", title = "Bacteroidetes class level") +
         theme_hc()
ggsave("diversity_barplot_bacteroidetes.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_class, Phyla == "Chloroflexi"), fill="Class", facet_grid = ~type,
         x="station", title = "Chloroflexi class level") +
         theme_hc()
ggsave("diversity_barplot_chloroflexi.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_family, Phyla == "Cyanobacteria"), fill="Family", facet_grid = ~type,
         x="station", title = "Cyanobacteria family level") +
         theme_hc()
ggsave("diversity_barplot_cyanobacteria.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_genus, Phyla == "Deinococcus-Thermus"), fill="Genus", facet_grid = ~type,
         x="station", title = "Deinococcus genus level") +
         theme_hc()
ggsave("diversity_barplot_deinococcus.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_family, Phyla == "Firmicutes"), fill="Family", facet_grid = ~type,
         x="station", title = "Firmicutes family level") +
         theme_hc()
ggsave("diversity_barplot_firmicutes.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_class, Phyla == "GAL15"), fill="Class", facet_grid = ~type,
         x="station", title = "Uncultured phyla GAL15 class level") +
         theme_hc()
ggsave("diversity_barplot_gal15.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_family, Phyla == "Nitrospirae"), fill="Family", facet_grid = ~type,
         x="station", title = "Nitrospirae family level") +
         theme_hc()
ggsave("diversity_barplot_nitrospirae.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_class, Phyla == "Patescibacteria"), fill="Class", facet_grid = ~type,
         x="station", title = "Patescibacteria class level") +
         theme_hc()
ggsave("diversity_barplot_patescibacteria.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_class, Phyla == "Planctomycetes"), fill="Class", facet_grid = ~type,
         x="station", title = "Planctomycetes class level") +
         theme_hc()
ggsave("diversity_barplot_planctomycetes.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_family, Phyla == "Rokubacteria"), fill="Family", facet_grid = ~type,
         x="station", title = "Rokubacteria Family level") +
         theme_hc()
ggsave("diversity_barplot_rokubacteria.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_family, Phyla == "Verrucomicrobia"), fill="Family", facet_grid = ~type,
         x="station", title = "Verrucomicrobia family level") +
         theme_hc()
ggsave("diversity_barplot_verrucomicrobia.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_class, Phyla == "Proteobacteria"), fill="Class", facet_grid = ~type,
         x="station", title = "Proteobacteria class level") +
         theme_hc()
ggsave("diversity_barplot_proteobacteria.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_order, Class == "Alphaproteobacteria"), fill="Order", facet_grid = ~type,
         x="station", title = "Alphaproteobacteria order level") +
         theme_hc()
ggsave("diversity_barplot_alphaproteobacteria.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_order, Class == "Deltaproteobacteria"), fill="Order", facet_grid = ~type,
         x="station", title = "Deltaaproteobacteria order level") +
         theme_hc()
ggsave("diversity_barplot_deltaproteobacteria.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_order, Class == "Gammaproteobacteria"), fill="Order", facet_grid = ~type,
         x="station", title = "Gammaproteobacteria order level") +
         theme_hc()
ggsave("diversity_barplot_gammaproteobacteria.svg", width=10, height=8)

plot_bar(subset_taxa(bac_ra_family, Order == "Betaproteobacteriales"), fill="Family", facet_grid = ~type,
         x="station", title = "Betaproteobactriales family level") +
         theme_hc()
ggsave("diversity_barplot_betaproteobactriales.svg", width=10, height=8)

# Trees can also be reduced usin tip_glom before plotting the tree
plot_tree(subset_taxa(bac_ra, Phyla == "Aquificae"), color = "province", label.tips = "Genus",
          size = "abundance", plot.margin = 0.5, ladderize = TRUE, title="Aquificae tree") + facet_grid(~type)
ggsave("tree_aquificae.svg", width=10, height=8)

# Trees can also be reduced usin tip_glom before plotting the tree
plot_tree(subset_taxa(bac_ra, Phyla == "Bacteria_unclassified"), color = "province", label.tips = "Genus",
          size = "abundance", plot.margin = 0.5, ladderize = TRUE, title="Unclassified Bacteria tree") + facet_grid(~type)
ggsave("tree_unclassified.svg", width=10, height=8)

# Trees can also be reduced usin tip_glom before plotting the tree
plot_tree(subset_taxa(bac_ra, Phyla == "Chloroflexi"), color = "province", label.tips = "Genus",
          size = "abundance", plot.margin = 0.5, ladderize = TRUE, title="Chloroflexi tree") + facet_grid(~type)
ggsave("tree_chloroflexi.svg", width=10, height=8)

# Trees can also be reduced usin tip_glom before plotting the tree
plot_tree(subset_taxa(bac_ra, Phyla == "Acidobacteria"), color = "province", label.tips = "Genus",
          size = "abundance", plot.margin = 0.5, ladderize = TRUE, title="Acidobacteria tree") + facet_grid(~type)
ggsave("tree_acidobacteria.svg", width=10, height=8)

# Trees can also be reduced usin tip_glom before plotting the tree
plot_tree(subset_taxa(bac_ra, Phyla == "Bacteroidetes"), color = "province", label.tips = "Genus",
          size = "abundance", plot.margin = 0.5, ladderize = TRUE, title="Bacteroidetes tree") + facet_grid(~type)
ggsave("tree_bacteroidetes.svg", width=10, height=8)

# Trees can also be reduced usin tip_glom before plotting the tree
plot_tree(subset_taxa(bac_ra, Phyla == "Cyanobacteria"), color = "province", label.tips = "Genus",
          size = "abundance", plot.margin = 0.5, ladderize = TRUE, title="Cyanobacteria tree") + facet_grid(~type)
ggsave("tree_cyanobacteria.svg", width=10, height=8)

# Trees can also be reduced usin tip_glom before plotting the tree
plot_tree(subset_taxa(bac_ra, Phyla == "Nitrospirae"), color = "province", label.tips = "Genus",
          size = "abundance", plot.margin = 0.5, ladderize = TRUE, title="Nitrospirae tree") + facet_grid(~type)
ggsave("tree_nitrospirae.svg", width=10, height=8)

# Trees can also be reduced usin tip_glom before plotting the tree
plot_tree(subset_taxa(bac_ra, Phyla == "Planctomycetes"), color = "province", label.tips = "Genus",
          size = "abundance", plot.margin = 0.5, ladderize = TRUE, title="Planctomycetes tree") + facet_grid(~type)
ggsave("tree_planctomycetes.svg", width=10, height=8)

# Trees can also be reduced usin tip_glom before plotting the tree
plot_tree(subset_taxa(bac_ra, Phyla == "Rokubacteria"), color = "province", label.tips = "Genus",
          size = "abundance", plot.margin = 0.5, ladderize = TRUE, title="Rokubacteria tree") + facet_grid(~type)
ggsave("tree_rokubacteria.svg", width=10, height=8)

# Trees can also be reduced usin tip_glom before plotting the tree
plot_tree(subset_taxa(bac_ra, Phyla == "Verrucomicrobia"), color = "province", label.tips = "Genus",
          size = "abundance", plot.margin = 0.5, ladderize = TRUE, title="Verrucomicrobia tree") + facet_grid(~type)
ggsave("tree_verrucomicrobia.svg", width=10, height=8)

# Trees can also be reduced usin tip_glom before plotting the tree
plot_tree(subset_taxa(bac_ra, Phyla == "Firmicutes"), color = "province", label.tips = "Genus",
          size = "abundance", plot.margin = 0.5, ladderize = TRUE, title="Firmicutes tree") + facet_grid(~type)
ggsave("tree_firmicutes.svg", width=10, height=8)

# Trees can also be reduced usin tip_glom before plotting the tree
plot_tree(subset_taxa(bac_ra, Class == "Alphaproteobacteria"), color = "province", label.tips = "Genus",
          size = "abundance", plot.margin = 0.5, ladderize = TRUE, title="Alphaproteobacteria tree") + facet_grid(~type)
ggsave("tree_alphaproteobacteria.svg", width=10, height=8)

# Trees can also be reduced usin tip_glom before plotting the tree
plot_tree(subset_taxa(bac_ra, Class == "Deltaproteobacteria"), color = "province", label.tips = "Genus",
          size = "abundance", plot.margin = 0.5, ladderize = TRUE, title="Deltaproteobacteria tree") + facet_grid(~type)
ggsave("tree_deltaproteobacteria.svg", width=10, height=8)

# Trees can also be reduced usin tip_glom before plotting the tree
plot_tree(subset_taxa(bac_ra, Class == "Gammaproteobacteria"), color = "province", label.tips = "Genus",
          size = "abundance", plot.margin = 0.5, ladderize = TRUE, title="Gammaproteobacteria tree") + facet_grid(~type)
ggsave("tree_gammaproteobacteria.svg", width=10, height=8)

# Trees can also be reduced usin tip_glom before plotting the tree
plot_tree(subset_taxa(bac_ra, Class == "Deltaproteobacteria"), color = "temp", label.tips = "Genus",
          size = "abundance", plot.margin = 0.5, ladderize = TRUE, title="Deltaproteobacteria tree") + facet_grid(~type)
ggsave("tree_deltaproteobacteria_temp.svg", width=10, height=8)

# Trees can also be reduced usin tip_glom before plotting the tree
plot_tree(subset_taxa(bac_ra, Class == "Gammaproteobacteria"), color = "temp", label.tips = "Genus",
          size = "abundance", plot.margin = 0.5, ladderize = TRUE, title="Gammaproteobacteria tree") + facet_grid(~type)
ggsave("tree_gammaproteobacteria_temp.svg", width=10, height=8)

# Trees can also be reduced usin tip_glom before plotting the tree
plot_tree(subset_taxa(bac_ra, Phyla == "Aquificae"), color = "temp", label.tips = "Genus",
          size = "abundance", plot.margin = 0.5, ladderize = TRUE, title="Aquificae tree") + facet_grid(~type)
ggsave("tree_aquificae_temp.svg", width=10, height=8)




#####################################################
## Starting the multivariate ecological analysis
## PCoA unifrac weighted with relative counts
bac_pcoa_w <- ordinate(bac_data, method = "PCoA", distance = "unifrac", weighted=T)
evals_w <- bac_pcoa_w$values$Eigenvalues

## PCoA unifrac unweighted with relative counts
bac_pcoa_un <- ordinate(bac_data, method = "PCoA", distance = "unifrac", weighted=F)
evals <- bac_pcoa_un$values$Eigenvalues

## nMDS with Jaccard and Bray-Curtis distance
bac_nmds_j <- ordinate(bac_data, method = "NMDS", distance = "jaccard", weighted=T, trymax=100)
bac_nmds_bc <- ordinate(bac_data, method = "NMDS", distance = "bray", weighted=T, trymax=100)

## Compare ordinations colored by temperature
grid.arrange(nrow = 2, ncol=2,
             plot_ordination(bac_data, bac_pcoa_w, type = "samples", color = "temp", shape="type", label="code", title="PCoA weighted Unifrac") +
                 theme(legend.position = "none") +
               coord_fixed(sqrt(evals_w[2] / evals_w[1])),
             plot_ordination(bac_data, bac_pcoa_un, type = "samples", color = "temp",  shape="type", label="code",title="PCoA weighted Unifrac") +
                 theme(legend.position = "none")+
               coord_fixed(sqrt(evals[2] / evals[1])),
             plot_ordination(bac_data, bac_nmds_j, color="temp", shape="type", label="code", title="nMDS Jaccard distances") +  theme(legend.position = "none"),
             plot_ordination(bac_data, bac_nmds_bc, color="temp", shape="type", label="code", title="nMDS Bray-Curtis distances") +   theme(legend.position = "none")
             )

## Compare ordinations by ph
grid.arrange(nrow = 2, ncol=2,
             plot_ordination(bac_data, bac_pcoa_w, type = "samples", color = "ph", shape="type", label="code", title="PCoA weighted Unifrac") +
                 theme(legend.position = "none") +
               coord_fixed(sqrt(evals_w[2] / evals_w[1])),
             plot_ordination(bac_data, bac_pcoa_un, type = "samples", color = "ph",  shape="type", label="code",title="PCoA weighted Unifrac") +
                 theme(legend.position = "none")+
               coord_fixed(sqrt(evals[2] / evals[1])),
             plot_ordination(bac_data, bac_nmds_j, color="ph", shape="type", label="code", title="nMDS Jaccard distances") +  theme(legend.position = "none"),
             plot_ordination(bac_data, bac_nmds_bc, color="ph", shape="type", label="code", title="nMDS Bray-Curtis distances") +   theme(legend.position = "none")
             )

## Compare ordinations by geological province
grid.arrange(nrow = 2, ncol=2,
             plot_ordination(bac_data, bac_pcoa_w, type = "samples", color = "geol_prov", shape="type", label="code", title="PCoA weighted Unifrac") +
                 theme(legend.position = "none") +
               coord_fixed(sqrt(evals_w[2] / evals_w[1])),
             plot_ordination(bac_data, bac_pcoa_un, type = "samples", color = "geol_prov",  shape="type", label="code",title="PCoA weighted Unifrac") +
                 theme(legend.position = "none")+
               coord_fixed(sqrt(evals[2] / evals[1])),
             plot_ordination(bac_data, bac_nmds_j, color="geol_prov", shape="type", label="code", title="nMDS Jaccard distances") +  theme(legend.position = "none"),
             plot_ordination(bac_data, bac_nmds_bc, color="geol_prov", shape="type", label="code", title="nMDS Bray-Curtis distances") +   theme(legend.position = "none")
             )

## Compare ordinations by cations groups as define in the ternary plot
grid.arrange(nrow = 2, ncol=2,
             plot_ordination(bac_data, bac_pcoa_w, type = "samples", color = "cations", shape="type", label="code", title="PCoA weighted Unifrac") +
                 theme(legend.position = "none") +
               coord_fixed(sqrt(evals_w[2] / evals_w[1])),
             plot_ordination(bac_data, bac_pcoa_un, type = "samples", color = "cations",  shape="type", label="code",title="PCoA weighted Unifrac") +
                 theme(legend.position = "none")+
               coord_fixed(sqrt(evals[2] / evals[1])),
             plot_ordination(bac_data, bac_nmds_j, color="cations", shape="type", label="code", title="nMDS Jaccard distances") +  theme(legend.position = "none"),
             plot_ordination(bac_data, bac_nmds_bc, color="province", shape="type", label="code", title="nMDS Bray-Curtis distances") +   theme(legend.position = "none")
             )

## Compare ordinations by anions as defined in the thernary plot
grid.arrange(nrow = 2, ncol=2,
             plot_ordination(bac_data, bac_pcoa_w, type = "samples", color = "anions", shape="type", label="code", title="PCoA weighted Unifrac") +
                 theme(legend.position = "none") +
               coord_fixed(sqrt(evals_w[2] / evals_w[1])),
             plot_ordination(bac_data, bac_pcoa_un, type = "samples", color = "anions",  shape="type", label="code",title="PCoA weighted Unifrac") +
                 theme(legend.position = "none")+
               coord_fixed(sqrt(evals[2] / evals[1])),
             plot_ordination(bac_data, bac_nmds_j, color="anions", shape="type", label="code", title="nMDS Jaccard distances") +  theme(legend.position = "none"),
             plot_ordination(bac_data, bac_nmds_bc, color="anions", shape="type", label="code", title="nMDS Bray-Curtis distances") +   theme(legend.position = "none")
             )

## Summary plot of ordinations
#svg("summary_ordinations.svg")
png("summary_ordinations.png", width=600, height=900)
grid.arrange(nrow = 3, ncol=2,
 plot_ordination(bac_data, bac_nmds_bc, color="temp", shape="type", label="code", title="nMDS Bray-Curtis Temp") + theme(legend.position = "none"),
 plot_ordination(bac_data, bac_nmds_bc, color="ph", shape="type", label="code", title="nMDS Bray-Curtis pH") + theme(legend.position = "none"),
 plot_ordination(bac_data, bac_nmds_bc, color="geol_prov", shape="type", label="code", title="nMDS Bray-Curtis Geology") + theme(legend.position = "none"),
 plot_ordination(bac_data, bac_nmds_bc, color="rocks", shape="type", label="code", title="nMDS Bray-Curtis solid Rocks") + theme(legend.position = "none"),
 plot_ordination(bac_data, bac_nmds_bc, color="fe_s", shape="type", label="code", title="nMDS Bray-Curtis solid [Fe]") + theme(legend.position = "none"),
 plot_ordination(bac_data, bac_nmds_bc, color="dic_d13", shape="type", label="code", title="nMDS Bray-Curtis solid [DIC]") + theme(legend.position = "none")
  )
dev.off()

## Poas Lake is not a hot spring, but a acidic volcanic lake, therefore it will be removed from downstream analyses
# Removing Poas Lake sample
bac_ra_nopl <- subset_samples(bac_ra, sample_names(bac_ra) != "DCO_LLO_Bv4v5..PLS_PL170224")
bac_data_nopl <- subset_samples(bac_data, sample_names(bac_data) != "DCO_LLO_Bv4v5..PLS_PL170224")
bac_ra_nopl <- filter_taxa(bac_ra_nopl, function(x) sum(x) > 0, TRUE) # After removing samples filter the taxa left with zero global abundance
bac_data_nopl <- filter_taxa(bac_data_nopl, function(x) sum(x) > 0, TRUE) # After removing samples filter the taxa left with zero global abundance
bac_ra_nopl # get stats on the dataset

## PCoA unifrac weighted with relative counts
bac_pl_pcoa_w <- ordinate(bac_data_nopl, method = "PCoA", distance = "unifrac", weighted=T)
evals_pl_w <- bac_pl_pcoa_w$values$Eigenvalues

## PCoA unifrac unweighted with relative counts
bac_pl_pcoa_un <- ordinate(bac_data_nopl, method = "PCoA", distance = "unifrac", weighted=F)
evals_pl <- bac_pl_pcoa_un$values$Eigenvalues

## nMDS with Jaccard and Bray-Curtis distance
bac_pl_nmds_j <- ordinate(bac_data_nopl, method = "NMDS", distance = "jaccard", weighted=T, trymax=100)
bac_pl_nmds_bc <- ordinate(bac_data_nopl, method = "NMDS", distance = "bray", weighted=T, trymax=100)

## Compare ordinations colored by temperature with Poas removed
grid.arrange(nrow = 2, ncol=2,
             plot_ordination(bac_data_nopl, bac_pl_pcoa_w, type = "samples", color = "temp", shape="type", label="code", title="PCoA weighted Unifrac") +
                 theme(legend.position = "none") +
               coord_fixed(sqrt(evals_pl_w[2] / evals_pl_w[1])),
             plot_ordination(bac_data_nopl, bac_pl_pcoa_un, type = "samples", color = "temp",  shape="type", label="code",title="PCoA weighted Unifrac") +
                 theme(legend.position = "none")+
               coord_fixed(sqrt(evals_pl[2] / evals_pl[1])),
             plot_ordination(bac_data_nopl, bac_pl_nmds_j, color="temp", shape="type", label="code", title="nMDS Jaccard distances") +  theme(legend.position = "none"),
             plot_ordination(bac_data, bac_nmds_j, color="temp", shape="type", label="code", title="nMDS Jaccard distances no POAS") +  theme(legend.position = "none")
             )

## Adonis analysis on the beta diversity
### Test for factors difference in controlling beta diversity
adonis(distance(bac_data, method="wunifrac") ~ type*plate*province*region*geol_prov*rocks*volcano*anions*cations, data = data.frame(sample_data(bac_data)), perm = 999)
adonis(distance(bac_data, method="unifrac") ~ type*plate*province*region*geol_prov*rocks*volcano*anions*cations, data = data.frame(sample_data(bac_data)), perm = 999)
adonis(distance(bac_data, method="jaccard") ~ type*plate*province*region*geol_prov*rocks*volcano*anions*cations, data = data.frame(sample_data(bac_data)), perm = 999)

### Test for factors difference in controlling beta diversity removing PL
adonis(distance(bac_data_nopl, method="wunifrac") ~ type*plate*province*region*geol_prov*rocks*volcano*anions*cations, data = data.frame(sample_data(bac_data_nopl)), perm = 999)
adonis(distance(bac_data_nopl, method="unifrac") ~ type*plate*province*region*geol_prov*rocks*volcano*anions*cations, data = data.frame(sample_data(bac_data_nopl)), perm = 999)
adonis(distance(bac_data_nopl, method="jaccard") ~ type*plate*province*region*geol_prov*rocks*volcano*anions*cations, data = data.frame(sample_data(bac_data_nopl)), perm = 999)

# Phyloseq uses a vegan wrapper to make the ordinations, so the results of the ordinations are
# identical although stored in two different objects in vegan and phyloseq
# Using custom scripts for converting the phyloseq into a vegan friendly data.frame
bac_data.v<-psotu2veg(bac_data_nopl) # custom function to export phyloseq objects to vegan
bac.v_nmds_j<-metaMDS(bac_data.v, methods="jaccard", trymax=100) #nMDS with jaccard distances
plot(bac.v_nmds_j, display = "sites", main = "nMDS Jaccard")

# Partition the environmental data into meaninful groups for the linear vector fitting
print(names(data.frame(sample_data(bac_data_nopl))))

ambient <- data.frame(sample_data(bac_data_nopl))[,18:26]
isotopes <- data.frame(sample_data(bac_data_nopl))[,27:36]
gas <- data.frame(sample_data(bac_data_nopl))[,37:47]
hydrocarbons <- data.frame(sample_data(bac_data_nopl))[,49:52]
om <- data.frame(sample_data(bac_data_nopl))[,53:58]
geochem_met <- data.frame(sample_data(bac_data_nopl))[,59:68]
geochem_an <- data.frame(sample_data(bac_data_nopl))[,69:82]
geochem_s <- data.frame(sample_data(bac_data_nopl))[,83:91]

# Vector fitting
envfit(bac.v_nmds_j, ambient, perm = 9999, na.rm = T)

# Vector fitting
envfit(bac.v_nmds_j, isotopes, perm = 9999, na.rm = T)

# Vector fitting
envfit(bac.v_nmds_j, gas, perm = 9999, na.rm = T)

# Vector fitting
envfit(bac.v_nmds_j, hydrocarbons, perm = 9999, na.rm = T)

# Vector fitting
envfit(bac.v_nmds_j, om, perm = 9999, na.rm = T)

# Vector fitting
envfit(bac.v_nmds_j, geochem_met, perm = 9999, na.rm = T)

# Vector fitting
envfit(bac.v_nmds_j, geochem_an, perm = 9999, na.rm = T)

# Vector fitting
envfit(bac.v_nmds_j, geochem_s, perm = 9999, na.rm = T)

# Testing identified correlations using the Pearson moment correlation
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$temp)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$ph)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$dic)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$dic_d13)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$doc_d13)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$tn)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$tn_d15)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$cha)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$al_ac)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$cr_ac)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$fe_ac)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$ni_ac)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$zn_ac)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$cl)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$po4)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$ca)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$al_s)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$cr_s)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$ni_s)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$cu_s)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$ch4)
cor.test(scores(bac.v_nmds_j)[,1], data.frame(sample_data(bac_data_nopl))$k)


## NMDS2
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$temp)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$ph)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$dic)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$dic_d13)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$doc_d13)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$tn)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$tn_d15)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$cha)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$al_ac)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$cr_ac)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$fe_ac)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$ni_ac)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$zn_ac)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$cl)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$po4)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$ca)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$al_s)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$cr_s)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$ni_s)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$cu_s)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$ch4)
cor.test(scores(bac.v_nmds_j)[,2], data.frame(sample_data(bac_data_nopl))$k)

# Plot of vector fitting correlations against nMDS1
#png("vector_fitting_correlations_nmds1.png", width=700, height=700)
#svg("vector_fitting_correlations_nmds1.svg")
par(mfrow=c(4,4), mar=c(4,4,0.5,0.5))
plot(scores(bac.v_nmds_j)[,1],data.frame(sample_data(bac_data_nopl))$temp, xlab="nMDS1", ylab="temp")
plot(scores(bac.v_nmds_j)[,1],data.frame(sample_data(bac_data_nopl))$ph, xlab="nMDS1", ylab="pH")
plot(scores(bac.v_nmds_j)[,1],data.frame(sample_data(bac_data_nopl))$dic, xlab="nMDS1", ylab="DIC")
plot(scores(bac.v_nmds_j)[,1],data.frame(sample_data(bac_data_nopl))$dic_d13, xlab="nMDS1", ylab="DIC d13C")
plot(scores(bac.v_nmds_j)[,1],data.frame(sample_data(bac_data_nopl))$doc_d13, xlab="nMDS1", ylab="DOC d13C")
plot(scores(bac.v_nmds_j)[,1],data.frame(sample_data(bac_data_nopl))$al_ac, xlab="nMDS1", ylab="Al_aq")
plot(scores(bac.v_nmds_j)[,1],data.frame(sample_data(bac_data_nopl))$cr_ac, xlab="nMDS1", ylab="Cr_aq")
plot(scores(bac.v_nmds_j)[,1],data.frame(sample_data(bac_data_nopl))$fe_ac, xlab="nMDS1", ylab="Fe_aq")
plot(scores(bac.v_nmds_j)[,1],data.frame(sample_data(bac_data_nopl))$ni_ac, xlab="nMDS1", ylab="Ni_aq")
plot(scores(bac.v_nmds_j)[,1],data.frame(sample_data(bac_data_nopl))$zn_ac, xlab="nMDS1", ylab="zn_aq")
plot(scores(bac.v_nmds_j)[,1],data.frame(sample_data(bac_data_nopl))$cl, xlab="nMDS1", ylab="Cl")
plot(scores(bac.v_nmds_j)[,1],data.frame(sample_data(bac_data_nopl))$po4, xlab="nMDS1", ylab="PO4")
plot(scores(bac.v_nmds_j)[,1],data.frame(sample_data(bac_data_nopl))$ca, xlab="nMDS1", ylab="Ca")
plot(scores(bac.v_nmds_j)[,1],data.frame(sample_data(bac_data_nopl))$cr_s, xlab="nMDS1", ylab="Cr_s")
plot(scores(bac.v_nmds_j)[,1],data.frame(sample_data(bac_data_nopl))$ni_s, xlab="nMDS1", ylab="Ni_s")
plot(scores(bac.v_nmds_j)[,1],data.frame(sample_data(bac_data_nopl))$k, xlab="nMDS1", ylab="K")
#dev.off()

# Plot of vector fitting correlations againts nMDS2
#png("vector_fitting_correlations_nmds2.png", width=700, height=700)
#svg("vector_fitting_correlations_nmds2.svg")
par(mfrow=c(3,2), mar=c(4,4,0.5,0.5))
plot(scores(bac.v_nmds_j)[,2],data.frame(sample_data(bac_data_nopl))$tn, xlab="nMDS2", ylab="TN")
plot(scores(bac.v_nmds_j)[,2],data.frame(sample_data(bac_data_nopl))$tn_d15, xlab="nMDS2", ylab="TN d15N")
plot(scores(bac.v_nmds_j)[,2],data.frame(sample_data(bac_data_nopl))$cr_ac, xlab="nMDS2", ylab="Cr_aq")
plot(scores(bac.v_nmds_j)[,2],data.frame(sample_data(bac_data_nopl))$ca, xlab="nMDS2", ylab="Ca")
plot(scores(bac.v_nmds_j)[,2],data.frame(sample_data(bac_data_nopl))$al_s, xlab="nMDS2", ylab="Al_s")
plot(scores(bac.v_nmds_j)[,2],data.frame(sample_data(bac_data_nopl))$cu_s, xlab="nMDS2", ylab="Cu_s")
#dev.off()




################################################################################
## Prepare the data for the co-occurrence analysis
## Low abundance and low prevalence ASV need to be removed

# Prevalence filtering on datatset with PL removed
bac_data.pf <- prune_taxa(taxa_sums(bac_data_nopl) > 20, bac_data_nopl) # remove taxa with less than X reads across samples
bac_data.pf <- filter_taxa(bac_data.pf, function(x) sum(x>0) >= 2.8, TRUE) # Filter based on presence in at least 3 samples
bac_data.pf  <- prune_taxa(taxa_sums(bac_data.pf) > 0, bac_data.pf) # Remove ASV with 0 presence to clean up
# Inspecting ASVs and read count numbers after the filtering
bac_data_pf
readcount(bac_data.pf)

# Barplot phylum level diversity before and after the last filtering step to see if global apttern of the major phyla have changed
bac_data_f.ra = transform_sample_counts(bac_data.pf, function(x){x / sum(x)})
bac_data_f.ra_phyla = tax_glom(bac_data_f.ra, "Phyla", NArm = TRUE)
plot_bar(bac_data_f.ra_phyla, fill="Phyla", x="station", facet_grid=~type) +
theme_hc()
ggsave("diversity_phyla_postpf.svg", width=10, height=8)


## Building the co-occurrence network object
# Extract OTU table and Taxonomy tables from the prevelence filtered dataset
bac_form <- microbiomeutilities::format_to_besthit(bac_data.pf)
bac.otu <- t(otu_table(bac_form)@.Data)
bac.tax <- as.data.frame(tax_table(bac_form)@.Data)

## Co-occurrence matrix construction
# Retaining edges with Spearman correlation > 0.65
bac.cor <- cor(bac.otu, method="spearman")
bac.cor[bac.cor < 0.7] = 0
bac.cor.ig <- graph.adjacency(bac.cor, mode='undirected', add.rownames = TRUE, weighted = TRUE)
bac.cor.ig <- igraph::simplify(bac.cor.ig)

# Set Phyla coloring
bac.cor.phyla <- map_levels(colnames(bac.otu), from = "best_hit", to = "Phylum", tax_table(bac_form))

# Set up the colors for the phyla
c47 <- c("#000000","#000080","#00008B","#0000CD","#0000FF","#006400","#008000","#008080","#008B8B","#00BFFF","#00CED1",
"#00FA9A","#00FF00","#00FF7F","#00FFFF","#00FFFF","#191970","#1E90FF","#20B2AA","#228B22","#2E8B57","#2F4F4F","#32CD32","#3CB371","#40E0D0","#4169E1","#4682B4",
"#483D8B","#48D1CC","#4B0082","#556B2F","#5F9EA0","#6495ED","#66CDAA","#696969","#6A5ACD","#6B8E23","#708090","#778899","#7B68EE","#7CFC00","#7FFF00","#7FFFD4",
"#800000","#800080","#808000","#808080","#87CEEB","#87CEFA","#8A2BE2","#8B0000","#8B008B","#8B4513","#8FBC8F","#90EE90","#9370DB","#9400D3","#98FB98","#9932CC",
"#9ACD32","#A0522D","#A52A2A","#A9A9A9","#ADD8E6","#ADFF2F","#AFEEEE","#B0C4DE","#B0E0E6","#B22222","#BA55D3","#BC8F8F","#BDB76B","#C0C0C0","#C71585","#CD5C5C",
"#CD853F","#D2691E","#D2B48C","#D3D3D3","#D8BFD8","#DA70D6","#DAA520","#DB7093","#DC143C","#DCDCDC","#DDA0DD","#DEB887","#E0FFFF","#E6E6FA","#E9967A","#EE82EE",
"#EEE8AA","#F08080","#F0E68C","#F0F8FF","#F0FFF0","#F0FFFF","#F4A460","#F5DEB3","#F5F5DC","#F5F5F5","#F5FFFA","#F8F8FF","#FA8072","#FAEBD7","#FAF0E6","#FAFAD2",
"#FDF5E6","#FF0000","#FF00FF","#FF00FF","#FF1493","#FF4500","#FF6347","#FF69B4","#FF7F50","#FF8C00","#FFA500","#FFB6C1","#FFC0CB","#FFD700","#FFDAB9",
"#FFDEAD","#FFE4B5","#FFE4C4","#FFE4E1","#FFEBCD","#FFEFD5","#FFF0F5","#FFF5EE","#FFF8DC","#FFFACD","#FFFAF0","#FFFAFA","#FFFF00","#FFFFE0","#FFFFF0","#FFFFFF"
)

my_color=c47[as.numeric(as.factor(bac.cor.phyla))]

# Plot the computed netwrok using the FR force-directed layout
#l <- layout_with_kk(bac.cor.ig)
svg("network_kk.svg", height=8, width=8)
#png("network_kk.png", height=1200, width=1200)
plot(bac.cor.ig, layout=l, vertex.color=my_color, vertex.size = 4, vertex.label = NA, rescale=T)
legend("left", legend=levels(as.factor(bac.cor.phyla)) , col = c47 , bty = "n", pch=20 , pt.cex = 1, cex = 0.5, text.col="black", horiz = FALSE, inset = c(0.1, 0.1))
dev.off()

# Network statistics
message("Edge and vertex count:")
ecount(bac.cor.ig)
vcount(bac.cor.ig)
message("Assortativity degree:")
assortativity_degree(bac.cor.ig, directed=F)





###############################################################################
## Modularity analysis to identify cliques of ASV in the networks
## Testing different clustering functions
# Identify cliques using greedy clustering methods
bac.cor.greedy <- cluster_fast_greedy(bac.cor.ig)
modularity(bac.cor.greedy)
V(bac.cor.ig)$color=bac.cor.greedy$membership

# Identify cliques using random walks
bac.cor.walk <- cluster_walktrap(bac.cor.ig)
modularity(bac.cor.walk)
V(bac.cor.ig)$color=bac.cor.walk$membership

# Identify cliques using label propagation
bac.cor.label <- cluster_label_prop(bac.cor.ig)
modularity(bac.cor.label)
V(bac.cor.ig)$color=bac.cor.label$membership

# Identify cliques using Louvain
bac.cor.louvain <- cluster_louvain(bac.cor.ig)
modularity(bac.cor.louvain)
V(bac.cor.ig)$color=bac.cor.louvain$membership

# Check cliques dimensions
sizes(bac.cor.greedy)
sizes(bac.cor.walk)
sizes(bac.cor.label)
sizes(bac.cor.louvain)

# Plot the Cliques
#svg("network_cliques.svg", height=8, width=8)
#png("network_cliques.png", height=1200, width=1200)
plot(bac.cor.ig, layout=l, col = bac.cor.louvain, vertex.size = 3, vertex.label = NA, rescale=T)
#dev.off()

# Remove unconnected nodes
unconnected_nodes <- which(degree(bac.cor.ig)==0)
bac.cor.ig <- delete.vertices(bac.cor.ig, unconnected_nodes)
l <- l[-unconnected_nodes,]


## Prepare the files for visualization in the Gephi software
# Make edges (note, columns must be titled "Source", and "Target") to import in Gephi
edges <- get.edgelist(bac.cor.ig)
colnames(edges)<-c("Source","Target")
write.csv(edges, file="edges.csv", row.names=F)

# Make node table to import in Gephi
nodes <- get.vertex.attribute(bac.cor.ig)
nodes <- nodes$name
nodes <- cbind(nodes, bac.cor.phyla, bac.cor.phyla, bac.cor.louvain$membership)
colnames(nodes)<-c("Id","Label", "Phyla", "Clique")
write.csv(nodes, file="nodes.csv", row.names=F)

### Extracting the list of ASV for the main cliques identified by the Louvain function
V(bac.cor.ig)[bac.cor.louvain$membership == 5] -> bac.cor.group5
V(bac.cor.ig)[bac.cor.louvain$membership == 17] -> bac.cor.group17
V(bac.cor.ig)[bac.cor.louvain$membership == 19] -> bac.cor.group19
V(bac.cor.ig)[bac.cor.louvain$membership == 28] -> bac.cor.group28
V(bac.cor.ig)[bac.cor.louvain$membership == 35] -> bac.cor.group35
V(bac.cor.ig)[bac.cor.louvain$membership == 36] -> bac.cor.group36
V(bac.cor.ig)[bac.cor.louvain$membership == 37] -> bac.cor.group37
V(bac.cor.ig)[bac.cor.louvain$membership == 40] -> bac.cor.group40
V(bac.cor.ig)[bac.cor.louvain$membership == 42] -> bac.cor.group42
V(bac.cor.ig)[bac.cor.louvain$membership == 44] -> bac.cor.group44

## Writing the file to a CSV for parsing outside R
## Edit the file to remove the first column and row and to parse out the taxonomy identifier after the : and the OTU- start of each row
write.csv(names(bac.cor.group5), "group5_names_species.csv")
write.csv(names(bac.cor.group17), "group17_names_species.csv")
write.csv(names(bac.cor.group19), "group19_names_species.csv")
write.csv(names(bac.cor.group28), "group28_names_species.csv")
write.csv(names(bac.cor.group35), "group35_names_species.csv")
write.csv(names(bac.cor.group36), "group36_names_species.csv")
write.csv(names(bac.cor.group37), "group37_names_species.csv")
write.csv(names(bac.cor.group40), "group40_names_species.csv")
write.csv(names(bac.cor.group42), "group42_names_species.csv")
write.csv(names(bac.cor.group44), "group44_names_species.csv")

# Read back the file if you have parsed it outside R
group5_names <- as.matrix(read.csv("group5_names_species.csv", header =F))
group17_names <- as.matrix(read.csv("group17_names_species.csv", header =F))
group19_names <- as.matrix(read.csv("group19_names_species.csv", header =F))
group28_names <- as.matrix(read.csv("group28_names_species.csv", header =F))
group35_names <- as.matrix(read.csv("group35_names_species.csv", header =F))
group36_names <- as.matrix(read.csv("group36_names_species.csv", header =F))
group37_names <- as.matrix(read.csv("group37_names_species.csv", header =F))
group40_names <- as.matrix(read.csv("group40_names_species.csv", header =F))
group42_names <- as.matrix(read.csv("group42_names_species.csv", header =F))
group44_names <- as.matrix(read.csv("group44_names_species.csv", header =F))

bac_g5 <- subset_taxa(bac_data.pf, rownames(tax_table(bac_data.pf)) %in% as.vector(group5_names))
bac_g17 <- subset_taxa(bac_data.pf, rownames(tax_table(bac_data.pf)) %in% as.vector(group17_names))
bac_g19 <- subset_taxa(bac_data.pf, rownames(tax_table(bac_data.pf)) %in% as.vector(group19_names))
bac_g28 <- subset_taxa(bac_data.pf, rownames(tax_table(bac_data.pf)) %in% as.vector(group28_names))
bac_g35 <- subset_taxa(bac_data.pf, rownames(tax_table(bac_data.pf)) %in% as.vector(group35_names))
bac_g36 <- subset_taxa(bac_data.pf, rownames(tax_table(bac_data.pf)) %in% as.vector(group36_names))
bac_g37 <- subset_taxa(bac_data.pf, rownames(tax_table(bac_data.pf)) %in% as.vector(group37_names))
bac_g40 <- subset_taxa(bac_data.pf, rownames(tax_table(bac_data.pf)) %in% as.vector(group40_names))
bac_g42 <- subset_taxa(bac_data.pf, rownames(tax_table(bac_data.pf)) %in% as.vector(group42_names))
bac_g44 <- subset_taxa(bac_data.pf, rownames(tax_table(bac_data.pf)) %in% as.vector(group44_names))


###############################################################################
## Test each clique abundance against environmental predictors using Person, Spearman and scatterplots
# Make dataframe for clique abundance
sample_data_frame <- data.frame(sample_data(bac_data.pf))
sample_data_frame <- sample_data_frame[c(-1, -3, -18, -20, -24, -25, -26, -46, -48)] # remove unwanted variables
print(names(sample_data_frame))

## Starting with cluster 17 (Clique 1 in the paper) and progressing up to Clique 10
## Only variables with a p<0.01 are printed to screen
# Scatterplots
par(mfrow=c(3,3))
for (i in 1:length(sample_data_frame)) {
    plot(sample_data_frame[,i], log(colSums(otu_table(bac_g17))), xlab=colnames(sample_data_frame)[i])
}

message("Test with Pearson correlation:")
# Pearson
for (i in 16:length(sample_data_frame)) {
    a <- cor.test(sample_data_frame[,i], colSums(otu_table(bac_g17)))
       if (a$p.value<0.01) {
           print(paste(i,colnames(sample_data_frame)[i],a$estimate, a$parameter, a$p.value))
       }
}
message("Test with Spearman correlation:")
# Spearman
for (i in 16:length(sample_data_frame)) {
    a <- suppressWarnings(cor.test(sample_data_frame[,i], colSums(otu_table(bac_g17)), method = "spearman"))
       if (a$p.value<0.01) {
           print(paste(i,colnames(sample_data_frame)[i],a$estimate, a$parameter, a$p.value))
       }
}

#Clique 2
message("Test with Pearson correlation:")
for (i in 16:length(sample_data_frame)) {
    a <- cor.test(sample_data_frame[,i], colSums(otu_table(bac_g36)))
       if (a$p.value<0.01) {
           print(paste(i,colnames(sample_data_frame)[i],a$estimate, a$parameter, a$p.value))
       }
}
message("Test with Spearman correlation:")
# Spearman
for (i in 16:length(sample_data_frame)) {
    a <- suppressWarnings(cor.test(sample_data_frame[,i], colSums(otu_table(bac_g36)), method = "spearman"))
       if (a$p.value<0.01) {
           print(paste(i,colnames(sample_data_frame)[i],a$estimate, a$parameter, a$p.value))
       }
}
# Scatterplots
par(mfrow=c(3,3))
for (i in 1:length(sample_data_frame)) {
    plot(sample_data_frame[,i], log(colSums(otu_table(bac_g36))), xlab=colnames(sample_data_frame)[i])
}

# Clique 3
message("Test with Pearson correlation:")
for (i in 16:length(sample_data_frame)) {
    a <- cor.test(sample_data_frame[,i], colSums(otu_table(bac_g19)))
       if (a$p.value<0.05) {
           print(paste(i,colnames(sample_data_frame)[i],a$estimate, a$parameter, a$p.value))
       }
}
message("Test with Spearman correlation:")
# Spearman
for (i in 16:length(sample_data_frame)) {
    a <- suppressWarnings(cor.test(sample_data_frame[,i], colSums(otu_table(bac_g19)), method = "spearman"))
       if (a$p.value<0.01) {
           print(paste(i,colnames(sample_data_frame)[i],a$estimate, a$parameter, a$p.value))
       }
}
# Scatterplots
par(mfrow=c(3,3))
for (i in 1:length(sample_data_frame)) {
    plot(sample_data_frame[,i], log(colSums(otu_table(bac_g19))), xlab=colnames(sample_data_frame)[i])
}

#Clique 4
# Scatterplots
par(mfrow=c(3,3))
for (i in 1:length(sample_data_frame)) {
    plot(sample_data_frame[,i], log(colSums(otu_table(bac_g5))), xlab=colnames(sample_data_frame)[i])
}
message("Test with Pearson correlation:")
for (i in 16:length(sample_data_frame)) {
    a <- cor.test(sample_data_frame[,i], colSums(otu_table(bac_g5)))
       if (a$p.value<0.05) {
           print(paste(i,colnames(sample_data_frame)[i],a$estimate, a$parameter, a$p.value))
       }
}
message("Test with Spearman correlation:")
# Spearman
for (i in 16:length(sample_data_frame)) {
    a <- suppressWarnings(cor.test(sample_data_frame[,i], colSums(otu_table(bac_g5)), method = "spearman"))
       if (a$p.value<0.01) {
           print(paste(i,colnames(sample_data_frame)[i],a$estimate, a$parameter, a$p.value))
       }
}

#Clique 5
message("Test with Pearson correlation:")
for (i in 16:length(sample_data_frame)) {
    a <- cor.test(sample_data_frame[,i], colSums(otu_table(bac_g44)))
       if (a$p.value<0.05) {
           print(paste(i,colnames(sample_data_frame)[i],a$estimate, a$parameter, a$p.value))
       }
}
message("Test with Spearman correlation:")
# Spearman
for (i in 16:length(sample_data_frame)) {
    a <- suppressWarnings(cor.test(sample_data_frame[,i], colSums(otu_table(bac_g44)), method = "spearman"))
       if (a$p.value<0.01) {
           print(paste(i,colnames(sample_data_frame)[i],a$estimate, a$parameter, a$p.value))
       }
}
# Scatterplots
par(mfrow=c(3,3))
for (i in 1:length(sample_data_frame)) {
    plot(sample_data_frame[,i], log(colSums(otu_table(bac_g44))), xlab=colnames(sample_data_frame)[i])
}

#Clique 6
message("Test with Pearson correlation:")
for (i in 16:length(sample_data_frame)) {
    a <- cor.test(sample_data_frame[,i], colSums(otu_table(bac_g37)))
       if (a$p.value<0.05) {
           print(paste(i,colnames(sample_data_frame)[i],a$estimate, a$parameter, a$p.value))
       }
}
message("Test with Spearman correlation:")
# Spearman
for (i in 16:length(sample_data_frame)) {
    a <- suppressWarnings(cor.test(sample_data_frame[,i], colSums(otu_table(bac_g37)), method = "spearman"))
       if (a$p.value<0.01) {
           print(paste(i,colnames(sample_data_frame)[i],a$estimate, a$parameter, a$p.value))
       }
}
# Scatterplots
par(mfrow=c(3,3))
for (i in 1:length(sample_data_frame)) {
    plot(sample_data_frame[,i], log(colSums(otu_table(bac_g37))), xlab=colnames(sample_data_frame)[i])
}

# Clique 7
message("Test with Pearson correlation:")
for (i in 16:length(sample_data_frame)) {
    a <- cor.test(sample_data_frame[,i], colSums(otu_table(bac_g42)))
       if (a$p.value<0.05) {
           print(paste(i,colnames(sample_data_frame)[i],a$estimate, a$parameter, a$p.value))
       }
}
message("Test with Spearman correlation:")
# Spearman
for (i in 16:length(sample_data_frame)) {
    a <- suppressWarnings(cor.test(sample_data_frame[,i], colSums(otu_table(bac_g42)), method = "spearman"))
       if (a$p.value<0.01) {
           print(paste(i,colnames(sample_data_frame)[i],a$estimate, a$parameter, a$p.value))
       }
}
# Scatterplots
par(mfrow=c(3,3))
for (i in 1:length(sample_data_frame)) {
    plot(sample_data_frame[,i], log(colSums(otu_table(bac_g42))), xlab=colnames(sample_data_frame)[i])
}

#Clique 8
message("Test with Pearson correlation:")
for (i in 16:length(sample_data_frame)) {
    a <- cor.test(sample_data_frame[,i], colSums(otu_table(bac_g28)))
       if (a$p.value<0.05) {
           print(paste(i,colnames(sample_data_frame)[i],a$estimate, a$parameter, a$p.value))
       }
}
message("Test with Spearman correlation:")
# Spearman
for (i in 16:length(sample_data_frame)) {
    a <- suppressWarnings(cor.test(sample_data_frame[,i], colSums(otu_table(bac_g28)), method = "spearman"))
       if (a$p.value<0.01) {
           print(paste(i,colnames(sample_data_frame)[i],a$estimate, a$parameter, a$p.value))
       }
}
# Scatterplots
par(mfrow=c(3,3))
for (i in 1:length(sample_data_frame)) {
    plot(sample_data_frame[,i], log(colSums(otu_table(bac_g28))), xlab=colnames(sample_data_frame)[i])
}

#Clique 9
message("Test with Pearson correlation:")
for (i in 16:length(sample_data_frame)) {
    a <- cor.test(sample_data_frame[,i], colSums(otu_table(bac_g40)))
       if (a$p.value<0.01) {
           print(paste(i,colnames(sample_data_frame)[i],a$estimate, a$parameter, a$p.value))
       }
}
message("Test with Spearman correlation:")
# Spearman
for (i in 16:length(sample_data_frame)) {
    a <- suppressWarnings(cor.test(sample_data_frame[,i], colSums(otu_table(bac_g40)), method = "spearman"))
       if (a$p.value<0.01) {
           print(paste(i,colnames(sample_data_frame)[i],a$estimate, a$parameter, a$p.value))
       }
}
# Scatterplots
par(mfrow=c(3,3))
for (i in 1:length(sample_data_frame)) {
    plot(sample_data_frame[,i], log(colSums(otu_table(bac_g40))), xlab=colnames(sample_data_frame)[i])
}

#Clique 10
message("Test with Pearson correlation:")
for (i in 16:length(sample_data_frame)) {
    a <- cor.test(sample_data_frame[,i], colSums(otu_table(bac_g35)))
       if (a$p.value<0.05) {
           print(paste(i,colnames(sample_data_frame)[i],a$estimate, a$parameter, a$p.value))
       }
}
message("Test with Spearman correlation:")
# Spearman
for (i in 16:length(sample_data_frame)) {
    a <- suppressWarnings(cor.test(sample_data_frame[,i], colSums(otu_table(bac_g35)), method = "spearman"))
       if (a$p.value<0.01) {
           print(paste(i,colnames(sample_data_frame)[i],a$estimate, a$parameter, a$p.value))
       }
}
# Scatterplots
par(mfrow=c(3,3))
for (i in 1:length(sample_data_frame)) {
    plot(sample_data_frame[,i], log(colSums(otu_table(bac_g35))), xlab=colnames(sample_data_frame)[i])
}

################################################################################
## Random Forests approach to selection of important variable predictors on the clique identified in the co-occurrence
## network analysis

# Extract environmental predictor and remove non numeric variables (i.e. factors) and unwanted variables (air temperature)
bac_data.pf_pred <- data.frame(sample_data(bac_data.pf))

# Contruct the vetors for each clique abundance response variable and normalize them as Z-scores
clique1 <- scale(data.frame(colSums(otu_table(bac_g17))), center = TRUE, scale = TRUE)
clique2 <- scale(data.frame(colSums(otu_table(bac_g36))), center = TRUE, scale = TRUE)
clique3 <- scale(data.frame(colSums(otu_table(bac_g19))), center = TRUE, scale = TRUE)
clique4 <- scale(data.frame(colSums(otu_table(bac_g5))), center = TRUE, scale = TRUE)
clique5 <- scale(data.frame(colSums(otu_table(bac_g44))), center = TRUE, scale = TRUE)
clique6 <- scale(data.frame(colSums(otu_table(bac_g37))), center = TRUE, scale = TRUE)
clique7 <- scale(data.frame(colSums(otu_table(bac_g42))), center = TRUE, scale = TRUE)
clique8 <- scale(data.frame(colSums(otu_table(bac_g28))), center = TRUE, scale = TRUE)
clique9 <- scale(data.frame(colSums(otu_table(bac_g40))), center = TRUE, scale = TRUE)
clique10 <- scale(data.frame(colSums(otu_table(bac_g35))), center = TRUE, scale = TRUE)

#Check if the predictors contain NAs
if (grep("NA", bac_data.pf_pred)!=0){
  print("There are NAs! Need to impute variables")
} else {
  print("No NAs in the predictors, good to go directly to the Random Forests")
} # The predictors dataset c ontains NAs

# Impute missing data using a iterative random forest approach with the missForest package
bac_data.pf_pred2 <- missForest(bac_data.pf_pred, maxiter=10, variablewise = TRUE, verbose = TRUE, ntree=5000)
bac_data.pf_pred2$OOBerror
bac_data.pf_pred2$ximp # dataframe with imputed missing values to be used in downstrewam analysis

# Run the the Random Forests model for variable selection using VSURF
clique1_RF <- VSURF(x = data.frame(bac_data.pf_pred2$ximp), y = clique1[,1], parallel=T, ntree=2000, mtry=23)
summary(clique1_RF)
plot(clique1_RF)
colnames(bac_data.pf_pred2$ximp)[clique1_RF$varselect.thres]
colnames(bac_data.pf_pred2$ximp)[clique1_RF$varselect.interp]

# Repeat for each cliques
clique2_RF <- VSURF(x = data.frame(bac_data.pf_pred2$ximp), y = clique2[,1], parallel=T, ntree=2000, mtry=23)
plot(clique2_RF)
colnames(bac_data.pf_pred2$ximp)[clique2_RF$varselect.thres]
colnames(bac_data.pf_pred2$ximp)[clique2_RF$varselect.interp]

clique3_RF <- VSURF(x = data.frame(bac_data.pf_pred2$ximp), y = clique3[,1], parallel=T, ntree=2000, mtry=23)
plot(clique3_RF)
colnames(bac_data.pf_pred2$ximp)[clique3_RF$varselect.thres]
colnames(bac_data.pf_pred2$ximp)[clique3_RF$varselect.interp]

clique4_RF <- VSURF(x = data.frame(bac_data.pf_pred2$ximp), y = clique4[,1], parallel=T, ntree=2000, mtry=23)
plot(clique4_RF)
colnames(bac_data.pf_pred2$ximp)[clique4_RF$varselect.thres]
colnames(bac_data.pf_pred2$ximp)[clique4_RF$varselect.interp]

clique5_RF <- VSURF(x = data.frame(bac_data.pf_pred2$ximp), y = clique5[,1], parallel=T, ntree=2000, mtry=23)
plot(clique5_RF)
colnames(bac_data.pf_pred2$ximp)[clique5_RF$varselect.thres]
colnames(bac_data.pf_pred2$ximp)[clique5_RF$varselect.interp]

clique6_RF <- VSURF(x = data.frame(bac_data.pf_pred2$ximp), y = clique6[,1], parallel=T, ntree=2000, mtry=23)
plot(clique6_RF)
colnames(bac_data.pf_pred2$ximp)[clique6_RF$varselect.thres]
colnames(bac_data.pf_pred2$ximp)[clique6_RF$varselect.interp]

clique7_RF <- VSURF(x = data.frame(bac_data.pf_pred2$ximp), y = clique7[,1], parallel=T, ntree=2000, mtry=23)
plot(clique7_RF)
colnames(bac_data.pf_pred2$ximp)[clique7_RF$varselect.thres]
colnames(bac_data.pf_pred2$ximp)[clique7_RF$varselect.interp]

clique8_RF <- VSURF(x = data.frame(bac_data.pf_pred2$ximp), y = clique8[,1], parallel=T, ntree=2000, mtry=23)
plot(clique8_RF)
colnames(bac_data.pf_pred2$ximp)[clique8_RF$varselect.thres]
colnames(bac_data.pf_pred2$ximp)[clique8_RF$varselect.interp]

clique9_RF <- VSURF(x = data.frame(bac_data.pf_pred2$ximp), y = clique9[,1], parallel=T, ntree=2000, mtry=23)
plot(clique9_RF)
colnames(bac_data.pf_pred2$ximp)[clique9_RF$varselect.thres]
colnames(bac_data.pf_pred2$ximp)[clique9_RF$varselect.interp]

clique10_RF <- VSURF(x = data.frame(bac_data.pf_pred2$ximp), y = clique10[,1], parallel=T, ntree=2000, mtry=23)
plot(clique10_RF)
colnames(bac_data.pf_pred2$ximp)[clique10_RF$varselect.thres]
colnames(bac_data.pf_pred2$ximp)[clique10_RF$varselect.interp]

save.image() # Save all this good work!

## Figures have been edited and/or assembled in Inkscape

################################################
## END OF THE ANALYSIS
################################################
