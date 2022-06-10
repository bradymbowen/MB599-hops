
## Used to instal phyloseq
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("phyloseq")

# load required packages
library(phyloseq)
library(dplyr)
library(magrittr)
library(Hmisc)
library(randomForest)
library(ANCOMBC)
library(microbiome)
library(ape)
library(ggpubr)
library(ggfortify)
library(gridExtra)

# Read in phyloseq object
psdata <- readRDS("./mb599data.Rds")

# To view each dataframe within the phyloseq object
# View(as.data.frame(otu_table(psdata)))
# View(as.data.frame(tax_table(psdata)))
# View(as.matrix.data.frame(sample_data(psdata))


# Factorize and set the levels on the metadata
week_levels <- c("1", "2", "3", "4", "5")
sample_data(psdata)$Week = factor(week_levels)


ids <- sample_data(psdata)$Participant.ID
participant_levels <- c("101", "102", "103", "104", "105", "106", "107", 
                        "109", "110", "111", "112", "113", "114", "115", 
                        "116", "117", "118", "119", "120", "121", "122", 
                        "123", "124" ,"125", "126", "127", "128", "129", 
                        "130", "131")
sample_data(psdata)$Participant.ID <-  factor(ids, levels = participant_levels)


group_levels <- c("treatment", "control")
group <- get_variable(psdata, "Group")
sample_data(psdata)$Group = factor(group, levels = group_levels)


# Add Metabotypes for DXN and 8PN
sample_data(psdata)$metabo_dxn <- ifelse(sample_data(psdata)$DXN > 2292, "high", "low")
sample_data(psdata)$metabo_8pn <- ifelse(sample_data(psdata)$X8PN > 1910, "high", "low")
# Facotize DXN and 8PN metabotypes
metabotypes <- c("high", "low")
sample_data(psdata)$metabo_dxn <- factor(sample_data(psdata)$metabo_dxn, levels = metabotypes)
sample_data(psdata)$metabo_8pn <- factor(sample_data(psdata)$metabo_8pn, levels = metabotypes)


# Agglomerate to the genus level
ps_genera <- psdata %>% tax_glom(taxrank = "Genus")
# Remove taxa not seen more than 3 times in at least 20% of the samples
ps_counts <- ps_genera %>% filter_taxa(function(x) sum(x > 3) > (0.2*length(x)), TRUE)
# Convert from counts to relative abundance
ps_relab <- ps_counts %>% transform_sample_counts(function(x) x / sum(x))
# Filter out low abundance (>1e-5) taxa
ps <- ps_relab %>% filter_taxa(function(x) mean(x) > 1e-5, TRUE)




#Extract relevant sub-data:
ps_treat <- ps %>% subset_samples(Group == "treament")
ps_wk_1_2 <- ps %>% subset_samples(Week %in% c("1", "2"))
weeks_after_base <- ps %>% subset_samples(Week %in% c("2", "3", "4", "5"))




############################# PCoA Analysis ###################################




# PCoA analysis



# Extract relevant sub-data
pswk1 <- ps %>% subset_samples(Week == '1')
pswk1C <- pswk1 %>% subset_samples(Group == 'control')
pswk2 <- ps %>% subset_samples(Week == '2')
pswk3 <- ps %>% subset_samples(Week == '3')
pswk3C <- pswk3 %>% subset_samples(Group == 'control')
pswk3T <- pswk3 %>% subset_samples(Group == 'treatment')
pswk4 <- ps %>% subset_samples(Week == '4')
pswk5 <- ps %>% subset_samples(Week == '5')

# Calculate eigenvalues by Bray Curtis distances
ord1 <- ordinate(pswk1, method = 'PCoA', distance = 'bray')
ord2 <- ordinate(pswk2, method = 'PCoA', distance = 'bray')
ord3 <- ordinate(pswk3, method = 'PCoA', distance = 'bray')
ord3C <- ordinate(pswk3C, method = 'PCoA', distance = 'bray')
ord3T <- ordinate(pswk3T, method = 'PCoA', distance = 'bray')
ord4 <- ordinate(pswk4, method = 'PCoA', distance = 'bray')
ord5 <- ordinate(pswk5, method = 'PCoA', distance = 'bray')

# Plot PCoA
pcoawk1 <- plot_ordination(pswk1, ord1, color = 'Group') +
  geom_point(aes(),
             size = 3) +
  ggtitle('Beta Diversity of Fecal Samples in Week 0')

pcoawk2 <- plot_ordination(pswk2, ord2, color = 'Group') +
  geom_point(aes(),
             size = 3) +
  ggtitle('Beta Diversity of Fecal Samples in Week 2')

pcoawk3 <- plot_ordination(pswk3, ord3, color = 'Group') +
  geom_point(aes(),
             size = 3) +
  ggtitle('Beta Diversity of Fecal Samples in Week 4')

pcoawk3C <- plot_ordination(pswk3C, ord3C, color = 'Group') +
  geom_point(aes(),
             size = 3) +
  ggtitle('Beta Diversity of Fecal Samples in Week 4, Control Subsample')

pcoawk3T <- plot_ordination(pswk3T, ord3T, color = 'Group') +
  geom_point(aes(),
             size = 3) +
  ggtitle('Beta Diversity of Fecal Samples in Week 4, Treatment Subsample')

pcoawk4 <- plot_ordination(pswk4, ord4, color = 'Group') +
  geom_point(aes(),
             size = 3) +
  ggtitle('Beta Diversity of Fecal Samples in Week 6')

pcoawk5 <- plot_ordination(pswk5, ord5, color = 'Group') +
  geom_point(aes(),
             size = 3) +
  ggtitle('Beta Diversity of Fecal Samples in Week 8')

# Comparison of the components of the top two axes across the 5 measurements
ordTop = cbind(ord1$vectors[,1],ord2$vectors[,1],ord3$vectors[,1],ord4$vectors[,1],ord5$vectors[,1])
colnames(ordTop) = c("Week 1", "Week 2", "Week 3", "Week 4", "Week 5")

barplot(t(ordTop),
        main="Differences in the Components of the Axis with First Largest Eigenvalue",
        legend=colnames(ordTop),
        beside=TRUE,
        las=2)

ordSecond = cbind(ord1$vectors[,2],ord2$vectors[,2],ord3$vectors[,2],ord4$vectors[,2],ord5$vectors[,2])
colnames(ordSecond) = c("Week 1", "Week 2", "Week 3", "Week 4", "Week 5")

barplot(t(ordSecond),
        main="Differences in the Components of the Axis with Second Largest Eigenvalue",
        legend=colnames(ordSecond),
        beside=TRUE,
        las=2,
        args.legend = list(x = 'bottomright'))

# Comparison of the components of the largest eigenvalue of the two separated groups of control and treatment in week 3

# Work In Progress
# ord3TopCT = cbind(ord3C$vectors[,1],ord3T$vectors[,1])
# colnames(ord3TopCT) = c("Control", "Treatment")
# 
# barplot(t(ord3TopCT),
#         main="Control vs Treatment PCoA, Vector with the Largest Eigenvalue, Week 3",
#         legend=colnames(ord3TopCT),
#         beside=TRUE,
#         las=2,
#         args.legend = list(x = 'bottomright'))









############################# PCA Analysis ###################################


# PCA analysis


# Extract relevant sub-data (need data generated in PCoA)
matrixwk1 <- cbind(sample_data(pswk1)[,6:10],otu_table(pswk1))
matrixwk2 <- cbind(sample_data(pswk2)[,6:10],otu_table(pswk2))
matrixwk3 <- cbind(sample_data(pswk3)[,6:10],otu_table(pswk3))
matrixwk4 <- cbind(sample_data(pswk4)[,6:10],otu_table(pswk4))
matrixwk5 <- cbind(sample_data(pswk5)[,6:10],otu_table(pswk5))
# matrixwk1C <- cbind(sample_data(pswk1C)[,6:10],otu_table(pswk1C))
# matrixwk3C <- cbind(sample_data(pswk3C)[,6:10],otu_table(pswk3C))


# Calculate eigenvectors in respect to ASV's and X6PN
pcawk1 <- prcomp(matrixwk1,scale = TRUE)
pcawk2 <- prcomp(matrixwk2,scale = TRUE)
pcawk3 <- prcomp(matrixwk3,scale = TRUE)
pcawk4 <- prcomp(matrixwk4,scale = TRUE)
pcawk5 <- prcomp(matrixwk5,scale = TRUE)
# pcawk1C <- prcomp(matrixwk1C,scale = TRUE)
# pcawk3C <- prcomp(matrixwk3C,scale = TRUE)


# Calculate variance captured by each eigenvector
var1=pcawk1$sdev^2/sum(pcawk1$sdev^2)
var2=pcawk2$sdev^2/sum(pcawk2$sdev^2)
var3=pcawk3$sdev^2/sum(pcawk3$sdev^2)
var4=pcawk4$sdev^2/sum(pcawk4$sdev^2)
var5=pcawk5$sdev^2/sum(pcawk5$sdev^2)
# var1C=pcawk1C$sdev^2/sum(pcawk1C$sdev^2)
# var3C=pcawk3C$sdev^2/sum(pcawk3C$sdev^2)

# Plotting variance drop-off for each week (none are that good)
pcawk1plot <- qplot(c(1:27), var1, main="Week 0", xlab = "Axis number", ylab = "Variance explained")
pcawk2plot <- qplot(c(1:27), var2, main="Week 2", xlab = "Axis number", ylab = "Variance explained")
pcawk3plot <- qplot(c(1:27), var3, main="Week 4", xlab = "Axis number", ylab = "Variance explained")
pcawk4plot <- qplot(c(1:27), var4, main="Week 6", xlab = "Axis number", ylab = "Variance explained")
pcawk5plot <- qplot(c(1:27), var5, main="Week 8", xlab = "Axis number", ylab = "Variance explained")
grid.arrange(pcawk2plot,pcawk3plot,pcawk4plot,pcawk5plot,nrow=2
             ,top=textGrob("Variance explained by each principle component axis",gp=gpar(fontsize=18)))
# qplot(c(1:length(var1C)), var1C)
# qplot(c(1:length(var3C)), var3C)

# Following code was an attempt to PCA on the control and then add in the treatment patients to the plot,
# it doesn't work for some reason, instead giving vastly different values
#
# testplotdata1 = as.matrix(matrixwk1) %*% pcawk1C$rotation
# plotwk1 <- ggplot(testplotdata1,aes(x=PC1, y=PC2)) + geom_point()
# plotwk1 + labs(colour = "Group")
# 
# testplotdata3 = as.matrix(matrixwk3) %*% pcawk3C$rotation
# plotwk3 <- ggplot(testplotdata,aes(x=PC10, y=PC2)) + geom_point()
# plotwk3 + labs(colour = "Group")

# Plotting top two eigenvector axes with data on control vs treatment
# (maybe other axes would be able to separate these two groups more?)
autoplot(pcawk1,data=sample_data(pswk1), colour = "Group")
ggplot(cbind(pcawk1$x,sample_data(pswk1)[,3]),aes(x=PC8, y=PC12, color=Group)) + geom_point()
autoplot(pcawk2,data=sample_data(pswk2), colour = "Group")
autoplot(pcawk3,data=sample_data(pswk3), colour = "Group")
ggplot(cbind(pcawk3$x,sample_data(pswk3)[,3]),aes(x=PC7, y=PC11, color=Group)) + geom_point()
#autoplot(pcawk3C,data=sample_data(pswk3C), colour = "Group")
autoplot(pcawk4,data=sample_data(pswk4), colour = "Group")
autoplot(pcawk5,data=sample_data(pswk5), colour = "Group")