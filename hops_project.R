
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
library(microbiome)
library(ape)
library(ggpubr)
library(ggfortify)

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





###################### Alpha and Beta Diversity Analysis#################################




# Alpha Diversity Metrics

adiv <- psdata %>%
  estimate_richness(measures = c('Observed', 'Shannon', 'Simpson')) %>% 
  cbind(., data.frame(sample_data(psdata))) %>% 
  mutate(wkc = ifelse(Week == 1, "Baseline",
                      ifelse(Week == 2, "Week 2", 
                             ifelse(Week == 3, 'Week 4', 
                                    ifelse(Week == 4, 'Week 6', 'Week 8'))))) %>% 
  rename_with(tolower) %>% 
  dplyr::select('observed', 'shannon', 'simpson', 'wkc', 'week', 'participant.id', 'group')



# FIGURE: 
# Plot alpha-diversity metrics by week:

adivplot <- ggplot(adiv, aes(x = group, y = shannon, color = group)) +
  geom_point(size = 1) +
  theme_bw() +
  facet_wrap(~wkc, scales = 'free', ncol = 5)+
  theme(axis.text.x = element_text(angle = 315), 
        axis.title.x = element_blank(),
        legend.title = element_blank()) +
  ylab("Shannon") +
  ggtitle('Alpha Diversity Metrics')




# Significance testing of alpha-diversity metrics: 

adivstats <- adiv %>% 
  pivot_longer(c('observed', 'shannon', 'simpson'), 
               names_to = 'measure', values_to = 'value') %>% 
  group_by(week, measure) %>% 
  nest() %>% 
  mutate(wrs = map(data, function(x) wilcox.test(value ~ group, data = x, paired = FALSE))) %>% 
  mutate(p_val = map_dbl(wrs, function(x) x$p.value))





# Beta Diversity Metrics 




# PCoA of treatment and placebo groups

# calculate ordications by Bray Curtis distance in treatment and controls 
group_ord <- ordinate(ps, method = "PCoA", distance = "bray")
# plot coor 1 and 2 for pcoa (color = treatment / control)
pcoa_group <- plot_ordination(ps, group_ord, color = 'Group') +
  facet_wrap('Week')




#PCoA of treatment group

# subset samples for treatment group excluding baseline measurement
pstreat_new <- ps %>% 
  subset_samples(Group %in% 'treatment')

# calculate ordination by Bray Curtis distance
treat_ord <- ordinate(pstreat_new, method = "PCoA", distance = 'bray')

# FIGURE: plot PCoA of treatment group by week 
pcoatre <- plot_ordination(pstreat_new, treat_ord, color = "Participant.ID") +
  facet_wrap(~Week)



# PCoA of treatment: component 1 by DXN concentration
pcoadf <- plot_ordination(pstreat_new, treat_ord, color = "Participant.ID", justDF = TRUE)

# FIGURE: axis 1 by DXN conc
ggplot(pcoadf, aes(x = Axis.1, y = DXN, color = Participant.ID)) +
  geom_point()+
  facet_wrap(~Week)







########################### Random Forest Analysis ######################################



#create data partition

set.seed(1)
index_train <- createDataPartition(psdata@sam_data$Group, p = 0.7)[[1]]
train <- psdata@otu_table[index_train, ]
test <- psdata@otu_table[-index_train, ]

#spliting the phyloseq objects
pstrain <- phyloseq(otu_table(train, taxa_are_rows = FALSE), psdata@sam_data[index_train, ])
pstest <- phyloseq(otu_table(test, taxa_are_rows = FALSE), psdata@sam_data[-index_train, ])




#Find optimum mtry

set.seed(1)
datatrain = data.frame(train)
datatrain$sample_group = pstrain@sam_data$Group
control <- trainControl(method='repeatedcv', 
                        number=3, 
                        repeats=3,
                        allowParallel = F)

tunegrid <- expand.grid(.mtry=c(1:20)) 
rf <- train(sample_group ~., 
            data= datatrain, 
            method='rf', 
            metric='Accuracy', 
            tuneGrid=tunegrid, 
            trControl=control)
print(rf)



#Testing model performance                               
}
mtry_best = as.numeric(rf$bestTune)
model = randomForest(train, y = as.factor(pstrain@sam_data$Group), mtry = mtry_best)

model 


#Performance on test set                              

preds = predict(model, test)
print(paste("Accuracy: ", sum(preds == as.factor(pstest@sam_data$Group)) / nsamples(pstest)))

#Visualize on test dataset
ord <- ordinate(pstest, method = "PCoA", distance = 'euclidean')
pstest@sam_data$rf_predictions = predict(model, pstest@otu_table)
plot_ordination(pstest, ord, 'samples', color = 'Group', shape = 'rf_predictions') + geom_point(size = 4) +
  ggtitle("Microbiome classification between groups based on random forest algoritm")



#obtain the important taxa

par(mfrow=c(1,2))
model = randomForest(pstest@otu_table, y = as.factor(pstest@sam_data$Group), mtry = mtry_best)
varImpPlot(model, type = 2)
model = randomForest(pstest@otu_table, y = as.factor(pstest@sam_data$DXN), mtry = mtry_best)
varImpPlot(model, type = 2)




#Important variables

imp_list <- list()
for(i in 1:20){
  model = randomForest(pstest@otu_table, y = as.factor(pstest@sam_data$Group), mtry = mtry_best)
  imp_list[i] <- varImp(model)
}

imp_df <- do.call(rbind.data.frame, imp_list)
colnames(imp_df) <- colnames(train)
colMeans(imp_df)
barplot(sort(colMeans(imp_df)), horiz = T, las = 1, xlab = "Mean variable taxa")          
                               





########################## Correlation Analysis #########################################



# Spearman's correlation agglomerated to genus level



# Agglomerate by Genus
psgen <- pstreat %>% 
  tax_glom(taxrank = "Genus") %>% 
  psmelt()
# Manipulate dataframe
corrgen <- psgen %>% 
  rename('ASV' = 'OTU') %>% 
  rename_with(tolower) %>% 
  select(asv, sample, abundance, participant.id, week, ixn, x8pn, dxn, genus, family) %>% 
  pivot_longer(c(ixn, x8pn, dxn), names_to = 'metab', values_to = 'conc')
# Find Spearmans correlation with Benjamini-Hochberg correction for multiple testing
cbg <- corrgen %>% 
  group_by(metab, genus, week) %>% 
  nest() %>% 
  mutate(spr = map(data, function(x) cor.test(x$abundance, x$conc, method = 'spearman'))) %>% 
  mutate(rho = map_dbl(spr, function(x) x$estimate)) %>% 
  mutate(pvalue = map_dbl(spr, function(x) x$p.value)) %>% 
  mutate(padj = p.adjust(pvalue, method = 'BH'))
# Summarize results
cbg2 <- cbg %>% 
  filter(padj <= 0.05 & rho >= 0.68 & week %in% c('2', '3', '4', '5')) %>% 
  group_by(week, metab, genus) %>% 
  arrange(., .by_group = TRUE) %>% 
  select(-c(data, spr))




# FIGURE: Genus x Metabolite

# Change metabolite names and week names
cbg_alt <- cbg %>% 
  mutate(metab = ifelse(metab == 'ixn', 'IXN', 
                        ifelse(metab == 'x8pn', '8PN', 'DXN'))) %>% 
  mutate(week = ifelse(week == 1, 'Baseline', 
                       ifelse(week == 2, 'Week 2', 
                              ifelse(week == 3, 'Week 4', 
                                     ifelse(week == 4, 'Week 6', 'Week 8')))))


# plot heatmap
f3 <- ggplot(cbg_alt, aes(x = genus, y = metab, fill = rho)) +
  geom_raster() +
  #coord_fixed() +
  scale_fill_gradient2(high = 'blue', low = 'red', limits = c(-1,1)) +
  cowplot::theme_cowplot() +
  #facet_wrap(~week, ncol = 1)
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5),
        axis.text.y = element_text(size = 10, vjust = 0.8),
        axis.title.x = element_text(size = 18, vjust = 0),
        axis.title.y = element_text(size = 18, vjust = 1)) +
  xlab('Genus') +
  ylab('Metabolite') +
  ggtitle('Genus vs Metabolite')





# Spearmans correlation agglomerated to Family level 


# Agglomerate by Family
psfam <- pstreat %>% 
  tax_glom(taxrank = "Family") %>% 
  psmelt()

# Manipulate dataframe
corrfam <- psfam %>%   
  rename_with(tolower) %>% 
  dplyr::rename('asv' = 'otu') %>% 
  select(asv, sample, abundance, participant.id, week, family, ixn, x8pn, dxn) %>% 
  pivot_longer(c(ixn, x8pn, dxn), names_to = 'metab', values_to = 'conc')

# Find Spearman's correlation with Benjamini-Hochberg correction for multiple testing
cyb <- corrfam %>%
  rename_with(tolower) %>%
  group_by(metab, family, week) %>%
  nest() %>%
  #mutate(spr = map_dbl(data, sprmfun))
  mutate(spr = map(data, function(x) cor.test(x$abundance, x$conc, method = 'spearman'))) %>%
  mutate(rho = map_dbl(spr, function(x) x$estimate)) %>% 
  mutate(pvalue = map_dbl(spr, function(x) x$p.value)) %>% 
  ungroup() %>%
  mutate(padj = p.adjust(pvalue, method = "BH"))
# Filter and summarize results
cyb2 <- cyb %>% 
  filter(padj <= 0.1 & rho >= 0.5 ) %>% 
  group_by(week, metab, family) %>% 
  arrange(., .by_group = TRUE) %>% 
  select(-c(data, spr))




# FIGURE: family x metabolite


# change metabolite names and week names
cyb_alt <- cyb %>% 
  mutate(metab = ifelse(metab == 'ixn', 'IXN', 
                        ifelse(metab == 'x8pn', '8PN', 'DXN'))) %>% 
  mutate(week = ifelse(week == 1, 'Baseline', 
                       ifelse(week == 2, 'Week 2', 
                              ifelse(week == 3, 'Week 4', 
                                     ifelse(week == 4, 'Week 6', 'Week 8')))))

# plot figure
f2 <- ggplot(cyb_alt, aes(x = family, y = metab, fill = rho)) +
  geom_raster() +
  coord_fixed() +
  scale_fill_gradient2(high = 'blue', low = 'red', limits = c(-1,1)) +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5),
        axis.text.y = element_text(size = 10, vjust = 0.8),
        axis.title.x = element_text(size = 18, vjust = 0),
        axis.title.y = element_text(size = 18, vjust = 1)) +
  facet_wrap(~week, ncol = 1) +
  #theme_bw() +
  xlab('Family') +
  ylab('Metabolite') +
  ggtitle('Family vs Metabolite')






############################# PCoA Analysis ###################################




# PCoA analysis



# Extract relevant sub-data
pswk1 <- ps %>% subset_samples(Week == '1')
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
matrixwk1 <- cbind(sample_data(pswk1)[,7],otu_table(pswk1))
matrixwk2 <- cbind(sample_data(pswk2)[,7],otu_table(pswk2))
matrixwk3 <- cbind(sample_data(pswk3)[,7],otu_table(pswk3))
matrixwk4 <- cbind(sample_data(pswk4)[,7],otu_table(pswk4))
matrixwk5 <- cbind(sample_data(pswk5)[,7],otu_table(pswk5))


# Calculate eigenvectors in respect to ASV's and X6PN
pcawk1 <- prcomp(matrixwk1,scale = TRUE)
pcawk2 <- prcomp(matrixwk2,scale = TRUE)
pcawk3 <- prcomp(matrixwk3,scale = TRUE)
pcawk4 <- prcomp(matrixwk4,scale = TRUE)
pcawk5 <- prcomp(matrixwk5,scale = TRUE)


# Calculate variance captured by each eigenvector
var1=pcawk1$sdev^2/sum(pcawk1$sdev^2)
var2=pcawk2$sdev^2/sum(pcawk2$sdev^2)
var3=pcawk3$sdev^2/sum(pcawk3$sdev^2)
var4=pcawk4$sdev^2/sum(pcawk4$sdev^2)
var5=pcawk5$sdev^2/sum(pcawk5$sdev^2)

# Plotting variance drop-off for each week (none are that good)
qplot(c(1:27), var1)
qplot(c(1:27), var2)
qplot(c(1:27), var3)
qplot(c(1:27), var4)
qplot(c(1:27), var5)

# Plotting top two eigenvector axes with data on control vs treatment
# (maybe other axes would be able to separate these two groups more?)
autoplot(pcawk1,data=sample_data(pswk1), colour = "Group")
autoplot(pcawk2,data=sample_data(pswk1), colour = "Group")
autoplot(pcawk3,data=sample_data(pswk1), colour = "Group")
autoplot(pcawk4,data=sample_data(pswk1), colour = "Group")
autoplot(pcawk5,data=sample_data(pswk1), colour = "Group")



########################## Differential Abundance Analysis #############################



## ANCOMBC: Differential Abundance Analysis between 'high' and 'low' DXN producers

# Subset relevant phyloseq data with absolute abundance 
treat0counts <- wk0_counts %>% 
  subset_samples(Group == 'treatment')
treat2counts <- wk2_counts %>% 
  subset_samples(Group == 'treatment')
treat4counts <- wk4_counts %>% 
  subset_samples(Group == 'treatment')
treat6counts <- wk6_counts %>% 
  subset_samples(Group == 'treatment')
treat8counts <- wk8_counts %>% 
  subset_samples(Group == 'treatment')


# Creat column names for results
col_names <- c('beta', 'se', 'w', 'p_val', 'q_val', 'diff_abn')
wk2names <- map_chr(col_names, function(x) paste0('wk2_dxn_low_', x))
wk4names <- map_chr(col_names, function(x) paste0('wk4_dxn_low_', x))
wk6names <- map_chr(col_names, function(x) paste0('wk6_dxn_low_', x))
wk8names <- map_chr(col_names, function(x) paste0('wk8_dxn_low_', x))


## ANCOMBC week 2
out_dxn2 <- ancombc(treat2counts, formula = 'metabo_dxn', 
                    p_adj_method = 'holm', zero_cut = 0.90, lib_cut = 1000, 
                    group = 'metabo_dxn', struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)

res_dxn2 <- as.data.frame(out_dxn2$res)
# filter sig results and add taxa names
filres2 <- res_dxn2 %>% filter(metabo_dxnlow.5 == TRUE) %>% 
  dplyr::rename(dxn_low_beta = metabo_dxnlow, dxn_low_se = metabo_dxnlow.1,
                dxn_low_w = metabo_dxnlow.2, dxn_low_p_val = metabo_dxnlow.3,
                dxn_low_q_val = metabo_dxnlow.4, dxn_low_diff_abn = metabo_dxnlow.5) %>% 
  mutate(week = 2) %>% 
  cbind(., data.frame(tax_table(ps)[row.names(.), ]))

## ANCOMBC week 4
out_dxn4 <- ancombc(treat4counts, formula = 'metabo_dxn', 
                    p_adj_method = 'holm', zero_cut = 0.90, lib_cut = 1000, 
                    group = 'metabo_dxn', struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)

res_dxn4 <- as.data.frame(out_dxn4$res)
# filter sig results and add taxa names
filres4 <- res_dxn4 %>% filter(metabo_dxnlow.5 == TRUE) %>% 
  dplyr::rename(dxn_low_beta = metabo_dxnlow, dxn_low_se = metabo_dxnlow.1,
                dxn_low_w = metabo_dxnlow.2, dxn_low_p_val = metabo_dxnlow.3,
                dxn_low_q_val = metabo_dxnlow.4, dxn_low_diff_abn = metabo_dxnlow.5) %>%
  mutate(week = 4) %>% 
  cbind(., data.frame(tax_table(ps)[row.names(.), ]))

## ANCOMBC week 6
out_dxn6 <- ancombc(treat6counts, formula = 'metabo_dxn', 
                    p_adj_method = 'holm', zero_cut = 0.90, lib_cut = 1000, 
                    group = 'metabo_dxn', struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)

res_dxn6 <- as.data.frame(out_dxn6$res)
# filter sig results and add taxa names
filres6 <- res_dxn6 %>% filter(metabo_dxnlow.5 == TRUE) %>%
  dplyr::rename(dxn_low_beta = metabo_dxnlow, dxn_low_se = metabo_dxnlow.1,
                dxn_low_w = metabo_dxnlow.2, dxn_low_p_val = metabo_dxnlow.3,
                dxn_low_q_val = metabo_dxnlow.4, dxn_low_diff_abn = metabo_dxnlow.5) %>%
  mutate(week = 6) %>% 
  cbind(., data.frame(tax_table(ps)[row.names(.), ]))

## ANCOMBC week 8
out_dxn8 <- ancombc(treat8counts, formula = 'metabo_dxn', 
                    p_adj_method = 'holm', zero_cut = 0.90, lib_cut = 1000, 
                    group = 'metabo_dxn', struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)

res_dxn8 <- as.data.frame(out_dxn8$res)
# filter sig results and add taxa names
filres8 <- res_dxn8 %>% filter(metabo_dxnlow.5 == TRUE) %>% 
  dplyr::rename(dxn_low_beta = metabo_dxnlow, dxn_low_se = metabo_dxnlow.1,
                dxn_low_w = metabo_dxnlow.2, dxn_low_p_val = metabo_dxnlow.3,
                dxn_low_q_val = metabo_dxnlow.4, dxn_low_diff_abn = metabo_dxnlow.5) %>%
  mutate(week = 8) %>% 
  cbind(., data.frame(tax_table(ps)[row.names(.), ]))




# Create dataframe for ASVs at time points

# compile all differentially abundant ASV names
names <- unique(c(row.names(filres2), row.names(filres4), 
                  row.names(filres6), row.names(filres8)))
# create dataframe of ASVs by timepoint
wk2 <- as.numeric(names %in% row.names(filres2))
wk4 <- as.numeric(names %in% row.names(filres4))
wk6 <- as.numeric(names %in% row.names(filres6))
wk8 <- as.numeric(names %in% row.names(filres8))
df_dxn <- cbind(wk2, wk4, wk6, wk8)
df_dxn <- df_dxn %>% data.frame() %>% 
  rowwise() %>% 
  mutate(total = sum(wk2, wk4, wk6, wk8))
# add ASV as rownames
rownames(df_dxn) <- names

# ASVs in at least 3 time points
T3 <- df_dxn %>%
  rownames_to_column %>% 
  filter(total >= 3)

# combine results from time points
# keep ASVs greater in 'high' producers and filter for occurance in at least 3 time points
# RESULTS 
dxn_res <- rbind(filres2, filres4, filres6, filres8) %>% 
  filter(dxn_low_beta < 0) %>% 
  filter(rownames(.) %in% T3$rowname)






## ANCOMBC: Differential Abundance Analysis between 'high' and 'low' 8PN producers

# Subset relevant phyloseq data with absolute abundance 
treat0counts <- wk0_counts %>% 
  subset_samples(Group == 'treatment')
treat2counts <- wk2_counts %>% 
  subset_samples(Group == 'treatment')
treat4counts <- wk4_counts %>% 
  subset_samples(Group == 'treatment')
treat6counts <- wk6_counts %>% 
  subset_samples(Group == 'treatment')
treat8counts <- wk8_counts %>% 
  subset_samples(Group == 'treatment')


# Creat column names for results
col_names <- c('beta', 'se', 'w', 'p_val', 'q_val', 'diff_abn')
wk2names <- map_chr(col_names, function(x) paste0('wk2_dxn_low_', x))
wk4names <- map_chr(col_names, function(x) paste0('wk4_dxn_low_', x))
wk6names <- map_chr(col_names, function(x) paste0('wk6_dxn_low_', x))
wk8names <- map_chr(col_names, function(x) paste0('wk8_dxn_low_', x))


## ANCOMBC week 2
out_8pn2 <- ancombc(treat2counts, formula = 'metabo_8pn', 
                    p_adj_method = 'holm', zero_cut = 0.90, lib_cut = 1000, 
                    group = 'metabo_8pn', struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)

res_8pn2 <- as.data.frame(out_8pn2$res)
# filter sig results and add taxa names
res2 <- res_8pn2 %>% filter(metabo_8pnlow.5 == TRUE) %>% 
  dplyr::rename(x8pn_low_beta = metabo_8pnlow, x8pn_low_se = metabo_8pnlow.1,
                x8pn_low_w = metabo_8pnlow.2, x8pn_low_p_val = metabo_8pnlow.3,
                x8pn_low_q_val = metabo_8pnlow.4, x8pn_low_diff_abn = metabo_8pnlow.5) %>% 
  mutate(week = 2) %>% 
  cbind(., data.frame(tax_table(ps)[row.names(.), ]))

## ANCOMBC week 4
out_8pn4 <- ancombc(treat4counts, formula = 'metabo_8pn', 
                    p_adj_method = 'holm', zero_cut = 0.90, lib_cut = 1000, 
                    group = 'metabo_8pn', struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)

res_8pn4 <- as.data.frame(out_8pn4$res)
# filter sig results and add taxa names
res4 <- res_8pn4 %>% filter(metabo_8pnlow.5 == TRUE) %>% 
  dplyr::rename(x8pn_low_beta = metabo_8pnlow, x8pn_low_se = metabo_8pnlow.1,
                x8pn_low_w = metabo_8pnlow.2, x8pn_low_p_val = metabo_8pnlow.3,
                x8pn_low_q_val = metabo_8pnlow.4, x8pn_low_diff_abn = metabo_8pnlow.5) %>%
  mutate(week = 4) %>% 
  cbind(., data.frame(tax_table(ps)[row.names(.), ]))

## ANCOMBC week 6
out_8pn6 <- ancombc(treat6counts, formula = 'metabo_8pn', 
                    p_adj_method = 'holm', zero_cut = 0.90, lib_cut = 1000, 
                    group = 'metabo_8pn', struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)

res_8pn6 <- as.data.frame(out_8pn6$res)
# filter sig results and add taxa names
res6 <- res_8pn6 %>% filter(metabo_8pnlow.5 == TRUE) %>%
  dplyr::rename(x8pn_low_beta = metabo_8pnlow, x8pn_low_se = metabo_8pnlow.1,
                x8pn_low_w = metabo_8pnlow.2, x8pn_low_p_val = metabo_8pnlow.3,
                x8pn_low_q_val = metabo_8pnlow.4, x8pn_low_diff_abn = metabo_8pnlow.5) %>%
  mutate(week = 6) %>% 
  cbind(., data.frame(tax_table(ps)[row.names(.), ]))

## ANCOMBC week 8
out_8pn8 <- ancombc(treat8counts, formula = 'metabo_8pn', 
                    p_adj_method = 'holm', zero_cut = 0.90, lib_cut = 1000, 
                    group = 'metabo_8pn', struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                    max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)

res_8pn8 <- as.data.frame(out_8pn8$res)
# filter sig results and add taxa names
res8 <- res_8pn8 %>% filter(metabo_8pnlow.5 == TRUE) %>% 
  dplyr::rename(x8pn_low_beta = metabo_8pnlow, x8pn_low_se = metabo_8pnlow.1,
                x8pn_low_w = metabo_8pnlow.2, x8pn_low_p_val = metabo_8pnlow.3,
                x8pn_low_q_val = metabo_8pnlow.4, x8pn_low_diff_abn = metabo_8pnlow.5) %>%
  mutate(week = 8) %>% 
  cbind(., data.frame(tax_table(ps)[row.names(.), ]))




# Create dataframe for ASVs at time points

# compile all differentially abundant ASV names
names_8pn <- unique(c(row.names(res2), row.names(res4), 
                      row.names(res6), row.names(res8)))
# create dataframe of ASVs by timepoint
wk2_8pn <- as.numeric(names_8pn %in% row.names(res2))
wk4_8pn <- as.numeric(names_8pn %in% row.names(res4))
wk6_8pn <- as.numeric(names_8pn %in% row.names(res6))
wk8_8pn <- as.numeric(names_8pn %in% row.names(res8))
df_8pn <- cbind(wk2_8pn, wk4_8pn, wk6_8pn, wk8_8pn)

df_8pn <- df_8pn %>% data.frame() %>% 
  rowwise() %>% 
  mutate(total = sum(wk2_8pn, wk4_8pn, wk6_8pn, wk8_8pn))
# add ASV as rownames
rownames(df_8pn) <- names_8pn

# ASVs in at least 3 time points
T4 <- df_8pn %>%
  rownames_to_column %>% 
  filter(total >= 3)

# combine results from time points
# keep ASVs greater in 'high' producers and filter for occurance in at least 3 time points
# RESULTS
x8pn_res <- rbind(res2, res4, res6, res8) %>% 
  filter(x8pn_low_beta < 0) %>% 
  filter(rownames(.) %in% T4$rowname)









################################## Other Figures #########################################



# FIGURE: Relative Abundance (family) at Baseline

# Make dataframe
fadf <- ps %>% 
  subset_samples(Week == 1) %>% 
  tax_glom('Family') %>% 
  psmelt() %>% 
  rename(ASV = OTU) %>% 
  rename(Subject = Participant.ID) %>% 
  arrange(Family, Group)

# Plot

f1 <- ggplot(fadf, aes(x=Subject, y=Abundance, fill=Family)) +
  geom_bar(stat = 'identity', position = 'fill') +
  facet_wrap(~Group, scales = 'free_x') +
  scale_x_discrete() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5),
        axis.text.y = element_text(size = 10, vjust = 0.8),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.5, 'cm'),
        axis.title.x = element_text(size = 15, vjust = 0),
        axis.title.y = element_text(size = 15, vjust = 1)) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_percent()) +
  labs(x = 'Participant ID',
       y = 'Microbial Abundance') +
  ggtitle('Relative Abundance at Baseline')



# FIGURE: Family Relative Abundance at Week 8 (study close)

# Make dataframe
fadf2 <- ps %>% 
  subset_samples(Week == 5) %>% 
  tax_glom('Family') %>% 
  psmelt() %>% 
  rename(ASV = OTU) %>% 
  rename(Subject = Participant.ID) %>% 
  arrange(Family, Group)

# Plot values
f6 <- ggplot(fadf2, aes(x=Subject, y=Abundance, fill=Family)) +
  geom_bar(stat = 'identity', position = 'fill') +
  facet_wrap(~Group, scales = 'free_x') +
  scale_x_discrete() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5),
        axis.text.y = element_text(size = 10, vjust = 0.8),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        axis.title.x = element_text(size = 15, vjust = 0),
        axis.title.y = element_text(size = 15, vjust = 1)) +
  scale_y_continuous(expand = c(0,0), labels = scales::label_percent()) +
  labs(x = 'Subject ID',
       y = 'Microbial Abundance') +
  ggtitle('Relative Abundance at Week 8')







# FIGURE: Metabolite concentration by week

# manipulate data 

newfam <- psfam %>%   
  pivot_longer(c("DXN", "X8PN", "IXN"), names_to = "metabolites", values_to = "int") %>%
  rename_with(tolower) %>%
  #rename_with(~gsub('.', '_', .x)) %>%
  mutate(wkc = ifelse(week == 1, 0,
                      ifelse(week == 2, 2,
                             ifelse(week == 3, 4,
                                    ifelse(week == 4, 6, 8))))) %>%
  modify_at('wkc', factor) %>% 
  mutate(metabolites = ifelse(metabolites == 'X8PN', '8PN', metabolites)) %>% 
  modify_at('metabolites', factor) %>% 
  mutate(conc = int / 1000) %>% 
  rename('Participant' = 'participant.id')

# line plots by metabolite

f4 <- ggplot(newfam, aes(x = wkc, y = conc, color = Participant)) +
  geom_point() +
  facet_wrap(~metabolites, scales = 'free_y') +
  geom_line(aes(group = Participant)) +
  theme(axis.text.x = element_text(size = 15, vjust = 1),
        axis.text.y = element_text(size = 15, vjust = 1),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.5, 'cm')) +
  labs(y = 'Concentration (µg/g dry wt)',
       x = 'Weeks') +
  ggtitle('Metabolite Concentration by Week')




