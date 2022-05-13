
## Used to instal phyloseq
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("phyloseq")

# load required packages
library(phyloseq)
library(dplyr)
library(magrittr)

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





# random forest

sample_names <- row.names(otu_table(ps))

library(randomForest)

rf_psnew <- psmelt(ps) %>% 
  dplyr::select(OTU, Sample, Abundance, Group, DXN, X8PN) %>% 
  pivot_wider(names_from = OTU, values_from = Abundance) %>% 
  arrange(factor(Sample, levels = sample_names))


## I have not tested random forest but this code might work...

rftreat <- randomForest(Group ~ ., data = rf_psnew, ntree = 100)








# spearmans correlation agglomerated to Family level 

psfam <- pstreat %>% 
  tax_glom(taxrank = "Family") %>% 
  psmelt()

corrfam <- psfam %>%   
  rename_with(tolower) %>% 
  dplyr::rename('asv' = 'otu') %>% 
  select(asv, sample, abundance, participant.id, week, family, ixn, x8pn, dxn) %>% 
  pivot_longer(c(ixn, x8pn, dxn), names_to = 'metab', values_to = 'conc')

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

cyb_adj <- cyb %>% 
  filter(padj <= 0.05)
## family = Enterobacteriaceae with 8PN [rho = 0.877, adj pvalue = 0.0193]

# fitler by week = 3 (which is actually week 4)

cyb3 <- cyb %>%
  filter(week == 3)

ggplot(cyb3, aes(x = family, y = metabolite, fill = rho)) +
  geom_raster() +
  coord_fixed() +
  scale_fill_gradient(high = 'green', low = 'red')+ #, limits = c(-1,1)) +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab('Family') +
  ylab('Metabolite') +
  ggtitle('Family vs Metabolite')




# Spearman's correlation agglomerated to genus level

psgen <- pstreat %>% 
  tax_glom(taxrank = "Genus") %>% 
  psmelt()

corrgen <- psgen %>% 
  rename('ASV' = 'OTU') %>% 
  rename_with(tolower) %>% 
  select(asv, sample, abundance, participant.id, week, ixn, x8pn, dxn, genus, family) %>% 
  pivot_longer(c(ixn, x8pn, dxn), names_to = 'metab', values_to = 'conc')

cbg <- corrgen %>% 
  group_by(metab, genus, week) %>% 
  nest() %>% 
  mutate(spr = map(data, function(x) cor.test(x$abundance, x$conc, method = 'spearman'))) %>% 
  mutate(rho = map_dbl(spr, function(x) x$estimate)) %>% 
  mutate(pvalue = map_dbl(spr, function(x) x$p.value)) %>% 
  mutate(padj = p.adjust(pvalue, method = 'BH'))

cbg2 <- cbg %>% 
  filter(padj <= 0.05 & rho >= 0.7)

cbg3 <- cbg %>% 
  filter(rho > 0)

# filtering for week 3

cbg3 <- cbg %>% 
  filter(week == 3)

ggplot(cbg3, aes(x = genus, y = metab, fill = rho)) +
  geom_raster() +
  coord_fixed() +
  scale_fill_gradient(high = 'green', low = 'red', limits = c(-1,1)) +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab('Genus') +
  ylab('Metabolite') +
  ggtitle('Genus vs Metabolite')





# Figure: Metabolite levels by week

sample_data(psfam) 

newfam <- as.data.frame(otu_table(psfam)) %>%
  cbind(as.data.frame(sample_data(psfam))) %>%
  pivot_longer(starts_with('ASV'), names_to = 'ASV', values_to = "relab") %>% 
  pivot_longer(c("DXN", "X8PN", "IXN"), names_to = "metabolites", values_to = "int") %>%
  rename_with(tolower) %>%
  #rename_with(~gsub('.', '_', .x)) %>%
  mutate(wkc = ifelse(week == 1, 0,
                      ifelse(week == 2, 2,
                             ifelse(week == 3, 4,
                                    ifelse(week == 4, 6, 8))))) %>%
  modify_at('wkc', factor)
# group_by(Week, metabolites, ASV) %>%
# nest() %>%
# mutate(pval = map_dbl(data, function(x) cor.test(x$relab, x$int)$p.value))

metab <- ggplot(newfam, aes(x = wkc, y = int, color = participant.id)) +
  geom_point() +
  geom_path(aes(group = participant.id)) +
  facet_wrap(~metabolites, scales = 'free_y') +
  labs(x = Week)








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

