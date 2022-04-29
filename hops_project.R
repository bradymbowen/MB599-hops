
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



# Agglomerate to the family level
ps_family <- psdata %>% tax_glom(taxrank = "Family")


# Spearman's Correlation of relative abundance at each time point (total = 5) 


# Spearman's correlation between change in relative abundance and change in concentration of metabolite


















