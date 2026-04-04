#loading libraries
library(tidyverse)
library(phyloseq)
library(vegan)

#loading ps
ps = readRDS('../datasets/ps_filtered.rds') # TODO check file path

# Subset hindgut samples 
hindgut = subset_samples(ps, sample_type == "hindgut")

### Remove obvious outliers ###
# 1/ 13414.100.hg.R1.fastq.gz has MDS1 of -8542.9979 when included
# 2/ 13414.34.hg.R1.fastq.gz has MDS2 of -3.4208732669 when included

# mds_data_rm_outliers = mds_data %>% 
#   filter(Row.names != "13414.100.hg.R1.fastq.gz") %>%
#   filter(Row.names != "13414.34.hg.R1.fastq.gz")

to_remove <- c("13414.100.hg.R1.fastq.gz", "13414.34.hg.R1.fastq.gz")
# ps <- subset_samples(ps, !(SampleID %in% to_remove))
ps <- prune_samples(!(sample_names(ps) %in% to_remove), ps)



# rarefying data
psrare = hindgut %>% rarefy_even_depth(sample.size = 2525, rngseed = 1)
# sample depth selected from alpha diversity


## Calculating metrics (bray-curtis) ##
ps_bray = phyloseq::distance(psrare, method = "bray") #distance object
View(as.matrix(ps_bray)) # for human interpretation


# MDS scaling
set.seed(1)
mds = metaMDS(ps_bray) 

# Extracting data
mds_data = mds$points %>% as.data.frame %>%  
  merge(sample_data(psrare), by='row.names', sort=F) %>%
  
  # Sort data by swim performance group (Fast, Moderate, Slow)
  mutate(speed_category = ifelse(swim_performance %in% c('flow refuging', 'burrowing'), 'slow',
                            ifelse(swim_performance %in% c('generalist'), 'moderate', 'fast')), .before = swim_performance)



# Plotting
p = mds_data_rm_outliers %>%
  ggplot(aes(MDS1,MDS2,color = speed_category)) + 
  geom_point() + #scatterplot 
  stat_ellipse() + # 95% confidence interval around median
  labs(color = "Swim Performance") +
  theme_classic(base_size=18) +
  scale_color_manual(values = c("fast" = "#8E7CA6", 
                                "moderate" = "#E07A9A", 
                                "slow" = "#fcba65")) #colour palette
p

# # log scale x for human interpretability
# p_logx = mds_data_rm_outliers %>%
#   ggplot(aes(MDS1,MDS2,color = speed_category)) + 
#   geom_point() + #scatterplot 
#   stat_ellipse() + # 95% confidence interval around median
#   labs(color = "Speed") +
#   theme_classic(base_size=18) +
#   scale_x_log10() +
#   scale_color_manual(values = c("fast" = "#8E7CA6", 
#                                 "moderate" = "#E07A9A", 
#                                 "slow" = "#fcba65")) #colour palette
# p_logx



### PERMANOVA
metadata = sample_data(mds_data) %>% data.frame() #%>% na.omit() 

# colSums(is.na(metadata)) 


# Running adonis on data to get Single variable
stats_univar = adonis2(ps_bray ~ speed_category, data = metadata, na.action = na.omit)
# p-value = 0.012 
# R2 = 0.0681346
# our model is not the best at explaining variance in our diversity 0.06 



### Comparing Two Speeds ###
# 1/ Compare fast and slow
# 2/ Compare fast and moderate
# 3/ Compare moderate and slow

to_compare_fast = c('accelerator', 'crusier sprinter', 'manoeuvrer')
to_compare_mod = c('generalist')
to_compare_slow = c('flow refuging', 'burrowing')


# list of samples belonging to specified speeds
# 1/ fast & slow
samples_to_keep_fast_slow = metadata %>% 
  filter(swim_performance %in% c(to_compare_fast, to_compare_slow)) %>% 
  select(Row.names) %>% 
  pull() #to turn output to char vector
# Subsetting distance matrix to only include fast & slow
ps_bray_sub_fast_slow = as.matrix(ps_bray)[samples_to_keep_fast_slow,samples_to_keep_fast_slow] %>% as.dist()
# Running PERMANOVA on fast and slow
stats_fast_slow = adonis2(ps_bray_sub_fast_slow ~ swim_performance, #metadata has speed_category, ps doesn't
                          data = metadata %>% filter(swim_performance %in% c(to_compare_fast, to_compare_slow)))
  # fast and slow are not significant 
  # p=0.771
  # R2=0.2291487 model explains it fairly well

  
# 2/ fast & moderate
samples_to_keep_fast_mod = metadata %>% 
  filter(swim_performance %in% c(to_compare_fast, to_compare_mod)) %>% 
  select(Row.names) %>% 
  pull() 
# Subsetting distance matrix to only include fast & moderate
ps_bray_sub_fast_mod = as.matrix(ps_bray)[samples_to_keep_fast_mod,samples_to_keep_fast_mod] %>% as.dist()
# Running PERMANOVA on fast and moderate
stats_fast_mod = adonis2(ps_bray_sub_fast_mod ~ swim_performance, #metadata has speed_category, ps doesn't
                          data = metadata %>% filter(swim_performance %in% c(to_compare_fast, to_compare_mod)))
# fast and moderate are not significant 
# p=0.08
# R2=0.09559119 model doesn't explains this well

  
# 3/ moderate & slow
samples_to_keep_mod_slow = metadata %>% 
  filter(swim_performance %in% c(to_compare_mod, to_compare_slow)) %>% 
  select(Row.names) %>% 
  pull() 
# Subsetting distance matrix to only include moderate & slow
ps_bray_sub_mod_slow = as.matrix(ps_bray)[samples_to_keep_mod_slow,samples_to_keep_mod_slow] %>% as.dist()
# Running PERMANOVA on moderate and slow
stats_mod_slow = adonis2(ps_bray_sub_mod_slow ~ swim_performance, #metadata has speed_category, ps doesn't
                         data = metadata %>% filter(swim_performance %in% c(to_compare_mod, to_compare_slow)))
# moderate and slow are not significant 
# p=0.346
# R2=0.08274644 model doesn't explains this well



# #temp = as.matrix(ps_bray)
# #setdiff(samples_to_keep, rownames(as.matrix(ps_bray)))
# 
# # Error in as.matrix(ps_bray)[samples_to_keep, samples_to_keep] : 
# #   subscript out of bounds
# # Missing Samples: The samples_to_keep vector contains names that are not present in rownames(as.matrix(ps_bray)).
# # Check for missing samples: setdiff(samples_to_keep, rownames(as.matrix(ps_bray))).


# Combining stats into a table

stats = bind_rows('Univariate' = stats_univar %>% as.data.frame() %>% 
                    rownames_to_column('Variable'),
                  'Fast v. Slow' = stats_fast_slow %>% as.data.frame() %>% 
                    rownames_to_column('Variable'),
                  'Fast v. Moderate' = stats_fast_mod %>% as.data.frame() %>% 
                    rownames_to_column('Variable'),
                  'Moderate v. Slow' = stats_mod_slow %>% as.data.frame() %>% 
                    rownames_to_column('Variable'),
                  .id = 'Model') # First column is Model, either Uni or Multi



stats_clean = stats %>% rename(Pval = `Pr(>F)`) %>% #renaming column to Pval
  filter(!is.na(Pval)) # Removing irrelevant lines

stats_clean # swim_performance seems to be significant



### Saving Output
# Save plot
ggsave('../plots/Beta Diversity.jpeg',
       plot = p,
       height= 4, width = 6)

ggsave('../plots/Beta Diversity Logx.jpeg',
       plot = p_logx,
       height= 4, width = 6)


# Save statistics as an Excel file
writexl::write_xlsx(list('Beta Diversity' = stats_clean),
                    '..//Beta Diversity Stats.xlsx')

