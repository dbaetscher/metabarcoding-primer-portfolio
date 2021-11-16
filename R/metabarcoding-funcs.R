library(tidyverse)
library(readr)
library(microDecon)


## DATA DECONTAMINATION

# use the R package microdecon to decontaminate sequence data using extraction blanks
# must be done per-locus because contamination varies tremendously across loci

# select the locus and the appropriate samples to be decontaminated with EXB 1 and 2
loc_blank_decon <- function(consensus_df, loc, num_samples, num_blanks){ 
  exb12_decon.nscoi <- consensus_df %>%
    filter(locus == loc) %>%
    filter(group %in% c("full_reference", "vouchered_reference") | group == "blank" & sample %in% c("EXB1", "EXB2")) %>%
    select(seq, sample, count, taxon) %>%
    spread(sample, value = count)
  
  # more the taxa column to the end
  exb12_format.nscoi <-exb12_decon.nscoi %>%
    select(1, 3:length(exb12_decon.nscoi), 2)
  
  # replace NAs with 0
  exb12_format.nscoi[is.na(exb12_format.nscoi)] = 0
  
  # turns out it doesn't work on tibbles and needs to be a data frame.
  exb12_df.nscoi <- as.data.frame(exb12_format.nscoi)
  
  # decontaminate
  tmp2 <- decon(exb12_df.nscoi, numb.blanks = num_blanks, numb.ind = num_samples, taxa = T, runs = 2)
  
  print(tmp2$reads.removed)
  print(tmp2$OTUs.removed)
  
  # remove contamination
  exb12_deconed.nscoi <- remove.cont(exb12_df.nscoi, numb.blanks = num_blanks, taxa = T, runs = 2)
  
  # reformat
  slim.df <- exb12_deconed.nscoi %>%
    select(-Mean.blank) %>%
    select(seq, taxon, 2:(length(exb12_deconed.nscoi)-1))
  
  exb12_decon_df_nscoi_long <- slim.df %>%
    pivot_longer(cols = 3:length(slim.df), names_to = "sample", values_to = "count")
  
  print(exb12_decon_df_nscoi_long)
}

## format for microdecon

microdecon_format <- function(df, loc){
  tmp.df <- df %>%
    filter(locus == loc) %>%
    select(seq, sample, count, taxon) %>%
    spread(sample, value = count)
  
  # more the taxa column to the end
  tmp1 <- tmp.df %>%
    select(1, 3:length(tmp.df), 2)
  
  # replace NAs with 0
  tmp1[is.na(tmp1)] = 0
  
  # turns out it doesn't work on tibbles and needs to be a data frame.
  decon.df.formatted <- as.data.frame(tmp1)
  
}

## reformat to long-tibble from microdecon
# reformat
reformat_from_microdecon <- function(decontaminated_df){
  slim.df <- decontaminated_df %>%
  select(seq, taxon, 2:(length(decontaminated_df)-1))

  slim.df %>%
  pivot_longer(cols = 3:length(slim.df), names_to = "sample", values_to = "count")

}


#
#
#

### Species occupancy detection models


# to remove false positives
# based on Bayesian implementation from Ryan Kelly's GitHub

### SODM for aquafeeds - full dataset
sodm.by.locus <- function(feature_table, site, n_replicates, loc){
  # select the "site" replicates and the locus
  tmp.df <- feature_table %>%
    filter(stri_detect(sample, regex = site)) %>% # comment this to test the ASV-centric df
    filter(locus == loc)
  
  # format the dataframe appropriately
  wide.loc.frame <- tmp.df %>%
    mutate(count = 1) %>% # change counts to presence/absence
    select(seq, sample, count) %>%
    group_by(seq, sample, count) %>%
    unique() %>% # eliminate duplicate entries for the same ASV and sample (these seem to be a thing for aquaF2)
    pivot_wider(names_from = sample, values_from = count)
  
  # Replace NA with 0, since technically, these are absences.
  wide.loc.frame[is.na(wide.loc.frame)] <- 0                
  
  # format must be dataframe, not tibble for converting rownames properly
  wide.loc.frame <- data.frame(wide.loc.frame)
  
  # helper function for maintaining rownames in matrix format
  matrix.please <- function(x) {
    m<-as.matrix(x[,-1])
    rownames(m)<-x[,1]
    m
  }
  
  # convert df to matrix
  wide.loc.matrix <- matrix.please(wide.loc.frame)
  
  # then call the Stan code to work on the dataset we created above. 
  # reformatting our data from above to match the Stan code's inputs
  testData <<- 
    data.frame(K = n_replicates,   #trials per species (row)
               N = rowSums(wide.loc.matrix),  #detections per species
               z = ifelse(rowSums(wide.loc.matrix) > 0, 1, 0)  #was it ever detected at this site? (integer that helps estimate psi)
    )
  
  #return(testData)
  
  
  mySOMmodel <- stan(file = "Stan_SOM_demo.stan",
                     data = list(
                       S = nrow(testData),
                       K = testData$K,
                       N = testData$N,
                       z = ifelse(testData$N > 0, 1, 0)
                     ),
                     chains = 4,   #number of chains
                     iter = 5000   #number of iterations per chain
  )
  
  # print a plot of that output so that I can easily determine the cutoff-threshold/breakpoint
  # tmp.plot <- tidy(mySOMmodel) %>%
  #   filter(!term %in% c("psi", "p11", "p10")) %>%
  #   arrange(desc(estimate)) %>%
  #   ggplot(aes(x = reorder(term, estimate), y = estimate)) +
  #   geom_point() +
  #   theme(
  #     axis.text.x = element_blank()
  #   )
  # 
  # print(tmp.plot)
  
  # tidy(mySOMmodel) %>% # this is redundant, but less informative than the CSV output with the ASV information
  #   arrange(desc(estimate)) %>%
  #   write_csv(paste0("csv_outputs/sodm_tax_only_",loc,"_",site,".csv"))
  
  p1 <- mcmc_areas(mySOMmodel, 
                   pars = c("psi", "p11", "p10"))
  p2 <- mcmc_intervals(mySOMmodel, 
                       regex_pars = "Occupancy_prob")
  
  list(plot(p1), plot(p2))
  
  # must have run modelText prior to this working
  write.table(modelText, "Stan_SOM_allPossibilities.stan", row.names = F, quote = F, col.names = F)
  SOMmodel_0_10_detections <- stan(file = "Stan_SOM_allPossibilities.stan", 
                                   data = list(
                                     S = nrow(testData),
                                     K = testData$K,
                                     N = testData$N,
                                     z = ifelse(testData$N > 0, 1, 0)
                                   ), 
                                   chains = 4,   #number of chains
                                   iter = 100000   #number of iterations per chain
  )
  
  # print the output of that plot
  # plot.occu <- mcmc_intervals(SOMmodel_0_10_detections, 
  #                             regex_pars = "Occupancy_prob") +
  #   xlab("Probability of Occupancy") +
  #   ylab("Number of detections (from 0 to 10 out of 10)")
  # 
  # ggsave(paste0("pdf_outputs/SODM/occupancy_tax_only_",loc,"_",site,"_5kiter.pdf"))
  
  
  # also need to save the output in a form that connects the ASVs to the occupancy estimates
  # grab the name of the ASVs
  asvs <- data.frame(wide.loc.frame[,1]) %>%
    rename(seq = wide.loc.frame...1.)
  
  # combine the ASV names with the sodm output
  asv.sodm.output <- tidy(mySOMmodel) %>%
    filter(!term %in% c("psi", "p11", "p10")) %>%
    bind_cols(., asvs) 
  
  # write that to a CSV for now 
  ## NOTE: this output directory was edited for the downsampling!!
  asv.sodm.output %>%
    #write_csv(paste0("csv_outputs/sodm/",loc,"_",site,"_ASV_SODM_100kiter.csv"))
    write_csv(paste0("/Users/dianabaetscher/Documents/git-repos/metabarcoding-primer-portfolio/downsampled_loci/csv_outputs/sodm/",loc,"_",site,"_ASV_SODM_100kiter.csv"))
  
}



### SODM by aquafeed (all loci together)
# this works because there are no redundant ASVs across the different loci
# so locus identity should be irrelevant

# again, based on the fct above and Ryan Kelly's demo.
# @param feature_table
# @param site = sample, "AF002"
# @param n_replicates = 4, e.g.


sodm.by.feed <- function(feature_table, site, n_replicates, loc){
  # select the "site" replicates
  tmp.df <- feature_table %>%
    filter(stri_detect(sample, regex = site)) # comment this to test the ASV-centric df
  
  # report how many replicates are present
  # n_reps <- tmp.df %>%
  #   ungroup() %>%
  #   select(sample) %>%
  #   unique() %>%
  #   tally()
  # 
  # if(n_reps$n < 4) warning(`num reps is < 4`)
  #   
  
  # format the dataframe appropriately
  wide.loc.frame <- tmp.df %>%
    mutate(count = 1) %>% # change counts to presence/absence
    select(seq, sample, count) %>%
    group_by(seq, sample, count) %>%
    unique() %>% # eliminate duplicate entries for the same ASV and sample (these seem to be a thing for aquaF2)
    pivot_wider(names_from = sample, values_from = count)
  
  # Replace NA with 0, since technically, these are absences.
  wide.loc.frame[is.na(wide.loc.frame)] <- 0                
  
  # format must be dataframe, not tibble for converting rownames properly
  wide.loc.frame <- data.frame(wide.loc.frame)
  
  # helper function for maintaining rownames in matrix format
  matrix.please <- function(x) {
    m<-as.matrix(x[,-1])
    rownames(m)<-x[,1]
    m
  }
  
  # convert df to matrix
  wide.loc.matrix <- matrix.please(wide.loc.frame)
  
  # then call the Stan code to work on the dataset we created above. 
  # reformatting our data from above to match the Stan code's inputs
  testData <<- 
    data.frame(K = n_replicates,   #trials per species (row)
               N = rowSums(wide.loc.matrix),  #detections per species
               z = ifelse(rowSums(wide.loc.matrix) > 0, 1, 0)  #was it ever detected at this site? (integer that helps estimate psi)
    )
  
  #return(testData)
  
  
  mySOMmodel <- stan(file = "Stan_SOM_demo.stan",
                     data = list(
                       S = nrow(testData),
                       K = testData$K,
                       N = testData$N,
                       z = ifelse(testData$N > 0, 1, 0)
                     ),
                     chains = 4,   #number of chains
                     iter = 5000   #number of iterations per chain
  )
  
  # print a plot of that output so that I can easily determine the cutoff-threshold/breakpoint
  # tmp.plot <- tidy(mySOMmodel) %>%
  #   filter(!term %in% c("psi", "p11", "p10")) %>%
  #   arrange(desc(estimate)) %>%
  #   ggplot(aes(x = reorder(term, estimate), y = estimate)) +
  #   geom_point() +
  #   theme(
  #     axis.text.x = element_blank()
  #   )
  # 
  # print(tmp.plot)
  
  # tidy(mySOMmodel) %>% # this is redundant, but less informative than the CSV output with the ASV information
  #   arrange(desc(estimate)) %>%
  #   write_csv(paste0("csv_outputs/sodm_tax_only_",loc,"_",site,".csv"))
  
  # p1 <- mcmc_areas(mySOMmodel, 
  #                  pars = c("psi", "p11", "p10"))
  # p2 <- mcmc_intervals(mySOMmodel, 
  #                      regex_pars = "Occupancy_prob")
  # 
  # list(plot(p1), plot(p2))
  # 
  # must have run modelText prior to this working
  write.table(modelText, "Stan_SOM_allPossibilities.stan", row.names = F, quote = F, col.names = F)
  SOMmodel_0_10_detections <- stan(file = "Stan_SOM_allPossibilities.stan", 
                                   data = list(
                                     S = nrow(testData),
                                     K = testData$K,
                                     N = testData$N,
                                     z = ifelse(testData$N > 0, 1, 0)
                                   ), 
                                   chains = 4,   #number of chains
                                   iter = 5000   #number of iterations per chain
  )
  
  # print the output of that plot
  # plot.occu <- mcmc_intervals(SOMmodel_0_10_detections, 
  #                             regex_pars = "Occupancy_prob") +
  #   xlab("Probability of Occupancy") +
  #   ylab("Number of detections (from 0 to 10 out of 10)")
  # 
  # ggsave(paste0("pdf_outputs/SODM/occupancy_tax_only_",site,"_100kiter.pdf"))
  # 
  
  # also need to save the output in a form that connects the ASVs to the occupancy estimates
  
  # grab the name of the ASVs
  asvs <- data.frame(wide.loc.frame[,1]) %>%
    rename(seq = wide.loc.frame...1.)
  
  # combine the ASV names with the sodm output
  asv.sodm.output <- tidy(mySOMmodel) %>%
    filter(!term %in% c("psi", "p11", "p10")) %>%
    bind_cols(., asvs) 
  
  # write that to a CSV for now
  asv.sodm.output %>%
    write_csv(paste0("csv_outputs/SODM/",loc,"/",site,"_",loc,"_ASV_SODM_5Kiter.csv"))
  
}


### False positive: true positive ratio

# generate the ratio of false-to-true positives
# at each taxonomic level
# because of syntax, I have two variables for taxonomic level:
# taxon_level should have "", e.g., "genus"
# while level should NOT have "", e.g., genus


false_true_ratio <- function(reference_df, feature_taxonomy_df, ASV_filtered_df, taxon_level, level){
  # generate the ratio for the original, unfiltered dataframe
  tmp_f_t_ratio <- reference_df %>%
    full_join(., feature_taxonomy_df, by = c(taxon_level = "taxon")) %>% 
    filter(taxonomic_level == taxon_level) %>% # need to match the taxonomic level to the appropriate degree of specificity from the joining step
    mutate(true_positive = ifelse(is.na(control), "false", "true")) %>%
    select(locus, true_positive, level) %>%
    unique() %>%
    group_by(locus) %>%
    summarise(true_pos = sum(true_positive=="true"),
              false_pos = sum(true_positive=="false")) %>%
    mutate(f_t_ratio = ifelse((true_pos == 0 & false_pos > 0), 1, false_pos/true_pos))
  
  
  # generate the ratio for the filtered data
  tmp_sodm80 <- reference_df %>%
    full_join(., ASV_filtered_df, by = c(taxon_level = "taxon")) %>% 
    filter(taxonomic_level == taxon_level) %>% # need to match the taxonomic level to the appropriate degree of specificity from the joining step
    mutate(true_positive = ifelse(is.na(control), "false", "true")) %>%
    select(locus, true_positive, level) %>%
    unique() %>%
    group_by(locus) %>%
    summarise(true_pos = sum(true_positive=="true"),
              false_pos = sum(true_positive=="false")) %>%
    mutate(f_t_ratio_80_sodm = false_pos/true_pos)
  
  
  # Join the two and quickly plot them up to take a look
  tmp_sodm80 %>%
    left_join(., tmp_f_t_ratio, by = "locus") %>%
    mutate(ratio_diff = f_t_ratio_80_sodm-f_t_ratio)

}



### Dissimilarity assessment using Bray-Curtis in VEGAN


## Dissimilarity analyses wrapped up for all loci
# @param sodm_filtered_df = vrp_sodm_filtered_spp98
# @param loc = "mifish"
# @param sample = "VRP"

bray_nmds_complete <- function(sodm_filtered_df, loc, sample){
  
  loc.sample.df <- sodm_filtered_df %>%
    filter(str_detect(sample, sample)) %>%
    filter(locus == loc) %>% # select the locus for this iteration
    select(seq, sample, count) %>%
    unique() %>% # ensure there are no duplicates because they will cause the pivot to fail later
    group_by(seq, sample) %>%
    summarise(count = sum(count))

### community-by-species matrix
# To run the Bray-Curtis and NMDS, we will use the function metaMDS. The function requires a community-by-species matrix.

  # reformat the dataframe - wider
  locus.comm.wide.df <- loc.sample.df %>%
    pivot_wider(names_from = seq, values_from = count)
  
  # replace all NAs with 0
  locus.comm.wide.df[is.na(locus.comm.wide.df)] <- 0
  locus.comm.matrix <- column_to_rownames(locus.comm.wide.df)
  
  # split apart the sample name in the data frame to get the strata
  community.df.sep <- locus.comm.wide.df %>%
    separate(sample, into = c("reference", "pool", "pcr_rep"), remove = FALSE)
  
  # Calculating relative abundance and creating new dataframe with relative abundance data
  relative.abund <- decostand(locus.comm.matrix, method = "total") # standardization method = total, which probably makes the most sense because it sums over the rows (and not the columns - which would be among replicates... which is what we're curious to calculate the dissimilarity among)
  
  loc.matrix.prop <- as.matrix(relative.abund)
  
  # jaccard similarity (presence/absence)
  jaccard.dist <- vegdist(relative.abund, method = "jaccard")
  
  # generate an easily-filtered dataframe with any comparisons that are > 0.49 dissimilar
  # Creating easy to view matrix
  jaccard_distmat2 <- as.matrix(jaccard.dist, labels = T)
  
    # as.data.frame.table(jaccard_distmat2) %>%
    # rename(rep1 = Var1, rep2 = Var2, Jaccard_dist = Freq) %>%
    # filter(Jaccard_dist > 0.49) %>% 
    # write_csv(paste0("csv_outputs/jaccard_dissimilar/",loc,"_",sample,"_jaccard_outliers.csv"))
  
  # Permanova for jaccard 
  set.seed(129)
  jaccard.pool.perm <- adonis2(jaccard.dist~pool, data = community.df.sep, permutations = 999, method = "jaccard")
  
  # Bray-Curtis dissimilarity
  # Calculate distance matrix
  bray.distmat <- vegdist(relative.abund, method = "bray")
  
  # generate an easily-filtered dataframe with any comparisons that are > 0.49 dissimilar
  # Creating easy to view matrix
  bray_distmat2 <- as.matrix(bray.distmat, labels = T)
  
  bray_df <- as.data.frame.table(bray_distmat2) %>%
    rename(rep1 = Var1, rep2 = Var2, BrayCurtis_dist = Freq) %>%
    filter(BrayCurtis_dist > 0.49)
  
 bray_df %>% write_csv(paste0("csv_outputs/bray_dissimilar/",loc,"_",sample,"_bray_outliers.csv"))
  
  # Permanova for Bray-Curtis distance
  # I can ask whether all PCR replicates are dissimilar
  set.seed(129)
  bray.perm <- adonis2(bray.distmat~sample, data = community.df.sep, permutations = 999, method = "bray")
  bray.perm
  
  # I can ask whether the pools of DNA extractions are dissimilar
  set.seed(129)
  bray.perm.pool <- adonis2(bray.distmat~pool, data = community.df.sep, permutations = 999, method = "bray")
  bray.perm.pool
  
  ## Multivariate dispersion
  # Does the amount of (within group) variance differ between groups?
  dispersion<-betadisper(bray.distmat, group=community.df.sep$pool)
  dispersion.bray <- permutest(dispersion)
  
  disp.plot <- plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
  
  ## NMDS
  # using the vegan `metaMDS` function and the relative abundance dataframe rather than the distance matrix
  fishMDS<-metaMDS(relative.abund, distance="bray", k=2, trymax=50, autotransform=TRUE) ##k is the number of dimensions
  #fishMDS ##metaMDS takes eaither a distance matrix or your community matrix (then requires method for 'distance=')
  
  nmds.stress <- stressplot(fishMDS)
  
  # plot the NMDS output in ggplot
  ## pull points from MDS
  NMDS1 <- fishMDS$points[,1] ##also found using scores(fishMDS)
  NMDS2 <- fishMDS$points[,2]
  fish.plot<-cbind(community.df.sep, NMDS1, NMDS2)
  
  # plot ordination
  nmds.plot <- ggplot(fish.plot, aes(NMDS1, NMDS2, color=pool)) +
    geom_point(size = 2) +
    theme_minimal() +
    geom_text(data=fish.plot,aes(NMDS1, NMDS2, label=sample), position = position_jitter(0.005)) +
    annotate("text", x=max(NMDS1), y=min(NMDS2), hjust=0.9, vjust=0.5, label=paste('Stress =',round(fishMDS$stress,3))) +
    ggtitle(paste0(loc," ",sample," NMDS"))
  # stat_ellipse(type='t',size =1) ##draws 95% confidence interval ellipses
  
  # ggsave(plot = nmds.plot, filename = paste0("pdf_outputs/nmds_plots/",loc,"_",sample,"_NMDS.pdf"))

  
  ## outputs
  output.list <- list(dispersion.bray,
    bray.perm.pool,
    bray.perm,
    #bray_df,
    jaccard.pool.perm)
  
  
  # list.save(output.list, paste0("rds_outputs/dissimilarity/",loc,"_",sample,"_dissimilarity.rds"))
  
  #return(output.list)
  
  # plots
  all_plots <- list(disp.plot,
       nmds.stress,
       nmds.plot)
  

  return(all_plots)
  
}



## simple bray-curtis (downsampling, VRP)
simple.bray <- function(sodm_filtered_df, loc, sample){
  
  loc.sample.df <- sodm_filtered_df %>%
    filter(str_detect(sample, sample)) %>%
    filter(locus == loc) %>% # select the locus for this iteration
    select(seq, sample, count) %>%
    unique() %>% # ensure there are no duplicates because they will cause the pivot to fail later
    group_by(seq, sample) %>%
    summarise(count = sum(count))
  
  ### community-by-species matrix
  # To run the Bray-Curtis and NMDS, we will use the function metaMDS. The function requires a community-by-species matrix.
  
  # reformat the dataframe - wider
  locus.comm.wide.df <- loc.sample.df %>%
    pivot_wider(names_from = seq, values_from = count)
  
  # replace all NAs with 0
  locus.comm.wide.df[is.na(locus.comm.wide.df)] <- 0
  locus.comm.matrix <- column_to_rownames(locus.comm.wide.df)
  
  # split apart the sample name in the data frame to get the strata
  community.df.sep <- locus.comm.wide.df %>%
    separate(sample, into = c("reference", "pool", "pcr_rep"), remove = FALSE)
  
  # Calculating relative abundance and creating new dataframe with relative abundance data
  relative.abund <- decostand(locus.comm.matrix, method = "total") # standardization method = total, which probably makes the most sense because it sums over the rows (and not the columns - which would be among replicates... which is what we're curious to calculate the dissimilarity among)
  
  loc.matrix.prop <- as.matrix(relative.abund)
  
  # jaccard similarity (presence/absence)
  jaccard.dist <- vegdist(relative.abund, method = "jaccard")
  
  # generate an easily-filtered dataframe with any comparisons that are > 0.49 dissimilar
  # Creating easy to view matrix
  jaccard_distmat2 <- as.matrix(jaccard.dist, labels = T)
  
  # as.data.frame.table(jaccard_distmat2) %>%
  # rename(rep1 = Var1, rep2 = Var2, Jaccard_dist = Freq) %>%
  # filter(Jaccard_dist > 0.49) %>% 
  # write_csv(paste0("csv_outputs/jaccard_dissimilar/",loc,"_",sample,"_jaccard_outliers.csv"))
  
  # Permanova for jaccard 
  set.seed(129)
  jaccard.pool.perm <- adonis2(jaccard.dist~pool, data = community.df.sep, permutations = 999, method = "jaccard")
  
  # Bray-Curtis dissimilarity
  # Calculate distance matrix
  bray.distmat <- vegdist(relative.abund, method = "bray")
  
  # generate an easily-filtered dataframe with any comparisons that are > 0.49 dissimilar
  # Creating easy to view matrix
  bray_distmat2 <- as.matrix(bray.distmat, labels = T)
  
  bray_df <- as.data.frame.table(bray_distmat2) %>%
    rename(rep1 = Var1, rep2 = Var2, BrayCurtis_dist = Freq) %>%
    filter(BrayCurtis_dist > 0.49)
  
  bray_df %>% write_csv(paste0("csv_outputs/bray_dissimilar/",loc,"_",sample,"_bray_outliers.csv"))
}


## New function for mock feed ASVs
# @param sodm_filtered_df = tibble with ASVs filtered by occupancy modeling removed
# @param loc = locus (primer set)
# @param site = sample (that has replicates to compare)
bray_nmds_mock_feed <- function(sodm_filtered_df, loc, site){
  # first, make a dir if it doesn't exist
  dir.create("csv_outputs/bray_dissimilar/mockfeeds/")
  
  # filter the dataframe for the appropriate sample/locus
  loc.sample.df <- sodm_filtered_df %>%
    filter(str_detect(sample, site)) %>%
    filter(locus == loc) %>% # select the locus for this iteration
    select(seq, sample, count) %>%
    unique() %>% # ensure there are no duplicates because they will cause the pivot to fail later
    group_by(seq, sample) %>%
    summarise(count = sum(count))
  
  ### community-by-species matrix
  # To run the Bray-Curtis and NMDS, we will use the function metaMDS. The function requires a community-by-species matrix.
  
  # reformat the dataframe - wider
  locus.comm.wide.df <- loc.sample.df %>%
    pivot_wider(names_from = seq, values_from = count)
  
  # replace all NAs with 0
  locus.comm.wide.df[is.na(locus.comm.wide.df)] <- 0
  locus.comm.matrix <- column_to_rownames(locus.comm.wide.df)
  
  # split apart the sample name in the data frame to get the strata
  community.df.sep <- locus.comm.wide.df %>%
    separate(sample, into = c("filler", "perc_fishmeal", "pool", "pcr_rep"), remove = FALSE) # use this version for the feeds
    #separate(sample, into = c("reference", "pool", "pcr_rep"), remove = FALSE) # use this one for MEP and MFP
  
  # Calculating relative abundance and creating new dataframe with relative abundance data
  relative.abund <- decostand(locus.comm.matrix, method = "total") # standardization method = total, which probably makes the most sense because it sums over the rows (and not the columns - which would be among replicates... which is what we're curious to calculate the dissimilarity among)
  
  loc.matrix.prop <- as.matrix(relative.abund)
  
 
  # Bray-Curtis dissimilarity
  # Calculate distance matrix
  bray.distmat <- vegdist(relative.abund, method = "bray")
  
  # generate an easily-filtered dataframe with any comparisons that are > 0.49 dissimilar
  # Creating easy to view matrix
  bray_distmat2 <- as.matrix(bray.distmat, labels = T)
  
  bray_df <- as.data.frame.table(bray_distmat2) %>%
    rename(rep1 = Var1, rep2 = Var2, BrayCurtis_dist = Freq) %>%
    filter(BrayCurtis_dist > 0.49)
  
  bray_df %>% write_csv(paste0("csv_outputs/bray_dissimilar/mockfeeds/",loc,"_",site,"_bray_outliers.csv"))
  
  # Permanova for Bray-Curtis distance
  # I can ask whether all PCR replicates are dissimilar
  set.seed(129)
  bray.perm <- adonis2(bray.distmat~sample, data = community.df.sep, permutations = 999, method = "bray")
  bray.perm
  
  # I can ask whether the pools of DNA extractions are dissimilar
  # set.seed(129)
  # bray.perm.pool <- adonis2(bray.distmat~pool, data = community.df.sep, permutations = 999, method = "bray")
  # bray.perm.pool
  
  ## Multivariate dispersion
  # Does the amount of (within group) variance differ between groups?
  # dispersion<-betadisper(bray.distmat, group=community.df.sep$pool)
  # dispersion.bray <- permutest(dispersion)
  # 
  # disp.plot <- plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
  # 
  ## NMDS
  # using the vegan `metaMDS` function and the relative abundance dataframe rather than the distance matrix
  fishMDS<-metaMDS(relative.abund, distance="bray", k=2, trymax=50, autotransform=TRUE) ##k is the number of dimensions
  #fishMDS ##metaMDS takes eaither a distance matrix or your community matrix (then requires method for 'distance=')
  
  nmds.stress <- stressplot(fishMDS)
  
  # plot the NMDS output in ggplot
  ## pull points from MDS
  NMDS1 <- fishMDS$points[,1] ##also found using scores(fishMDS)
  NMDS2 <- fishMDS$points[,2]
  fish.plot<-cbind(community.df.sep, NMDS1, NMDS2)
  
  # plot ordination
  nmds.plot <- ggplot(fish.plot, aes(NMDS1, NMDS2, color=pool)) +
    geom_point(size = 2) +
    theme_minimal() +
    geom_text(data=fish.plot,aes(NMDS1, NMDS2, label=sample), position = position_jitter(0.005)) +
    annotate("text", x=max(NMDS1), y=min(NMDS2), hjust=0.9, vjust=0.5, label=paste('Stress =',round(fishMDS$stress,3))) +
    ggtitle(paste0(loc," ",site," NMDS"))
  # stat_ellipse(type='t',size =1) ##draws 95% confidence interval ellipses
  

  
  ## outputs
  output.list <- list(
                      #dispersion.bray,
                     # bray.perm.pool,
                      bray.perm)
                    
  
  # plots
  all_plots <- list(
                    #disp.plot,
                    nmds.stress,
                    nmds.plot)
  
  
  return(all_plots)
  
}

