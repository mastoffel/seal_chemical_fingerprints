# subsets and identification
library(vegan)
library(ggplot2)
require(dplyr)
require(magrittr)
library(vegan)
library(reshape2)

# similarity percentage analysis based on bray-curtis decomposition

# mother-offspring similarity, best substances ---------------------------------

# two most important compounds defining the similarity for each mother offspring
# pair and their average contribution over all pairs

mp_simp <- read.csv(".\\files\\simper_mp_results.csv", colClasses = c("character", "numeric"))

# mother offspring similarity

# overall
anosim(dat = scent[mp_simp$comp], grouping = factors$family, 
       distance = "bray", permutations = 1000)

# within colony 1 (Special study beach)
anosim(dat = scent[factors$colony == 1, mp_simp$comp], 
       grouping = factors[factors$colony == 1, ]$family, 
       distance = "bray", permutations = 1000)

# within colony 2 (Freshwater beach)
anosim(dat = scent[factors$colony == 2, mp_simp$comp], 
       grouping = factors[factors$colony == 2, ]$family, 
       distance = "bray", permutations = 1000)


# colony dissimilarity, best substances ----------------------------------------
simp_colony <- simper(scent, factors$colony)

# getting 15 best substances and their contribution to colony dissimilarity
simp_colony_names <- rownames(summary(simp_colony, ordered = TRUE)[[1]])[1:15]
contribution <- summary(simp_colony, ordered = TRUE)[[1]]$contr[1:15]

# indices of colony substances (58,62,68,74,86,89,90,98,106,107,110,164,181,189,211)
ind_col <- paste(which(names(scent)%in%simp_colony_names), collapse = ",")

# connect to data frame and compute contribution in percent
col_simp <- data.frame(comp = simp_colony_names, contrib = contribution*100, stringsAsFactors = FALSE)
col_simp

# overall (number of permutations is 1000 instead of 10,000 in the paper)
anosim(dat = scent[col_simp$comp], grouping = factors$colony, 
       distance = "bray", permutations = 1000)

# # just mothers
# anosim(dat = scent[factors$age == 1, ], grouping = factors$colony, 
#        distance = "bray", permutations = 1000)
# 
# # just pups
# anosim(dat = scent[factors$age == 2, ], grouping = factors$colony, 
#        distance = "bray", permutations = 1000)



################ BIO-ENV bootstrap procedure ###################################
#### run seperately on multicore server #####
#### aim: resampling test for finding the substances associated with genetic 
#### relatedness. Basic assumption: Each variable will be tested in many different 
#### environments (individuals, other variables), which will prevent spurious 
#### correlations, as the really important substances will occur in best subsets
#### in many different constellations. (see methods section)

# parallel computing using 40 cores, takes some days nevertheless and is just
# shown here.

library(vegan)
library(stringr)
library(dplyr)
library(snow)
library(snowfall)
source("bio.env.R")

# number of cores
ncores <- 2
# subset
scent_mum <- filter(scent, factors$age == 1)
relate_mum <-  relatedness[factors$age == 1, factors$age == 1]

# initialise results vector
all_best <- vector()

# initialise cluster
sfInit(parallel=TRUE, cpus=ncores, type="SOCK")

# export libraries and main function to all cores
sfSource("bio.env.R")
sfLibrary(vegan)
sfLibrary(stringr)
sfLibrary(dplyr)

bootstrap <- function(iter_comp) { # main resampling function
        
        for (i in 1:500) {
                # sample 20 out of 41 mothers, indices
                ind_obs <- sort(sample(1:41, size = 20, replace = F))
                # subset relate_mum and scent_mum
                reltemp <- 1-as.dist(relate_mum[ind_obs, ind_obs])
                abundtemp <- scent_mum[ind_obs, ]
                
                for (i in iter_comp) {
                        # sample 10 compounds
                        index_comps <- sort(sample(1:213, size = 10, replace = F))
                        abundtemp_sub <- abundtemp[, index_comps]
                        # get vector with 0 for null-column and 1 for non-null column
                        nullcomps <- apply(abundtemp_sub, 2, function(x) sum(x>0))
                        abundtemp_sub <- subset(abundtemp_sub, 
                                                subset = c(rep(TRUE,nrow(abundtemp_sub))), 
                                                select = (nullcomps >= 2))
                        # new iteration if too less substances left
                        if (ncol(abundtemp_sub) <= 2) next
                        
                        # main function: bio.env finds subset that mostly correlates
                        # with relatedness
                        results <- bio.env(reltemp, abundtemp_sub, 
                                           var.dist.method = "bray", 
                                           scale.fix = F, scale.var = F)
                        
                        mods <- results$best.model.vars
                        best <- unlist(str_split(mods, ","))
                        all_best <- append(all_best, best)
                        # write(best, file = "best.txt", append = TRUE, sep = " ")
                }
        }
        return(all_best)
}
# export objects
sfExportAll(except = NULL, debug = FALSE)
sfClusterEval(ls())

# create list of 500 iterations for all cores
vals <- list()
for (i in 1:ncores) {
        vals[[i]] <- 1:500
}
# run analysis
# best is a list of all best subsets
best <- sfLapply(vals, bootstrap)
# stop cluster
sfStop()
# bring all results 
results <- unlist(best)

############################## END #############################################

# analysing results from bio-env bootstrap results -----------------------------

# substance occurences are counted in sorted in a table
best_mums <- read.csv("bootstrap_mums.csv",row.names=1)

# subset
scent_mum <- filter(scent, factors$age == 1)
relate_mum <-  1-relatedness[factors$age == 1, factors$age == 1]

# get vectors of best substance names
sub_names_mums <- row.names(best_mums)

statm <- vector()
sigm <- vector()

# compute mantelR for an increasing set of best substances
for (i in 2:100) {
        bc_dist <- vegdist(scent_mum[, sub_names_mums[1:i]], method = "bray")
        mod <- mantel(relate_mum, bc_dist, na.rm = T, method = "spearman")
        statm <- append(statm, mod$statistic)
        sigm <- append(sigm, mod$sig)
}

stat_df <- data.frame(num_comps = 1:length(statm), mantelR = statm)

library(grid)
# simple plot
ggplot(stat_df, aes(x = num_comps, y = mantelR)) +
        stat_smooth(se = FALSE, span = 0.11, size = 1.5, method = "loess") +
        geom_point(colour = "black", size = 3, alpha = 0.4) +
        theme_minimal(base_size = 26) +
        theme(strip.text.x = element_text(vjust=1,size = 18),
              axis.title.x = element_text(vjust= -2 ,size = 28),
              axis.title.y = element_text(vjust=3,size = 28),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              plot.margin = (unit(c(.5, .5, 2, 2), "cm"))) +
        # scale_x_continuous(breaks=c(seq(from = 0.8, to = 1.20, by = 0.1))) +
        #geom_text(aes(0.85,80, label="(a) r = 0.34, p = 0.027"),size=4) +
        xlab("cumulative substances from bootstrap") +
        ylab("mantelR") 



# indices of the 10 best compounds associated with relatedness -----------------
comp_ind_m <- c(36,52,86,88,96,103,110,203,206,207)

# bray curtis similarity matrix based on this 10 compounds
scent_bc <- 1-(as.matrix(vegdist(as.matrix(scent[factors$age == 1, comp_ind_m])),
                      method = "bray"))
# relatedness matrix
rel_m <- relatedness[1:41, 1:41]

# mantel test for association between both
vegan::mantel(rel_m, scent_bc, method = "spearman", permutation = 1000, na.rm = TRUE)



