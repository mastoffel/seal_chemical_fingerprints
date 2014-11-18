# subsets and identification
library(vegan)

# similarity percentage analysis based on bray-curtis decomposition

# colony dissimilarity, best substances ----------------------------------------
simp_colony <- simper(scent, factors$colony)

# getting 15 best substances and their contribution to colony dissimilarity
simp_colony_names <- rownames(summary(simp_colony, ordered = TRUE)[[1]])[1:15]
contribution <- summary(simp_colony, ordered = TRUE)[[1]]$contr[1:15]

# mother-offspring similarity, best substances ---------------------------------
simp_family <- simper(scent, factors$family)

# contribution of each substance to within-pair similarity not available in vegan
# thus extracted from primer
source("extract_simper.R")
simp <- extract_simper((paste("C:\\Users\\Martin\\Studium\\",
                              "projects\\sealscent\\data_files\\",
                              "Rdata\\csv_files\\",
                              "allsimper.csv", sep = ""))) 

toptwo <- sapply(simp, function(x) x$comp[1:2])
vars <- sapply(simp, function(x) x$var[1:length(var)])

library(plyr)
library(dplyr)

# all best simper substances plus explained variances in one data.frame
allsimp <- ldply(simp, data.frame)

allsimp$comp <- as.factor(allsimp$comp)
allsimp$var <- as.numeric(allsimp$comp)

simpsum <- allsimp %>%
        group_by(comp) %>%
        summarise(
                meanvar = mean(var, na.rm = TRUE),
                meansd = sd(var, na.rm = TRUE))


hist(table(toptwo), breaks=c(1:41))

comps <- sort(table(toptwo), decreasing = TRUE)

# getting top compounds
topcomp <- names(sort(table(toptwo), decreasing = TRUE))
top_mp_ind <- which(names(scent) %in% topcomp)
