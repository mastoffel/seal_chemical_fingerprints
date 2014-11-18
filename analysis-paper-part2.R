# factor analysis --------------------------------------------------------------
library(HDMD)

# factor analysis with 4 factors, promax rotation ------------------------------
scent_fa <- factor.pa.ginv(scent, nfactors = 4, 
                           prerotate = T,rotate = "promax", 
                           scores = T, m = 3)

fa_scores <- as.data.frame(scent_fa$scores)


# glms for association between heterozygosity and factor scores in mothers -----
library(minmodelr)

# backwards deletion testing 

# bind heterozygosity and the factor scores in one data.frame and subset mothers
het_df <- cbind(heterozygosity, fa_scores)[factors$age == 1, ]
het_model <- glm(SH ~., data=het_df)
summary(het_model)

# model reduction via deletion testing
het_reduced <- MinMod(het_df)

# extract data frame
het_reduced_df <- het_reduced[[1]]

# extract reduced model
het_reduced_mod <- het_reduced[[2]]

# deviance explained from reduced model
dev_expl <- (het_reduced_mod$null.deviance - het_reduced_mod$deviance) / het_reduced_mod$null.deviance
summary(het_reduced_mod)

# deletion test for both variables in the reduced model
del_test <- DelTestVar(het_reduced_df)

# sum of factors as variable
het_df$F1F2 <- het_df$F1 + het_df$F2
del_test <- DelTestVar(as.data.frame(cbind(het_df$SH, het_df$F1F2)))
summary(lm(SH ~ F1F2, data = het_df))


# glms for association between relatedness and factor scores in mothers --------
source("get_pairdiff.R")

# get 4 matrices with pairwise differences in factors
fa_diff_mums <- get_pairdiff(relatedness[factors$age == 1, factors$age == 1],
                             fa_scores[factors$age == 1, ], df = F)

# assign pairwise difference factor matrices to names
for (i in seq_along(1:4)) {
        assign(paste("f", i, "_diff", sep=""), fa_diff_mums[, i+1])
}

#using mantel in ecodist
library(ecodist)
rel_dist <- as.dist(relatedness[factors$age == 1, factors$age == 1])
mantel(rel_dist ~ f1_diff + f2_diff + f3_diff + f4_diff, mrank = T, nperm = 1000)
mantel(rel_dist ~ f2_diff + f1_diff + f3_diff + f4_diff, mrank = T)
mantel(rel_dist ~ f3_diff + f2_diff + f1_diff + f4_diff, mrank = T)
mantel(rel_dist ~ f4_diff + f3_diff + f2_diff + f1_diff, mrank = T)

# glms for association between relatedness and factor scores in offspring ------
fa_diff_pups <- get_pairdiff(relatedness[factors$age == 2, factors$age == 2],
                             fa_scores[factors$age == 2, ], df = F)

# assign pairwise difference factor matrices to names
for (i in seq_along(1:4)) {
        assign(paste("f", i, "_diff", sep=""), fa_diff_pups[, i+1])
}

rel_dist <- as.dist(relatedness[factors$age == 2, factors$age == 2])
mantel(rel_dist ~ f1_diff + f2_diff + f3_diff + f4_diff, mrank = T)
mantel(rel_dist ~ f2_diff + f1_diff + f3_diff + f4_diff, mrank = T)
mantel(rel_dist ~ f3_diff + f2_diff + f1_diff + f4_diff, mrank = T)
mantel(rel_dist ~ f4_diff + f3_diff + f2_diff + f1_diff, mrank = T)

# colony differences in factor scores-------------------------------------------

# create data frame
col_df <- cbind(factors["colony"], fa_scores)

# reduce model by deletion testing
col_reduced <- MinMod(col_df)

# deletion test single variable
del_test <- DelTestVar(col_reduced[[1]])

#### relatedness correlation increase with marker number ####






