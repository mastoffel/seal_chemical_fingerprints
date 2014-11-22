# factor analysis --------------------------------------------------------------
library(HDMD)
library(minmodelr)
source("get_pairdiff.R")
# factor analysis with 4 factors, promax rotation ------------------------------
scent_fa <- factor.pa.ginv(scent, nfactors = 4, 
                           prerotate = T,rotate = "promax", 
                           scores = T, m = 3)
fa_scores <- as.data.frame(scent_fa$scores)

# screeplot, 4 factors left to the "scree"
plot(scent_fa$values[1:8], type="b", ylab = "eigenvalue", xlab = "factor", 
     main = "Screeplot")

# distribution of factor scores
df <- cbind(fa_scores, factors["colony"])
df$colony <- as.factor(df$colony)

for (i in c(1,2,4)) {
plot <- ggplot(df, aes_string(x = paste("F", i, sep = ""), fill = "colony")) +
        geom_density(alpha=0.8, size=0.5, aes(fill = colony),adjust=1.5) +
        scale_fill_manual(values = c("blue","red")) +
        guides(fill=guide_legend(title=NULL)) +
        theme_minimal(base_size = 20) + 
        theme(legend.position="none") +
        scale_x_continuous(breaks = c(seq(from = -1, to = 6, by = 1))) + 
        scale_y_continuous(breaks = c(seq(from = 0, to = 1.4, by = 0.4))) +  
        #scale_y_continuous(breaks = c(seq(from = 0, to = 1, by = 0.2))) +
        xlab(paste("Factor", i, sep = " ")) +
        ylab("Density")

assign(paste("f", i, "_plot", sep = ""), plot)
}

# using multiplot function from cookbook-r.com for plotting multiple ggplots
source("multiplot.R")
multiplot(f1_plot, f2_plot, f4_plot, cols = 2)

# glms for association between heterozygosity and factor scores in mothers -----

# bind heterozygosity and the factor scores in one data.frame and subset mothers
het_df <- cbind(heterozygosity, fa_scores)[factors$age == 1, ]
het_model <- glm(SH ~., data=het_df)
summary(het_model)

# model simplification via deletion testing. See ?MinMod
het_reduced <- MinMod(het_df)

# extract data frame
het_reduced_df <- het_reduced[[1]]

# extract reduced model
het_reduced_mod <- het_reduced[[2]]

# deletion testing for both variables in the reduced model. See ?DelTestVar
table <- DelTestVar(het_reduced_df)

# deviance explained by the reduced model
dev_expl <- (het_reduced_mod$null.deviance - het_reduced_mod$deviance) / het_reduced_mod$null.deviance

summary(het_reduced_mod)

# sum of factors as variable
het_df$F1F2 <- het_df$F1 + het_df$F2
table <- DelTestVar(as.data.frame(cbind(het_df$SH, het_df$F1F2)))
summary(lm(SH ~ F1F2, data = het_df))

# glms for association between relatedness and factor scores in mothers --------

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
col_reduced_df <- col_reduced[[1]]

# dev explained
dev_expl <- (col_reduced_df$null.deviance - col_reduced_df$deviance) / col_reduced_df$null.deviance

# deletion test single variable
table <- DelTestVar(col_reduced[[1]])








