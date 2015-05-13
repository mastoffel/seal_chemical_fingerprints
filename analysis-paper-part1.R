# Fur seal odour encodes colony membership, mother-offspring similarity, relatedness and genetic quality
# Stoffel, M.A., Caspers, B.A., Forcada, J., Giannakara, A., Baier, M.C., 
# Eberhart-Phillips, L.J. , Müller, C. & Hoffman, J.I.

# Any questions? -martin.adam.stoffel@gmail.com

library(vegan)
library(dplyr)
library(Rhh)
library(ggplot2)
library(grid)

# two self-written packages are included in the analysis. 
# Both can be downloaded from my github repo with the following code:
# Note: package "devtools" has to be installed

# library(devtools)
# install_github("mastoffel/minmodelr")
# install_github("mastoffel/g2")

library(minmodelr)
library(g2)

# loading data with observations in rows and compounds in columns---------------
scent_raw <- as.data.frame(t(read.csv(".\\files\\scent_raw.csv", row.names = 1)))

# loading colony, family and age (mothers vs. pups) factors 
factors <- read.csv(".\\files\\factors.csv",row.names=1) 
head(factors)

# pre-treatment ----------------------------------------------------------------

# standardising observation by total, such that every observation
# adds up to 100 % across compounds 
scent_stand <- as.data.frame(t(apply(scent_raw, 1, 
                                     function(x) (x/sum(x)) * 100)))
                             
# Log(x+1) transformation and analysis-ready data frame scent
scent <- log(scent_stand + 1)

# descriptives -----------------------------------------------------------------

# average number of compounds in all profiles
num_comp <- apply(scent, 1, function(x) length(x[x>0]))
mean(num_comp)

# number of compounds in mothers vs. pups
t.test(num_comp[factors$age == 1], num_comp[factors$age == 2])


# colony differences -----------------------------------------------------------

# overall (number of permutations is 1000 instead of 10,000 in the paper)
vegan::anosim(dat = scent, grouping = factors$colony, 
       distance = "bray", permutations = 1000)

# just mothers
vegan::anosim(dat = scent[factors$age == 1, ], grouping = factors$colony, 
       distance = "bray", permutations = 1000)

# just pups
vegan::anosim(dat = scent[factors$age == 2, ], grouping = factors$colony, 
       distance = "bray", permutations = 1000)

# genetic differentiation tested in the programs Arlequin and Structure


# mother-offspring similarity ---------------------------------------------------

# overall
vegan::anosim(dat = scent, grouping = factors$family, 
       distance = "bray", permutations = 1000)

# within colony 1 (Special study beach)
vegan::anosim(dat = scent[factors$colony == 1, ], 
       grouping = factors[factors$colony == 1, ]$family, 
       distance = "bray", permutations = 1000)

# within colony 2 (Freshwater beach)
vegan::anosim(dat = scent[factors$colony == 2, ], 
       grouping = factors[factors$colony == 2, ]$family, 
       distance = "bray", permutations = 1000)

# olfactory similarity vs. geographic distance on special study beach ----------
coord  <- read.csv(".\\files\\coordinates_beach1.csv", row.names=1) 

# converting coordinates to a pairwise euclidian distance matrix
dist_mat <- as.matrix(dist(coord, method = "euclidian"))

# pairwise bray curtis similarity of all individuals from beach 1
scent_bc <- as.matrix(vegdist(as.matrix(scent[factors$colony == 1, ])),
                      method = "bray")

# geographic distance vs. olfactory similarity in mothers
geo_mum <- dist_mat[1:20, 1:20]
scent_mum <- scent_bc[1:20, 1:20]
vegan::mantel(geo_mum, scent_mum, method = "spearman")

# # geographic distance vs. olfactory similarity in pups
geo_pup <- dist_mat[21:40, 21:40]
scent_pup <- scent_bc[21:40, 21:40]
vegan::mantel(geo_pup, scent_pup, method = "spearman")


# Relatedness and overall olfactory similarity ---------------------------------

# load pairwise relatedness based on 41 markers
relatedness <- as.matrix(read.csv(".\\files\\relatedness.csv",row.names=1))

# pairwise bray curtis similarity of all individuals
scent_bc <- 1-(as.matrix(vegdist(as.matrix(scent)), method = "bray"))

# mantel test between relatedness and bray curtis similarity of all individuals
vegan::mantel(relatedness, scent_bc, method = "spearman", permutation = 1000)

# mothers: mantel test between relatedness and bray curtis similarity 
vegan::mantel(relatedness[factors$age == 1, factors$age == 1], 
              scent_bc[factors$age == 1, factors$age == 1],
              method = "spearman", permutation = 1000)

# offspring: mantel test between relatedness and bray curtis similarity 
vegan::mantel(relatedness[factors$age == 2, factors$age == 2], 
              scent_bc[factors$age == 2, factors$age == 2],
              method = "spearman", permutation = 1000)


# Heterozygosity and number of compounds in odour profiles----------------------

# load standardised multilocus heterozygosity (sMLH) based on 41 markers
heterozygosity <- as.matrix(read.csv(".\\files\\heterozygosity.csv", row.names=1))

# number of compounds per individual
num_comp <- as.vector(apply(scent, 1, function(x) length(x[x>0])))

# mothers
het_mum <- heterozygosity[factors$age == 1]
num_comp_mum <- num_comp[factors$age==1]
het_mum_df <- data.frame("het" = het_mum, "comps" = num_comp_mum)

summary(glm(het_mum ~ num_comp_mum))
table <- DelTestVar(het_mum_df)


# pups
het_pup <- heterozygosity[factors$age == 2]
num_comp_pup <- num_comp[factors$age==2]
het_pup_df <- data.frame("het" = het_pup, "comps" = num_comp_pup)

summary(glm(het_pup ~ num_comp_pup))
table <- DelTestVar(het_pup_df)

# resampling plot  -------------------------------------------------------------

# heterozygosity package
library(Rhh)

# loading genotype data
genotypes <- read.table(".\\files\\genotypes.txt", row.names = 1)

# probably put that function in a seperate package
resample_loci <- function(genotypes, y, num_iter = 1000, subset_rows = 1:41) {
# "genotypes" is a matrix with ID´s as row.names and the genotypes in the others (one locus in two columns)
# Y is a vector against which to correlate heterozygosity
# num_iter is the number of resampling per added locus
# subset_rows specifies a subset of the data to compute the resampling test from

        # subset mothers
        y <- y[subset_rows]
        # calculate number of loci
        num_loci <- ncol(genotypes)/2
        # initialize results df
        results <- data.frame(matrix(nrow = num_iter, ncol = num_loci))
        # create index from which to sample loci from genotypes
        ind <- seq(from = 1, to = ncol(genotypes) - 1, by = 2)
        
        for (i in seq(from = 1, to = num_loci)) {
                for (k in seq_along(1:num_iter)) {
                        # create genotype matrix with subset of loci (defined by loop)
                        # loci to take
                        loci_ind <- as.list(sample(ind, i, replace = FALSE))
                        # subset genotypes
                        gen_sub <- lapply(loci_ind, function(x) genotypes[x:(x+1)])
                        gen_sub <- do.call("cbind", gen_sub)
                        # heterozygosity data mothers
                        het <- Rhh::sh(gen_sub)[subset_rows, ]
                        # fill in results
                        results[k, i] <- cor(het,y)
                }
        }   
results
}

# resampling with just 10 iterations instead of 1000 for computational purpose
resample_mums <- resample_loci(genotypes, num_comp_mum, num_iter = 1000, subset_rows = 1:41)

# summarising results: mean and se for the heterozygosity - number of compounds
# correlation, while heterozygosity is estimated by an increasing number of loci
sum_results <- function(resampling_output) {
        mean.cor <- apply(resampling_output,2,mean, na.rm=T)
        sd.cor <- apply(resampling_output,2,sd, na.rm=T)
        se.cor <- sd.cor/(sqrt(sd.cor))
        sum_results <- data.frame(locnum = 1:ncol(resampling_output), 
                                  cormean = mean.cor, corsd = sd.cor, corse = se.cor)
}

results_mums <- sum_results(resample_mums)

# plotting
ggplot(results_mums, aes(x = locnum, y = cormean)) +
        geom_errorbar(aes(ymin = cormean-corse, ymax = cormean+corse),
                      width=0.8, alpha=0.7, size = 0.8) +
        geom_point(size = 3, shape = 16) +
        geom_line(size = 0.8) +
        theme_minimal(base_size = 20) +
        theme(axis.title.x = element_text(vjust= -2 ,size = 22),
              axis.title.y = element_text(vjust=3,size = 22),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              plot.margin = (unit(c(.5, .5, 2, 2), "cm"))) +
        #geom_hline(yintercept=0.305) +
        ylab(expression(r [mean+-se])) +
        xlab("number of loci") +
        ggtitle("Correlation between olfactory complexity and heterozygosity \nestimated by for an increasing number of loci")
        
 
# estimating g2 from an increasing number of loci-------------------------------

resample_g2 <- function(genotypes, num_iter = 1000) {
        # "genotypes" is a matrix with ID´s as row.names and the genotypes in the others (one locus in two columns)
        # Y is a vector against which to correlate heterozygosity
        # num_iter is the number of resampling per added locus
        # subset_rows specifies a subset of the data to compute the resampling test from
        
        # calculate number of loci
        num_loci <- ncol(genotypes)/2
        # initialize results df
        results <- data.frame(matrix(nrow = num_iter, ncol = num_loci))
        # create index from which to sample loci from genotypes
        ind <- seq(from = 1, to = ncol(genotypes) - 1, by = 2)
        # starts with 2 loci, which is the minimum to compute g2
        for (i in seq(from = 2, to = num_loci)) {
                
                for (k in seq_along(1:num_iter)) {
                        # create genotype matrix with subset of loci (defined by loop)
                        # loci to take
                        loci_ind <- as.list(sample(ind, i, replace = FALSE))
                        # subset genotypes
                        gen_sub <- lapply(loci_ind, function(x) genotypes[x:(x+1)])
                        gen_sub <- do.call("cbind", gen_sub)
                        # get ID´s in the first column for calculating heterozygosity
                        g_val <- g2(gen_sub)
                        
                        results[k, i-1] <- g_val
                }
        }
        results      
}

# compute the data frame with all results
resample_g2_results <- resample_g2(genotypes, num_iter = 10)

# sum up the estimated g2 for increasing number of loci 
results_g2 <- sum_results(resample_g2_results)

# plotting
ggplot(results_g2, aes(x = locnum, y = cormean)) +
        geom_errorbar(aes(ymin = cormean-corse, ymax = cormean+corse),
                      width=0.8, alpha=0.7, size = 0.8) +
        geom_point(size = 3, shape = 16) +
        geom_line(size = 0.8) +
        theme_minimal(base_size = 20) +
        theme(axis.title.x = element_text(vjust= -2 ,size = 22),
              axis.title.y = element_text(vjust=3,size = 22),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              plot.margin = (unit(c(.5, .5, 2, 2), "cm"))) +
        #geom_hline(yintercept=0.305) +
        ylab(expression(r [mean+-se])) +
        xlab("number of loci") +
        ggtitle("g2 estimated by for an increasing number of loci")
