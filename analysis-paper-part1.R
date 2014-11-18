# analysis script as appendix for publication

# loading data with observations in rows and compounds in columns---------------
scent_raw <- as.data.frame(t(read.csv("scent_raw.csv", row.names = 1)))

# loading colony and family factor
factors <- read.csv("factors.csv",row.names=1) 
names(factors)  <- c("colony", "family")

# adding factor for age (mothers = 1, pups = 2)
factors$age <- c(rep(1,41), rep(2,41))

# pre-treatment ----------------------------------------------------------------

# standardising observation by total, such that every observation
# adds up to 100 % across compounds 
scent_stand <- as.data.frame(t(apply(scent_raw, 1, 
                                     function(x) (x/sum(x)) * 100)))
                             
# Log(x+1) transformation
scent <- log(scent_stand + 1)


# descriptives -----------------------------------------------------------------

# average number of compounds in all profiles
num_comp <- apply(scent, 1, function(x) length(x[x>0]))
mean(num_comp)

# number of compounds in mothers vs. pups
t.test(num_comp[factors$age == 1], num_comp[factors$age == 2])


# colony differences -----------------------------------------------------------

library(vegan)

# overall (number of permutations is 1000 instead of 10,000 in the paper)
anosim(dat = scent, grouping = factors$colony, 
       distance = "bray", permutations = 1000)

# just mothers
anosim(dat = scent[factors$age == 1, ], grouping = factors$colony, 
       distance = "bray", permutations = 1000)

# just pups
anosim(dat = scent[factors$age == 2, ], grouping = factors$colony, 
       distance = "bray", permutations = 1000)

# genetic differentiation tested in the programs Arlequin and Structure


# mother-offspring similarity ---------------------------------------------------

# overall
anosim(dat = scent, grouping = factors$family, 
       distance = "bray", permutations = 1000)

# within colony 1 (Special study beach)
anosim(dat = scent[factors$colony == 1, ], 
       grouping = factors[factors$colony == 1, ]$family, 
       distance = "bray", permutations = 1000)

# within colony 2 (Freshwater beach)
anosim(dat = scent[factors$colony == 2, ], 
       grouping = factors[factors$colony == 2, ]$family, 
       distance = "bray", permutations = 1000)

# olfactory similarity vs. geographic distance on special study beach ----------
library(vegan)
coord  <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                   "projects\\sealscent\\data_files\\",
                                   "Rdata\\csv_files\\",
                                   "coordinates_beach1.csv", sep = ""),
                                   row.names=1) 

# convert coordinates to pairwise euclidian distance matrix
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


# Relatedness and overall olfactory similarity ------------------------------------

# load pairwise relatedness based on 41 markers
relatedness <- as.matrix(read.csv("relatednessnew.csv",row.names=1))

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
heterozygosity <- as.matrix(read.csv("heterozygosity.csv", row.names=1))

# number of compounds per individual
num_comp <- as.vector(apply(scent, 1, function(x) length(x[x>0])))

# mothers
het_mum <- heterozygosity[factors$age == 1]
num_comp_mum <- num_comp[factors$age==1]

summary(glm(het_mum ~ num_comp_mum))

# pups
het_pup <- heterozygosity[factors$age == 2]
num_comp_pup <- num_comp[factors$age==2]

summary(glm(het_pup ~ num_comp_pup))


# resampling plot 1-------------------------------------------------------------

# g2 and plot-------------------------------------------------------------------







