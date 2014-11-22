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
        # get ID´s in the first column for calculating heterozygosity
        gen <- cbind(row.names(genotypes)[subset_rows], gen_sub)
        # heterozygosity data mothers
        het <- mlhWS(gen, na.string = "NA", n.digits = 4)[subset_rows, ]
        # fill in results
        results[k, i] <- cor(as.vector(het$SH),y)
        }
}

results

}