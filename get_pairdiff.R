get_pairdiff <- function(relate, scores, df=F) {
# creates data.frame with  
# 1) pairwise mums relatedness
# 2) pairwise differences in factor scores
        
# input should be: relatedness data frame (lower triangular), 
# data frame with factor scores in columns
# if df=TRUE, get_pairdiff will return a list of score-difference
# dataframes (for each component/factor) with pairwise pc-differences. 

# make sure to have data.frames
relate <- as.data.frame(relate)
scores <- as.data.frame(scores)

# copy similarity matrix and clear
score_mat <- relate
score_mat[, ] <- NA

# get vector of pairwise-rownames
allnames <- vector()
for (i in 1:ncol(relate)) {
        for (k in 1:nrow(relate)) {
                nametemp <- paste(names(relate)[i], 
                                  row.names(relate)[k], sep = "")
                allnames <- append(allnames, nametemp)
        }
}

# roll out as vector
relate_vec <- unlist(relate)
# label the rows
names(relate_vec) <- allnames
# delete na´s
relate_vec <- relate_vec[!is.na(relate_vec)]
# get new row-names vector
pairnamessub <- names(relate_vec)
# create raw data frame 
fac_diff_all <- data.frame("relatedness"= relate_vec)

# construct similarity matrix out of pairwise differences in factors
names <- rownames(relate)

row.names(scores) <- names

fac_diff_mats <- list()

for (z in 1:ncol(scores)) {
        for (i in names) {
                for (k in names) {
                        if (!(is.na(relate[i,k]))) {
                                diff_fac <- abs(scores[i,z] - scores[k,z])
                                score_mat[i,k] <- diff_fac
                        }
                }
        }
        
        # create list of data frames, containing difference matrices per Factor
        fac_diff_mats <- c(fac_diff_mats, list(score_mat))
        
        # turn into vector
        factor_diff <- as.vector(as.matrix(score_mat))
        factor_diff <- factor_diff[!is.na(factor_diff)]
        fac_diff_all <- cbind(fac_diff_all, factor_diff)
        
}

## check argument for what to return
if (df == T) {
        names(fac_diff_mats) <- names(scores)
        return(fac_diff_mats)
} else if (df == F) {
        names(fac_diff_all) <- c("relatedness",names(scores))
        row.names(fac_diff_all) <- pairnamessub
        return(fac_diff_all)
}

        
}