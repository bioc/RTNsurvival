#' Semi-supervised hierarchical clustering
#' 
#' Semi-supervised hierarchical clustering by chosen groups with hclust.
#'
#' @param data a data.frame to be clustered by rows
#' @param groups a list of vectors. If we unlist(groups), all elements must be
#'   present in the rownames of data. Each vector in the list will be treated as
#'   a separate group for the hierarchical clustering, and rejoined in order at
#'   the end.
#' @param dist_method a distance computation method. Must be one of "euclidean", 
#' "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman"
#' @param dist_p the power of the Minkowski distance, if chosen dist_method is "minkowski"
#' @param hclust_method an agglomeration method. Should be a method supported by
#' hclust, one of:  "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), 
#' "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' 
#' @return hclust_semisupervised returns a list. The first element of the list
#' is the data, reordered so that the merged hclust object will work. The second
#' element is the result of the semi-supervised hierarchical clustering.
#' @importFrom stats as.dendrogram dist as.dist as.hclust cor dist hclust
#' @export

hclust_semisupervised <- function(data, groups, dist_method = "euclidean",
                      dist_p = 2, hclust_method = "complete") {
    #-- Group checks
    if (!all(unlist(groups) %in% rownames(data))) {
        stop("vectors in `groups` must contain rownames of `data`")
    }
    if(anyDuplicated(unlist(groups))){
        stop("`groups` can't have elements duplicated within or in different 
             groups")
    }
    
    #-- dist_method checks
    alldists <- c("euclidean", "maximum", "manhattan", "canberra", "binary", 
                  "minkowski", "pearson", "spearman")
    if(!(dist_method %in% alldists)) {
        stop("`dist_method` must be one of: `euclidean`, `maximum`, `manhattan`, 
             `canberra`, `binary`, `minkowski`, `pearson`, `spearman` or `kendall`.")
    }
    
    #-- Get groups with 1 member
    g_size <- sapply(groups, length)
    if (any(g_size == 1)) {
        # s_groups <- groups[g_size == 1]
        groups <- groups[g_size != 1]
    }
    
    #-- Make distance matrices
    if (dist_method %in% c("pearson", "spearman", "kendall")) {
        distlist <- lapply(groups, function(group) {
            as.dist(1 - cor(t(data[group,]), method = dist_method, 
                            use="pairwise.complete.obs"))
        })
    } else {
        distlist <- lapply(groups, function (group) {
            dist(data[group,], method = dist_method, p = dist_p)
        })
    }
    #-- Use hclust
    hclist <- lapply(distlist, hclust, method = hclust_method)
    
    hc <- .merge_hclust(hclist)
    
    # #-- Join groups with one element
    # if(exists("s_groups")) {
    #     s_groups <- unlist(s_groups)
    #     s_hc <- hclust(dist(data[s_groups,]))
    #     hc <- .merge_hclust(list(hc, s_hc))
    #     data_reordered <- data[c(unlist(groups), s_groups),]
    # } else {
    #     data_reordered <- data[unlist(groups),]
    # }
    
    data_reordered <- data[unlist(groups),]

    return(list(data = data_reordered, 
                hclust = hc))
}

.merge_hclust <- function(hclist) {
    #-- Check
    if(!is.list(hclist)) {
        stop("`hclist` must be a list.")
    }
    if(!all(sapply(hclist, class) == "hclust")){
        stop("All objects in `hclist` must be `hclust-class`")
    }
    
    #-- Merge
    d <- hclist[[1]]
    d$height <- d$height/max(d$height)
    d <- as.dendrogram(d)
    for (i in 2:length(hclist)) {
        dd <- hclist[[i]]
        dd$height <- dd$height/max(dd$height)
        d <- merge(d, as.dendrogram(dd))
    }
    as.hclust(d)
}
