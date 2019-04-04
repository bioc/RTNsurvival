#' Subgroup Regulon Difference for TNS-class objects
#' 
#' This regulon evaluates differences between regulon activity of subgroups of
#' samples, given a grouping variable. It performs Wilcoxon-Mann-Whitney 
#' (2 subgroups) or Kruskal-Wallis (3+ subgroups) Rank Sum Tests to check whether
#' the activity scores of a given regulon are different between subgroups of samples.
#' 
#' @param tns A A \linkS4class{TNS} object.
#' @param subgroup a character vector. It must be the name of a column in the
#' survivalData featuring the grouping information as a categorical variable.
#' @param regs An optional string vector specifying regulons to use for the analysis.
#' @param pValueCutoff a single numeric value specifying the cutoff for 
#' p-values considered significant.
#' @param pAdjustMethod a single character value specifying the p-value 
#' adjustment method to be used (see 'p.adjust' for details).
#' @param verbose a logical value specifying whether to display messages and
#' progress bar.
#' 
#' @return A TNS-class object with the results of the subgroup regulon difference
#' added to the results slot. To recover the results, use tnsGet(tns, "regulonDifference")
#' 
#' @examples 
#' # load survival data
#' data(survival.data)
#' # load TNI-object
#' data(stni, package = "RTN")
#' 
#' # create TNS object
#' stns <- tni2tnsPreprocess(stni, survivalData = survival.data,
#'                           keycovar = c('Grade','Age'), time = 1, event = 2)
#' stns <- tnsGSEA2(stns)
#' 
#' # run subgroup regulon enrichment analysis
#' stns <- tnsSRD(stns, "ER+")
#' 
#' @docType methods
#' @rdname tnsSRD-methods
#' @aliases tnsSRD
#' @importFrom stats kruskal.test wilcox.test
#' @importFrom utils capture.output
#' @importFrom  dunn.test dunn.test
#' @export
setMethod("tnsSRD", "TNS", 
          function(tns, subgroup, pValueCutoff = 0.05, pAdjustMethod = "BH", 
                   regs = NULL, verbose = TRUE) {
    #-- Basic check + get survData
    .tns.checks(tns, type = "Activity")
    survData <- tnsGet(tns, "survivalData")

    #-- Checks
    .tns.checks(subgroup, survData, type = "subgroup")
    .tns.checks(regs, type = "regs")
    .tns.checks(pValueCutoff, type = "pValueCutoff")
    .tns.checks(pAdjustMethod, type = "pAdjustMethod")
    .tns.checks(verbose, type = "verbose")

    #-- Get data
    regact <- tnsGet(tns, "regulonActivity")$dif
    all_regels <- tni.get(tnsGet(tns, "TNI"), "regulatoryElements")
    
    #-- regs check
    if (!is.null(regs)) {
        if(all(regs %in% names(all_regels))) {
            regact <- regact[,regs] 
        } else {
            stop("All `regs` must be in the the regulatoryElements of the TNI-object.")
        }
    }
    
    #-- subgroup checks
    if (!is.numeric(subgroup) && !is.character(subgroup)) {
        stop("`subgroup` must be numeric or character")
    }
    if (is.numeric(subgroup)) {
        if (subgroup > ncol(survData) || subgroup < 0) {
            stop("`subgroup` must be > 0 and < number of features in the column annotation")
        }
    } else {
        if(!(subgroup %in% colnames(survData))) {
            stop("`subgroup` doesn't correspond to a column in the column annotation")
        }
    }
    gvec <- survData[,subgroup]
    if(!any(duplicated(gvec))) {
        stop("`subgroup` column doesn't contain useful information to divide the samples into subgroups")
    }
    
    grouping <- split(rownames(survData), survData[,subgroup])
    gfac <- as.factor(survData[,subgroup])
    
    #-- Create Results table
    
    if (length(grouping) == 2) {
        res_list <- lapply(colnames(regact), .regulonMW, regact, grouping)
        res_tb <- data.table::rbindlist(res_list)
        res_tb$pAdjusted <- p.adjust(res_tb$pValue, method = pAdjustMethod)
        
        #-- Reorder
        res_tb <- res_tb[order(res_tb$pAdjusted),]
    } else {
        if(verbose) {
            pb <- txtProgressBar(min = 0, max = ncol(regact), style = 3)
        } else {
            pb <- NULL
        }
        res_list <- lapply(colnames(regact), .regulonKW, regact, gfac, 
                           pAdjustMethod, pb, verbose)
        res_tb <- data.table::rbindlist(res_list)
        res_tb$KW.pValue <- p.adjust(res_tb$KW.pValue, method = pAdjustMethod)
        res_tb <- res_tb[order(res_tb$KW.pValue), ]
        
    }

    #-- Add to tns
    tns@results$subgroupDifference[[subgroup]] <- res_tb
    tns@para$srd$pValueCutoff[[subgroup]] <- pValueCutoff
    
    return(tns)
})

.regulonMW <- function(reg, regact, grouping) {
    regAct1 <- regact[grouping[[1]],reg]
    regAct2 <- regact[grouping[[2]],reg]
    mw_res <- wilcox.test(regAct1, regAct2, paired = FALSE)
    res <- list(
        Regulon = reg,
        pValue = mw_res$p.value
    )
    return(res)
}

.regulonKW <- function(reg, regact, grouping, pAdjustMethod, pb, verbose) {
    
    #-- For one regulon
    reg_dES <- split(regact[,reg], grouping)
    
    kw_res <- kruskal.test(reg_dES)
    capture.output(dunn_res <- dunn.test(reg_dES, method = pAdjustMethod))
    dunn_tb <- as.data.frame(dunn_res)
    
    groups <- names(reg_dES)
    res <- lapply(1:length(groups), function(i) {
        idx <- grep(i, dunn_tb$comparisons)
        pval <- max(dunn_tb[idx,"P.adjusted"])
        return(pval)
    })
    names(res) <- paste0(groups, ".specific.pValue")
    res$Regulon <- reg
    res$KW.pValue <- kw_res$p.value
    
    n <- length(groups)
    res <- res[c(n+1,n+2,1:n)]
    
    if(verbose) setTxtProgressBar(pb, which(colnames(regact) %in% reg))
    
    return(res)
}