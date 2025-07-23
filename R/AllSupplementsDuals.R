
##------------------------------------------------------------------------------
.survstatsDuals <- function(regulonActivity, survData, regs, excludeMid){
  
  #--- survData
  survData <- survData[,c("time","event")]
  survData <- survData[complete.cases(survData),]
  
  #--- tabstatus
  tabstatus <- regulonActivity$status[rownames(survData), regs]
  tabstatus <- data.frame(tabstatus, check.names = FALSE)
  idx <- rowSums(tabstatus==0)>0
  if(excludeMid){
    tabstatus <- tabstatus[!idx,]
    sections <- 1:4
    names(sections) <- c("-|-","-|+","+|-","+|+")
  } else {
    tabstatus[idx,] <- 0
    sections <- 1:5
    names(sections) <- c("0|0","-|-","-|+","+|-","+|+")
  }
  survData <- survData[rownames(tabstatus),]
  
  #--- regstatusChar
  tp1 <- as.character(tabstatus[[regs[1]]])
  tp1[tp1=="1"]  <- "+"
  tp1[tp1=="-1"] <- "-"
  tp2 <- as.character(tabstatus[[regs[2]]])
  tp2[tp2=="1"] <- "+"
  tp2[tp2=="-1"] <- "-"
  regstatusChar <- paste(tp1,tp2, sep = "|")
  names(regstatusChar) <- rownames(tabstatus)
  
  #--- regstatusNum
  regstatusNum <- sections[regstatusChar]
  names(regstatusNum) <- names(regstatusChar)
  nclass <- length(unique(regstatusNum))
  ddt <- survData[names(regstatusNum), ]
  ddt$class <- regstatusNum
  
  #---log-rank test
  survtb <- c(ChiSquare=NA, Pvalue=NA)
  survdf <- NA
  survft <- NA
  if(nclass > 1){
    survft <- survfit(Surv(time, event) ~ class, data = ddt)
    survdf <- survdiff(Surv(time, event) ~ class, data = ddt)
    pval <- 1 - pchisq(survdf$chisq, length(survdf$n) - 1)
    survtb[] <- c(survdf$chisq,pval)
  }
  res <- list(kmTable=survtb, survdiff=survdf, survfit=survft)
  return(res)
}

##------------------------------------------------------------------------------
.survplotDuals <- function(model, regulonActivity, survData, regs, endpoint,
                           excludeMid, ylab, xlab, colorPalette){
  
  #-- set endpoint
  survData$event[survData$time > endpoint] <- 0
  survData$time[survData$time > endpoint] <- endpoint
  
  #--- survData
  survData <- survData[,c("time","event")]
  survData <- survData[complete.cases(survData),]
  
  #--- tabstatus
  tabstatus <- regulonActivity$status[rownames(survData), regs]
  tabstatus <- data.frame(tabstatus)
  idx <- rowSums(tabstatus==0)>0
  if(excludeMid){
    tabstatus <- tabstatus[!idx,]
    sections <- 1:4
    names(sections) <- c("-|-","-|+","+|-","+|+")
  } else {
    tabstatus[idx,] <- 0
    sections <- 1:5
    names(sections) <- c("0|0","-|-","-|+","+|-","+|+")
  }
  survData <- survData[rownames(tabstatus),]
  
  #--- regstatusChar
  tp1 <- as.character(tabstatus[[regs[1]]])
  tp1[tp1=="1"]  <- "+"
  tp1[tp1=="-1"] <- "-"
  tp2 <- as.character(tabstatus[[regs[2]]])
  tp2[tp2=="1"] <- "+"
  tp2[tp2=="-1"] <- "-"
  regstatusChar <- paste(tp1,tp2, sep = "|")
  names(regstatusChar) <- rownames(tabstatus)
  
  #--- regstatusNum
  regstatusNum <- sections[regstatusChar]
  names(regstatusNum) <- names(regstatusChar)
  
  #-- get colors
  if (.is_singleString(colorPalette)){
    if (colorPalette == "reds"){
      cols <- pal1(4)
    } else if (colorPalette == "blues"){
      cols <- pal2(4)
    } else if (colorPalette %in% c("redblue","bluered")){
      cols <- pal3(4)
    }
    if(colorPalette!="redblue") 
      cols <- rev(cols)
  } else {
    cols <- colorPalette
  }
  if(excludeMid){
    cols <- cols[-3]
  } else {
    cols <- cols[c(3,1,2,4,5)]
  }
  names(cols) <- paste0("class=", sections)
  #--- adjusting graphical parameters
  op <- par(no.readonly=TRUE)
  par(mar = c(4, 4, 2, 2), mgp = c(2, 0.4, 0), cex=0.66)
  if(endpoint/3==round(endpoint/3)){
    length.out=4
  } else {
    length.out=5
  }
  labs <- as.integer(seq(0, endpoint, length.out = length.out))
  if (!endpoint %in% labs) labs <- pretty(c(0, endpoint))
  #-- survival analysis
  ddt <- survData[names(regstatusNum), ]
  ddt$class <- regstatusNum
  # res1 <- survfit(Surv(time, event) ~ class, data = ddt)
  res1 <- model$survfit
  plot(res1, col = cols[names(res1$strata)], lwd = 1.8, axes = FALSE, cex = 0.5,
       mark.time = TRUE, ylab = "", xlab = "", xlim = range(labs))
  title(ylab = ylab, adj = 0.5, cex.lab = 1.2, mgp = c(2.2, 0.4, 0))
  title(xlab = xlab, adj = 0.5, cex.lab = 1.2, mgp = c(1.6, 0.4, 0))
  axis(1, at = labs, labels = labs, tcl = -0.2, las = 1, lwd = 1.8, cex.axis = 1.2)
  axis(2, tcl = -0.2, las = 2, lwd = 1.8, cex.axis = 1.2)
  #--- get log-rank pval
  # res2 <- survdiff(Surv(time, event) ~ class, data = ddt)
  # pval <- 1 - pchisq(res2$chisq, length(res2$n) - 1)
  res2 <- model$survdiff
  pval <- model$pAdjustInteraction
  #---legends
  legs <- names(sections)
  names(legs) <- paste0("class=", sections)
  strata <- res2$n
  events <- res2$obs
  names(events) <- names(strata)
  strata <- strata[names(legs)]
  events <- events[names(legs)]
  names(events) <- names(strata) <- names(legs)
  strata[is.na(strata)] <- 0
  events[is.na(events)] <- 0
  if(!excludeMid)legs[1] <- "undetermined"
  legs <- paste(legs, " : ", strata, " (", events,")", sep = "")
  par(xpd=TRUE)
  legend("bottomleft", legend = rev(legs), col = rev(cols), bty = "n", pch = 15, 
    title = paste(paste(regs, collapse = " | "), "\nregulon status",sep=""),
    inset = c(0.01,0), cex = 0.8, pt.cex = 1, title.adj = 0, adj = 0)
  pval <- paste("Logrank P: ", format(pval, digits = 3, scientific = TRUE))
  legend("topright", cex = 1, legend = pval, bty = "n", inset = c(0,-0.05))
  par(op)
}
.namesCorrect <- function(regs) {
  xregs <- gsub("-|\\+|\\.|:|\\*|,|;", "_", regs)
  xregs <- gsub("\\s", "", xregs)
  xregs
}

##------------------------------------------------------------------------------
.getSurvplotCols <- function(regulonActivity, regs, excludeMid, colorPalette){
  
  #--- tabstatus
  tabstatus <- regulonActivity$status[, regs]
  tabstatus <- data.frame(tabstatus)
  idx <- rowSums(tabstatus==0)>0
  if(excludeMid){
    tabstatus <- tabstatus[!idx,]
    sections <- 1:4
    names(sections) <- c("-|-","-|+","+|-","+|+")
  } else {
    tabstatus[idx,] <- 0
    sections <- 1:5
    names(sections) <- c("0|0","-|-","-|+","+|-","+|+")
  }
  
  #--- regstatusChar
  tp1 <- as.character(tabstatus[[regs[1]]])
  tp1[tp1=="1"]  <- "+"
  tp1[tp1=="-1"] <- "-"
  tp2 <- as.character(tabstatus[[regs[2]]])
  tp2[tp2=="1"] <- "+"
  tp2[tp2=="-1"] <- "-"
  regstatusChar <- paste(tp1,tp2, sep = "|")
  names(regstatusChar) <- rownames(tabstatus)
  
  #--- regstatusNum
  regstatusNum <- sections[regstatusChar]
  names(regstatusNum) <- names(regstatusChar)
  
  #-- get colors
  if (.is_singleString(colorPalette)){
    if (colorPalette == "reds"){
      cols <- pal1(4)
    } else if (colorPalette == "blues"){
      cols <- pal2(4)
    } else if (colorPalette %in% c("redblue","bluered")){
      cols <- pal3(4)
    }
    if(colorPalette!="redblue") 
      cols <- rev(cols)
  } else {
    cols <- colorPalette
  }
  if(excludeMid){
    cols <- cols[-3]
  } else {
    cols <- cols[c(3,1,2,4,5)]
  }
  regstatusCol <- cols[regstatusNum]
  names(regstatusCol) <- names(regstatusNum)
  res <- list(numb=regstatusNum, char=regstatusChar, 
              cols=regstatusCol)
  return(res)
}

#-------------------------------------------------------------------------------
# IMPORTANT: This function automatically recodes preventive to risk factors 
# in order to assess interaction on additive scale. For background information,
# please refer to Knol et al. 2011 (https://doi.org/10.1007/s10654-011-9554-9).
# NOTE: This function is expecting a 'coxph' object with exact two factors and
# one interaction term, for example, generated from a call similar to:
# coxph(formula = Surv(time, event) ~ X1 * X2, data = dataset), in which
# 'X1' and 'X2' are numerical standardized variables (i.e. 'X1' and 'X2' are
# on a comparable scale and centered around zero).
.interaction.stats <- function(model, method="additive", conf.level = 0.95){
  if(!is(model,"coxph")) 
    stop("'model' must be a coxph object.", call. = FALSE)
  icoef <- c(1, 2, 3)
  coefficients <- summary(model)$coefficients
  bt1 <- as.numeric(coefficients[icoef[1], 1])
  bt2 <- as.numeric(coefficients[icoef[2], 1])
  bt3 <- as.numeric(coefficients[icoef[3], 1])
  bt1.se <- coefficients[icoef[1], 3]
  bt2.se <- coefficients[icoef[2], 3]
  bt3.se <- coefficients[icoef[3], 3]
  pvals <- coefficients[, 5]
  if(method=="additive"){
    cov.mat <- vcov(model)
    var12 <- cov.mat[icoef[1], icoef[2]]
    var13 <- cov.mat[icoef[1], icoef[3]]
    var23 <- cov.mat[icoef[2], icoef[3]]
    recode <- 1
    if(bt1<0 || bt2<0){
      if(sign(bt1)==sign(bt2)){
        bt1 <- -bt1; bt2 <- -bt2
        var13 <- -var13
        var23 <- -var23
      } else if(bt1<0){
        bt1 <- -bt1; bt3 <- -bt3
        var12 <- -var12
        var23 <- -var23
        recode <- -1
      } else {
        bt2 <- -bt2; bt3 <- -bt3
        var12 <- -var12
        var13 <- -var13
        recode <- -1
      }
    }
    p <- 1 - ((1 - conf.level)/2)
    z <- qnorm(p, mean = 0, sd = 1)
    rrab <- exp(bt1 + bt2 + bt3)
    reri <- rrab - exp(bt1) - exp(bt2) + 1
    a1 <- rrab - exp(bt1)
    a2 <- rrab - exp(bt2)
    a3 <- rrab
    reri.var <- (a1^2 * bt1.se^2) + (a2^2 * bt2.se^2) + (a3^2 * bt3.se^2) + 
      (2 * a1 * a2 * var12) + (2 * a1 * a3 * var13) + (2 * a2 * a3 * var23)
    reri.se <- sqrt(reri.var)
    reri.l <- reri - (z * reri.se)
    reri.u <- reri + (z * reri.se)
    reri.p <- pnorm(reri/reri.se, lower.tail = FALSE)
    res <- c(stat = reri*recode, lower = reri.l, upper = reri.u, pval=reri.p)
  } else {
    p.mult <- as.numeric(pvals[icoef[3]])
    mult <- as.numeric(exp(bt3))
    mult.ci <- confint(object = model, parm = icoef[3]) #suppressMessages
    mult.l <- as.numeric(exp(mult.ci[1]))
    mult.u <- as.numeric(exp(mult.ci[2]))
    res <- c(stat = mult, lower = mult.l, upper = mult.u, pval=p.mult)
  }
  return(res)
}

##------------------------------------------------------------------------------
# remove 'strong' dependencies using phi coefficient
.removeDependencies <- function(dualtb, regulonActivity, th=0.5){
  status <- regulonActivity$status
  status[status==0] <- NA
  phi <- sapply(seq_len(nrow(dualtb)), function(i){
    rg <- dualtb[i, ]
    tb <- table(status[,rg$reg1], status[,rg$reg2])
    # chisq.test(tb)$p.value
    if(prod(dim(tb)) == 4){
      res <- .phi(tb)
    } else {
      res <- 0
    }
    res
  })
  dualtb <- dualtb[phi < th, ]
  return(dualtb)
}
.phi <- function(tb) {
  r.sum <- rowSums(tb)
  c.sum <- colSums(tb)
  total <- sum(r.sum)
  r.sum <- r.sum/total
  c.sum <- c.sum/total
  v <- prod(r.sum, c.sum)
  res <- (tb[1, 1]/total - c.sum[1] * r.sum[1])/sqrt(v)
  res <- ifelse(is.na(res), 0, res)
  res <- as.numeric(res)
  res <- abs(res)
  return(res)
}

