# run masigpro, create design matrix on the way
require(dplyr)
require(maSigPro)
# get kinetic genes based on 2 conseq days thres criterion
# needs already lognormalized data frame

df_ids <- read_csv("data/annotation/annot_genes_to_id.csv")

get_kin_genes <- function(df, thres = 1.0){
  df <- as.data.frame(df)
  
  # filter genes that are above or below fold-change thres at two conseq. time points
  # shift data frame and compare arrays
  df1 <- df[ ,1:(ncol(df)-1)]
  df2 <- df[ ,2:ncol(df)]
  
  # threshold criterion
  df1 <- abs(df1) >= thres
  df2 <- abs(df2) >= thres
  df3 <- df1 & df2
  df4 <- rowSums(df3)
  # there should be two consecutive TRUEs so rowsums should be greater or equal to 2
  conseq_tp <- 2
  df4 <- df4 >= conseq_tp
  
  df$kinetic <- FALSE
  df$kinetic[df4] <- TRUE
  return(df)
}

get_masig_genes <- function(tfit, thres){
  allSigs <- get.siggenes(tfit, rsq = thres, vars = "all")
  kinetic <- allSigs$summary
  kinetic_genes <- convert_id(kinetic)
  return(kinetic_genes)
}

# take a list of ids and convert to genes
convert_id <- function(ids, annot_df = df_ids){
  genes <- df_ids$gene[df_ids$id %in% ids]
  return(genes)
}



make_des_mat <- function(df, cell, keep_timepoints = NULL){
  
  # cell should be a char vec that includes all cell types to be included separated 
  # e.g. c("Th1", "Th2"), the first element is used as control
  designInfo <- df[,1:3]
  designInfo <- as.data.frame(designInfo)
  
  # kick out timepoints that sould be ommitted
  if(length(keep_timepoints) != 0){
    designInfo <- designInfo[designInfo[,"timePoint"] %in% keep_timepoints, ]
  }  
  
  # add information about replicates
  colnames <- c("ThMix","Th1","Th2", "Th0")
  for(i in colnames){
    # add new column with title and fill with zeros
    designInfo[i] <- 0
    # add 1 if coltitle pops up in first column
    designInfo[which(designInfo[,1] == i), i] <- 1
  }
  
  # kick out celltype column now
  designInfo <- designInfo[,-1]
  # reorder if to put control group in 3rd column
  designInfo <- designInfo %>% select("timePoint", "replicate", everything())
  
  # keep only cells indicated by cell vector
  grepcells <- paste(cell, collapse = "|")
  designInfo <- designInfo[grep(grepcells, rownames(designInfo)), c("timePoint","replicate",cell)]
  # fill replicate column with sequence 1,2,.. but omit places with time point zero
  # this need to be done at last step, otherwise numbering will come out wrong
  sequence <- seq(1:length(designInfo[,1]))
  designInfo[,"replicate"] = 1
  # indexing vector for sequence (2 is correct)
  j <- 2
  for (i in seq_along(sequence)){
    if (designInfo[i,"timePoint"] != 0){
      designInfo[i,"replicate"] <- j
      j <- j+1
    } else {
      # for a commong starting point at t=0 set all replicate information to 1
      designInfo[i, 2:ncol(designInfo)] <- 1
    }
    
  }  
  
  designInfo <- as.matrix(designInfo)
  designInfo <- apply(designInfo,1:2,as.numeric)
  
  return(designInfo)
}



get_degs <- function(df, fdata, cell, n_degree = 4, alpha = 0.05, Q = 0.05){
  
  # get design info for cell type
  if (cell %in% c("Th1", "Th2", "ThMix", "Th1_Th2_ThMix", "Th1_Th2")){
    use_th0 = F
    toRemove <- grep("Th0",colnames(df))
    
    # remove th0 columns
    df <- df[,-toRemove]
    designInfo <- make_des_mat(fdata, cell, use_th0 = use_th0)  
  } else {
    use_th0 = T
    designInfo <- make_des_mat(fdata, cell, use_th0 = use_th0)  
  }
  
  # reduce df and design to appropriate cell type
  if ((cell!="Th1_Th2_ThMix") & (cell!="Th1_Th2_Th0") & (cell!= "Th1_Th2")){
    # subset expr data
    keep <- grep(cell, colnames(df))
    df <- df[,keep]
  }
  
  # run masigpro
  designMatrix <- make.design.matrix(designInfo, degree=n_degree)
  
  myFit <- p.vector(df, designMatrix, Q = Q, MT.adjust = "BH")
  return(myFit)
}  


# this function is same as pvector but returns pval for all (not set to NA)
p.vector_modified <- function (data, design, Q = 0.05, MT.adjust = "BH", min.obs = 6, 
                               counts = FALSE, family = NULL, theta = 10, epsilon = 1e-05, 
                               item = "gene") 
{
  if (is.data.frame(design) || is.matrix(design)) {
    dis <- design
    groups.vector = NULL
    edesign = NULL
  }
  else if (is.list(design)) {
    dis <- as.data.frame(design$dis)
    groups.vector <- design$groups.vector
    edesign <- design$edesign
  }
  if (is.null(family)) {
    if (!counts) {
      family = gaussian()
    }
    if (counts) {
      family = negative.binomial(theta)
    }
  }
  dat <- as.matrix(data)
  dat <- dat[, as.character(rownames(dis))]
  G <- nrow(dat)
  count.na <- function(x) (length(x) - length(x[is.na(x)]))
  dat <- dat[apply(dat, 1, count.na) >= min.obs, ]
  sumatot <- apply(dat, 1, sum)
  counts0 <- which(sumatot == 0)
  if (length(counts0) > 0) {
    dat <- dat[-counts0, ]
  }
  g <- dim(dat)[1]
  n <- dim(dat)[2]
  p <- dim(dis)[2]
  p.vector <- vector(mode = "numeric", length = g)
  for (i in 1:g) {
    y <- as.numeric(dat[i, ])
    div <- c(1:round(g/100)) * 100
    if (is.element(i, div)) 
      print(paste(c("fitting ", item, i, "out of", g), 
                  collapse = " "))
    model.glm <- glm(y ~ ., data = dis, family = family, 
                     epsilon = epsilon)
    if (model.glm$null.deviance == 0) {
      p.vector[i] = 1
    }
    else {
      model.glm.0 <- glm(y ~ 1, family = family, epsilon = epsilon)
      if (family$family == "gaussian") {
        test <- anova(model.glm.0, model.glm, test = "F")
        if (is.na(test[6][2, 1])) {
          p.vector[i] = 1
        }
        else {
          p.vector[i] = test[6][2, 1]
        }
      }
      else {
        test <- anova(model.glm.0, model.glm, test = "Chisq")
        if (is.na(test[5][2, 1])) {
          p.vector[i] = 1
        }
        else {
          p.vector[i] = test[5][2, 1]
        }
      }
    }
  }
  p.adjusted <- p.adjust(p.vector, method = MT.adjust, n = length(p.vector))
  #genes.selected <- rownames(dat)[which(p.adjusted <= Q)]
  # changing this should give me all pvalues
  genes.selected <- rownames(dat)
  FDR <- sort(p.vector)[length(genes.selected)]
  SELEC <- as.matrix(as.data.frame(dat)[genes.selected, ])
  if (nrow(SELEC) == 0) 
    print("no significant genes")
  p.vector <- as.matrix(p.vector)
  rownames(p.vector) <- rownames(dat)
  colnames(p.vector) <- c("p.value")
  output <- list(SELEC, p.vector, p.adjusted, G, g, FDR, nrow(SELEC), 
                 dis, dat, min.obs, Q, groups.vector, edesign, family)
  names(output) <- c("SELEC", "p.vector", "p.adjusted", "G", 
                     "g", "FDR", "i", "dis", "dat", "min.obs", "Q", "groups.vector", 
                     "edesign", "family")
  output
}


# T.fit but done write NANS, therefore, change alpha valua within code
T.fit_modified <- function (data, design = data$dis, step.method = "backward", 
                            min.obs = data$min.obs, alfa = data$Q, nvar.correction = FALSE, 
                            family = gaussian(), epsilon = 1e-05, item = "gene") 
{
  if (is.list(data)) {
    dat <- as.matrix(data$SELEC)
    dat <- rbind(c(rep(1, ncol(dat))), dat)
    groups.vector <- data$groups.vector
    groups.vector <- c(groups.vector[nchar(groups.vector) == 
                                       min(nchar(groups.vector))][1], groups.vector)
    edesign <- data$edesign
    G <- data$g
    family <- data$family
  }
  else {
    G <- nrow(data)
    data <- rbind(c(rep(1, ncol(data))), data)
    dat <- as.matrix(data)
    count.na <- function(x) (length(x) - length(x[is.na(x)]))
    dat <- dat[apply(dat, 1, count.na) >= min.obs, ]
    groups.vector = NULL
    edesign = NULL
  }
  dis <- as.data.frame(design)
  dat <- dat[, as.character(rownames(dis))]
  g <- (dim(dat)[1] - 1)
  n <- dim(dat)[2]
  p <- dim(dis)[2]
  vars.in <- colnames(dis)
  sol <- coefficients <- group.coeffs <- t.score <- sig.profiles <- NULL
  influ.info <- matrix(NA, nrow = nrow(dis), ncol = 1)
  rownames(influ.info) <- rownames(dis)
  if (nvar.correction) 
    alfa <- alfa/ncol(dis)
  for (i in 2:(g + 1)) {
    y <- as.numeric(dat[i, ])
    name <- rownames(dat)[i]
    if (step.method == "backward") {
      reg <- stepback(y = y, d = dis, alfa = alfa, family = family, 
                      epsilon = epsilon)
    }
    else if (step.method == "forward") {
      reg <- stepfor(y = y, d = dis, alfa = alfa, family = family, 
                     epsilon = epsilon)
    }
    else if (step.method == "two.ways.backward") {
      reg <- two.ways.stepback(y = y, d = dis, alfa = alfa, 
                               family = family, epsilon = epsilon)
    }
    else if (step.method == "two.ways.forward") {
      reg <- two.ways.stepfor(y = y, d = dis, alfa = alfa, 
                              family = family, epsilon = epsilon)
    }
    else stop("stepwise method must be one of backward, forward, two.ways.backward, two.ways.forward")
    div <- c(1:round(g/100)) * 100
    if (is.element(i, div)) 
      print(paste(c("fitting ", item, i, "out of", g), 
                  collapse = " "))
    lmf <- glm(y ~ ., data = as.data.frame(dis), family = family, 
               epsilon = epsilon)
    result <- summary(lmf)
    novar <- vars.in[!is.element(vars.in, names(result$coefficients[, 
                                                                    4]))]
    influ <- influence.measures(reg)$is.inf
    influ <- influ[, c(ncol(influ) - 3, ncol(influ) - 1)]
    influ1 <- which(apply(influ, 1, all))
    if (length(influ1) != 0) {
      paste.names <- function(a) {
        paste(names(a)[a], collapse = "/")
      }
      match <- match(rownames(dis), rownames(influ))
      influ <- as.data.frame(apply(influ, 1, paste.names))
      influ.info <- cbind(influ.info, influ[match, ])
      colnames(influ.info)[ncol(influ.info)] <- name
    }
    result <- summary(reg)

    if ((!(result$aic == -Inf) & !is.na(result$aic) & family$family == 
         "gaussian") | family$family != "gaussian") {
      k <- i
      model.glm.0 <- glm(y ~ 1, family = family, epsilon = epsilon)
      if (family$family == "gaussian") {
        test <- anova(model.glm.0, reg, test = "F")
        p.value = test[6][2, 1]
      }
      else {
        test <- anova(model.glm.0, reg, test = "Chisq")
        p.value = test[5][2, 1]
      }
      bondad <- (reg$null.deviance - reg$deviance)/reg$null.deviance
      if (bondad < 0) {
        bondad = 0
      }
      beta.coeff <- result$coefficients[, 1]
      beta.p.valor <- result$coefficients[, 4]
      coeff <- rep(0, (length(vars.in) + 1))
      if (length(novar) != 0) {
        for (m in 1:length(novar)) {
          coeff[position(dis, novar[m]) + 1] <- NA
        }
      }
      p.valor <- t <- as.numeric(rep(NA, (length(vars.in) + 
                                            1)))
      # here I changed alfa to 1
      print(result$coefficients)

      if (result$coefficients[, 4][rownames(result$coefficients) == 
                                   "(Intercept)"] < 1) {
        coeff[1] <- result$coefficients[, 1][rownames(result$coefficients) == 
                                               "(Intercept)"]
        p.valor[1] <- result$coefficients[, 4][rownames(result$coefficients) == 
                                                 "(Intercept)"]
        t[1] <- result$coefficients[, 3][rownames(result$coefficients) == 
                                           "(Intercept)"]
      }
      for (j in 2:length(coeff)) {
        # removing the is.element condition here
        if (is.element(vars.in[j - 1], rownames(result$coefficients))) {
          coeff[j] <- result$coefficients[, 1][rownames(result$coefficients) == 
                                                 vars.in[j - 1]]
          p.valor[j] <- result$coefficients[, 4][rownames(result$coefficients) == 
                                                   vars.in[j - 1]]
          t[j] <- result$coefficients[, 3][rownames(result$coefficients) == 
                                             vars.in[j - 1]]
        }
      }
      if (!all(is.na(p.valor))) {
        sol <- rbind(sol, as.numeric(c(p.value, bondad, 
                                       p.valor)))
        coefficients <- rbind(coefficients, coeff)
        t.score <- rbind(t.score, t)
        sig.profiles <- rbind(sig.profiles, y)
        h <- nrow(sol)
        rownames(sol)[h] <- name
        rownames(coefficients)[h] <- name
        rownames(t.score)[h] <- name
        rownames(sig.profiles)[h] <- name
      }
    }
  }
  if (!is.null(sol)) {
    sol <- as.data.frame(sol)
    coefficients <- as.data.frame(coefficients)
    coeffic <- coefficients
    t.score <- as.data.frame(t.score)
    sig.profiles <- as.data.frame(sig.profiles)
    colnames(sol) <- c("p-value", "R-squared", "p.valor_beta0", 
                       paste("p.valor_", vars.in, sep = ""))
    colnames(coefficients) <- c("beta0", paste("beta", vars.in, 
                                               sep = ""))
    colnames(t.score) <- c("t.score_beta0", paste("t.score_", 
                                                  vars.in, sep = ""))
    colnames(sig.profiles) <- colnames(dat)
    if (!is.null(groups.vector) & !is.null(edesign)) {
      groups <- colnames(edesign)[3:ncol(edesign)]
      degree <- (length(groups.vector)/length(groups)) - 
        1
      for (w in 1:nrow(coefficients)) {
        A <- NULL
        col.names <- NULL
        for (l in 1:length(groups)) {
          B <- reg.coeffs(coefficients = coefficients[w, 
          ], groups.vector = groups.vector, group = groups[l])
          cols <- paste(rep(groups[l], each = length(B)), 
                        paste("beta", c(0:(length(B) - 1)), sep = ""), 
                        sep = "_")
          A <- c(A, B)
          col.names <- c(col.names, cols)
        }
        group.coeffs <- (rbind(group.coeffs, A))
      }
      colnames(group.coeffs) <- col.names
      rownames(group.coeffs) <- rownames(coefficients)
    }
  }
  if (ncol(influ.info) > 2) {
    print(paste("Influence:", ncol(influ.info) - 1, "genes with influential data at slot influ.info. Model validation for these genes is recommended"))
  }
  influ.info <- influ.info[, -1]
  output <- list(sol, sig.profiles, coefficients, as.data.frame(group.coeffs), 
                 t.score, vars.in, G, g, dat, dis, step.method, groups.vector, 
                 edesign, influ.info)
  names(output) <- c("sol", "sig.profiles", "coefficients", 
                     "group.coeffs", "t.score", "variables", "G", "g", "dat", 
                     "dis", "step.method", "groups.vector", "edesign", "influ.info")
  output
}
