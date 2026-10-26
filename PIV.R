# rm(list = ls())
# data <- data.frame(
#   Y = rnorm(100),
#   Z1 = factor(sample(0:1,100,TRUE)),
#   Z2 = factor(sample(0:1,100,TRUE)),
#   Z3 = rnorm(100),
#   Z4 = rnorm(100),
#   X1  = factor(sample(1:3,100,TRUE)),
#   X2  = rnorm(100)
# )
# make_formula(data=data,Y_names='Y',Z_names=c('Z1','Z2','Z3','Z4'), X_names=NULL,m=0)
# make_formula(data=data,Y_names='Y',Z_names=c('Z1','Z2','Z3','Z4'), X_names=NULL,m=1)
# make_formula(data=data,Y_names='Y',Z_names=c('Z1','Z2','Z3','Z4'), X_names=c('L1'),m=0)
# make_formula(data=data,Y_names='Y',Z_names=c('Z1','Z2','Z3','Z4'), X_names=c('L1'),m=1)
# make_formula(data=data,Y_names='Y',Z_names=c('Z3'), X_names=NULL,m=0)
# make_formula(data=data,Y_names='Y',Z_names=c('Z3'), X_names=NULL,m=1)
# make_formula(data=data,Y_names='Y',Z_names=c('Z3'), X_names=c('L1'),m=0)
# make_formula(data=data,Y_names='Y',Z_names=c('Z3'), X_names=c('L1'),m=1)
# summary(gam(make_formula(data=data,Y_names='Y',Z_names=c('Z1','Z2','Z3','Z4'), X_names=NULL,m=3),data=data))


library(dplyr);library(tibble);library(mgcv);library(parallel)
library(grf);library(randomForest)

G_estimation <- function(formula, data, invalid=0,smoothness=5){
  parse_iv_formula <- function(formula){
    
    fml <- deparse(formula)
    fml <- gsub(" ", "", fml)
    fml <- sub("^formula\\((.*)\\)$", "\\1", fml)
    
    # 3. split by ~
    parts <- strsplit(fml, "~")[[1]]
    Y_part <- parts[1]; A_part <- parts[2]; XZ_part <- parts[3]
    
    # 4. split X | Z
    XZ_split <- strsplit(XZ_part, "\\|")[[1]]
    
    X_part <- XZ_split[1]
    Z_part <- XZ_split[2]
    
    # 5. helper：split "+"
    split_plus <- function(x){
      if (x == "") return(character(0))
      unlist(strsplit(x, "\\+"))
    }
    
    list(
      Y_names = split_plus(Y_part),
      A_names = split_plus(A_part),
      X_names = split_plus(X_part),
      Z_names = split_plus(Z_part)
    )
  }
  make_formula <- function(data,Y_names,Z_names,s,X_names=NULL, m=NULL,type=1) {
    auto_smooth_term <- function(vars, data, discrete_threshold = 10){
      is_discrete <- function(x){
        is.factor(x) || is.character(x) || length(unique(x)) <= discrete_threshold
      }
      
      D <- vars[sapply(vars, function(v) is_discrete(data[[v]]))]
      C <- setdiff(vars, D)
      if(length(D)==0 && length(C)==0){
        return('1')
      }
      
      # all the variables are discrete
      if(length(C) == 0){
        return(sprintf("interaction(%s)", paste(D, collapse = ",")))
      }
      # return(sprintf("s(%s)", paste(c(C,D), collapse = ",")))
      
      # all the variables are continuous
      if(length(D) == 0){
        return(sprintf("s(%s, k=%d)", paste(C, collapse = ","),s**length(C)))
      }
      # some continuous, some discrete
      return(c(
        sprintf("s(%s, by = interaction(%s, k=%d))",paste(C, collapse = ","),paste(D, collapse = ","),s**length(C)),
        sprintf("interaction(%s)",paste(D, collapse = ","))
      ))
    }
    
    # Y_names = 'Y'; A_names = 'A'; X_names = NULL; Z_names = c('Z1','Z2','Z3','Z4'); m=0
    if(is.null(m)){m <- length(Z_names)} # fit a saturated model
    if(m>length(Z_names)){stop('too large invalid value')}
    if(m==0 && is.null(X_names)){
      if(type==1){
        return(formula(paste(Y_names,'~ 1')))
      }else{
        return(paste0(Y_names,' ~ 1'))
      }
    }
    terms = matrix(combn(Z_names, m, function(idx){c(idx, X_names)}),nrow = m+length(X_names))
    Terms <- c()
    for(i in 1:ncol(terms)){Terms = c(Terms, auto_smooth_term(terms[,i],data))}
    Terms <- unique(Terms)
    if(type==1){
      res <- as.formula(paste(Y_names,"~", paste(Terms, collapse = " + ")))
      return(res)
    }else{
      return(paste0(Y_names,' ~ ',Terms))
    }
  }
  
  # 1. parse formula
  formula = as.formula(formula)
  fo <- parse_iv_formula(formula)
  
  A_names = fo$A_names
  Y_names = fo$Y_names
  # print(Y_names)
  if(length(fo$X_names)==1 && fo$X_names=='.'){
    X_names=NULL
  } else {
    X_names=fo$X_names
  }
  
  if(length(fo$Z_names)==1 && fo$Z_names=='.'){
    Z_names=NULL
  } else if(length(fo$Z_names)==1 && fo$Z_names=='..'){
    Z_names <- colnames(data)[-which(colnames(data) %in% c(A_names,Y_names,X_names))]
  } else {
    Z_names=fo$Z_names
  }
  
  Y <- data %>% select(all_of(Y_names))
  Z <- data %>% select(all_of(Z_names))
  A <- data %>% select(all_of(A_names))
  formula_A <- make_formula(data,A_names,Z_names=Z_names,X_names = X_names,m = invalid,s=smoothness)
  formula_Y <- make_formula(data,Y_names,Z_names=Z_names,X_names = X_names,m = invalid,s=smoothness)
  formula_Z <- list()
  for (j in 1:length(Z_names)) {
    formula_Z[[j]] <- make_formula(data,Z_names[j],Z_names=NULL,X_names = X_names,s=smoothness)  
  }
  # print(formula_A);print(formula_Y);print(formula_Z);return()
  Z_center <- Z
  A_center <- residuals(gam(formula_A, family = gaussian(), data = data))
  Y_center <- residuals(gam(formula_Y, family = gaussian(), data = data))
  for (j in 1:length(Z_names)) {
    Z_center[,j] <- residuals(gam(formula_Z[[j]], family = gaussian(), data = data))
  }
  restricted_product_matrix <- function(X, k = 0) {
    n <- nrow(X); p <- ncol(X)
    subsets <- list()
    for (size in k:p) {
      subsets <- c(subsets, combn(1:p, size, simplify = FALSE))
    }
    res <- sapply(subsets, function(idx) {
      apply(X[, idx, drop = FALSE], 1, prod)
    })
    if (is.null(dim(res))) {
      res <- matrix(res, ncol = 1)
    }
    colnames(res) <- sapply(subsets, function(idx) {
      paste(colnames(X)[idx], collapse = "_")
    })
    return(res)
  }
  
  Z_mul <- restricted_product_matrix(Z_center,k=invalid+1)
  Z_names_mul <- colnames(Z_mul)
  formula_ivreg <- formula(paste0(
    'Y ~ ',paste('A',collapse = ' + '),
    ' | ',paste(Z_names_mul,collapse = ' + ')
  ))
  # colnames(cbind(Y=Y_center,A=A_center))
  fit.ivreg <- AER::ivreg(formula_ivreg, data = data.frame(cbind(Y=Y_center,A=A_center,Z_mul)))
  # fit.ivreg <- AER::ivreg(formula_ivreg, data = data.frame(cbind(Y=Y,A=A,Z_mul)))
  return(summary(fit.ivreg)$coefficients['A',c('Estimate','Std. Error')])
}



# =====PIV_mgcv=====
PIV_mgcv <- function(formula, data, invalid=0, K=5, print_state=0, M = NULL, smoothness = 5){
  parse_iv_formula <- function(formula){
    
    fml <- deparse(formula)
    fml <- gsub(" ", "", fml)
    fml <- sub("^formula\\((.*)\\)$", "\\1", fml)
    
    # 3. split by ~
    parts <- strsplit(fml, "~")[[1]]
    Y_part <- parts[1]; A_part <- parts[2]; XZ_part <- parts[3]
    
    # 4. split X | Z
    XZ_split <- strsplit(XZ_part, "\\|")[[1]]
    
    X_part <- XZ_split[1]
    Z_part <- XZ_split[2]
    
    # 5. helper：split "+"
    split_plus <- function(x){
      if (x == "") return(character(0))
      unlist(strsplit(x, "\\+"))
    }
    
    list(
      Y_names = split_plus(Y_part),
      A_names = split_plus(A_part),
      X_names = split_plus(X_part),
      Z_names = split_plus(Z_part)
    )
  }
  make_formula <- function(
    data,Y_names,Z_names,s, X_names=NULL, 
    m=NULL,type=1
  ) {
    auto_smooth_term <- function(vars, data, discrete_threshold = 10){
      is_discrete <- function(x){
        is.factor(x) || is.character(x) || length(unique(x)) <= discrete_threshold
      }
      
      D <- vars[sapply(vars, function(v) is_discrete(data[[v]]))]
      C <- setdiff(vars, D)
      if(length(D)==0 && length(C)==0){
        return('1')
      }
      
      # all the variables are discrete
      if(length(C) == 0){
        return(sprintf("interaction(%s)", paste(D, collapse = ",")))
      }
      # return(sprintf("s(%s)", paste(c(C,D), collapse = ",")))
      
      # all the variables are continuous
      if(length(D) == 0){
        return(sprintf("s(%s,k=%d)", paste(C, collapse = ","), s**length(C)))
      }
      # some continuous, some discrete
      return(c(
        sprintf("s(%s, by = interaction(%s),k=%d)",paste(C, collapse = ","),paste(D, collapse = ","),s**length(C)),
        sprintf("interaction(%s)",paste(D, collapse = ","))
      ))
    }
    
    # Y_names = 'Y'; A_names = 'A'; X_names = NULL; Z_names = c('Z1','Z2','Z3','Z4'); m=0
    if(is.null(m)){m <- length(Z_names)} # fit a saturated model
    if(m>length(Z_names)){stop('too large invalid value')}
    if(m==0 && is.null(X_names)){
      if(type==1){
        return(formula(paste(Y_names,'~ 1')))
      }else{
        return(paste0(Y_names,' ~ 1'))
      }
    }
    terms = matrix(combn(Z_names, m, function(idx){c(idx, X_names)}),nrow = m+length(X_names))
    Terms <- c()
    for(i in 1:ncol(terms)){Terms = c(Terms, auto_smooth_term(terms[,i],data))}
    Terms <- unique(Terms)
    if(type==1){
      res <- as.formula(paste(Y_names,"~", paste(Terms, collapse = " + ")))
      return(res)
    }else{
      return(paste0(Y_names,' ~ ',Terms))
    }
  }
  
  # 1. parse formula
  formula = as.formula(formula)
  fo <- parse_iv_formula(formula)
  
  A_names = fo$A_names
  Y_names = fo$Y_names
  
  if(length(fo$X_names)==1 && fo$X_names=='.'){
    X_names=NULL
  } else {
    X_names=fo$X_names
  }
  
  if(length(fo$Z_names)==1 && fo$Z_names=='.'){
    Z_names=NULL
  } else if(length(fo$Z_names)==1 && fo$Z_names=='..'){
    Z_names <- colnames(data)[-which(colnames(data) %in% c(A_names,Y_names,X_names))]
  } else {
    Z_names=fo$Z_names
  }
  
  if(is.null(M)){M<-length(Z_names)}
  if(M>length(Z_names)){M<-length(Z_names)}
  # 2. formula construction
  fom_A_Z    <- make_formula(data,Y_names=A_names,Z_names=Z_names,X_names=X_names,m=M,s=smoothness)
  fom_Y_Z    <- make_formula(data,Y_names=Y_names,Z_names=Z_names,X_names=X_names,m=M,s=smoothness)
  fom_A <- list(); fom_Y <- list(); fom_pseudo <- list()
  for (j in 1:length(invalid)) {
    fom_A[[j]] <- make_formula(data,Y_names=A_names,Z_names=Z_names,X_names=X_names,m=invalid[j],s=smoothness)
    fom_Y[[j]] <- make_formula(data,Y_names=Y_names,Z_names=Z_names,X_names=X_names,m=invalid[j],s=smoothness)
    fom_pseudo[[j]]<- make_formula(
      data,Y_names="pseudo",Z_names=Z_names,X_names=X_names,
      m=0,s = smoothness*2
    )
  }
  
  
  
  
  if(print_state){
    lapply(list(fom_A_Z,fom_A,fom_Y_Z,fom_Y,fom_pseudo),
           function(f) cat(deparse(f), "\n"))
    if(print_state==2) return()
  }
  
  
  # 3. data preprocessing
  n <- nrow(data); mydat <- data
  # if(type == 'treat'){
  #   mydat <- mydat %>% mutate(Y = Y * A)
  # } else if(type == 'control'){
  #   mydat <- mydat %>% mutate(Y = Y * (A - 1))
  # }
  
  # 4. split
  cross_fit_split_index <- function(n, K){
    folds <- sample(rep(1:K, length.out=n))
    lapply(1:K, function(k){
      list(
        train_idx = which(folds != k),
        test_idx  = which(folds == k)
      )
    })
  }
  if(K>1){
    split_indices <- cross_fit_split_index(n, K)  
  }else{
    split_indices <- list(list(train_idx = seq(n),test_idx  = seq(n)))
  }
  
  
  # 5. parallel cluster
  cl <- makeCluster(min(K, detectCores() - 1))
  
  clusterExport(cl,
                varlist = c("mydat","split_indices", "invalid", "Y_names","A_names",
                            "fom_A_Z","fom_A","fom_Y_Z","fom_Y","fom_pseudo"),
                envir = environment()
  )
  clusterEvalQ(cl, library(mgcv))
  clusterEvalQ(cl, library(dplyr)); k <- 1
  
  # 6. parallel K-fold
  results <- parLapply(cl, 1:K, function(k){
    
    train_index <- split_indices[[k]]$train_idx
    test_index  <- split_indices[[k]]$test_idx
    train <- mydat[train_index,]; test <- mydat[test_index,]
    
    
    pi_model   <- gam(fom_A_Z, data=train, family=gaussian())
    xi_model   <- gam(fom_Y_Z, data=train, family=gaussian())
    pi_test <- predict(pi_model, test, type="response")
    xi_test <- predict(xi_model, test, type="response")
    pi_train <- predict(pi_model, train, type="response")
    xi_train <- predict(xi_model, train, type="response")
    
    delta_model <- list(); eta_model <- list()
    delta_test <- list(); eta_test <- list()
    kappa_model <- list(); mu_model <- list()
    kappa_test <- list(); mu_test <- list()
    delta_train <- list(); eta_train <- list()
    est_aug <- list(); est_ipw <-list()
    est_aug_Y <- list(); est_aug_A <- list(); est_tau <- list(); j<-1
    for (j in 1:length(invalid)) {
      
      delta_model[[j]]<- gam(fom_A[[j]], data=train, family=gaussian())
      eta_model[[j]]  <- gam(fom_Y[[j]], data=train, family=gaussian())
      delta_test[[j]] <- predict(delta_model[[j]], test, type="response")
      eta_test[[j]]   <- predict(eta_model[[j]], test, type="response")
      delta_train[[j]]<- delta_model[[j]]$fitted.values
      eta_train[[j]]  <- eta_model[[j]]$fitted.values
      
      Y_test <- test[Y_names][[1]]; A_test <- test[A_names][[1]]; Y_train <- train[Y_names][[1]]; A_train <- train[A_names][[1]]
      pseudo_kappa <- (A_train - delta_train[[j]])^2 - (A_train - pi_train)^2
      kappa_model[[j]]<- gam(fom_pseudo[[j]],data = train %>% mutate(pseudo = pseudo_kappa),family = gaussian(link="log"))
      kappa_test[[j]] <- predict(kappa_model[[j]], test, type="response")
      # kappa_test[[j]] <- pmax(kappa_test[[j]], quantile(kappa_test[[j]], 0.025))
      
      pseudo_mu <- ((A_train-pi_train)*xi_train - (A_train-delta_train[[j]])*eta_train[[j]] + (pi_train-delta_train[[j]])*Y_train)
      mu_model[[j]]  <- gam(fom_pseudo[[j]],data = train %>% mutate(pseudo = pseudo_mu), family = gaussian())
      mu_test[[j]] <- predict(mu_model[[j]], test, type="response")/kappa_test[[j]]
      
      est_aug[[j]] <- (pi_test - delta_test[[j]])/kappa_test[[j]] * Y_test +
        mu_test[[j]]/kappa_test[[j]] * (
          kappa_test[[j]] + (A_test - pi_test)^2 - (A_test - delta_test[[j]])^2
        ) +
        1/kappa_test[[j]] * (
          xi_test * (A_test - pi_test) -
            eta_test[[j]] * (A_test - delta_test[[j]])
        )
      est_ipw[[j]] <- Y_test * (pi_test - delta_test[[j]])/kappa_test[[j]]
      
      est_aug_Y[[j]] <- (pi_test-delta_test[[j]])*(Y_test-eta_test[[j]]) + 
        (A_test-pi_test)*(xi_test-eta_test[[j]])
      est_aug_A[[j]] <- (pi_test-delta_test[[j]])*(A_test-delta_test[[j]]) + 
        (A_test-pi_test)*(pi_test-delta_test[[j]])
      est_tau[[j]] <- (A_test-delta_test[[j]])^2 - (A_test-pi_test)^2
      # est_tau[[j]] <- -est_aug_A[[j]] 
    }
    list(
      test_index = test_index, kappa_test = kappa_test,
      est_aug = est_aug, est_or = mu_test, est_ipw = est_ipw, 
      est_aug_Y = est_aug_Y, est_aug_A = est_aug_A, est_tau = est_tau
    )
  })
  stopCluster(cl)
  # 7. aggregate
  EST_AUG <- matrix(0,n,length(invalid))
  EST_OR  <- matrix(0,n,length(invalid))
  EST_IPW <- matrix(0,n,length(invalid))
  kappa_test <- matrix(0,n,length(invalid))
  
  EST_AUG_Y <- matrix(0,n,length(invalid))
  EST_AUG_A <- matrix(0,n,length(invalid))
  EST_TAU <- matrix(0,n,length(invalid)); res <- results[[1]]
  
  for(res in results){
    idx <- res$test_index; j <- 1
    for(j in 1:length(invalid)){
      EST_AUG[idx,j] <- res$est_aug[[j]]
      EST_OR[idx,j]  <- res$est_or[[j]]
      EST_IPW[idx,j] <- res$est_ipw[[j]]
      kappa_test[idx,j] <- res$kappa_test[[j]]
      
      EST_AUG_Y[idx,j] <- res$est_aug_Y[[j]]
      EST_AUG_A[idx,j] <- res$est_aug_A[[j]]
      EST_TAU[idx,j] <- res$est_tau[[j]]
    }
  }
  if(print_state){
    par(mfrow=c(2,2));for(j in 1:length(invalid)){plot(data$X,EST_OR[,j],'p',cex=0.4,col='gray')}
    seq(-0.5,0.5,0.01) %>% lines(., 4 * sin(.*6) + 2,col='red');par(mfrow=c(1,1))
  }
  
  # 8. output
  est_aug2 <- apply(EST_AUG_Y, 2, mean)/apply(EST_AUG_A, 2, mean)
  est_YA <- apply(EST_AUG_Y-sweep(EST_AUG_A, 2, est_aug2, "*"), 2, sd)/abs(apply(EST_TAU, 2, mean))/sqrt(n)
  res <- list(
    num_of_crossing_folds=K, invalid=invalid,
    res=list(
      est_aug = apply(EST_AUG, 2, mean),
      # est_or  = apply(EST_OR, 2, mean),
      # est_ipw = apply(EST_IPW, 2, mean),
      sd_aug  = apply(EST_AUG, 2, sd)/sqrt(n),
      # sd_or   = apply(EST_OR, 2, sd)/sqrt(n),
      # sd_ipw  = apply(EST_IPW, 2, sd)/sqrt(n),
      est_aug2 = est_aug2,
      sd_aug2  = est_YA
    ), 
    formula=list(
      fom=formula,
      fom_A_Z=fom_A_Z,
      fom_A=fom_A,
      fom_Y_Z=fom_Y_Z,
      fom_Y=fom_Y
    )
  )
  class(res) <- 'PIV'
  return(res)
}





summary.PIV <- function(object, ...) {
  cat("\n================================================================\n")
  cat("Fortified Instrumental Variable (PIV) Estimator Summary\n")
  cat("================================================================\n\n")
  
  
  cat("Formula: ", deparse(object$formula$fom),"\n")
  # cross-fitting info
  cat("Cross-fitting folds: ", object$num_of_crossing_folds,"\n")
  
  # invalid IV info
  cat("Projected dimension: ", object$invalid,"\n")
  
  
  
  cat("----------------------------------------------------------------\n")
  cat("Estimation results:\n")
  cat("----------------------------------------------------------------\n")
  
  res <- object$res
  
  cat(sprintf("AIPW: %.6f (SE = %.6f)\n",
              res$est_aug, res$sd_aug))
  
  cat(sprintf("OR:   %.6f (SE = %.6f)\n",
              res$est_or, res$sd_or))
  
  cat(sprintf("IPW:  %.6f (SE = %.6f)\n",
              res$est_ipw, res$sd_ipw))
  
  cat("\n----------------------------------------------------------------\n")
  cat("Formulas used:\n")
  cat("----------------------------------------------------------------\n")
  
  invisible(lapply(
    names(object$formula), 
    function(nm){if(nm!='fom'){cat(nm,": ",deparse(object$formula[[nm]]),"\n")}}
  ))
  cat("\n================================================================\n")
  
  invisible(object)
}

print.PIV <- function(object, type=1,...){
  if(type==1){
    print(cbind(data.frame(invalid = object$invalid),data.frame(object$res)))
  }else{
    print(unlist(object$res))  
  }
}




# The default ATE = 2 reproduces Table "concrete DGPs" exactly.
sim_PIV <- function(n, case, P = 3, ATE = 2){
  if(length(case) != 1L || is.na(case) || !case %in% 1:7){
    stop("case must be one of 1, 2, 3, 4, 5, 6, or 7")
  }
  if(length(P) != 1L || is.na(P) || P != 3){
    stop("The seven specified DGPs are defined for exactly P = 3 instruments")
  }
  if(length(n) != 1L || is.na(n) || n < 1 || n != as.integer(n)){
    stop("n must be a positive integer")
  }
  
  # Unmeasured and measured confounders.
  U <- runif(n, min = -1, max = 1)
  X <- runif(n, min = -0.5, max = 0.5)
  
  # Independent disturbances used to generate Z and A.
  epsilon_z  <- matrix(runif(n * P, min = -1, max = 1),
                       nrow = n, ncol = P)
  
  # phi_1(s) = phi_2(s) = 2 I(s >= 0) - 1 and phi_3(s) = 0.6 s.
  phi_binary <- function(s) 2 * as.numeric(s >= 0) - 1
  Z_mat <- matrix(NA_real_, nrow = n, ncol = P)
  if(case %in% c(1, 2)){
    Z_mat[, 1] <- phi_binary(X + 6 * U)
  }else{
    Z_mat[, 1] <- phi_binary(X + 6 * epsilon_z[, 1])
  }
  Z_mat[, 2] <- phi_binary(X + 6 * epsilon_z[, 2])
  
  if(case == 7){
    Z_mat[, 3] <- phi_binary(X + 6 * epsilon_z[, 3])
  }else{
    Z_mat[, 3] <- 0.25 * (X + 6 * epsilon_z[, 3])
  }
  colnames(Z_mat) <- paste0("Z", seq_len(P))
  
  # Treatment mechanism.
  if(case %in% c(3, 4)){
    H <- rbinom(n, size = 1, prob = 0.7 + 0.3 * plogis(-2*X))
    V <- rowSums(Z_mat) + runif(n, min = -1, max = 1)
    A <- H * as.numeric(V >= 0) +
      (1 - H) * as.numeric(2 * U + Z_mat[, 1] >= 0)
  }else{
    A <- as.numeric((5 * plogis(5 * X) + 1) * rowSums(Z_mat) + U >= 0)
  }
  
  # Baseline potential outcome Y(0).
  if(case == 1){
    Y0 <- 0.05 * rowSums(Z_mat) * (0.5 + X^2) + U
  }else if(case == 2){
    Y0 <- 0.05 * rowSums(Z_mat) * (0.5 + X^2) + U
  }else if(case == 3){
    Y0 <- 0.05 * rowSums(Z_mat) * (0.5 + X^2) + U
  }else if(case == 4){
    # Y0 <- 0.15 * Z_mat[, 1] * (0.5 + X^2) + U
    Y0 <- 0.1 * Z_mat[, 1] * (0.5 + X^2) + U
  }else if(case %in% c(5,7)){
    Y0 <- 0.05 * rowSums(Z_mat) + U + X
  }else{
    Y0 <- 0.1 * Z_mat[, 1] + U + X
  }
  
  # Individual treatment effect Y(1) - Y(0).
  if(case == 1){
    Y10 <- - 4 * sin(6 * X) + ATE
  }else if(case == 2){
    Y10 <- - 4 * X + ATE
  }else if(case == 3){
    Y10 <- 2 * U - 4 * sin(6 * X) + ATE
  }else if(case == 4){
    # Y10 <- 2 * U - 4 * X + ATE
    Y10 <- 2 * U - 2 * X + ATE
  }else{
    Y10 <- rep(ATE, n)
  }
  
  tibble(
    X = X,
    Z1 = Z_mat[, 1],
    Z2 = Z_mat[, 2],
    Z3 = Z_mat[, 3],
    A = A, Y = Y0 + A * Y10
  )
}

