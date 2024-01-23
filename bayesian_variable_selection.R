library(MASS)
library(faraway)
library(glmnet)
library(sparsevb)
library(dplyr)
data(UScrime) # respon : y
data(diabetes, package="faraway") # respon : glyhb
head(UScrime)
head(diabetes)
dim(diabetes)

####### 변수선택 간단한 예제(전진선택법) ####
dat <- scale(UScrime)
MSE_LM <- c()
MSE_FW <- c()
for(i in 1:100){
  res.ind <- 16
  x <- dat
  train <- sample(1:nrow(x), 2*nrow(x)/3)
  dat_train <- dat[train,]
  dat_test <- dat[(-train),]
  dat_train_x <- dat_train[,-res.ind]
  dat_train_y <- dat_train[,res.ind]
  dat_test_x <- dat_test[,-res.ind]
  dat_test_y <- dat_test[,res.ind]
  null.model <- lm(y ~ 1, data=as.data.frame(dat_train))
  full.model <- lm(y ~ ., data=as.data.frame(dat_train))
  FW.model <- step(direction="forward", object=null.model,
                   scope=list(lower=null.model, upper=full.model))
  pred_y_LM <- predict(full.model, as.data.frame(dat_test_x))
  pred_y_FW <- predict(FW.model, as.data.frame(dat_test_x))
  y <- dat_test_y
  n <- nrow(dat_test)
  MSE_LM[i] <- mean((pred_y_LM-y)^2)
  MSE_FW[i] <- mean((pred_y_FW-y)^2)
}
mean(MSE_LM)
mean(MSE_FW)

###### simulation ####
N <- 100
p <- 15
name <- 'X1'
for(i in 2:p){
  name <- c(name,paste0('X',i))
}

#### mean((\beta - \hat{\beta})^2) ####
MSE_LM_beta <- c()
MSE_FW_beta <- c()
MSE_BE_beta <- c()
MSE_SW_beta <- c()
MSE_lasso_beta <- c()
MSE_VB_beta <- c()

#### 변수선택 오류율 ####
MSE_FW_01 <- c()
MSE_BE_01 <- c()
MSE_SW_01 <- c()
MSE_lasso_01 <- c()
MSE_VB_01 <- c()

for(i in 1:100){
  X <- matrix(rnorm(N*p), nrow = N, ncol = p)
  colnames(X) <- name
  true_betas <- c(c(0.1, 0.2, 0.3, 0.4, 0.5),
                  rep(0, 10))
  
  y <- rnorm(N, X %*% true_betas, sd = 1)
  dat <- data.frame(X,y)
  
  null.model <- lm(y ~ 1, data=dat)
  full.model <- lm(y ~ ., data=dat)
  FW.model <- step(direction="forward", object=null.model, scope=list(lower=null.model, upper=full.model))
  FW.model$coefficients
  BE.model <- step(direction="backward", object=full.model, scope=list(lower=null.model, upper=full.model))
  BE.model$coefficients
  SW.model <- step(direction="both", object=null.model, scope=list(lower=null.model, upper=full.model))
  SW.model$coefficients
  
  true_betas2 <- t(data.frame(true_betas,row.names = name))
  true_betas2 <- cbind(colnames(true_betas2),t(true_betas2))
  colnames(true_betas2) <- c('num','value')
  
  LM <- cbind(names(full.model$coefficients[-1]),full.model$coefficients[-1])
  colnames(LM) <- c('num','value')
  LM <- merge(true_betas2,LM,by='num',all.x=T) %>% arrange(num)
  LM[is.na(LM)] <- 0
  
  FW <- cbind(names(FW.model$coefficients[-1]),FW.model$coefficients[-1])
  colnames(FW) <- c('num','value')
  FW <- merge(true_betas2,FW,by='num',all.x=T) %>% arrange(num)
  FW[is.na(FW)] <- 0
  
  BE <- cbind(names(BE.model$coefficients[-1]),BE.model$coefficients[-1])
  colnames(BE) <- c('num','value')
  BE <- merge(true_betas2,BE,by='num',all.x=T) %>% arrange(num)
  BE[is.na(BE)] <- 0
  
  SW <- cbind(names(SW.model$coefficients[-1]),SW.model$coefficients[-1])
  colnames(SW) <- c('num','value')
  SW <- merge(true_betas2,SW,by='num',all.x=T) %>% arrange(num)
  SW[is.na(SW)] <- 0
  
  pred_y_LM <- predict(full.model, as.data.frame(X))
  pred_y_FW <- predict(FW.model, as.data.frame(X))
  pred_y_BE <- predict(BE.model, as.data.frame(X))
  pred_y_SW <- predict(SW.model, as.data.frame(X))
  
  MSE_LM_beta[i] <- mean((as.numeric(LM[,2])-as.numeric(LM[,3]))^2)
  MSE_FW_beta[i] <- mean((as.numeric(FW[,2])-as.numeric(FW[,3]))^2)
  MSE_BE_beta[i] <- mean((as.numeric(BE[,2])-as.numeric(BE[,3]))^2)
  MSE_SW_beta[i] <- mean((as.numeric(SW[,2])-as.numeric(SW[,3]))^2)
  
  MSE_FW_01[i] <- mean(((as.numeric(FW[,2]) == 0)-(as.numeric(FW[,3]) == 0))^2)
  MSE_BE_01[i] <- mean(((as.numeric(BE[,2]) == 0)-(as.numeric(BE[,3]) == 0))^2)
  MSE_SW_01[i] <- mean(((as.numeric(SW[,2]) == 0)-(as.numeric(SW[,3]) == 0))^2)
  
  grid <- 2^seq(10, -2, length=100)
  
  cv.out <- cv.glmnet(X, y, alpha=1)
  bestlam <- cv.out$lambda.min
  out <- glmnet(X, y, alpha=1, lambda=grid)
  lasso.coef <- predict(out, type="coefficients", s=bestlam)[2:(ncol(X)+1),]
  
  MSE_lasso_beta[i] <- mean((true_betas - lasso.coef)^2)
  MSE_lasso_01[i] <- mean(((true_betas == 0)-(lasso.coef == 0))^2)
  
  svbfit <- svb.fit(
    X,
    y,
    family = "linear",
    slab = "laplace",
    mu=rep(0,ncol(X)),
    sigma = rep(1, ncol(X)),
    gamma=rep(0.5,ncol(X)),
    alpha=2,
    beta=2,
    prior_scale = 1,
    intercept = FALSE,
    max_iter = 1000,
    tol = 1e-05
  )
  MSE_VB_beta[i] <- mean((true_betas - svbfit$mu*(svbfit$gamma > 0.5))^2)
  MSE_VB_01[i] <- mean(((true_betas == 0)-(svbfit$gamma < 0.5))^2)
}
MSE_LM_beta %>% mean
MSE_FW_beta %>% mean
MSE_BE_beta %>% mean
MSE_SW_beta %>% mean
MSE_lasso_beta %>% mean
MSE_VB_beta %>% mean

MSE_FW_01 %>% mean
MSE_BE_01 %>% mean
MSE_SW_01 %>% mean
MSE_lasso_01 %>% mean
MSE_VB_01 %>% mean

####### UScrime data ####
iter=100
dat <- scale(UScrime)
for(i in 1:iter){
  res.ind <- 16
  x <- dat
  y <- dat_test_y
  n <- nrow(dat_test)
  train <- sample(1:nrow(x), 2*nrow(x)/3)
  test <- (-train)
  
  # FW, BE, SW
  dat_train <- dat[train,]
  dat_test <- dat[test,]
  dat_train_x <- dat_train[,-res.ind]
  dat_train_y <- dat_train[,res.ind]
  dat_test_x <- dat_test[,-res.ind]
  dat_test_y <- dat_test[,res.ind]
  null.model <- lm(y ~ 1, data=as.data.frame(dat_train))
  full.model <- lm(y ~ ., data=as.data.frame(dat_train))
  FW.model <- step(direction="forward", object=null.model, scope=list(lower=null.model, upper=full.model))
  BE.model <- step(direction="backward", object=full.model, scope=list(lower=null.model, upper=full.model))
  SW.model <- step(direction="both", object=null.model, scope=list(lower=null.model, upper=full.model))
  pred_y_LM <- predict(full.model, as.data.frame(dat_test_x))
  pred_y_FW <- predict(FW.model, as.data.frame(dat_test_x))
  pred_y_BE <- predict(BE.model, as.data.frame(dat_test_x))
  pred_y_SW <- predict(SW.model, as.data.frame(dat_test_x))
  MSE_LM[i] <- mean((pred_y_LM-y)^2)
  MSE_FW[i] <- mean((pred_y_FW-y)^2)
  MSE_BE[i] <- mean((pred_y_BE-y)^2)
  MSE_SW[i] <- mean((pred_y_SW-y)^2)
  
  # LASSO
  x <- model.matrix(y ~.,as.data.frame(dat))[,-1]
  y <- dat[,res.ind]
  grid <- 2^seq(10, -2, length=100)
  y.test <- y[test]
  cv.out <- cv.glmnet(x[train ,], y[train], alpha=1)
  bestlam <- cv.out$lambda.min
  lasso.pred <- predict(lasso.mod, s=bestlam, newx=x[test,])
  out <- glmnet(x, y, alpha=1, lambda=grid)
  lasso.coef <- predict(out, type="coefficients", s=bestlam)[1:(ncol(x)+1),]
  MSE_lasso[i] <- mean((lasso.pred-y.test)^2)
  
  # VB
  svbfit <- svb.fit(
    dat_train_x,
    dat_train_y,
    family = "linear",
    slab = "laplace",
    mu=rep(0,ncol(dat_train_x)),
    sigma = rep(1, ncol(dat_train_x)),
    gamma=rep(0.5,ncol(dat_train_x)),
    alpha=2,
    beta=2,
    prior_scale = 1,
    intercept = FALSE,
    max_iter = 1000,
    tol = 1e-05
  )
  beta_vb <- (svbfit$gamma > 0.5) * svbfit$mu
  MSE_VB[i] <- mean((dat_test_x %*% beta_vb - dat_test_y)^2)
}
MSE_LM %>% mean
MSE_FW %>% mean
MSE_BE %>% mean
MSE_SW %>% mean
MSE_lasso %>% mean
MSE_VB %>% mean

###### diabetes data ######
dat <- na.omit(diabetes[,-1])
for(i in 1:100){
  res.ind <- 5
  x <- dat
  y <- dat_test_y
  n <- nrow(dat_test)
  train <- sample(1:nrow(x), 2*nrow(x)/3)
  test <- (-train)
  # FW, BE, SW
  dat_train <- dat[train,]
  dat_test <- dat[test,]
  dat_train_x <- dat_train[,-res.ind]
  dat_train_y <- dat_train[,res.ind]
  dat_test_x <- dat_test[,-res.ind]
  dat_test_y <- dat_test[,res.ind]
  null.model <- lm(glyhb ~ 1, data=as.data.frame(dat_train))
  full.model <- lm(glyhb ~ ., data=as.data.frame(dat_train))
  FW.model <- step(direction="forward", object=null.model, scope=list(lower=null.model, upper=full.model))
  BE.model <- step(direction="backward", object=full.model, scope=list(lower=null.model, upper=full.model))
  SW.model <- step(direction="both", object=null.model, scope=list(lower=null.model, upper=full.model))
  pred_y_LM <- predict(full.model, as.data.frame(dat_test_x))
  pred_y_FW <- predict(FW.model, as.data.frame(dat_test_x))
  pred_y_BE <- predict(BE.model, as.data.frame(dat_test_x))
  pred_y_SW <- predict(SW.model, as.data.frame(dat_test_x))
  
  MSE_LM[i] <- mean((pred_y_LM-y)^2)
  MSE_FW[i] <- mean((pred_y_FW-y)^2)
  MSE_BE[i] <- mean((pred_y_BE-y)^2)
  MSE_SW[i] <- mean((pred_y_SW-y)^2)
  
  # LASSO
  x <- model.matrix(glyhb ~.,as.data.frame(dat))[,-1]
  y <- dat[,res.ind]
  grid <- 2^seq(10, -2, length=100)
  y.test <- y[test]
  cv.out <- cv.glmnet(x[train ,], y[train], alpha=1)
  bestlam <- cv.out$lambda.min
  lasso.pred <- predict(lasso.mod, s=bestlam, newx=x[test,])
  out <- glmnet(x, y, alpha=1, lambda=grid)
  lasso.coef <- predict(out, type="coefficients", s=bestlam)[1:(ncol(x)+1),]
  MSE_lasso[i] <- mean((lasso.pred-y.test)^2)
  
  # VB
  svbfit <- svb.fit(
    matrix(unlist(dat_train_x),nrow=nrow(dat_train_x)),
    dat_train_y,
    family = "linear",
    slab = "laplace",
    mu=rep(0,ncol(dat_train_x)),
    sigma = rep(1, ncol(dat_train_x)),
    gamma=rep(0.5,ncol(dat_train_x)),
    alpha=2,
    beta=2,
    prior_scale = 1,
    intercept = FALSE,
    max_iter = 1000,
    tol = 1e-05
  )
  beta_vb <- (svbfit$gamma > 0.5) * svbfit$mu
  MSE_VB[i] <- mean((matrix(unlist(dat_test_x),nrow=nrow(dat_test_x)) %*% beta_vb - dat_test_y)^2)
}
MSE_LM %>% mean
MSE_FW %>% mean
MSE_BE %>% mean
MSE_SW %>% mean
MSE_lasso %>% mean
MSE_VB %>% mean
