# An example
# coefficient generating functions used in the simulation
Coef.gen<- function(s, h,q=30, size.A0, M, sig.beta,sig.delta1, sig.delta2, p, exact=T){
  beta0 <- c(rep(sig.beta,s), rep(0, p - s))
  W <- rep.col(beta0,  M) # ten prior estimates
  W[1,]<-W[1,]-2*sig.beta
  for(k in 1:M){
    if(k <= size.A0){
      if(exact){
        samp0<- sample(1:p, h, replace=F)
        W[samp0,k] <-W[samp0,k] + rep(-sig.delta1, h)
      }else{
        W[1:100,k] <-W[1:100,k] + rnorm(100, 0, h/100)
      }
    }else{
      if(exact){
        samp1 <- sample(1:p, q, replace = F)
        W[samp1,k] <- W[samp1,k] + rep(-sig.delta2,q)
      }else{
        W[1:100,k] <-W[1:100,k] + rnorm(100, 0, q/100)
      }
    }
  }

  return(list(W=W, beta0=beta0))
}

# Prepare test data for translasso
prep.data <- function(p, M, n0, size.A0, l1) {
  s = 16
  sig.beta = 0.3

  n.vec <- c(n0, rep(100, M))
  Sig.X <- diag(1, p)

  h=6
  A0 = 1:size.A0
  beta0<-
    coef.all <-Coef.gen( s, h = h, q = 2*s, size.A0 = size.A0,  M = M,   sig.beta = sig.beta,
                         sig.delta1 = sig.beta, sig.delta2 = sig.beta+0.2, p = p, exact=F)
  B <- cbind(coef.all$beta0, coef.all$W)
  beta0 <- coef.all$beta0

  ###generate the data###
  X <- NULL
  y <- NULL
  for (k in 1:(M + 1)) {
    X <- rbind(X, SimDesign::rmvnorm(n.vec[k], rep(0, p), Sig.X))
    ind.k <- ind.set(n.vec, k)
    y <- c(y, X[ind.k, ] %*% B[, k] + rnorm (n.vec[k], 0, 1))
  }

  return(list(X=X, y=y, beta0=beta0, p=p, M=M, n0=n0, size.A0=size.A0, l1=l1, n.vec=n.vec))
}

test_that("oracle translasso alogrithm works", {
  set.seed(123)
  data <- prep.data(500, 20, 150, 12, T)
  list2env(data, .GlobalEnv)

  otl <- las.kA(X, y, A0 = 1:size.A0, n.vec = n.vec, l1=l1)
  mse.val <- mse.fun(as.numeric(otl$beta.kA), beta0)$est.err

  expect_lt(mse.val, 0.15)

  # A method for comparison: it has the same pipeline of the Trans-Lasso
  # but with sparsity index R_k=\|w^{(k)}-\beta\|_1 and a naive aggregation (empirical risk minimization)
  prop.sp.re1 <- Trans.lasso.sp(X, y, n.vec, I.til = 1:50, l1 = l1)
  prop.sp.re2 <- Trans.lasso.sp(X, y, n.vec, I.til = 101:n.vec[1], l1=l1)
  if(size.A0 > 0 & size.A0< M){
    Rank.re.sp <- (sum(prop.sp.re1$rank.pi[1:size.A0]<=size.A0) +
                     sum(prop.sp.re2$rank.pi[1:size.A0]<=size.A0))/2/size.A0
  }else{ Rank.re.sp <-1 }
  beta.sp <- (prop.sp.re1$beta.sp + prop.sp.re2$beta.sp) / 2
  mse.val <- mse.fun(beta.sp, beta0)$est.err

  expect_lt(mse.val, 0.5)

  #print("OTL")
  #print(otl)
  #print("SP")
  #print(prop.sp.re1)
  expect_lt(max(abs(otl$beta.kA - prop.sp.re1$beta.sp)), 0.2)
  expect_lt(max(abs(otl$beta.kA - prop.sp.re2$beta.sp)), 0.2)

  rm(X, y, beta0, p, M, n0, size.A0, l1, n.vec, envir = .GlobalEnv)
})

test_that("translasso algorithm works", {
  set.seed(123)
  data <- prep.data(500, 20, 150, 12, T)
  list2env(data, .GlobalEnv)

  prop.re1 <- Trans.lasso(X, y, n.vec, I.til = 1:50, l1 = l1)
  prop.re2 <- Trans.lasso(X, y, n.vec, I.til = 101:n.vec[1], l1=l1)
  if(size.A0 > 0 & size.A0< M){ # Rank.re characterizes the performance of the sparsity index Rk
    Rank.re<- (sum(prop.re1$rank.pi[1:size.A0]<=size.A0) +
                 sum(prop.re2$rank.pi[1:size.A0]<=size.A0))/2/size.A0
  }else{ Rank.re <- 1 }
  beta.prop <- (prop.re1$beta.hat + prop.re2$beta.hat) / 2
  mse.val = mse.fun(beta.prop, beta0)$est.err

  expect_lt(mse.val, 0.5)

  # A method for comparison: it is the same as Trans-Lasso except
  # that the bias correction step (step 2 of Oracle Trans-Lasso) is omitted
  beta.pool<-(prop.re1$beta.pool+prop.re2$beta.pool)/2
  mse.val <- mse.fun(beta.pool, beta0)$est.err

  expect_lt(mse.val, 0.5)

  # Naive TransLasso: simply assumes A0=1:K
  otl <- las.kA(X, y, A0 = 1:M, n.vec = n.vec, l1=l1) # naive TransLasso
  mse.val <- mse.fun(as.numeric(otl$beta.kA), beta0)$est.err

  expect_lt(mse.val, 0.5)



  rm(X, y, beta0, p, M, n0, size.A0, l1, n.vec, envir = .GlobalEnv)
})

