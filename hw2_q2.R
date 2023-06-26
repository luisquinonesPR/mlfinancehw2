################
### HW2 - Q2 ###
################

setwd(".")

data = read.csv("datframe_big.csv")

library(xts)
library(quantmod)
library(CVXR)
library(mvtnorm)
library(PerformanceAnalytics)

T = nrow(data)

X = data[, 7:20]
F_FamaFrench = data[, 3:5]
SentIndx = data[, 21:34]
dates = data[,1]
dates = as.matrix(dates)
origin_date <- as.Date("2015-01-05")

T_trn = round(0.5*T)
X_trn = X[1:T_trn, ]
X_trn = as.matrix(X_trn) 
X_tst = X[(T_trn+1):T, ]
X_tst = as.matrix(X_tst) 
F_FamaFrench_trn = F_FamaFrench[1:T_trn, ]
F_FamaFrench_trn = as.matrix(F_FamaFrench_trn)
F_FamaFrench_tst = F_FamaFrench[(T_trn+1):T, ]
F_FamaFrench_tst = as.matrix(F_FamaFrench_tst)
SentIndx_trn = SentIndx[1:T_trn,]
SentIndx_trn = as.matrix(SentIndx_trn)
SentIndx_tst <-SentIndx[(T_trn+1):T,]
SentIndx_tst = as.matrix(SentIndx_tst)



# Fama-French 3-factor model
F_ <- cbind(ones = 1, F_FamaFrench_trn)
Gamma <- t(solve(t(F_) %*% F_, t(F_) %*% X_trn))
colnames(Gamma) <- c("alpha", "beta1", "beta2", "beta3")
alpha <- Gamma[, 1]
B <- Gamma[, 2:4]
#E <- xts(t(t(X_trn) - Gamma %*% t(F_)), order.by = as.POSIXct(dates,origin = origin_date))
E <-t(t(X_trn) - Gamma %*% t(F_))
PsiFF <- (1/(T_trn-4)) * t(E) %*% E
Sigma_FamaFrench <- B %*% cov(F_FamaFrench_trn) %*% t(B) + diag(diag(PsiFF))

mu_hat = colMeans(X_trn)

portfolioMaxReturnRobustEllipsoid <- function(mu_hat, S, kappa = 0.1) {
  S12 <- chol(S)  # t(S12) %*% S12 = Sigma
  w <- Variable(length(mu_hat))
  prob <- Problem(Maximize( t(w) %*% mu_hat - kappa*p_norm(S12 %*% w,p=2) ), 
                  constraints = list(w >= 0, sum(w) == 1))
  result <- solve(prob)
  return(as.vector(result$getValue(w)))
}


check = portfolioMaxReturnRobustEllipsoid(mu_hat, Sigma_FamaFrench)

print(check)

df <- as.data.frame(check)

library(grid)
library(gridExtra)

filename <- "output.png"

png(filename)
grid.arrange(tableGrob(df))
dev.off()
