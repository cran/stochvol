### R code from vignette source 'article.Rnw'

###################################################
### code chunk number 1: article.Rnw:1-4
###################################################
# This should be set to FALSE in case you want to redo the sampling as well.
# Default is TRUE to save CRAN some time.
usePreCalcResults <- TRUE


###################################################
### code chunk number 2: usd1
###################################################
set.seed(123)
library(stochvol)
data(exrates)
ret <- logret(exrates$USD, demean=TRUE)
par(mfrow=c(2, 1), mar = c(1.9, 1.9, 1.9, .5), mgp = c(2, .6, 0))
plot(exrates$date, exrates$USD, type = 'l', main = "Price of 1 EUR in USD")
plot(exrates$date[-1], ret, type = 'l', main = "Demeaned log returns")


###################################################
### code chunk number 3: usd2
###################################################
sim <- svsim(500, mu = -9, phi = 0.99, sigma = 0.1)
par(mfrow=c(2, 1))
plot(sim)


###################################################
### code chunk number 4: article.Rnw:213-215 (eval = FALSE)
###################################################
## res <- svsample(ret, priormu = c(-10, 1), priorphi = c(20, 1.1),
##                 priorsigma = .1)


###################################################
### code chunk number 5: article.Rnw:218-219
###################################################
res <- svsample(ret, priormu = c(-10, 1), priorphi = c(20, 1.1), priorsigma = .1)


###################################################
### code chunk number 6: article.Rnw:241-242
###################################################
summary(res, showlatent = FALSE)


###################################################
### code chunk number 7: usd3
###################################################
volplot(res, forecast = 100, dates = exrates$date[-1])


###################################################
### code chunk number 8: usd4
###################################################
res <- updatesummary(res, quantiles = c(.01, .1, .5, .9, .99))
volplot(res, forecast = 100, dates = exrates$date[-1])


###################################################
### code chunk number 9: usd5
###################################################
par(mfrow = c(3, 1))
paratraceplot(res)


###################################################
### code chunk number 10: usd6
###################################################
par(mfrow = c(1, 3))
paradensplot(res)


###################################################
### code chunk number 11: usd7
###################################################
plot(res)


###################################################
### code chunk number 12: usd8
###################################################
myresid <- resid(res)
plot(myresid, ret)


###################################################
### code chunk number 13: article.Rnw:374-380
###################################################
set.seed(123456)
T <- 1000
beta.true <- c(.1, .5)
sigma.true <- 0.01
X <- matrix(c(rep(1, T), rnorm(T, sd = sigma.true)), nrow = T)
y <- rnorm(T, X %*% beta.true, sigma.true)


###################################################
### code chunk number 14: article.Rnw:384-390
###################################################
draws <- 5000
burnin <- 100
b0 <- matrix(c(0, 0), nrow = ncol(X))
B0inv <- diag(c(10^-10, 10^-10))
c0 <- 0.001
C0 <- 0.001


###################################################
### code chunk number 15: article.Rnw:393-397
###################################################
p <- ncol(X)
preCov <- solve(crossprod(X) + B0inv)
preMean <- preCov %*% (crossprod(X, y) + B0inv %*% b0)
preDf <- c0 + T/2 + p/2


###################################################
### code chunk number 16: article.Rnw:400-403
###################################################
draws1 <- matrix(NA_real_, nrow = draws, ncol = p + 1)
colnames(draws1) <- c(paste("beta", 0:(p-1), sep='_'), "sigma")
sigma2draw <- 1


###################################################
### code chunk number 17: article.Rnw:406-413 (eval = FALSE)
###################################################
## for (i in -(burnin-1):draws) {
##  betadraw <- as.numeric(mvtnorm::rmvnorm(1, preMean, sigma2draw*preCov))
##  tmp <- C0 + .5*(crossprod(y - X%*%betadraw) +
##  		 crossprod((betadraw - b0), B0inv) %*% (betadraw - b0))
##  sigma2draw <- 1/rgamma(1, preDf, rate = tmp)
##  if (i > 0) draws1[i,] <- c(betadraw, sqrt(sigma2draw))
## }


###################################################
### code chunk number 18: article.Rnw:415-427
###################################################
#if (file.exists('~/tmprary/article_draws1.RData')) {
# load('~/tmprary/article_draws1.RData')
#} else {
for (i in -(burnin-1):draws) {
 betadraw <- as.numeric(mvtnorm::rmvnorm(1, preMean, sigma2draw*preCov))
 tmp <- C0 + .5*(crossprod(y - X%*%betadraw) +
 		 crossprod((betadraw - b0), B0inv) %*% (betadraw - b0))
 sigma2draw <- 1/rgamma(1, preDf, rate = tmp)
 if (i > 0) draws1[i,] <- c(betadraw, sqrt(sigma2draw))
}
# save(draws1, file = '~/tmprary/article_draws1.RData')
#}


###################################################
### code chunk number 19: article.Rnw:430-431
###################################################
colMeans(draws1)


###################################################
### code chunk number 20: homo
###################################################
par(mar = c(3.1, 1.8, 1.9, .5), mgp = c(1.8, .6, 0))
plot(coda::mcmc(draws1))


###################################################
### code chunk number 21: article.Rnw:453-458
###################################################
mu.true <- log(sigma.true^2)
phi.true <- 0.97
volvol.true <- 0.3
simresid <- svsim(T, mu = mu.true, phi = phi.true, sigma = volvol.true)
y <- X %*% beta.true + simresid$y


###################################################
### code chunk number 22: article.Rnw:461-467
###################################################
draws <- 50000
burnin <- 1000
thinning <- 10
priormu <- c(-10, 2)
priorphi <- c(20, 1.5)
priorsigma <- 1


###################################################
### code chunk number 23: article.Rnw:470-475
###################################################
draws2 <- matrix(NA_real_, nrow = floor(draws/thinning), ncol = 3 + T + p)
colnames(draws2) <- c("mu", "phi", "sigma", paste("beta", 0:(p-1), sep='_'),
                      paste("h", 1:T, sep='_'))
betadraw <- c(0, 0)
svdraw <- list(para = c(mu = -10, phi = .9, sigma = .2), latent = rep(-10, T))


###################################################
### code chunk number 24: article.Rnw:478-501 (eval = FALSE)
###################################################
## for (i in -(burnin-1):draws) {
## 
##  # draw latent volatilities and AR-parameters:
##  ytilde <- y - X %*% betadraw
##  svdraw <- .svsample(ytilde, startpara=para(svdraw),
##                      startlatent=latent(svdraw), priormu=priormu,
##                      priorphi=priorphi, priorsigma=priorsigma)
## 
##  # draw the betas:
##  normalizer <- as.numeric(exp(-latent(svdraw)/2))
##  Xnew <- X * normalizer
##  ynew <- y * normalizer
##  Sigma <- solve(crossprod(Xnew) + B0inv)
##  mu <- Sigma %*% (crossprod(Xnew, ynew) + B0inv %*% b0)
##  betadraw <- as.numeric(mvtnorm::rmvnorm(1, mu, Sigma))
## 
##  # store the results:
##  if (i > 0 & i %% thinning == 0) { 
##   draws2[i/thinning, 1:3] <- para(svdraw)
##   draws2[i/thinning, 4:5] <- betadraw
##   draws2[i/thinning, 6:(T+5)] <- latent(svdraw)
##  }
## }


###################################################
### code chunk number 25: article.Rnw:504-533
###################################################
if (usePreCalcResults) {
 load('vignette_sampling_draws2.RData')
} else {
for (i in -(burnin-1):draws) {

 # draw latent volatilities and AR-parameters:
 ytilde <- y - X %*% betadraw
 svdraw <- .svsample(ytilde, startpara=para(svdraw),
                     startlatent=latent(svdraw), priormu=priormu,
                     priorphi=priorphi, priorsigma=priorsigma)

 # draw the betas:
 normalizer <- as.numeric(exp(-latent(svdraw)/2))
 Xnew <- X * normalizer
 ynew <- y * normalizer
 Sigma <- solve(crossprod(Xnew) + B0inv)
 mu <- Sigma %*% (crossprod(Xnew, ynew) + B0inv %*% b0)
 betadraw <- as.numeric(mvtnorm::rmvnorm(1, mu, Sigma))

 # store the results:
 if (i > 0 & i %% thinning == 0) { 
  draws2[i/thinning, 1:3] <- para(svdraw)
  draws2[i/thinning, 4:5] <- betadraw
  draws2[i/thinning, 6:(T+5)] <- latent(svdraw)
 }
}
 draws2selection <- draws2[,4:8]
 #save(draws2selection, file = 'vignette_sampling_draws2.RData')
}


###################################################
### code chunk number 26: article.Rnw:536-537 (eval = FALSE)
###################################################
## colMeans(draws2[,4:8])


###################################################
### code chunk number 27: article.Rnw:540-541
###################################################
colMeans(draws2selection)


###################################################
### code chunk number 28: article.Rnw:544-546 (eval = FALSE)
###################################################
## par(mar = c(3.1, 1.8, 1.9, .5), mgp = c(1.8, .6, 0))
## plot(coda::mcmc(draws2[,4:7]))


###################################################
### code chunk number 29: hetero
###################################################
par(mar = c(3.1, 1.8, 1.9, .5), mgp = c(1.8, .6, 0))
plot(coda::mcmc(draws2selection[,1:4]))


###################################################
### code chunk number 30: scatter
###################################################
x <- exrates$USD[-length(exrates$USD)]
y <- exrates$USD[-1]
X <- matrix(c(rep(1, length(x)), x), nrow = length(x))
par(mfrow=c(1,1), mar = c(2.9, 2.9, 2.7, .5), mgp = c(1.8,.6,0), tcl = -.4)
plot(x,y, xlab=expression(p[t]), ylab=expression(p[t+1]),
    main="Scatterplot of lagged daily raw prices of 1 EUR in USD",
    col="#00000022", pch=16, cex=1)
abline(0,1)


###################################################
### code chunk number 31: article.Rnw:609-730
###################################################
if (usePreCalcResults) {
 load('vignette_sampling_realworld.RData')
} else {
set.seed(34567890)

#configuration parameters
draws <- 100000
thinning <- 100
burnin <- 10000
priormu <- c(-10, 1)
priorphi <- c(20, 1.5)
priorsigma <- .1
p <- 2

n <- length(x)

## design matrix:
X <- matrix(c(rep(1, n), x), nrow=n)

## prior for beta (or beta|sigma in homoskedastic regression):
b0 <- matrix(c(0, 0), nrow=ncol(X))
B0 <- diag(c(10^10, 10^10))

## prior for sigma^2 (only used in homoskedastic regression)
c0 <- 0.001
C0 <- 0.001

B0inv <- solve(B0)

## initialize some space to hold results:
realres <- vector('list', 2)

realres[[1]] <- matrix(NA_real_, nrow = floor(draws/thinning), ncol = 3 + n + p)
colnames(realres[[1]]) <- c("mu", "phi", "sigma", paste("beta", 0:(p-1), sep='_'),
			   paste("h", 1:n, sep='_'))

## some indicators:
paras <- 1:3
betas <- 3+(1:p)
latents <- 3+p+(1:n)

## starting values:
betadraw <- rep(.1, p)
svdraw <- list(para = c(mu = -10, phi = .8, sigma = .2),
	       latent = rep(-10, n))

## sampler:
tim <- system.time(
for (i in -(burnin-1):draws) {
 if (i%%1000 == 0) cat("Iteration", i, "done.\n")
 # draw latent volatilities and AR-parameters:
 ytilde <- y - X %*% betadraw
 svdraw <- .svsample(ytilde, startpara=para(svdraw),
		      startlatent=latent(svdraw), priormu=priormu,
		      priorphi=priorphi, priorsigma=priorsigma)

 # draw the betas:
 normalizer <- as.numeric(exp(-latent(svdraw)/2))
 Xnew <- X * normalizer
 ynew <- y * normalizer
 Sigma <- solve(crossprod(Xnew) + B0inv)
 mu <- Sigma %*% (crossprod(Xnew, ynew) + B0inv %*% b0)
 betadraw <- as.numeric(mvtnorm::rmvnorm(1, mu, Sigma))

 # store the results:
 if (i > 0 & i %% thinning == 0) { 
  realres[[1]][i/thinning, paras] <- para(svdraw)
  realres[[1]][i/thinning, latents] <- latent(svdraw)
  realres[[1]][i/thinning, betas] <- betadraw
 }
}
)[["elapsed"]]

## standard sampler:
## pre-calculate some values:
preCov <- solve(crossprod(X) + B0inv)
preMean <- preCov %*% (crossprod(X, y) + B0inv %*% b0)
preDf <- c0 + n/2 + p/2

## allocate space:
realres[[2]] <- matrix(as.numeric(NA), nrow=floor(draws/thinning), ncol=p+1)
colnames(realres[[2]]) <- c(paste("beta", 0:(p-1), sep='_'), "sigma")

## starting values:
betadraw <- rep(0, p)
sigma2draw <- var(y)

## sampler:
tim2 <- system.time(
for (i in -(burnin-1):draws) {
 if (i%%1000 == 0) cat("Iteration", i, "done.\n")
 # draw beta:
 betadraw <- as.numeric(mvtnorm::rmvnorm(1, preMean, sigma2draw*preCov))
 
 # draw sigma^2:
 tmp <- C0 + .5*(crossprod(y - X%*%betadraw) +
		crossprod((betadraw-b0), B0inv) %*% (betadraw-b0))
 sigma2draw <- 1/rgamma(1, preDf, rate=tmp)

 # store the results:
 if (i > 0 & i %% thinning == 0) { 
  realres[[2]][i/thinning, 1:p] <- betadraw
  realres[[2]][i/thinning, p+1] <- sqrt(sigma2draw)
 }
}
)[["elapsed"]]

#Standardized residuals for model checking:
#homoskedastic
predmeans2 <- tcrossprod(X, realres[[2]][,1:2])
standresidmeans2 <- rowMeans((y - predmeans2)/realres[[2]][,3])

#heteroskedastic
predmeans <- tcrossprod(X, realres[[1]][,4:5])
standresidmeans <- rowMeans((y - predmeans)/exp(t(realres[[1]][,6:ncol(realres[[1]])])/2))

realresselection <- vector('list', 2)
realresselection[[1]] <- realres[[1]][,c('beta_0', 'beta_1')]
realresselection[[2]] <- realres[[2]][,c('beta_0', 'beta_1')]
#save(realresselection, standresidmeans, standresidmeans2, file = 'vignette_sampling_realworld.RData')
}


###################################################
### code chunk number 32: betapost
###################################################
smootherfactor <- 1.8
par(mar = c(2.9, 2.9, 2.7, .5), mgp = c(1.8,.6,0), tcl = -.4)
layout(matrix(c(1,2,3,3), byrow=TRUE, nrow=2))
plot(density(realresselection[[1]][,"beta_0"], bw="SJ", adjust=smootherfactor),
    main=bquote(paste("p(",beta[0],"|y)", sep="")),
    xlab="", ylab="", xlim=c(-.0025, .0045))
lines(density(realresselection[[2]][,"beta_0"], bw="SJ", adjust=smootherfactor), lty=2, col=2)
#points(ols$coefficients[1],0, col=2)
#lines(confint(ols)[1,], rep(0,2), col=2)
legend("topleft", c("SV", "homosked."), col=1:2, lty=1:2)

plot(density(realresselection[[1]][,"beta_1"], bw="SJ", adjust=smootherfactor),
    main=bquote(paste("p(",beta[1],"|y)", sep="")),
    xlab="", ylab="", xlim=c(0.9965, 1.0022))
lines(density(realresselection[[2]][,"beta_1"], bw="SJ", adjust=smootherfactor), lty=2, col=2)
#points(ols$coefficients[2],0, col=2)
#lines(confint(ols)[2,], rep(0,2), col=2)
legend("topright", c("SV", "homosked."), col=1:2, lty=1:2)

plotorder <- sample.int(2*nrow(realresselection[[1]]))
cols <- rep(c("#00000055", "#ff000055"), each=nrow(realresselection[[1]]))[plotorder]
pchs <- rep(1:2, each=nrow(realresselection[[1]]))[plotorder]
beta0 <- c(realresselection[[1]][,"beta_0"], realresselection[[2]][,"beta_0"])[plotorder]
beta1 <- c(realresselection[[1]][,"beta_1"], realresselection[[2]][,"beta_1"])[plotorder]

plot(beta0, beta1, col=cols, pch=pchs,
    xlab=bquote(paste("p(",beta[0],"|y)")),
    ylab=bquote(paste("p(",beta[1],"|y)")),
    main="Scatterplot of posterior draws")
legend("topright", c("SV", "homosked."), col=1:2, pch=1:2)


###################################################
### code chunk number 33: qqplot
###################################################
par(mfrow=c(2,2), mar = c(3.1, 3.3, 2.0, .5), mgp = c(1.7,.5,0),
    tcl = -.4)
plot(standresidmeans2, main="Residual scatterplot (homoskedastic errors)",
     xlab="Time", ylab="Standardized residuals")
qqplot(qnorm(ppoints(length(standresidmeans2)), mean=0, sd=1),
       standresidmeans2, main="Residual Q-Q plot (homoskedastic errors)",
       xlab = "Theoretical N(0,1)-quantiles", ylab = "Empirical Quantiles")
abline(0,1)
plot(standresidmeans, main="Residual scatterplot (SV errors)",
     xlab="Time", ylab="Standardized residuals")
qqplot(qnorm(ppoints(length(standresidmeans)), mean=0, sd=1),
       standresidmeans, main="Residual Q-Q plot (SV errors)",
       xlab = "Theoretical N(0,1)-quantiles", ylab = "Empirical Quantiles")
abline(0,1)


