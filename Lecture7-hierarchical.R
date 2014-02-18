
## ----chunk_options, echo=FALSE, message=FALSE----------------------------
opts_chunk$set(fig.width=6, fig.height=5, out.width='.8\\linewidth', fig.align='center', size='tiny')
####################################
# L6 - Bayesian hypothesis testing #
####################################
library(reshape2)
library(plyr)
library(ggplot2)
library(rjags)


## ----dawkins_data, eval=FALSE--------------------------------------------
## d = read.csv("dawkins.csv")
## d = rbind(d, data.frame(date=NA, opponent='Total', made=sum(d$made), attempts=sum(d$attempts)))
## d$a = 0.5 + d$made
## d$b = 0.5 + d$attempts-d$made
## d$lcl = qbeta(.025, d$a, d$b)
## d$ucl = qbeta(.975, d$a, d$b)
## d$game = 1:nrow(d)
## p = ggplot(d, aes(x=lcl, xend=ucl, y=game, yend=game,color=1+(opponent=="Total")))+
##   geom_segment(lwd=2)+ theme(legend.position="none")+labs(x=expression(theta))
## print(p)
## d = d[d$opponent!="Total",]


## ----dawkins_data_plot, echo=FALSE---------------------------------------
d = read.csv("dawkins.csv")
d = rbind(d, data.frame(date=NA, opponent='Total', made=sum(d$made), attempts=sum(d$attempts)))
d$a = 0.5 + d$made
d$b = 0.5 + d$attempts-d$made
d$lcl = qbeta(.025, d$a, d$b)
d$ucl = qbeta(.975, d$a, d$b)
d$game = 1:nrow(d)
p = ggplot(d, aes(x=lcl, xend=ucl, y=game, yend=game,color=1+(opponent=="Total")))+
  geom_segment(lwd=2)+ theme(legend.position="none")+labs(x=expression(theta))
print(p)
d = d[d$opponent!="Total",]


## ----proper_prior--------------------------------------------------------
n = 1e4
end = 100
size = rlnorm(n, 0, 1)
mean = runif(n, 0 , 1)
summary(size)
summary(mean)
alpha = size*mean
beta  = size*(1-mean)
summary(alpha)
summary(beta)


## ----proper_prior_plot, fig.width=8--------------------------------------
par(mfrow=c(2,3))
hist(size)
hist(mean)
plot(size, mean, log='x')
hist(alpha)
hist(beta)
plot(alpha,beta)


## ----jags, eval=FALSE----------------------------------------------------
## model = "
## model {
##   for (i in 1:length(y)) {
##     y[i] ~ dbin(theta[i], n[i])
##     theta[i] ~ dbeta(alpha, beta)
##   }
##   alpha <- size*mean
##   beta  <- size*(1-mean)
## 
##   mean ~ dbeta(1,1)
##   size ~ dlnorm(0,1)
## }
## "
## 
## dat = list(y=d$made, n=d$attempts)
## 
## m = jags.model(textConnection(model), data=dat)
## res = coda.samples(m, c("theta","mean","size"), 10000)
## plot(res[,c("mean","size")])


## ----jags_run, echo=FALSE, cache=TRUE, message=FALSE, fig.width=8--------
model = "
model {
  for (i in 1:length(y)) {
    y[i] ~ dbin(theta[i], n[i])
    theta[i] ~ dbeta(alpha, beta)
  }
  alpha <- size*mean
  beta  <- size*(1-mean)

  mean ~ dbeta(1,1)
  size ~ dlnorm(0,1)
}
"

dat = list(y=d$made, n=d$attempts)

m = jags.model(textConnection(model), data=dat)
res = coda.samples(m, c("theta","mean","size"), 10000)
plot(res[,c("mean","size")])


## ----quantiles, dependson="jags_run", eval=FALSE-------------------------
## tmp = data.frame(summary(res)$quantiles)
## d$model = "independent"
## new_d = d
## new_d$model = "hierarchical"
## new_d$lcl = tmp[-c(1:2),1]
## new_d$ucl = tmp[-c(1:2),5]
## combined = rbind(d, new_d)
## 
## e = 0.2
## p = ggplot(combined, aes(x=lcl, xend=ucl, y=game+e*(model=="hierarchical"), yend=game+e*(model=="hierarchical"), color=model))+
##   geom_segment(lwd=2, alpha=0.5) + labs(x=expression(theta), y="game")
## print(p)


## ----quantiles_plot, echo=FALSE------------------------------------------
tmp = data.frame(summary(res)$quantiles)
d$model = "independent"
new_d = d
new_d$model = "hierarchical"
new_d$lcl = tmp[-c(1:2),1]
new_d$ucl = tmp[-c(1:2),5]
combined = rbind(d, new_d)

e = 0.2
p = ggplot(combined, aes(x=lcl, xend=ucl, y=game+e*(model=="hierarchical"), yend=game+e*(model=="hierarchical"), color=model))+
  geom_segment(lwd=2, alpha=0.5) + labs(x=expression(theta), y="game")
print(p)


## ----vague_prior, fig.width=8--------------------------------------------
end = 2
size = exp(runif(n, -10^end, 10^end))
mean = 1/(1+1/exp(runif(n, -10^end, 10^end)))
plot(size,mean,log='x', main="Seemingly vague prior")


## ----jags2, eval=FALSE---------------------------------------------------
## model = "
## model {
##   for (i in 1:length(y)) {
##     y[i] ~ dbin(theta[i], n[i])
##     theta[i] ~ dbeta(alpha, beta)
##   }
##   alpha ~ dunif(0, 1000)
##   beta  ~ dunif(0, 1000)
## 
##   phi <- 5/2 * log(alpha+beta)
##   zero ~ dpois(phi)
## 
##   mean <- alpha/(alpha+beta)
##   size <- alpha+beta
## }
## "
## 
## dat = list(y=d$made, n=d$attempts, zero=0)
## 
## m = jags.model(textConnection(model), data=dat, n.adapt=2e4)
## res = coda.samples(m, c("theta","mean","size"), 2e4)
## plot(res[,c("mean","size")])


## ----jags2_run, echo=FALSE, cache=TRUE, message=FALSE, fig.width=8-------
model = "
model {
  for (i in 1:length(y)) {
    y[i] ~ dbin(theta[i], n[i])
    theta[i] ~ dbeta(alpha, beta)
  }
  alpha ~ dunif(0, 1000)
  beta  ~ dunif(0, 1000)

  phi <- 5/2 * log(alpha+beta)
  zero ~ dpois(phi)

  mean <- alpha/(alpha+beta)
  size <- alpha+beta
}
"

dat = list(y=d$made, n=d$attempts, zero=0)

m = jags.model(textConnection(model), data=dat, n.adapt=2e4)
res = coda.samples(m, c("theta","mean","size"), 2e4)
plot(res[,c("mean","size")])


## ----quantiles2, dependson="jags_run", eval=FALSE------------------------
## tmp = data.frame(summary(res)$quantiles)
## d$model = "independent"
## new_d = d
## new_d$model = "hierarchical"
## new_d$lcl = tmp[-c(1:2),1]
## new_d$ucl = tmp[-c(1:2),5]
## combined = rbind(d, new_d)
## 
## e = 0.2
## p = ggplot(combined, aes(x=lcl, xend=ucl, y=game+e*(model=="hierarchical"), yend=game+e*(model=="hierarchical"), color=model))+
##   geom_segment(lwd=2, alpha=0.5) + labs(x=expression(theta), y="game")
## print(p)


## ----quantiles2_plot, echo=FALSE-----------------------------------------
tmp = data.frame(summary(res)$quantiles)
d$model = "independent"
new_d = d
new_d$model = "hierarchical"
new_d$lcl = tmp[-c(1:2),1]
new_d$ucl = tmp[-c(1:2),5]
combined = rbind(d, new_d)

e = 0.2
p = ggplot(combined, aes(x=lcl, xend=ucl, y=game+e*(model=="hierarchical"), yend=game+e*(model=="hierarchical"), color=model))+
  geom_segment(lwd=2, alpha=0.5) + labs(x=expression(theta), y="game")
print(p)


