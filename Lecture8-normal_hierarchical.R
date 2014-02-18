
## ----chunk_options, echo=FALSE, message=FALSE----------------------------
opts_chunk$set(fig.width=6, fig.height=5, out.width='.8\\linewidth', fig.align='center', size='tiny')
####################################
# L8 - Normal hierarchical models  #
####################################
library(reshape2)
library(plyr)
library(ggplot2)
library(rjags)


## ----simulation----------------------------------------------------------
J = 10
n_per_group = 9
n = rep(n_per_group,J)
sigma = 1

set.seed(1)
df = data.frame(group = rep(1:J, each=n_per_group))
df$y1 = rnorm(nrow(df))
df$y2 = rnorm(nrow(df),df$group)


## ----data----------------------------------------------------------------
par(mfrow=c(1,2))
plot(df$group, df$y1, main="Sim1: data", xlab="group", ylab="data", pch=19, col=df$group,
     ylim=range(df$y1, df$y2))
plot(df$group, df$y2, main="Sim2: data", xlab="group", ylab="data", pch=19, col=df$group,
     ylim=range(df$y1, df$y2))


## ----exploratory_analysis------------------------------------------------
ddply(df, .(group), summarize, 
      n=length(y1), 
      mean_y1=mean(y1), sd_y1=sd(y1), 
      mean_y2=mean(y2), sd_y2=sd(y2))

df$group = factor(df$group)


## ----analysis, echo=FALSE------------------------------------------------
tau_log_posterior = function(tau, ybar, sigmaj2) 
{
  spt = sigmaj2+tau^2
  Vmu = 1/sum(1/spt)
  mu = sum(ybar/spt)*Vmu
  0.5*log(Vmu)+sum(-0.5*log(spt)-(ybar-mu)^2/(2*spt))
}

V_tau_log_posterior = Vectorize(tau_log_posterior,"tau")

rposterior = function(n_samples, y, gp) 
{
  ybar = by(y,gp,mean)
  n_groups = nlevels(gp)

  # Used throughout
  sigmaj2 = rep(sigma/n_per_group, n_groups)

  # Sample from tau|y
  half_width = 0.05
  tau_xx = seq(0,10,by=2*half_width)+half_width
  tau_log_post = V_tau_log_posterior(tau_xx, ybar, sigmaj2)
  tau_post = exp(tau_log_post)
  tau = tau_xx[sample(1:length(tau_xx), n_samples, replace=T, prob=exp(tau_log_post))]+
         runif(n_samples, -half_width, half_width)
  
  # Sample from mu|tau,y
  Vmu = muhat = rep(NA,n_samples)
  for (i in 1:n_samples)
  {
    spt = sigmaj2 + tau[i]^2
    Vmu[i] = 1/sum(1/spt)
    muhat[i] = sum(ybar/spt)*Vmu[i]
  }
  mu = rnorm(n_samples, muhat, sqrt(Vmu))

  # Sample from theta|mu,tau,y
  theta = matrix(NA, n_samples, n_groups)
  for (i in 1:n_samples)
  {
    tau2 = tau[i]^2
    Vjs = 1/(1/sigmaj2+1/tau2)
    thetahat = (ybar/sigmaj2+mu[i]/tau2)*Vjs
    theta[i,] = rnorm(n_groups, thetahat, sqrt(Vjs))
  }

  return(list(tau=tau, mu=mu, theta=theta))
}

res1 = rposterior(1e4, df$y1, df$group)
res2 = rposterior(1e4, df$y2, df$group)


## ----tau, fig.width=8----------------------------------------------------
par(mfrow=c(1,2))
hist(res1$tau, freq=F, xlim=range(res1$tau,res2$tau), main="Sim1: tau", ylim=c(0,5))
hist(res2$tau, freq=F, xlim=range(res1$tau,res2$tau), main="Sim2: tau", ylim=c(0,5))


## ----, fig.width=8-------------------------------------------------------
par(mfrow=c(1,2))
hist(res1$mu, freq=F, xlim=range(res1$mu,res2$mu), main="Sim1: mu", ylim=c(0,3))
hist(res2$mu, freq=F, xlim=range(res1$mu,res2$mu), main="Sim2: mu", ylim=c(0,3))


## ----theta, eval=FALSE---------------------------------------------------
## par(mfrow=c(1,2))
## plot(0,0, type="n", xlab="theta", ylab="p(theta|y)", xlim=range(res1$theta, res2$theta),
##      main="Sim1: thetas", ylim=c(0,2.5))
## for (i in 1:nlevels(df$group))
## {
##   lines(density(res1$theta[,i]), col=i, lty=i)
## }
## 
## plot(0,0, type="n", xlab="theta", ylab="p(theta|y)", xlim=range(res1$theta, res2$theta),
##      main="Sim2: thetas", ylim=c(0,2.5))
## for (i in 1:nlevels(df$group))
## {
##   lines(density(res2$theta[,i]), col=i, lty=i)
## }


## ----theta_plot, echo=FALSE----------------------------------------------
par(mfrow=c(1,2))
plot(0,0, type="n", xlab="theta", ylab="p(theta|y)", xlim=range(res1$theta, res2$theta),
     main="Sim1: thetas", ylim=c(0,2.5))
for (i in 1:nlevels(df$group)) 
{
  lines(density(res1$theta[,i]), col=i, lty=i)
}

plot(0,0, type="n", xlab="theta", ylab="p(theta|y)", xlim=range(res1$theta, res2$theta),
     main="Sim2: thetas", ylim=c(0,2.5))
for (i in 1:nlevels(df$group)) 
{
  lines(density(res2$theta[,i]), col=i, lty=i)
}


## ----jags, message=FALSE-------------------------------------------------
model = "
model {
  for (i in 1:length(y)) {
    y[i] ~ dnorm(theta[group[i]], 1/sigma^2)
  }

  for (i in 1:ngroups) {
    theta[i] ~ dnorm(mu, 1/tau^2)
  }

  mu    ~ dnorm(0,1e-5)
  sigma ~ dunif(0,1000)
  tau   ~ dunif(0,1000)
}
"
df$group = as.numeric(df$group)

dat = list(group=df$group, y=df$y1, ngroups=max(df$group))
m = jags.model(textConnection(model), dat, n.chains=3)
res1 = coda.samples(m, c("theta","sigma","mu","tau"), 1000)

dat = list(group=df$group, y=df$y2, ngroups=max(df$group))
m = jags.model(textConnection(model), dat, n.chains=3)
res2 = coda.samples(m, c("theta","sigma","mu","tau"), 1000)


