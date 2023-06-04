rm(list=ls())
par(mfrow=c(1,1))
library(MCMCpack)
library(readr)
library(MASS)
library(dplyr)
#五维正态分布
MLE <- function(data)
{
  dim = dim(data)
  myu_estimate = colMeans(data)
  dim_col = dim[2]
  dim_row = dim[1]
  dim_estimate = dim_col + 1
  for (i in 1:dim_col)
  {
    data[,i] = data[,i] - myu_estimate[i]
  }
  sigma_estimate = (t(data) %*% data)/dim_row
  #return(myu_estimate)
  estimate = matrix(NA,ncol = dim_col,nrow = dim_estimate)
  estimate[1:dim_col,] = sigma_estimate
  estimate[dim_estimate,] = myu_estimate
  return(estimate)
}

Gibbs <- function (myu, sigma,burnin,n,x0,by)
{
  dim_x0 = length(x0)
  dim_myu = length(myu)
  dim_sigma = dim(sigma)
  dim_sigma_col = dim_sigma[1]
  dim_sigma12 = dim_sigma - 1
  x_matrix = matrix(NA,ncol = dim_sigma_col,nrow = n)
  if (dim_myu != dim_sigma_col & dim_myu != dim_x0)
  {
    print("dim error")
  }
  else
  {
    for (i in 1:n)
    {
      for (j in 1:dim_sigma_col)
      {
        myu2 = matrix(myu[-j],nrow = dim_sigma12)
        myu1 = myu[j]
        sigma11 = sigma[j,j]
        sigma12 = matrix(sigma[j,-j],ncol = dim_sigma12)
        sigma21 = matrix(sigma[-j,j],nrow = dim_sigma12)
        sigma22 = sigma[-j,-j]
        a = matrix(x0[-j],nrow = dim_sigma12)
        myu_middle = myu1 + sigma12 %*% ginv(sigma22) %*% (a-myu2)
        sigma_middle = sigma11 - sigma12 %*% ginv(sigma22) %*% sigma21
        x0[j] = rnorm(1,mean = myu_middle,sd = sqrt(sigma_middle))
      }
      x_matrix[i,] = x0
    }
  }
  x_numbers = floor((n - burnin)/by)
  x_matrix_final = matrix(NA,ncol = dim_sigma_col,nrow = x_numbers)
  for (h in 1:x_numbers)
  {
    nn = burnin + h * by
    x_matrix_final[h,] = x_matrix[nn,] 
  }
  return(x_matrix_final)
}

myu = c(50,50,50,50,50);myu
sigma = matrix(c(5,0,0,0,0,
                 0,5,2,3,0,
                 0,2,5,0,1,
                 0,3,0,5,0,
                 0,0,1,0,5),ncol = 5);sigma
x0 = c(0,0,0,0,0)
Gibbs_sample_by = Gibbs(myu = myu,sigma = sigma,burnin = 1000,n = 5000,x0 = x0,by = 20)
MLE(Gibbs_sample_by)
sigma
myu
plot(Gibbs_sample_by[,1],type='l')
acf(Gibbs_sample_by[,1])
Gibbs_sample = Gibbs(myu = myu,sigma = sigma,burnin = 1,n = 5000,x0 = x0,by = 1)
#index检验
N=5000-1
av<-c(1:N)
for(i in 1:N)
{
  av[i]<-mean(Gibbs_sample[1:i,1])
}
plot(av,type='l')


#生成手写数字样本
data = read_csv("mnist_train.csv")
data %>%
  count(label)#count labels
data1 = subset(data,data$label==1)#select labels=1
data1 = subset(data1,select=-c(label))#delete label col,keep others
data1_matrix = as.matrix(data1)
dim(data1_matrix)
paraments = MLE(data1_matrix)
dim(paraments)
myu = paraments[785,]
sigma = paraments[1:784,]
x0 = rep(0,784)
Gibbs_sample = Gibbs(myu = myu, sigma = sigma,burnin = 100,n = 1000,x0 = x0,by = 10)
MLE_estimate = MLE(Gibbs_sample)
MLE_estimate[785,]


#Gibbs简单线性回归
#simulated data
x = -15:15
y <- 3 + 5 * x + rnorm(n=length(x),mean=0,sd=7)
plot(x,y)

#parameters setting
N = 31
mu0 = 0
tau0 = 1
mu1 = 0
tau1 = 1
alpha = 2
beta = 1

#posterior distribution sampling
#beta_0
sample_beta_0 = function(beta1,tau){
  
  mean = (tau0*mu0 + tau*sum(y - beta1*x))/(tau0+tau*N)
  sd = sqrt(1/(tau0+tau*N))
  
  rnorm(1,mean = mean, sd = sd)
  
}

#beta_1
sample_beta_1 = function(beta0,tau){
  
  mean = (tau1*mu1+tau*sum((y-beta0)*x))/(tau1+tau*sum(x^2))
  sd = sqrt(1/(tau1+tau*sum(x^2)))
  
  rnorm(1,mean = mean, sd = sd)
  
}


#tau
sample_tau = function(beta0,beta1){
  
  alpha = alpha+N/2
  beta = beta + sum(((y-beta0-beta1*x)^2)/2)
  
  rgamma(1,shape=alpha,rate = beta)
  
}

#iteration
beta0_r = c(0)  #three parameters are set 0,0,2
beta1_r = c(0)
tau_r = c(2)

n = 1e4

for(i in 1:n){
  
  beta0_r = c(beta0_r,sample_beta_0(beta1_r[i], tau_r[i]))
  
  beta1_r = c(beta1_r,sample_beta_1(beta0_r[i+1],tau_r[i]))
  
  tau_r = c(tau_r,sample_tau(beta0_r[i+1],beta1_r[i+1]))
  
}

tau_r = sqrt(1/tau_r)  #sigma^2=1/tau, get sigma
#plots
h = 5000:n

par(mfrow=c(2,3))

plot(beta0_r[h])
abline(h=3,col="red")

plot(beta1_r[h])
abline(h=5,col="red")

plot(tau_r[h])
abline(h=7,col="red")

hist(beta0_r[h])
hist(beta1_r[h])
hist(tau_r[h])
#Gibbs的结果
beta0 = mean(beta0_r[h])
beta1 = mean(beta1_r[h])
beta0;beta1
mean((y-beta0-beta1*x)^2)
mean(tau_r[h])#对tau的估计
#线性回归结果
model = lm(y~x)
mean(model$residuals^2)
model


#MH简单线性回归
#模拟数据
trueA <- 5
trueB <- 0
trueSd <- 10
sampleSize <- 31

# create independent x-values 
x <- (-(sampleSize-1)/2):((sampleSize-1)/2)
# create dependent values according to ax + b + N(0,sd)
y <-  trueA * x + trueB + rnorm(n=sampleSize,mean=0,sd=trueSd)

plot(x,y, main="Test Data")
#样本的似然函数
likelihood <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  
  pred = a*x + b
  singlelikelihoods = dnorm(y, mean = pred, sd = sd, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)   
}

# 参数的先验分布似然
prior <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  aprior = dunif(a, min=0, max=10, log = T)
  bprior = dnorm(b, sd = 5, log = T)
  sdprior = dunif(sd, min=0, max=30, log = T)
  return(aprior+bprior+sdprior)
}

#联合后验分布似然
posterior <- function(param){
  return (likelihood(param) + prior(param))
}


######## Metropolis 算法 ################
proposalfunction <- function(param){            #建议密度函数
  return(rnorm(3,mean = param, sd= c(0.1,0.5,0.3)))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,3))  #按列分开了参数
  chain[1,] = startvalue
  for (i in 1:iterations){
    
    proposal = proposalfunction(chain[i,])
    
    probab = exp(posterior(proposal) - posterior(chain[i,])) #前面取了对数，这里取指数
    
    if (runif(1) < probab){       #使用 mcmc 接受-拒绝样本，获得(beta_0,beta_1,sigma)的多维联合样本
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
    
  }
  return(chain)
}

startvalue = c(4,0,10)
chain = run_metropolis_MCMC(startvalue, 10000)

burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))


par(mfrow = c(2,3))

hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of a", xlab="True value = red line" )
abline(v = mean(chain[-(1:burnIn),1]))
abline(v = trueA, col="red" )

hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of b", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),2]))
abline(v = trueB, col="red" )

hist(chain[-(1:burnIn),3],nclass=30, main="Posterior of sd", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),3]) )
abline(v = trueSd, col="red" )

plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of a", )
abline(h = trueA, col="red" )

plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of b", )
abline(h = trueB, col="red" )

plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = "Chain values of sd", )
abline(h = trueSd, col="red" )


#MCMC logistic
#模拟数据
trueA <- 5
trueB <- 0
trueSd <- 10
sampleSize <- 31

# create independent x-values 
x <- (-(sampleSize-1)/2):((sampleSize-1)/2)
# create dependent values according to ax + b + N(0,sd)
y <-  trueA * x + trueB + rnorm(n=sampleSize,mean=0,sd=trueSd)

plot(x,y, main="Test Data")
#样本的似然函数
likelihood <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  
  pred = a*x + b
  singlelikelihoods = dnorm(y, mean = pred, sd = sd, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)   
}

# 参数的先验分布似然
prior <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  aprior = dunif(a, min=0, max=10, log = T)
  bprior = dnorm(b, sd = 5, log = T)
  sdprior = dunif(sd, min=0, max=30, log = T)
  return(aprior+bprior+sdprior)
}

#联合后验分布似然
posterior <- function(param){
  return (likelihood(param) + prior(param))
}


######## Metropolis 算法 ################
proposalfunction <- function(param){            #建议密度函数
  return(rnorm(3,mean = param, sd= c(0.1,0.5,0.3)))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,3))  #按列分开了参数
  chain[1,] = startvalue
  for (i in 1:iterations){
    
    proposal = proposalfunction(chain[i,])
    
    probab = exp(posterior(proposal) - posterior(chain[i,])) #前面取了对数，这里取指数
    
    if (runif(1) < probab){       #使用 mcmc 接受-拒绝样本，获得(beta_0,beta_1,sigma)的多维联合样本
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
    
  }
  return(chain)
}

startvalue = c(4,0,10)
chain = run_metropolis_MCMC(startvalue, 10000)

burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))


par(mfrow = c(2,3))

hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of a", xlab="True value = red line" )
abline(v = mean(chain[-(1:burnIn),1]))
abline(v = trueA, col="red" )

hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of b", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),2]))
abline(v = trueB, col="red" )

hist(chain[-(1:burnIn),3],nclass=30, main="Posterior of sd", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),3]) )
abline(v = trueSd, col="red" )

plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of a", )
abline(h = trueA, col="red" )

plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of b", )
abline(h = trueB, col="red" )

plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = "Chain values of sd", )
abline(h = trueSd, col="red" )



#MCMC package
data = read_csv("kickData.csv")
head(data,10)
posterior <- MCMClogit(kick_result~kick_distance + game_seconds_remaining + total_home_score + posteam_score, b0=0, B0=.001,
                       data=data)
plot(posterior)
summary(posterior)

model = glm(kick_result~kick_distance + game_seconds_remaining + total_home_score + posteam_score,data=data)
summary(model)