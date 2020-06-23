# Load Package
# Install the package first if the package is not installed
# install.packages("Rlab")
# install.packages("ggplot2")
# install.packages("reshape2")
# install.packages("latex2exp")
# install.packages("tcltk2")
library(Rlab)
library(ggplot2)
library(reshape2)
library(latex2exp)
library(tcltk2)


# Simulation Parameters
## number of iter
N = 500
## number of arms
A = 15


# Model Parameters
mu <- c(0.1,0.2,0.23,0.27,0.32,0.32,0.34,0.41,0.43,0.54,0.55,0.56,0.67,0.71,0.79)
sigma_2 = c(0.05,0.34,0.28,0.09,0.23,0.72,0.19,0.14,0.44,0.53,0.24,0.36,0.56,0.49,0.85)
rho = 0
## a list of rho will used in simulation
rho.list <- c(0,0.001,0.003,0.01,0.03,0.1,0.3,1,3,10,30,100,300,1000)
#c(0,0.001,0.01,0.1,0.3,1,3,10,100,1000)


#B = max(sigma_2)

# Temparary Variables
arm.ts.t <- 0
arm.ucb.t <- 0
#arm.ucb.our.t <- 0
arm.ts.total.t <- 0
arm.ts.variance.t <- 0
arm.ts.total.m.t <- 0

# Support Function
## These functions are used to generate vectorised samples
h <- function(x){
  (x-1-log(x))/2
}
gauss.sample <- function(para.vec){
  rnorm(para.vec[1],para.vec[2],para.vec[3])
}
gamma.sample <- function(para.vec) {
  rgamma(para.vec[1],para.vec[2],para.vec[3])
}
ngamma.mv <- function(para.vec) {
  inv.sig <- rgamma(1, para.vec[3],para.vec[4])
  mean <- rnorm(1, para.vec[1],sqrt(1/(para.vec[2]*inv.sig)))
  1/inv.sig - rho*mean
}

ngamma.mvm <- function(para.vec) {
  inv.sig <- rgamma(1, para.vec[3],para.vec[4])
  mean <- rnorm(1, para.vec[1],sqrt(1/para.vec[2]))
  1/inv.sig - rho*mean
}

# Arm Selection Function
ts.arm <- function(){
  temp.mat <- cbind(rep(1,A), mean.ts, sqrt(1/count.ts.t))
  sample.ts <- apply(temp.mat,1,gauss.sample)
  index.ts <- variance.ts - rho * sample.ts
  which.min(index.ts)
}

ucb.arm <- function(){
  index.ucb <- variance.ucb - rho * mean.ucb - b * sqrt(log(t)/count.ucb.t)
  which.min(index.ucb)
}


#ucb.our.arm <- function(){
#  index.ucb <- variance.our.ucb - rho * mean.our.ucb - sqrt(log(1+t*log(t)^2)/count.our.ucb.t)
#  which.min(index.ucb)
#}

ts.total.arm <- function(){
  temp.mat <- cbind(mean.ts.total, count.ts.total.t, alpha.total.t, beta.total.t)
  index.ts.total <- apply(temp.mat, 1, ngamma.mv)
  which.min(index.ts.total)
}

ts.total.m.arm <- function(){
  temp.mat <- cbind(mean.ts.total, count.ts.total.m.t,alpha.total.m.t,beta.total.m.t)
  index.ts.total.m <- apply(temp.mat, 1, ngamma.mvm)
  which.min(index.ts.total.m)
}

ts.variance.arm <- function(){
  temp.mat <- cbind(rep(1,A), alpha.t, beta.t)
  sample.ts <- apply(temp.mat,1,gamma.sample)
  index.ts <- 1/sample.ts - rho * mean.variance.ts
  which.min(index.ts)
}



# simulation
Time <- 30000
Time.sim <- c(100,200,300,1000,2000,3000,5000,8000,13000,21000,30000)
#1000,2000,3000,5000,8000,13000,21000,30000
x <- matrix(ncol = length(Time.sim),nrow = 0)
reward.ts <- data.frame(x)
reward.variance.ts <- data.frame(x)
reward.ucb <- data.frame(x)
#reward.our.ucb <- data.frame(x)
reward.ts.total <- data.frame(x)
reward.ts.total.m <- data.frame(x)
colnames(reward.ts) <- Time.sim 
colnames(reward.variance.ts) <- Time.sim 
colnames(reward.ucb) <- Time.sim 
colnames(reward.ts.total) <- Time.sim
colnames(reward.ts.total.m) <- Time.sim
#colnames(reward.our.ucb) <- Time.sim 
regret <- NULL

# generate two progress bar to indicate the progress of simulation
p.inner <- tkProgressBar("Inner"," Finished %",  0, 100)
p.outer <- tkProgressBar("Outer", " Finished %", 0, 100)

for (rho in rho.list) {
  info <- sprintf("Finished %d%%", round(which(rho == rho.list)*100/length(rho.list)))
  
  setTkProgressBar(p.outer, round(which(rho == rho.list)*100/length(rho.list)), sprintf("Outer progress (%s)", info), info)
  
  b = sqrt(3)*(2 + rho)/sqrt(2)
  for (iter in 1:N) {
    info <- sprintf("Finished %d%%", round(iter*100/N))
    
    setTkProgressBar(p.inner, iter*100/N, sprintf("Inner progress (%s)", info), info)
    #progress(iter, progress.bar = TRUE)
    reward.ts.t <- rep(0, Time)
    reward.variance.ts.t <- rep(0, Time)
    reward.ucb.t <- rep(0, Time)
    reward.ts.total.t <- rep(0, Time)
    reward.ts.total.m.t <- rep(0, Time)
    #  reward.our.ucb.t <- rep(0, Time)
    count.ts.t <- rep(1,A)
    count.variance.ts.t <- rep(1, A)
    count.ucb.t <- rep(1,A)
    count.ts.total.t <- rep(1,A)
    count.ts.total.m.t <- rep(1,A)
    #  count.our.ucb.t <- rep(1,A)
    mean.ts <- rep(0,A)
    mean.variance.ts <- rep(0,A)
    mean.ucb <- rep(0,A)
    mean.ts.total <- rep(0,A)
    mean.ts.total.m <- rep(0,A)
    #  mean.our.ucb <- rep(0,A)
    variance.ts <- rep(0,A)
    variance.variance.ts <- rep(0,A)
    variance.ucb <- rep(0,A)
    variance.ts.total <- rep(0,A)
    variance.ts.total.m <- rep(0,A)
    #  variance.our.ucb <- rep(0,A)
    alpha.t <- 1/2*rep(1,A)
    beta.t <- 1/2*rep(1,A)
    alpha.total.t <- 1/2*rep(1,A)
    beta.total.t <- 1/2*rep(1,A)
    alpha.total.m.t <- 1/2*rep(1,A)
    beta.total.m.t <- 1/2*rep(1,A)
    #### Pre Algo: first A period
    for (t in 1:A) {
      mean.ts[t] <- rnorm(1,mu[t],sqrt(sigma_2[t]))
      mean.variance.ts[t] <- rnorm(1,mu[t],sqrt(sigma_2[t]))
      mean.ucb[t] <- rnorm(1,mu[t],sqrt(sigma_2[t]))
      mean.ts.total[t] <- rnorm(1,mu[t],sqrt(sigma_2[t]))
      mean.ts.total.m[t] <- rnorm(1,mu[t],sqrt(sigma_2[t]))
      #    mean.our.ucb[t] <- rnorm(1,mu[t],sqrt(sigma_2[t]))
      reward.ts.t[t] <- mean.ts[t]
      reward.variance.ts.t[t] <- mean.variance.ts[t]
      reward.ucb.t[t] <- mean.ucb[t]
      reward.ts.total.t[t] <- mean.ts.total[t]
      reward.ts.total.m.t[t] <- mean.ts.total.m[t]
      #    reward.our.ucb.t[t] <- mean.our.ucb[t]
    }
    for (t in (A+1):Time) {
      # Compute the arm
      arm.ts.t <- ts.arm()
      arm.ts.variance.t <- ts.variance.arm()
      arm.ucb.t <- ucb.arm()
      arm.ts.total.t <- ts.total.arm()
      arm.ts.total.m.t <- ts.total.m.arm()
      #    arm.ucb.our.t <- ucb.our.arm()
      # Record reward
      reward.ts.t[t] <- rnorm(1, mu[arm.ts.t], sqrt(sigma_2[arm.ts.t]))
      reward.variance.ts.t[t] <- rnorm(1, mu[arm.ts.variance.t], sqrt(sigma_2[arm.ts.variance.t]))
      reward.ucb.t[t] <- rnorm(1, mu[arm.ucb.t], sqrt(sigma_2[arm.ucb.t]))
      reward.ts.total.t[t] <- rnorm(1, mu[arm.ts.total.t], sqrt(sigma_2[arm.ts.total.t]))
      reward.ts.total.m.t[t] <- rnorm(1, mu[arm.ts.total.m.t],sqrt(sigma_2[arm.ts.total.m.t]))
      #    reward.our.ucb.t[t] <- rnorm(1, mu[arm.ucb.our.t], sqrt(sigma_2[arm.ucb.our.t]))
      # Update para
      alpha.t[arm.ts.variance.t] <- alpha.t[arm.ts.variance.t] + 1/2
      beta.t[arm.ts.variance.t] <- beta.t[arm.ts.variance.t] + count.variance.ts.t[arm.ts.variance.t]/(count.variance.ts.t[arm.ts.variance.t] + 1) * (reward.variance.ts.t[t]-mean.variance.ts[arm.ts.variance.t])^2/2
      alpha.total.t[arm.ts.total.t] <- alpha.total.t[arm.ts.total.t] + 1/2
      beta.total.t[arm.ts.total.t] <- beta.total.t[arm.ts.total.t] + count.ts.total.t[arm.ts.total.t]/(count.ts.total.t[arm.ts.total.t] + 1) * (reward.ts.total.t[t]-mean.ts.total[arm.ts.total.t])^2/2
      alpha.total.m.t[arm.ts.total.m.t] <- alpha.total.m.t[arm.ts.total.m.t] + 1/2
      beta.total.m.t[arm.ts.total.m.t] <- beta.total.m.t[arm.ts.total.m.t] + count.ts.total.m.t[arm.ts.total.m.t] /(count.ts.total.m.t[arm.ts.total.m.t]+1) * (reward.ts.total.m.t[t] - mean.ts.total.m[arm.ts.total.m.t])^2/2
      ### Update mean-variance
      old.mean.ts <- mean.ts[arm.ts.t]
      old.mean.variance.ts <- mean.variance.ts[arm.ts.variance.t]
      old.mean.ucb <- mean.ucb[arm.ucb.t]
      old.mean.ts.total <- mean.ts.total[arm.ts.total.t]
      old.mean.ts.total.m <- mean.ts.total.m[arm.ts.total.m.t]
      #    old.mean.our.ucb <- mean.our.ucb[arm.ucb.our.t]
      mean.ts[arm.ts.t] <- (reward.ts.t[t] + count.ts.t[arm.ts.t]*mean.ts[arm.ts.t])/(count.ts.t[arm.ts.t] + 1)
      mean.variance.ts[arm.ts.variance.t] <- (reward.variance.ts.t[t] + count.variance.ts.t[arm.ts.variance.t]*mean.variance.ts[arm.ts.variance.t])/(count.variance.ts.t[arm.ts.variance.t] + 1)
      mean.ucb[arm.ucb.t] <- (reward.ucb.t[t] + count.ucb.t[arm.ucb.t]*mean.ucb[arm.ucb.t])/(count.ucb.t[arm.ucb.t] + 1)
      mean.ts.total[arm.ts.total.t] <- (reward.ts.total.t[t] + count.ts.total.t[arm.ts.total.t]*mean.ts.total[arm.ts.total.t])/(count.ts.total.t[arm.ts.total.t] + 1)
      mean.ts.total.m[arm.ts.total.m.t] <- (reward.ts.total.m.t[t] + count.ts.total.m.t[arm.ts.total.m.t]*mean.ts.total.m[arm.ts.total.m.t])/(count.ts.total.m.t[arm.ts.total.m.t]+1)
      #    mean.our.ucb[arm.ucb.our.t] <- (reward.our.ucb.t[t] + count.our.ucb.t[arm.ucb.our.t]*mean.our.ucb[arm.ucb.our.t])/(count.our.ucb.t[arm.ucb.our.t] + 1)
      # Update Variance
      variance.ts[arm.ts.t] <- (count.ts.t[arm.ts.t]*(variance.ts[arm.ts.t] + old.mean.ts^2) + reward.ts.t[t]^2 - (count.ts.t[arm.ts.t] + 1)*mean.ts[arm.ts.t]^2)/(count.ts.t[arm.ts.t] + 1)
      variance.variance.ts[arm.ts.variance.t] <- (count.variance.ts.t[arm.ts.variance.t]*(variance.variance.ts[arm.ts.variance.t] + old.mean.variance.ts^2) + reward.variance.ts.t[t]^2 - (count.variance.ts.t[arm.ts.variance.t] + 1)*mean.variance.ts[arm.ts.variance.t]^2)/(count.variance.ts.t[arm.ts.variance.t] + 1)
      variance.ucb[arm.ucb.t] <- (count.ucb.t[arm.ucb.t]*(variance.ucb[arm.ucb.t] + old.mean.ucb^2) + reward.ucb.t[t]^2 - (count.ucb.t[arm.ucb.t] + 1)*mean.ucb[arm.ucb.t]^2)/(count.ucb.t[arm.ucb.t] + 1)
      variance.ts.total[arm.ts.total.t] <- (count.ts.total.t[arm.ts.total.t]*(variance.ts.total[arm.ts.total.t] + old.mean.ts.total^2) + reward.ts.total.t[t]^2 - (count.ts.total.t[arm.ts.total.t] + 1)*mean.ts.total[arm.ts.total.t]^2)/(count.ts.total.t[arm.ts.total.t] + 1)
      variance.ts.total.m[arm.ts.total.m.t] <- (count.ts.total.m.t[arm.ts.total.m.t]*(variance.ts.total.m[arm.ts.total.m.t] + old.mean.ts.total.m^2)+ reward.ts.total.m.t[t]^2 - (count.ts.total.m.t[arm.ts.total.m.t]+1)*mean.ts.total.m[arm.ts.total.m.t]^2)/(count.ts.total.m.t[arm.ts.total.m.t]+1)
      #    variance.our.ucb[arm.ucb.our.t] <- (count.our.ucb.t[arm.ucb.our.t]*(variance.our.ucb[arm.ucb.our.t] + old.mean.our.ucb^2) + reward.our.ucb.t[t]^2 - (count.our.ucb.t[arm.ucb.our.t] + 1)*mean.our.ucb[arm.ucb.our.t]^2)/(count.our.ucb.t[arm.ucb.our.t] + 1)
      count.ts.t[arm.ts.t] = count.ts.t[arm.ts.t] + 1
      count.variance.ts.t[arm.ts.variance.t] = count.variance.ts.t[arm.ts.variance.t] + 1
      count.ucb.t[arm.ucb.t] = count.ucb.t[arm.ucb.t] + 1
      count.ts.total.t[arm.ts.total.t] = count.ts.total.t[arm.ts.total.t] + 1
      count.ts.total.m.t[arm.ts.total.m.t] = count.ts.total.m.t[arm.ts.total.m.t] + 1
      #    count.our.ucb.t[arm.ucb.our.t] = count.our.ucb.t[arm.ucb.our.t] + 1
    }
    for (k in 1:length(Time.sim)) {
      reward.ts[iter, k] <- var(reward.ts.t[1:Time.sim[k]]) - rho * mean(reward.ts.t[1:Time.sim[k]])
      reward.ucb[iter, k] <- var(reward.ucb.t[1:Time.sim[k]]) - rho * mean(reward.ucb.t[1:Time.sim[k]])
      reward.variance.ts[iter, k] <- var(reward.variance.ts.t[1:Time.sim[k]]) - rho * mean(reward.variance.ts.t[1:Time.sim[k]])
      reward.ts.total[iter, k] <- var(reward.ts.total.t[1:Time.sim[k]]) - rho * mean(reward.ts.total.t[1:Time.sim[k]])
      reward.ts.total.m[iter, k] <- var(reward.ts.total.m.t[1:Time.sim[k]]) - rho * mean(reward.ts.total.m.t[1:Time.sim[k]])
      #    reward.our.ucb[iter, k] <- var(reward.our.ucb.t[1:Time.sim[k]]) - rho * mean(reward.our.ucb.t[1:Time.sim[k]])
    }
  }
  opt <- min(sigma_2 - rho * mu)
  
  regret.rho <- data.frame(Time = Time.sim, 
                           TS.mean = apply(reward.ts,2,mean) - opt,
                           
                           TS.variance = apply(reward.variance.ts,2,mean) - opt,
                           
                           UCB = apply(reward.ucb,2,mean) - opt,
                           
                           TS.total = apply(reward.ts.total,2,mean) - opt,
                           
                           TS.total.m = apply(reward.ts.total.m,2,mean) - opt,
                           
                           TS.mean.sd = apply(reward.ts,2, sd),
                           TS.variance.sd = apply(reward.variance.ts,2, sd),
                           UCB.sd = apply(reward.ucb,2, sd),
                           TS.total.sd = apply(reward.ts.total,2, sd),
                           TS.total.m.sd = apply(reward.ts.total.m,2,sd))
  #,UCB.our = apply(reward.our.ucb,2,mean) - opt
  
  # Plot
  regret.data.long <- melt(regret.rho, id = "Time")
  colnames(regret.data.long) <- c("Time","Algorithms", "Regret")
  regret = rbind(regret, cbind(regret.data.long, rho))
}

# close the progress bar
close(p.inner)
close(p.outer)

regret.value <- regret[!grepl("sd",regret$Algorithms),]
regret.sd <- regret[grepl("sd",regret$Algorithms),]
regret.temp <- cbind(regret.value, low.conf = regret.value$Regret - regret.sd$Regret, up.conf = regret.value$Regret + regret.sd$Regret)

regret.temp$Algorithms <- as.character(regret.temp$Algorithms)
regret.temp$Algorithms[which(regret.temp$Algorithms=="TS.mean")] <- "MTS"
regret.temp$Algorithms[which(regret.temp$Algorithms=="TS.variance")] <- "VTS"
regret.temp$Algorithms[which(regret.temp$Algorithms=="TS.total")] <- "MVTS.o"
regret.temp$Algorithms[which(regret.temp$Algorithms=="UCB")] <- "MV-LCB"
regret.temp$Algorithms[which(regret.temp$Algorithms=="TS.total.m")] <- "MVTS"


# Generate the figure of regret v.s rho
regret.int <- regret.temp[which(regret.temp$Time == Time & regret.temp$Algorithms!="MVTS.o"),]

pd <- position_dodge(0.1)
ggplot(regret.int, aes(x = log(rho), y = log(Time*Regret), colour = Algorithms, linetype = Algorithms), ylab = "regret") + 
  labs(x = latex2exp("log($\\rho$)"), y = "log(Regret)")  + 
  geom_smooth(se = FALSE,size = 1.3) +
  #geom_point() +
  theme_bw() +
  theme(plot.caption = element_text(hjust=0.5, size=25), 
        legend.position = c(0.15, 0.78), 
        axis.text = element_text(size = 23),
        axis.title = element_text(vjust = 1, size = 25),
        legend.text = element_text(size = 23),
        legend.title = element_text(size = 25)) 

# Generate the figure of regret v.s. time
## choose a rho from rho.list first
# rho.list <- c(0,0.001,0.003,0.01,0.03,0.1,0.3,1,3,10,30,100,300,1000)
regret.rrr <- regret.temp[which(regret.temp$rho == rho.list[10] & regret.temp$Algorithms != "MVTS.o"),]

ggplot(regret.rrr, aes(x = Time, y = log(Time*Regret), colour = Algorithms, linetype = Algorithms), ylab = "regret") + 
  geom_line(size = 1.2)  + 
  geom_point(size = 5,aes( shape = Algorithms)) + 
  theme_bw() + ylim(c(5,15))+
  labs(x = "Time", y = "Log(Regret)") +   
  theme(plot.caption = element_text(hjust=0.5, size=25), 
        legend.position = c(0.15, 0.85), 
        axis.text = element_text(size = 23),
        axis.title = element_text(vjust = 1, size = 25),
        legend.text = element_text(size = 23),
        legend.title = element_text(size = 25)) 


