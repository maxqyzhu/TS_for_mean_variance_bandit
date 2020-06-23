# Load Package
library(Rlab)
library(ggplot2)
library(reshape2)
library(latex2exp)
library(tcltk2)


# Simulation Parameters
N = 500
A = 15


# Model Parameters
p <- c(0.1,0.2,0.23,0.27,0.32,0.32,0.34,0.41,0.43,0.54,0.55,0.56,0.67,0.71,0.79)
rho = 0
rho.list <- seq(0,1,length = 10)
#c(0,0.001,0.01,0.1,0.3,1,3,10,100,1000)


#B = max(sigma_2)

# Temparary Variables
arm.UCB.t <- 0
arm.MVTS.total.t <- 0

# Support Function
beta.sample <- function(para.vec){
  rbeta(para.vec[1],para.vec[2],para.vec[3])
}

# Arm Selection Function
MVTS.arm <- function(){
  temp.mat <- cbind(rep(1,A), alpha.t, beta.t)
  sample.MVTS <- apply(temp.mat,1,beta.sample)
  index.MVTS <- rho*sample.MVTS - sample.MVTS*(1-sample.MVTS)
  which.max(index.MVTS)
}

UCB.arm <- function(){
  mean.temp <- alpha.UCB.t/count.UCB.t
  index.UCB <- rho*mean.temp - mean.temp*(1-mean.temp) + b * sqrt(log(t)/count.UCB.t)
  which.max(index.UCB)
}

# simulation
Time <- 30000
Time.sim <- c(100,200,300,1000,2000,3000,5000,8000,13000,21000,30000)
#100,200,300,1000,2000,3000,5000,8000,13000,20000
x <- matrix(ncol = length(Time.sim),nrow = 0)
reward.MVTS <- data.frame(x)
reward.UCB <- data.frame(x)
colnames(reward.MVTS) <- Time.sim 
colnames(reward.UCB) <- Time.sim 
regret <- NULL
#regret.MVTS <- NULL
#regret.variance.MVTS <- NULL
#regret.UCB <- NULL
#regret.MVTS.total <- NULL
#regret.our.UCB <- NULL

p.inner <- tkProgressBar("Inner"," Finished %",  0, 100)
p.outer <- tkProgressBar("Outer", " Finished %", 0, 100)

for (rho in rho.list) {
  info <- sprintf("Finished %d%%", round(which(rho == rho.list)*100/length(rho.list)))
  
  setTkProgressBar(p.outer, round(which(rho == rho.list)*100/length(rho.list)), sprintf("Outer progress (%s)", info), info)
  
  b = sqrt(3)*(2 + rho)*sqrt(2)
  
  for (iter in 1:N) {
    info <- sprintf("Finished %d%%", round(iter*100/N))
    setTkProgressBar(p.inner, iter*100/N, sprintf("Inner progress (%s)", info), info)
    #progress(iter, progress.bar = TRUE)
    reward.MVTS.t <- rep(0, Time)
    reward.UCB.t <- rep(0, Time)
    #  reward.our.UCB.t <- rep(0, Time)
    count.MVTS.t <- rep(1,A)
    count.UCB.t <- rep(1,A)
    #  count.our.UCB.t <- rep(1,A)
    mean.MVTS <- rep(0,A)
    mean.UCB <- rep(0,A)
    alpha.t <- rep(1,A)
    beta.t <- rep(1,A)
    alpha.UCB.t <- rep(1,A)
    beta.UCB.t <- rep(1,A)
    #### Pre Algo
    for (t in 1:Time) {
      # Compute the arm
      arm.MVTS.t <- MVTS.arm()
      arm.UCB.t <- UCB.arm()
      # Record reward
      reward.MVTS.t[t] <- rbern(1,p[arm.MVTS.t])
      reward.UCB.t[t] <- rbern(1,p[arm.UCB.t])
      #    reward.our.UCB.t[t] <- rnorm(1, mu[arm.UCB.our.t], sqrt(sigma_2[arm.UCB.our.t]))
      # Update para
      alpha.t[arm.MVTS.t] <- alpha.t[arm.MVTS.t] + reward.MVTS.t[t] 
      beta.t[arm.MVTS.t] <- beta.t[arm.MVTS.t] + 1 - reward.MVTS.t[t] 
      alpha.UCB.t[arm.UCB.t] <- alpha.UCB.t[arm.UCB.t] + reward.UCB.t[t] 
      beta.UCB.t[arm.UCB.t] <- beta.UCB.t[arm.UCB.t] + 1 - reward.UCB.t[t] 
      ### Update mean-variance
      count.MVTS.t[arm.MVTS.t] = count.MVTS.t[arm.MVTS.t] + 1
      count.UCB.t[arm.UCB.t] = count.UCB.t[arm.UCB.t] + 1
    }
    for (k in 1:length(Time.sim)) {
      reward.MVTS[iter, k] <- rho * mean(reward.MVTS.t[1:Time.sim[k]]) - var(reward.MVTS.t[1:Time.sim[k]])
      reward.UCB[iter, k] <- rho * mean(reward.UCB.t[1:Time.sim[k]]) - var(reward.UCB.t[1:Time.sim[k]])
    }
  }
  opt <- max(rho*p-p*(1-p))
  
  regret.rho <- data.frame(Time = Time.sim, 
                           MVTS = opt - apply(reward.MVTS,2,mean),
                           
                           UCB = opt - apply(reward.UCB,2,mean),
                           
                           UCB.sd = apply(reward.UCB,2, sd),
                           MVTS.sd = apply(reward.MVTS,2, sd)
                           )

  # Plot
  regret.data.long <- melt(regret.rho, id = "Time")
  colnames(regret.data.long) <- c("Time","Algorithms", "Regret")
  regret = rbind(regret, cbind(regret.data.long, rho))
}

close(p.inner)
close(p.outer)

regret.value <- regret[!grepl("sd",regret$Algorithms),]
regret.sd <- regret[grepl("sd",regret$Algorithms),]
regret.temp <- cbind(regret.value, low.conf = regret.value$Regret - regret.sd$Regret, up.conf = regret.value$Regret + regret.sd$Regret)

regret.temp$Algorithms <- as.character(regret.temp$Algorithms)
regret.temp$Algorithms[which(regret.temp$Algorithms=="MVTS")] <- "BMVTS"
regret.temp$Algorithms[which(regret.temp$Algorithms=="UCB")] <- "LCB"


regret.int <- regret.temp[which(regret.temp$Time == Time),]
pd <- position_dodge(0.1)
ggplot(regret.int, aes(x = rho, y = log(Time*Regret)), ylab = "regret") + 
  labs(x = latex2exp("$\\rho$"), y = "log(Regret)")  + 
  geom_smooth(se = FALSE,size = 1.3, aes(linetype = Algorithms,colour = Algorithms)) +
  scale_color_manual(values = c('skyblue', 'red')) + 
  theme_bw() + ylim(5,15) + 
  theme(plot.caption = element_text(hjust=0.5, size=25), 
        legend.position = c(0.15, 0.78), 
        axis.text = element_text(size = 23),
        axis.title = element_text(vjust = 1, size = 25),
        legend.text = element_text(size = 23),
        legend.title = element_text(size = 25)) 

regret.rrr <- regret.temp[which(regret.temp$rho == rho.list[1]),]

ggplot(regret.rrr, aes(x = Time, y = log(Time*Regret), colour = Algorithms, linetype = Algorithms), ylab = "regret") + 
  geom_line(size = 1.2)  + 
  scale_fill_manual(values=c("red", "blue")) +
  geom_point(size = 5,aes(shape = Algorithms)) + 
  theme_bw() + ylim(c(1,10))+
  labs(x = "Time", y = "Log(Regret)") +   
  theme(plot.caption = element_text(hjust=0.5, size=25), 
        legend.position = c(0.1, 0.85), 
        axis.text = element_text(size = 23),
        axis.title = element_text(vjust = 1, size = 25),
        legend.text = element_text(size = 23),
        legend.title = element_text(size = 25)) 


