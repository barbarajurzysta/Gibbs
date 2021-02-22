# 2. Implement the Gibbs sampler sketched above for the Bayesian network in Figure 1 
# and draw 100 samples from the joint probability distribution P(R, C | S = T, W = T)
suppressMessages(library(graph))
suppressMessages(library(gRbase))
suppressMessages(library(Rgraphviz))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(coda))


p.CgRT <- 0.04 / (0.04 + 0.05)  # P(C=T | R=T,S=T,W=T)
p.CgRF <- 0.01 / (0.01 + 0.2)  # P(C=T | R=F,S=T,W=T)
p.RgCT <- 0.792/ (0.792 + 0.18)  # P(R=T | C=T,S=T,W=T)
p.RgCF <- 0.198 / (0.198 + 0.72)  # P(R=T | C=F,S=T,W=T)

mygibbs <- function(n, n_burn, step){
  samples <- matrix(0, 2, n)
  last <- matrix(1,2,1)  # initialise: C=T, R=T
  k <- 1  # first free column in samples matrix
  # sampling
  for (t in 1:(n_burn + n*step)) {
    i = sample(1:2, 1)
    r = runif(1)
    if (i==1){  # sampling C
      if (last[2]==0){  # R(t-1) = F
        if (r < p.CgRF){
          last[1] = 1
        }else{
          last[1] = 0
        }
      } else{  # R(t-1) = T
        if (r < p.CgRT){
          last[1] = 1
        }else{
          last[1] = 0
        }
      }
    }
    else {  # sampling R
      if (last[1]==0){  # C(t-1) = F
        if (r < p.RgCF){
          last[2] = 1
        }else{
          last[2] = 0
        }
      } else{  # C(t-1) = T
        if (r < p.RgCT){
          last[2] = 1
        }else{
          last[2] = 0
        }
      }
    }
    if ((t > n_burn)  & ((t-n_burn) %% step == 0)){
      samples[,k] = last
      k = k+1
    }
  }
  rownames(samples) <- c("C","R")
  return(samples)
}

# drawing 100 samples from the joint probability distribution P(R, C | S = T, W = T)
# burn-in=0, step=1
p.CR <- mygibbs(100, 0, 1)

# 3. Estimate the marginal probability of rain, given that the sprinkler is on 
# and the grass is wet P(R = T | S = T, W = T) from the 100 samples
(p.R <- sum(p.CR[2,])/ncol(p.CR))

samp1 <- mygibbs(50000,0,1)
samp2 <- mygibbs(50000,0,1)
freq1 <- matrix(0, 2, 50000)
freq2 <- matrix(0, 2, 50000)

freq1[,1] = samp1[,1]
freq2[,1] = samp2[,1]
for (t in 2:50000){
  freq1[,t] = freq1[,t-1] + samp1[,t]
  freq2[,t] = freq2[,t-1] + samp2[,t]
  freq1[,t-1] = freq1[,t-1]/(t-1)
  freq2[,t-1] = freq2[,t-1]/(t-1)
}
freq1[,50000] = freq1[,50000]/(50000)
freq2[,50000] = freq2[,50000]/(50000)

rownames(freq1) <- c("C1","R1")
rownames(freq2) <- c("C2","R2")
freqs <- cbind(as.data.frame(t(freq1)), as.data.frame(t(freq2)))

plot1 <- ggplot(freqs, aes(x=1:50000)) + 
  geom_line(aes(y = C1), color = "darkred") + 
  geom_line(aes(y = C2), color="steelblue") +
  labs(x = "t", y = "C=T frequency")
plot2 <- ggplot(freqs, aes(x=1:50000)) + 
  geom_line(aes(y = R1), color = "darkred") + 
  geom_line(aes(y = R2), color="steelblue") +
  labs(x = "t", y = "R=T frequency")
ggsave("convergence.pdf", grid.arrange(plot1, plot2))

n <- 10000
plot1 <- ggplot(freqs[1:n,], aes(x=1:n)) + 
  geom_line(aes(y = C1), color = "darkred") + 
  geom_line(aes(y = C2), color="steelblue") +
  labs(x = "t", y = "C=T frequency")
plot2 <- ggplot(freqs[1:n,], aes(x=1:n)) + 
  geom_line(aes(y = R1), color = "darkred") + 
  geom_line(aes(y = R2), color="steelblue") +
  labs(x = "t", y = "R=T frequency")
ggsave("convergence1.pdf", grid.arrange(plot1, plot2))


s1 <- as.mcmc(t(samp1))
s2 <- as.mcmc(t(samp2))
pdf("gelman.pdf", width = 9, height = 6)
gelman.plot(mcmc.list(s1, s2))
dev.off()

pdf("autocorr_R.pdf", width = 9, height = 6)
acorr <- acf(s1[,'R'], plot = FALSE)
plot(acorr, main = "R autocorrelation")
dev.off()

pdf("autocorr_C.pdf", width = 9, height = 6)
acorr <- acf(s1[,'C'], plot = FALSE)
plot(acorr, main = "C autocorrelation")
dev.off()

# Re-estimating
p.CR <- mygibbs(100, 10000, 15)
(p.R <- sum(p.CR[2,])/ncol(p.CR))