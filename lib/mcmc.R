#Setup
### Libraries etc.
library(gsl)
library(dplyr)
library(tidyr)
library(coda)
library(utils)
library(cowplot)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


### Functions

#confluent hypergeometric
f1<-function(a,b,z){ return(log(hyperg_1F1(a,b,z))) }

#prochhammer symbol
proch<-function(x,n){return(lgamma(x+n)-lgamma(x))}

#Proposal Distribution
normalProposal <- function(theta,std){
    theta.prime = rnorm(n=1,mean=theta,sd=std)
    ln.hastings.ratio = 0.0 # symmetrical so hastings ratio is 0. hastings ratio is prob of moving from A->B vs prob moving B->A which should just be ratio of posteriors but may not be depending on proposal function
    return(c(abs(theta.prime),ln.hastings.ratio))
}

scaleProposal <- function(theta,lambda) {
    m <- exp(lambda * runif(1,-0.5,0.5))
    theta.prime <- theta * m
    ln.hastings.ratio <- log(m) 
    return(c(theta.prime=theta.prime,ln.hastings.ratio=ln.hastings.ratio))
}
#likelihood function
like<-function(conditional,k,Ne,mu,nu,s,my_sfs){
    #Generate SFS for these values
    alpha=4*mu*Ne
    beta=4*nu*Ne
    gamma=4*s*Ne
    n=max(k)
    #unfolded SFS with fixed/absent, equation 16b from Charlesworh & Jain 2014 
    prob_SFS <- sapply(k,function(K){
        log(choose(n,K))+(f1(beta+K,alpha+beta+n,gamma)+proch(beta,K)+proch(alpha,n-K))-(f1(beta,alpha+beta,gamma)+proch(alpha+beta,n))
    })
    
    #normalize
    prob_SFS=prob_SFS-max(prob_SFS)
    prob_SFS=exp(prob_SFS)/sum(exp(prob_SFS))
    
    #Make conditional
    if(conditional==TRUE){
        cond_SFS <- sapply(1:length(prob_SFS), function(x) prob_SFS[x]/sum(prob_SFS[2:(length(prob_SFS)-1)]))  #divide by total p(seg)
        prob_SFS <- cond_SFS[-c(1,length(cond_SFS))]
    }
    return((dmultinom(my_sfs,prob=prob_SFS,log=TRUE)))
}


### Set some values 

conditional=FALSE #use only conditional likelihood
Ne=150000 #replace with estimated Ne from SNP data
ngen = 10000 # Set the number of generations.
sample.freq = 100 # Set the sample frequency.
l.samples = rep(NA,ngen/sample.freq) # Initialize a likelihood vector with length equal to the number of samples.
p.samples = vector("list", ngen/sample.freq)  # Initialize a prior list with length equal to the number of samples.
mu.samples=rep(NA,ngen/sample.freq) # Initialize a posterior vector for each param 
nu.samples = rep(NA,ngen/sample.freq)  #
s.samples = rep(NA,ngen/sample.freq)  #
acceptances = vector("list", ngen)  # List of when param changes were accepted.



# Data
#Either Fake, Assuming sample size 20 chromosomes (10 diploid dudes) with 10K SNPs, or real

fake=FALSE

if(fake==TRUE){
    snps=10000 # only variant sites, used for conditional model only
    sites=500000 # total number of variant and invariant sites, used for complete model (conditional==FALSE) only
    k=0:40
    n=max(k)
    fake.alpha=rexp(1,rates[1])*4*Ne
    fake.beta=rexp(1,rates[2])*4*Ne
    fake.gamma=rexp(1,rates[3])*4*Ne
    
    #neutral
    #my_sfs=(rmultinom(1,theta,(theta/1:(length(k)-2)))) 
    
    my_sfs <- sapply(k,function(K){
        log(choose(n,K))+(f1(fake.beta+K,fake.alpha+fake.beta+n,fake.gamma)+proch(fake.beta,K)+proch(fake.alpha,n-K))-(f1(fake.beta,fake.alpha+fake.beta,fake.gamma)+proch(fake.alpha+fake.beta,n))
    })
    my_sfs=my_sfs-max(my_sfs)
    my_sfs=exp(my_sfs)/sum(exp(my_sfs))
    
    if(conditional==TRUE){
        c_csfs <- sapply(1:length(my_sfs), function(x) my_sfs[x]/sum(my_sfs[2:(length(my_sfs)-1)]))  #divide by 
        my_sfs <- round(c_csfs[-c(1,length(c_csfs))]*snps)
    } else{
        my_sfs <- round(my_sfs*sites)
    }
} else{
    download.file("https://gist.githubusercontent.com/rossibarra/71d0d22bcb6a7c4a786fd99fdf42fcab/raw/47ecd73ec50a92258044618c322d2e83ea5370cb/sfsPC","PCsfs.csv")
    sfs_data<-read.table("PCsfs.csv",header=T)
    my_sfs=sfs_data$Freq
}
k=0:(length(my_sfs)-1)
n=max(k)
plot(my_sfs~k,pch=19,cex=2,ylab="counts",xlab="number of chromosomes",cex.lab=1.5)


# MCMCBC
### Initial Values

### Priors and proposal distribution. $\mu$ is mutation $A_2 \rightarrow A_1$, $\nu$ is $A_1 \rightarrow A_2$, $s$ is advantage of $A_2$ over $A_1$ under simple $1-s$, $1-s/2$, $1$ model for $A_1A_1,A_1A_2$ and $A_2A_2$ genotypes


rates=c(1E7,1E8,1E6) # rates for mu, nu, s (in that order)

sd=c(0.05,0.05,0.05) # lambda for scale proposal dist for mu, nu, s (in that order)
#sd=c(3E-9,7E-10,7E-8) # sd for proposal dist for mu, nu, s (in that order)
params<-rexp(3,rates) # initial values of mu,nu,s (in that order) from exponential

priors<-dexp(params,rates,log=TRUE) # Get the initial prior values of mu,nu,s (in that order) from exponential
l=like(conditional,c(0:(length(my_sfs)-1)),Ne,params[1],params[2],params[3],my_sfs) # initial likelihood


### Run the MCMC

textbar = txtProgressBar(style=ifelse(interactive(),1,3),width=50, file = stderr())

for(i in 1:ngen){ # For each generation...
    #choose which param
    params.prime = params
    random.param=sample(c(1:3),1)
    acceptances[[i]]=c(NA,NA,NA)
    
    # Propose a value based on the previous values.
    proposal = scaleProposal(params[random.param],sd[random.param])
    #proposal = normalProposal(params[random.param],sd[random.param])
    params.prime[random.param] = proposal[1]
    
    # Calculate the proposed likelihood.
    l.prime = like(conditional,k,Ne,params.prime[1],params.prime[2],params.prime[3],my_sfs) 
    
    # Calculate the proposed prior probability.
    priors.prime=dexp(params.prime,rates,log=TRUE) # from exponential
    
    # Calculate the acceptance probability.
    R = (l.prime-l)+(priors.prime[random.param]-priors[random.param]) + proposal[2]
    
    # If r < R, accept the proposed parameters.
    r = runif(1)
    acceptances[[i]][random.param]=0
    if(log(r) < R){ 
        acceptances[[i]][random.param]=1
        params[random.param] = params.prime[random.param] # Set the current value to the proposed value.
        l = l.prime # Set current likelihood to  proposed likelihood.
        priors = priors.prime # Set current prior probability to  proposed prior probability.
    }
    
    # Sample from posterior
    if(i %% sample.freq == 0){ 
        mu.samples[i/sample.freq] = params[1] # Save the current param values.
        nu.samples[i/sample.freq] = params[2]
        s.samples[i/sample.freq] = params[3]
        l.samples[i/sample.freq] = l # Save the current likelihood value.
        p.samples[[i/sample.freq]] = priors # Save the current prior value.
        setTxtProgressBar(textbar,(i/ngen)) # Progress bar.
    }
}


#Calculate acceptances to evaluate parameters of proposal distribution

#Acceptance rate
mu.acc=round(mean(sapply(acceptances, function (x) x[1]),na.rm=T),3)
nu.acc=round(mean(sapply(acceptances, function (x) x[2]),na.rm=T),3)
s.acc=round(mean(sapply(acceptances, function (x) x[3]),na.rm=T),3)


#Traces: toss first 25% as burnin

s.samples=s.samples[(length(s.samples)*0.25):length(s.samples)]
strace=ggplot(data=data.frame(s.samples),aes(y=s.samples,x=1:(length(s.samples))))+
    geom_line(color=cbPalette[4])+
    ylab("s")+
    xlab(paste("generations\n(",round(effectiveSize(10^9*s.samples))," eff. samples)\n(acc. rate ",s.acc,")",sep="")) +
    theme(axis.text=element_text(size=10),axis.title=element_text(size=10))
nu.samples=nu.samples[(length(nu.samples)*0.25):length(nu.samples)]
ntrace=ggplot(data=data.frame(nu.samples),aes(y=nu.samples,x=1:(length(nu.samples))))+
    geom_line(color=cbPalette[3])+
    ylab(expression(nu))+
    xlab(paste("generations\n(",round(effectiveSize(10^9*nu.samples))," eff. samples)\n(acc. rate ",nu.acc,")",sep="")) +
    theme(axis.text=element_text(size=10),axis.title=element_text(size=10))
mu.samples=mu.samples[(length(mu.samples)*0.25):length(mu.samples)]
mtrace=ggplot(data=data.frame(mu.samples),aes(y=mu.samples,x=1:(length(mu.samples))))+
    geom_line(color=cbPalette[2])+
    ylab(expression(mu))+
    xlab(paste("generations\n(",round(effectiveSize(10^9*mu.samples))," eff. samples)\n(acc. rate ",mu.acc,")",sep="")) +
    theme(axis.text=element_text(size=10),axis.title=element_text(size=10))
l.samples=l.samples[(length(l.samples)*0.25):length(l.samples)]
ltrace=ggplot(data=data.frame(l.samples),aes(y=l.samples,x=1:(length(l.samples))))+
    geom_line(color=cbPalette[1])+
    ylab("Likelihood")+
    xlab("generations") +
    theme(axis.text=element_text(size=10),axis.title=element_text(size=10) )



#Posterior vs Prior

#ALPHA
prior.mu=rexp(length(mu.samples),rates[1])
post.mu=mu.samples
mode.mu=density(post.mu)$x[which(density(post.mu)$y==max(density(post.mu)$y))]
muplot<-ggplot(data.frame(post.mu,prior.mu)) +
    geom_histogram(aes(post.mu),fill=cbPalette[2],bins=30) + 
    geom_histogram(aes(prior.mu),bins=30,alpha=0.2,fill=cbPalette[1])+
    scale_x_log10()+
    xlab(expression(mu))
if(fake==TRUE){ muplot=muplot+geom_vline(xintercept = fake.alpha/(4*Ne))} else{ muplot=muplot+geom_vline(xintercept = mean(mu.samples)) }

muplotzoom<-ggplot(data.frame(post.mu)) +
    geom_histogram(aes(post.mu),fill=cbPalette[2],bins=30)+
    xlab(expression(mu))+  
    theme(axis.text=element_text(size=6))
if(fake==TRUE){ muplotzoom=muplotzoom+geom_vline(xintercept = fake.alpha/(4*Ne))} else{ muplotzoom=muplotzoom+geom_vline(xintercept = mean(mu.samples)) }

#BETA
prior.nu=rexp(length(nu.samples),rates[2])
post.nu=nu.samples
mode.nu=density(post.nu)$x[which(density(post.nu)$y==max(density(post.nu)$y))]
nuplot<-ggplot(data.frame(post.nu,prior.nu)) +
    geom_histogram(aes(post.nu),fill=cbPalette[3],bins=30) + 
    geom_histogram(aes(prior.nu),bins=30,alpha=0.2,fill=cbPalette[1])+
    scale_x_log10()+
    xlab(expression(nu))
if(fake==TRUE){ nuplot=nuplot+geom_vline(xintercept = fake.beta/(4*Ne))} else{ 
    nuplot=nuplot+geom_vline(xintercept = mean(nu.samples)) }

nuplotzoom<-ggplot(data.frame(post.nu)) +
    geom_histogram(aes(post.nu),fill=cbPalette[3],bins=30)+
    xlab(expression(nu))+ 
    theme(axis.text=element_text(size=6))
if(fake==TRUE){ nuplotzoom=nuplotzoom+geom_vline(xintercept = fake.beta/(4*Ne))} else{ nuplotzoom=nuplotzoom+geom_vline(xintercept = mean(nu.samples)) }

#GAMMA
prior.s=rexp(length(s.samples),rates[3])
post.s=s.samples
mode.s=density(post.s)$x[which(density(post.s)$y==max(density(post.s)$y))]
splot<-ggplot(data.frame(post.s,prior.s)) + 
    geom_histogram(aes(post.s),fill=cbPalette[4],bins=30) + 
    geom_histogram(aes(prior.s),bins=30,alpha=0.2,fill=cbPalette[1])+
    scale_x_log10()+
    xlab("s")
if(fake==TRUE){ splot=splot+geom_vline(xintercept = fake.gamma/(4*Ne))} else{ splot=splot+geom_vline(xintercept = mean(s.samples)) }

splotzoom<-ggplot(data.frame(post.s)) + 
    geom_histogram(aes(post.s),fill=cbPalette[4],bins=30)+
    xlab("s")+
    theme(axis.text=element_text(size=6))
if(fake==TRUE){ splotzoom=splotzoom+geom_vline(xintercept = fake.gamma/(4*Ne))}  else{ splotzoom=splotzoom+geom_vline(xintercept = mean(s.samples)) }

#PLOT
plot_grid(mtrace,ntrace,strace,muplot,nuplot,splot,muplotzoom,nuplotzoom,splotzoom,ncol=3,rel_heights=c(1.5,1,1),align="v")
plot(ltrace)

#Plot observed data and estimate from mean and MAP:

#plot mean
plot(my_sfs~(c(0:max(k))),pch=19,cex=2,ylab="counts",xlab="number of chromosomes",cex.lab=1.5)

post_sfs=sapply(k,function(K){
    log(choose(n,K))+(f1(mean(post.nu)*4*Ne+K,mean(post.mu)*4*Ne+mean(post.nu)*4*Ne+n,mean(post.s)*4*Ne)+proch(mean(post.nu)*4*Ne,K)+proch(mean(post.mu)*4*Ne,n-K))-(f1(mean(post.nu)*4*Ne,mean(post.mu)*4*Ne+mean(post.nu)*4*Ne,mean(post.s)*4*Ne)+proch(mean(post.mu)*4*Ne+mean(post.nu)*4*Ne,n))})

post_sfs=post_sfs-max(post_sfs)
post_sfs=exp(post_sfs)/sum(exp(post_sfs))*sum(my_sfs)
points(post_sfs~c(0:max(k)),cex=1,col=cbPalette[2],pch=19)
legend("top",legend=c("observed","mean of posterior"),fill=c("black",cbPalette[2]))

#Plot maximum a posteriori (mode)
plot(my_sfs~(c(0:max(k))),pch=19,cex=2,ylab="counts",xlab="number of chromosomes",cex.lab=1.5)
post_sfs=sapply(k,function(K){
    log(choose(n,K))+(f1(mode.nu*4*Ne+K,mode.mu*4*Ne+mode.nu*4*Ne+n,mode.s*4*Ne)+proch(mode.nu*4*Ne,K)+proch(mode.mu*4*Ne,n-K))-(f1(mode.nu*4*Ne,mode.mu*4*Ne+mode.nu*4*Ne,mode.s*4*Ne)+proch(mode.mu*4*Ne+mode.nu*4*Ne,n))})

post_sfs=post_sfs-max(post_sfs)
post_sfs=exp(post_sfs)/sum(exp(post_sfs))*sum(my_sfs)
points(post_sfs~c(0:max(k)),cex=1,col=cbPalette[2],pch=19)
legend("top",legend=c("observed","mode of posterior"),fill=c("black",cbPalette[2]))

#WRITE TO FILE
write.table(file="./posts.txt",data.frame(post.mu,post.nu,post.s),quote=F,row.name=F,sep="\t")

