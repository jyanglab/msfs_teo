### http://rpubs.com/rossibarra/mcmcbc
### By JRI
### http://rpubs.com/rossibarra/179515

library(gsl) #Gnu scientific Library is a collection of numerical routines for scientific computing.
library(coda) #Output Analysis and Diagnostics for MCMC
library(utils)
library(tidyr)

#confluent hypergeometric function 2
f1 <- function(a,b,z){ return(log(hyperg_1F1(a,b,z))) }

#prochhammer symbol
proch <- function(x,n){return(lgamma(x+n)-lgamma(x))}

#Proposal Distribution
normalProposal <- function(theta,std){
    theta.prime = rnorm(n=1,mean=theta,sd=std)
    ln.hastings.ratio = 0.0 
    # symmetrical so hastings ratio is 0. hastings ratio is prob of moving 
    # from A->B vs prob moving B->A which should just be ratio of posteriors 
    # but may not be depending on proposal function
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
        log(choose(n,K))+(f1(beta+K,alpha+beta+n,gamma)+proch(beta,K)+proch(alpha,n-K))-
            (f1(beta,alpha+beta,gamma)+proch(alpha+beta,n))
    })
    #normalize
    prob_SFS=prob_SFS-max(prob_SFS)
    prob_SFS=exp(prob_SFS)/sum(exp(prob_SFS))
    #Make conditional
    if(conditional==TRUE){
        cond_SFS <- sapply(1:length(prob_SFS), function(x) prob_SFS[x]/sum(prob_SFS[2:(length(prob_SFS)-1)]))  
        #divide by total p(seg)
        prob_SFS <- cond_SFS[-c(1,length(cond_SFS))]
    }
    return((dmultinom(my_sfs,prob=prob_SFS,log=TRUE)))
}

MCMCBC <- function(my_sfs, rates, sd, k,
                   conditional=FALSE, Ne, ngen,
                   verbose=TRUE){
    ## Fix some stuff for our model
    
    #conditional=FALSE #use only conditional likelihood
    #Ne=150000 #replace with estimated Ne from SNP data
    #ngen = 100000 # Set the number of generations.
    sample.freq = 100 # Set the sample frequency.
    l.samples = rep(NA,ngen/sample.freq) # Initialize a likelihood vector with length equal to the number of samples.
    p.samples = vector("list", ngen/sample.freq)  # Initialize a prior list with length equal to the number of samples.
    mu.samples=rep(NA,ngen/sample.freq) # Initialize a posterior vector for each param 
    nu.samples = rep(NA,ngen/sample.freq)  #
    s.samples = rep(NA,ngen/sample.freq)  #
    acceptances = vector("list", ngen)  # List of when param changes were accepted.
    
    #rates=c(1E7,1E8,1E6) # rates for mu, nu, s (in that order)
    #sd=c(0.05,0.05,0.05) # lambda for scale proposal dist for mu, nu, s (in that order)
    
    ### Initial values
    params<-rexp(3,rates) # initial values of mu,nu,s (in that order) from exponential
    priors<-dexp(params,rates,log=TRUE) # Get the initial prior values of mu,nu,s (in that order) from exponential
    l=like(conditional,c(0:(length(my_sfs)-1)),Ne,params[1],params[2],params[3],my_sfs) # initial likelihood
    
    ### MCMC BC
    textbar = txtProgressBar(style=3, width=50, file = stderr())
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
            if(verbose) setTxtProgressBar(textbar,(i/ngen)) # Progress bar.
        }
    }
    
    #Calculate acceptances to evaluate parameters of proposal distribution
    
    #Acceptance rate
    mu.acc=round(mean(sapply(acceptances, function (x) x[1]),na.rm=T),3)
    nu.acc=round(mean(sapply(acceptances, function (x) x[2]),na.rm=T),3)
    s.acc=round(mean(sapply(acceptances, function (x) x[3]),na.rm=T),3)
    
    #Posterior vs Prior

    ### posterior SFS
    n <- max(k)
    #Plot maximum a posteriori (mode)
    
    ## ALPHA
    prior.mu=rexp(length(mu.samples[-c(1:(0.1*ngen/sample.freq))]),rates[1])
    post.mu=mu.samples[-c(1:(0.1*ngen/sample.freq))]
    mode.mu=density(post.mu)$x[which(density(post.mu)$y==max(density(post.mu)$y))]
    
    #BETA
    prior.nu=rexp(length(nu.samples[-c(1:(0.1*ngen/sample.freq))]),rates[2])
    post.nu=nu.samples[-c(1:(0.1*ngen/sample.freq))]
    mode.nu=density(post.nu)$x[which(density(post.nu)$y==max(density(post.nu)$y))]
    
    #GAMMA
    prior.s=rexp(length(s.samples[-c(1:(0.1*ngen/sample.freq))]),rates[3])
    post.s=s.samples[-c(1:(0.1*ngen/sample.freq))]
    mode.s=density(post.s)$x[which(density(post.s)$y==max(density(post.s)$y))]
    
    post_sfs=sapply(k,function(K){
        log(choose(n,K))+(f1(mode.nu*4*Ne+K,mode.mu*4*Ne+mode.nu*4*Ne+n,mode.s*4*Ne)+proch(mode.nu*4*Ne,K)+proch(mode.mu*4*Ne,n-K))-(f1(mode.nu*4*Ne,mode.mu*4*Ne+mode.nu*4*Ne,mode.s*4*Ne)+proch(mode.mu*4*Ne+mode.nu*4*Ne,n))})
    
    post_sfs=post_sfs-max(post_sfs)
    post_sfs=exp(post_sfs)/sum(exp(post_sfs))*sum(my_sfs)

    return(list(acc=list(mu=mu.acc, nu=nu.acc, s=s.acc),
                samples=list(s=s.samples, nu=nu.samples, mu=mu.samples, l=l.samples, p=p.samples),
                my_sfs=my_sfs, post_sfs=post_sfs))
    
}


