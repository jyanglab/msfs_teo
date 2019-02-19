
library(dplyr)
library(cowplot)

mplot <- function(res, burnin=0.25, rates=c(1E7,1E8,1E6), plot=TRUE){
    
    cbPalette=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    s.samples <- res[['samples']]$s
    nu.samples <- res[['samples']]$nu
    mu.samples <- res[['samples']]$mu
    l.samples <- res[['samples']]$l
    mu.acc <- res[['acc']]$mu
    nu.acc <- res[['acc']]$nu
    s.acc <- res[['acc']]$s
    
    s.samples=s.samples[(length(s.samples)*burnin+1):length(s.samples)]
    strace=ggplot(data=data.frame(s.samples),aes(y=s.samples,x=1:(length(s.samples))))+
        geom_line(color=cbPalette[4])+
        ylab("s")+
        xlab(paste("generations\n(",round(effectiveSize(10^9*s.samples))," eff. samples)\n(acc. rate ",s.acc,")",sep="")) +
        theme(axis.text=element_text(size=10),axis.title=element_text(size=18))
    
    nu.samples=nu.samples[(length(nu.samples)*burnin+1):length(nu.samples)]
    ntrace=ggplot(data=data.frame(nu.samples),aes(y=nu.samples,x=1:(length(nu.samples))))+
        geom_line(color=cbPalette[3])+
        ylab(expression(nu))+
        xlab(paste("generations\n(",round(effectiveSize(10^9*nu.samples))," eff. samples)\n(acc. rate ",nu.acc,")",sep="")) +
        theme(axis.text=element_text(size=10),axis.title=element_text(size=18))
    mu.samples=mu.samples[(length(mu.samples)*burnin+1):length(mu.samples)]
    mtrace=ggplot(data=data.frame(mu.samples),aes(y=mu.samples,x=1:(length(mu.samples))))+
        geom_line(color=cbPalette[2])+
        ylab(expression(mu))+
        xlab(paste("generations\n(",round(effectiveSize(10^9*mu.samples))," eff. samples)\n(acc. rate ",mu.acc,")",sep="")) +
        theme(axis.text=element_text(size=10),axis.title=element_text(size=18))
    l.samples=l.samples[(length(l.samples)*burnin+1):length(l.samples)]
    ltrace=ggplot(data=data.frame(l.samples),aes(y=l.samples,x=1:(length(l.samples))))+
        geom_line(color=cbPalette[1])+
        ylab("Likelihood")+
        xlab("generations") +
        theme(axis.text=element_text(size=10),axis.title=element_text(size=10) )
    
    #####
    fake=FALSE
    prior.mu=rexp(length(mu.samples),rates[1])
    post.mu=mu.samples
    mode.mu=density(post.mu)$x[which(density(post.mu)$y==max(density(post.mu)$y))]
    
    prior.nu=rexp(length(nu.samples),rates[2])
    post.nu=nu.samples
    mode.nu=density(post.nu)$x[which(density(post.nu)$y==max(density(post.nu)$y))]
    
    prior.s=rexp(length(s.samples),rates[3])
    post.s=s.samples
    mode.s=density(post.s)$x[which(density(post.s)$y==max(density(post.s)$y))]
    
    #ALPHA
    muplot<-ggplot(data.frame(post.mu,prior.mu)) +
        geom_histogram(aes(post.mu),fill=cbPalette[2],bins=30) + 
        geom_histogram(aes(prior.mu),bins=30,alpha=0.2,fill=cbPalette[1])+
        scale_x_log10() +
        ylab("") +xlab("") + theme(axis.text=element_text(size=10),axis.title=element_text(size=10))
        #+ xlab(expression(mu))
    if(fake){ muplot=muplot+geom_vline(xintercept = fake.alpha/(4*Ne)) 
      } else{ muplot=muplot+geom_vline(xintercept = mode.mu) }
    
    muplotzoom<-ggplot(data.frame(post.mu)) +
        geom_histogram(aes(post.mu),fill=cbPalette[2],bins=30)+
        xlab(expression(mu))+ ylab("") + 
      theme(axis.text=element_text(size=10),axis.title=element_text(size=18))
    if(fake==TRUE){ muplotzoom=muplotzoom+geom_vline(xintercept = fake.alpha/(4*Ne))
      } else{ muplotzoom=muplotzoom+geom_vline(xintercept = mode.mu) }
    
    #BETA
    nuplot<-ggplot(data.frame(post.nu,prior.nu)) +
        geom_histogram(aes(post.nu),fill=cbPalette[3],bins=30) + 
        geom_histogram(aes(prior.nu),bins=30,alpha=0.2,fill=cbPalette[1])+
        scale_x_log10()+ 
        ylab("") +xlab("") + theme(axis.text=element_text(size=10),axis.title=element_text(size=18))
        #xlab(expression(nu))
    if(fake==TRUE){ nuplot=nuplot+geom_vline(xintercept = fake.beta/(4*Ne)) 
      }else{ nuplot=nuplot+geom_vline(xintercept = mode.nu) }
    
    nuplotzoom<-ggplot(data.frame(post.nu)) +
        geom_histogram(aes(post.nu),fill=cbPalette[3],bins=30)+
        xlab(expression(nu))+ ylab("") +
      theme(axis.text=element_text(size=10),axis.title=element_text(size=18))
    if(fake==TRUE){ nuplotzoom=nuplotzoom+geom_vline(xintercept = fake.beta/(4*Ne)) 
      }else{ nuplotzoom=nuplotzoom+geom_vline(xintercept = mode.nu) }
    
    #GAMMA
    splot<-ggplot(data.frame(post.s,prior.s)) + 
        geom_histogram(aes(post.s),fill=cbPalette[4],bins=30) + 
        geom_histogram(aes(prior.s),bins=30,alpha=0.2,fill=cbPalette[1])+
        scale_x_log10()+ ylab("") + theme(axis.text=element_text(size=10),axis.title=element_text(size=10)) + xlab("")
    if(fake==TRUE){ splot=splot+geom_vline(xintercept = fake.gamma/(4*Ne))
      }else{ splot=splot+geom_vline(xintercept = mode.s) }
    
    splotzoom<-ggplot(data.frame(post.s)) + 
        geom_histogram(aes(post.s),fill=cbPalette[4],bins=30)+
        xlab("s")+ ylab("") +
      theme(axis.text=element_text(size=10),axis.title=element_text(size=18))
    if(fake==TRUE){ splotzoom=splotzoom+geom_vline(xintercept = fake.gamma/(4*Ne))  
      }else{ splotzoom=splotzoom+geom_vline(xintercept = mode.s) }
    
    #PLOT
    message(sprintf("posterior mu [ %s ], nu [ %s ] and s [ %s ]", mode.mu, mode.nu, mode.s))
    if(plot){
      print(plot_grid(mtrace,ntrace,strace,muplot,nuplot,splot,muplotzoom,nuplotzoom,splotzoom,
                      ncol=3,rel_heights=c(1.5,1,1), align="v"))
      return(c(mode.mu, mode.nu, mode.s))
    }else{
      p <- plot_grid(mtrace,ntrace,strace,muplot,nuplot,splot,muplotzoom,nuplotzoom,splotzoom,
                     ncol=3,rel_heights=c(1.5,1,1), align="v")
      return(p)
    }
}

sfsplot <- function(res, burnin=0.2,rates=c(1E8,1E8,1E8), sfsplot=NULL, Ne=150000, k=0:40, plot=TRUE){
    cbPalette=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    s.samples <- res[['samples']]$s
    nu.samples <- res[['samples']]$nu
    mu.samples <- res[['samples']]$mu
    l.samples <- res[['samples']]$l
    mu.acc <- res[['acc']]$mu
    nu.acc <- res[['acc']]$nu
    s.acc <- res[['acc']]$s
    n=max(k)
    
    s.samples=s.samples[(length(s.samples)*burnin+1):length(s.samples)]
    nu.samples=nu.samples[(length(nu.samples)*burnin+1):length(nu.samples)]
    mu.samples=mu.samples[(length(mu.samples)*burnin+1):length(mu.samples)]
    
    l.samples=l.samples[(length(l.samples)*burnin+1):length(l.samples)]
    #####
    fake=FALSE
    prior.mu=rexp(length(mu.samples),rates[1])
    post.mu=mu.samples
    mode.mu=density(post.mu)$x[which(density(post.mu)$y==max(density(post.mu)$y))]
    
    prior.nu=rexp(length(nu.samples),rates[2])
    post.nu=nu.samples
    mode.nu=density(post.nu)$x[which(density(post.nu)$y==max(density(post.nu)$y))]
    
    prior.s=rexp(length(s.samples),rates[3])
    post.s=s.samples
    mode.s=density(post.s)$x[which(density(post.s)$y==max(density(post.s)$y))]
    
    my_sfs <- res$my_sfs
    if(sfsplot == "plotmean"){
      if(plot){
        #plot mean
        plot(my_sfs~(c(0:max(k))), ylim=c(0, max(my_sfs)*1.3), pch=19,cex=2,ylab="counts",xlab="number of chromosomes",cex.lab=1.5)
        post_sfs=sapply(k,function(K){
          log(choose(n,K))+(f1(mean(post.nu)*4*Ne+K,mean(post.mu)*4*Ne+mean(post.nu)*4*Ne+n,mean(post.s)*4*Ne)+proch(mean(post.nu)*4*Ne,K)+proch(mean(post.mu)*4*Ne,n-K))-(f1(mean(post.nu)*4*Ne,mean(post.mu)*4*Ne+mean(post.nu)*4*Ne,mean(post.s)*4*Ne)+proch(mean(post.mu)*4*Ne+mean(post.nu)*4*Ne,n))})
        post_sfs=post_sfs-max(post_sfs)
        post_sfs=exp(post_sfs)/sum(exp(post_sfs))*sum(my_sfs)
        points(post_sfs~c(0:max(k)),cex=1,col=cbPalette[2],pch=19)
        legend("top",legend=c("observed","mean of posterior"),fill=c("black",cbPalette[2]))
      }else{
        #plot mean
        post_sfs=sapply(k,function(K){
          log(choose(n,K))+(f1(mean(post.nu)*4*Ne+K,mean(post.mu)*4*Ne+mean(post.nu)*4*Ne+n,mean(post.s)*4*Ne)+proch(mean(post.nu)*4*Ne,K)+proch(mean(post.mu)*4*Ne,n-K))-(f1(mean(post.nu)*4*Ne,mean(post.mu)*4*Ne+mean(post.nu)*4*Ne,mean(post.s)*4*Ne)+proch(mean(post.mu)*4*Ne+mean(post.nu)*4*Ne,n))})
        post_sfs=post_sfs-max(post_sfs)
        post_sfs=exp(post_sfs)/sum(exp(post_sfs))*sum(my_sfs)
        
        fsize=18
        df1 <- data.frame(alleles= 0:max(k), frq=my_sfs, type="Observed")
        df2 <- data.frame(alleles= 0:max(k), frq=post_sfs, type="Posterior")
        df <- rbind(df1, df2)
        p <- ggplot(df, aes(x=alleles, y=frq, fill=factor(type))) +
          #geom_point() +
          geom_bar(stat = "identity", position = "dodge") +
          #scale_size_manual(values=c(2.5, 1.5)) +
          labs(y="Frequency", x=NULL) + 
          theme(legend.position="none") +
          #scale_y_log10() +
          theme(#axis.text=element_text(size=fsize),
                axis.title=element_text(size=fsize, face="bold"),
                legend.title = element_text(size=fsize, face="bold"),
                legend.text = element_text(size=fsize))
        
        #plot(~(), ylim=c(0, max(my_sfs)*1.3), pch=19,cex=2,ylab="counts",xlab="number of chromosomes",cex.lab=1.5)
        
        #post_sfs=exp(post_sfs)/sum(exp(post_sfs))*sum(my_sfs)
        #points(post_sfs~c(0:max(k)),cex=1,col=cbPalette[2],pch=19)
        #legend("top",legend=c("observed","mean of posterior"),fill=c("black",cbPalette[2]))
        
        return(p)
      }
        
    }
    if(sfsplot == "plotmode"){
        #Plot maximum a posteriori (mode)
        plot(my_sfs~(c(0:max(k))), ylim=c(0, max(my_sfs)*1.3), pch=19,cex=2,ylab="counts",xlab="number of chromosomes",cex.lab=1.5)
        post_sfs=sapply(k,function(K){
            log(choose(n,K))+(f1(mode.nu*4*Ne+K,mode.mu*4*Ne+mode.nu*4*Ne+n,mode.s*4*Ne)+proch(mode.nu*4*Ne,K)+proch(mode.mu*4*Ne,n-K))-(f1(mode.nu*4*Ne,mode.mu*4*Ne+mode.nu*4*Ne,mode.s*4*Ne)+proch(mode.mu*4*Ne+mode.nu*4*Ne,n))})
        post_sfs=post_sfs-max(post_sfs)
        post_sfs=exp(post_sfs)/sum(exp(post_sfs))*sum(my_sfs)
        points(post_sfs~c(0:max(k)),cex=1,col=cbPalette[2],pch=19)
        legend("top",legend=c("observed","mode of posterior"),fill=c("black",cbPalette[2]))
    }
}