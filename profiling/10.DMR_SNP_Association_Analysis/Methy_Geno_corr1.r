library(data.table)
library(Ropt)
corr=function(start,end)
{
d_sweep=fread("DMR_SWEEP.txt",head=T,data.table = F)
d_methy=fread("DMR_methy.txt",head=T,data.table = F)
outf=qq("{start}_{end}_sig.txt")
g=fread("Teo_Lan_Maize_dp3miss05maf002.txt",head=T,data.table = F,na.strings = "NA")
res=NULL
for(i in start:end)
{cat(i,"\n")
 d1= d_sweep[i,]
 g1=g[g[,1]==d1[1,1] & g[,2]>=d1[1,6] & g[,2]<=d1[1,7],]
 if(nrow(g1)<1){next}
cont=unlist(strsplit(d1[1,9],"_"))[1]
ty=unlist(strsplit(d1[1,9],"_"))[4]
phe=d_methy[d_methy[,1]==cont & d_methy[,6]==ty & d_methy[,2]==d1[1,1] & d_methy[,3]==d1[1,2] & d_methy[,4]==d1[1,3],]
phe=t(phe)
dd=as.data.frame(cbind(colnames(d_methy)[-c(1:6)],phe[-c(1:6),1]))
colnames(dd)=c("ID","Methy")
###100 times bootstrap
sp=NULL
for(t in 1:100)
{
  p=NULL
  for(j in 1:nrow(g1))
  {
    g2=g1[j,]
    g3= as.data.frame(cbind(colnames(g2)[-c(1:4)],t(g2[1,-c(1:4)])))
    colnames(g3)=c("ID","Geno")
    sphe=dd[sample(1:51),]
    p_g=cbind(sphe,g3[,2])
    p_g=na.omit(p_g)
    p_va=cor.test(as.numeric(p_g[,2]),as.numeric(p_g[,3]))$p.value
    p=c(p,p_va)
  }
  p1=min(p,na.rm=T)
  sp=c(sp,p1)
  }
thr=quantile(sp,0.01,na.rm=T)
cat("Finish row ",i," threhold ",thr,"\n")
####################
p=NULL
for(k in 1:nrow(g1))
{
  g2=g1[k,]
 g3= as.data.frame(cbind(colnames(g2)[-c(1:4)],t(g2[1,-c(1:4)])))
 colnames(g3)=c("ID","Geno")
 p_g=merge(dd,g3,by="ID")
 p_g=na.omit(p_g)
 p_va=cor.test(as.numeric(p_g[,2]),as.numeric(p_g[,3]))$p.value
 p=c(p,p_va)
}
p_sig=p[p<=thr]
if(length(p_sig)>0)
{
  r=unlist(c(d1[1,],length(p_sig),min(p_sig)))
}else{
  r=unlist(c(d1[1,],0,"Non-sig"))
  }
res=rbind(res,r)
cat(r,"\n")
cat("Finshed now ",i,"\n")
}
colnames(res)[10:11]=c("sig_num","min_P-value")
write.table(res,file=outf,sep="\t",col.names = T,row.names = F,quote=F)
}
arg=getarg()
print(arg)

do.call(corr,arg)