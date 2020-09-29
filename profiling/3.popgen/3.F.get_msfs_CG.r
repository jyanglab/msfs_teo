library(data.table)
meth=fread("/common/jyanglab/shared/Gen_Xu/07-08-2019-mr_100bptile/CG_meth.txt",head=T,data.table=F)
d=meth[,1:3]
for(i in 4:ncol(meth))
{
  meth[which(meth[,i]<0.3),i]=0
  meth[which(meth[,i]>0.7),i]=2
  meth[which(meth[,i]>=0.3 & meth[,i]<=0.7),i]=1
  cat(i,"\n")
}
write.table(meth,file="CG_meth_012.txt",row.names = F,col.names = T,sep="\t",quote=F)
msfs=function(x)
{
  x=x[which(x>0)]
  f0=length(which(x==0))
  f1=length(which(x==1))
  f2=length(which(x==2))
  (f1+2*f2)/(2*length(x)) ###need to check the formula.
}
d1=meth[,-c(1:3)]
m=apply(d1, 1,msfs)
res=cbind(meth,m)
colnames(res[ncol(res)])="mSFS"
write.table(res,file="CG_meth_012_mSFS.txt",row.names = F,col.names = T,sep="\t",quote=F)

