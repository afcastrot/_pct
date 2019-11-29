
library(survey); library(nnet); library(multiwayvcov)
library(effects)
library(maptools); library(mvtnorm)
library(TraMineR); library(ade4); library(gdata); library(RColorBrewer)

# Formula to fit the multinomial model
f0<-ml5_fam ~ ed3_ses * car_mig + b10_dem + rel_dem + rac_dem
s0<-summary(m0<-multinom(f0, db, rweig, debug=F))
v0<-vcov(m0)

# Predicted probabilities
df<-data.frame(effect('ed3_ses * car_mig', mod=m0))
dp<-df[, sort(grep('prob.', colnames(df), value=T)[1:6])]
ds<-df[, sort(grep('prob.', colnames(df), value=T)[7:12])]
dc<-ds/dp
head(dp); head(ds)

dl<-df[, sort(grep('logit', colnames(df), value=T)[1:6])]
# exp(dl)/(1+exp(dl)) = dp
ds<-df[, sort(grep('logit', colnames(df), value=T)[7:12])]

dq<-NULL
for(i in 1:nrow(dp)){
  n<-1000
  dq<-rbind(dq, rmvnorm(n, mean=as.numeric(dl[i,]), diag(ds[i, ])^2))
}
dq<-exp(dq)/(1+exp(dq))
apply(dq, 2, mean)
apply(dp, 2, mean)
dim(dq)
center<-apply(dp, 2, mean)

acb<-dudi.pca(dp, scannf = FALSE, nf = ncol(dp), center=T)
acp<-dudi.pca(dq, scannf = FALSE, nf = ncol(dp), center=center)


setwd()
load('example.RData')
library(ade4); library(RColorBrewer)

inertia.dudi(acb)
inertia.dudi(acp)

lim<-range(acp$li[,1:2]); cf<-3.8
pch4<-c(1,2,15,16)
lwd_mig<-c(1,1,2,2,2,2)
lty_mig<-c(1,2,3,1,4,2)
#colmig<-c("black","gray30","#F98500","#F90000","#00C700","#056C9E")
colmig<-c("black","gray30",brewer.pal(10, 'Paired')[c(2,6,8,10)])
lt<-c('Never married','Delayed','Normative','Unstable',"Single parent","Early")

coo<-c(1,2)
par(mar=c(2,2,0.5,0.5), oma=c(4,1,1,1))
plot(acp$li[,coo], xlim=lim*1.02, ylim=lim*1.02, type='n', xaxt="n", yaxt="n")
mtext(side=1, line=0.5, "(-)  Intensity of family trajectories  (+)")
mtext(side=2, line=0.5, "(+)  Normative trajectories  (-)")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray92")
sds<-apply(acp$li[,coo], 2, sd)
abline(v=seq(-5*sds[1],5*sds[1],sds[1]), h=seq(-5*sds[2],5*sds[2],sds[2]), col='white')
abline(v=0, h=0, col='gray60', lty=2)
#points(acp$li, col="#44444410", pch=16, cex=1)
#pointLabel(acp$co[,coo]*cf, labels=lt, font=2, cex=1, col='gray70', pos=c(3,3,1,1,3,1))
points(acb$co[,coo]*cf, pch=4, col='gray30')
text(acb$co[,coo]*cf, labels=lt, font=2, cex=1, col='gray30', pos=c(3,3,1,3,3,1))

for(j in 1:6){
  lines(acb$li[df$car_mig==levels(df$car_mig)[j], coo][1:4,], 
        col=colmig[j], type='b', pch=pch4, lty=lty_mig[j], lwd=lwd_mig[j])
}
legend(x=-4, y=-5, legend=levels(db$car_mig), xpd=NA, ncol=3, lty=lty_mig,
       col=colmig, lwd=3, bty='n', title=expression(bold('Race/ethnicity, age at migration')))
legend(x=2, y=-5, legend=levels(db$ed3_ses), xpd=NA, ncol=3, pch=pch4, bty='n',
       title=expression(bold('Educational attainment')))


