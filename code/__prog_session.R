
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
dsp<-df[, sort(grep('prob.', colnames(df), value=T)[7:12])]
dc<-ds/dp
head(dp); head(ds)

dl<-df[, sort(grep('logit', colnames(df), value=T)[1:6])]
# exp(dl)/(1+exp(dl)) = dp
dsl<-df[, sort(grep('logit', colnames(df), value=T)[7:12])]
#dr<-df[, sort(grep('se.p', colnames(df), value=T)[1:6])]



dq<-NULL; dt<-NULL; colsim<-NULL; groups<-NULL
for(i in 1:nrow(dp)){
  n<-1500
  #p<-t(dp[i,])
  
  #variance<-p * (1-p)
  #covarian<-matrix(outer(p, p), ncol=6)
  #diag(covarian) <-  as.numeric(dsp[i,])
  #diag(covarian)<-variance
  #dt<-rbind(dt, rmvnorm(n, mean=as.numeric(dp[i,]), sigma=covarian))
  
  dq<-rbind(dq, mvtnorm::rmvnorm(n, mean=as.numeric(dl[i,]), sigma=diag(dsl[i,]^2)))
  colsim<-c(colsim, rep(as.vector(t(matrix(rep(colmig, 4), ncol=4)))[i], n))
  groups<-c(groups, rep(paste(df$ed3_ses[i], df$car_mig[i], sep='_'), n))
}

dq<-inv.logit(dq)
dq<-dq/apply(dq, 1, sum)
apply(dp, 2, mean)
apply(dq, 2, mean); summary(apply(dq, 1, sum))

dim(dq)
center<-apply(dp, 2, mean)

acb<-dudi.pca(dp, scannf = FALSE, nf = ncol(dp), center=T)
acp<-dudi.pca(dq, scannf = FALSE, nf = ncol(dp), center=center)


inertia.dudi(acb)
inertia.dudi(acp)

lim<-range(acp$li[,1:2]); cf<-3.8
pch4<-c(1,2,15,16)
lwd_mig<-c(1,1,2,2,2,2)
lty_mig<-c(1,2,3,1,4,2)
colmig<-c("black","gray30","#F98500","#F90000","#00C700","#056C9E")
#colmig<-c("black","gray30",brewer.pal(10, 'Paired')[c(2,6,8,10)])
lt<-c('Never married','Delayed','Normative','Unstable',"Single parent","Early")

coo<-c(1,2)
par(mar=c(2,2,0,0), oma=c(4,1,1,1))
plot(acp$li[,coo], xlim=lim*1.02, ylim=lim*1.02, type='n', xaxt="n", yaxt="n")
mtext(side=1, line=0.5, "First factorial axis, 55% total variance")
mtext(side=2, line=0.5, "Second factorial axis, 22% total variance")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray94")
sds<-apply(acp$li[,coo], 2, sd)
abline(v=seq(-5*sds[1],5*sds[1],sds[1]), h=seq(-5*sds[2],5*sds[2],sds[2]), col='white')
abline(v=0, h=0, col='gray60', lty=2)
#points(acp$li, col="#44444410", pch=16, cex=1)
#points(acp$li, col=alpha(colsim, alpha=0.15), pch=16, cex=1)
#pointLabel(acp$co[,coo]*cf, labels=lt, font=2, cex=1, col='gray70', pos=c(3,3,1,1,3,1))
points(acb$co[,coo]*cf, pch=4, col='gray30')

for(i in 1:16){
  g<-groups==unique(groups)[i]
  #points(acp$li[g, coo], col=alpha(colsim[g], alpha=0.15), pch=16, cex=1)
  ellipse(colMeans(acb$li[i, coo]), cov(acp$li[g, coo]), alpha=0.3, col=unique(colsim[g]), 
          lwd=1, npoints=100, lty=2)
}

for(j in 1:4){
  lines(acb$li[df$car_mig==levels(df$car_mig)[j], coo][1:4,], 
        col=alpha(colmig[j], alpha=.3), type='b', pch=pch4, lty=1, lwd=5)
  
  lines(acb$li[df$car_mig==levels(df$car_mig)[j], coo][1:4,], 
        col=alpha(colmig[j], alpha=1), type='b', pch=pch4, lty=lty_mig[j], lwd=lwd_mig[j])
}
text(acb$co[,coo]*cf, labels=lt, font=2, cex=1, col='gray30', pos=c(3,3,1,3,3,1))


legend(x=-4, y=-4.1, legend=levels(db$car_mig), xpd=NA, ncol=3, lty=lty_mig,
       col=colmig, lwd=3, bty='n', title=expression(bold('Race/ethnicity, age at migration')))
legend(x=2, y=-4.1, legend=levels(db$ed3_ses), xpd=NA, ncol=3, pch=pch4, bty='n',
       title=expression(bold('Educational attainment')))


