#histogram of prior and posterior
par=1
fact1 = 1000
fact2 = 100
prd = rnorm(10000,-0.11,cov_mat[par,par]*fact2)
pod = posterior_samples[40000:50000,par]-0.08
pr<-hist(prd,breaks=30)
po<-hist(pod,breaks=30)
xmin = min(min(prd),min(pod))
xmax = max(max(prd),max(pod))
xmin = -0.2
xmax = 0.05
plot(pr,xlim=c(xmin,xmax),ylim=c(0,1000),col="red",main=c("preference parameter for distance(Star)"),cex.lab=1.5, cex.main=2, cex.axis=1.5)
plot(po, add=T, col="blue")
legend('topleft',c('Prior','Posterior'),
       fill = c("red", "green", bty = 'n',
       border = NA)

#histogram of interaction factor
fact1 = 1
fact2 = 1
prd = posterior_samples[40000:50000,9]
pod = posterior_samples[40000:50000,10]
pr<-hist(prd,breaks=30)
po<-hist(pod,breaks=30)
xmin = min(min(prd),min(pod))
xmax = max(max(prd),max(pod))
xmin = -60
xmax = 60
plot(pr,xlim=c(xmin,xmax),ylim=c(0,1500),col="blue",main=c("Interaction factor"),cex.lab=1.5, cex.main=2, cex.axis=1.5)
plot(po, add=T, col="blue")
legend('topleft',c('Star','SkyTeam'),
       fill = c("red", "green", bty = 'n',
       border = NA)

#histogram of airport presence
par=9
fact1 = 1
fact2 = 100
prd = rnorm(10000,mean_vect[par]*fact1,cov_mat[par,par]*fact2)
pod = posterior_samples[40000:50000,par]
pr<-hist(prd,breaks=30)
po<-hist(pod,breaks=30)
xmin = min(min(prd),min(pod))
xmax = max(max(prd),max(pod))
xmin = 0
xmax = 0.4
plot(pr,xlim=c(xmin,xmax),ylim=c(0,1500),col="red",main=c("preference parameter for demand(Star)"),cex.lab=1.5, cex.main=2, cex.axis=1.5)
plot(po, add=T, col="blue")
legend('topleft',c('Prior','Posterior'),
       fill = c("red", "green", bty = 'n',
       border = NA)



#raw individual posterior plots
plot(posterior_samples[,1],main="posterior for cost (Star)",type='o',xlab='MCMC sample',ylab='poserior sample value')
plot(posterior_samples[,2],main="posterior for demand (Star)",type='o')
plot(posterior_samples[,3],main="Scatterplot for third theta component",type='o')
plot(posterior_samples[,4],main="Scatterplot for fourth theta component",type='o')
plot(posterior_samples[,5],main="Scatterplot for fifth theta component",type='o')
plot(posterior_samples[,6],main="Scatterplot for sixth theta component",type='o')

 plot(posterior_samples[,2],main="posterior for demand (Star)",type='o',xlab='MCMC sample',ylab='poserior sample value', cex.lab=1.5, cex.main=1.5, cex.axis=1.5)


#classifying equilibrium regions
burn = 15000;N = 20000;dp=2
B11 = mean(posterior_samples[burn:N,1])
B12 = mean(posterior_samples[burn:N,2])
if(dp>2) {B13 = mean(posterior_samples[burn:N,3])}
B21 = mean(posterior_samples[burn:N,(dp+1)])
B22 = mean(posterior_samples[burn:N,(dp+2)])
if(dp>2) {B23 = mean(posterior_samples[burn:N,(2*dp)])}
L1 = mean(posterior_samples[burn:N,(2*dp+1)])
L2 = mean(posterior_samples[burn:N,(2*dp+2)])
A1 = mean(posterior_samples[burn:N,(2*dp+3)])
A2 = mean(posterior_samples[burn:N,(2*dp+4)])
B = c(B11,B12,B21,B22,del1,del2,A1,A1)
if(dp>2) {B = c(B11,B12,B13,B21,B22,B23,del1,del2,A1,A1)}

R = length(x_11)
eps1 = rnorm(R)
eps2 = rnorm(R)
tie = runif(R)

eqm = matrix(data=NA,nrow=R,ncol=1)
region = matrix(data=NA,nrow=R,ncol=1)

for (r in 1:R)
    {
e1 = eps1[r]
e2 = eps2[r]
xB1 = -(x_11[r]*B11 + x_12[r]*B12)
xB2 = -(x_21[r]*B21 + x_22[r]*B22)
if(dp>2)
{
xB1 = -(x_11[r]*B11 + x_12[r]*B12 + x_13[r]*B13)
xB2 = -(x_21[r]*B21 + x_22[r]*B22 + x_23[r]*B23)
}
xBL1 = xB1-L1
xBL2 = xB2-L2
T[r] = (1+exp(A1+A2*ap[r]))^(-1)
if(e1<=xB1 & e2<=xB2) {region[r] = 1; eqm[r]=1}
if(e1>=xBL1 & e2>=xBL2) {region[r] = 4; eqm[r]=4}
if((e1>=xBL1 & e2<=xBL2)|(e1>=xB1 & e1<=xBL1 & e2<=xB2)) {region[r] = 2; eqm[r] = 2}
if((e1<=xBL1 & e2>=xBL2)|(e1<=xB1 & e2<=xBL2 & e2>=xB2)) {region[r] = 3; eqm[r] = 3}
if(e1<=xBL1 & e1>=xB1 & e2<=xBL2 & e2>=xB2) 
{
region[r] = 5
if(tie[r]<=T[r]) {eqm[r]=2}
if(tie[r]>T[r]) {eqm[r]=3}
}
}

as.data.frame(table(z_data))
as.data.frame(table(eqm))
as.data.frame(table(region))

#plotting all posteriors
layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
plot(posterior_samples[,1]/100,main="preference for distance (Star)",type='o',xlab='MCMC sample number',ylab='poserior sample value',cex.lab=1.5, cex.main=2, cex.axis=1.5)
plot(posterior_samples[,2],main="posterior for demand (Star)",type='o',xlab='MCMC sample',ylab='poserior sample value',cex.lab=1.5, cex.main=2, cex.axis=1.5)
plot(posterior_samples[,3],main="preference for cost (Star)",type='o',xlab='MCMC sample',ylab='poserior sample value',cex.lab=1.5, cex.main=2, cex.axis=1.5)
plot(posterior_samples[,4]/100,main="preference for distance (SkyTeam)",type='o',xlab='MCMC sample number',ylab='poserior sample value',cex.lab=1.5, cex.main=2, cex.axis=1.5)
plot(posterior_samples[,5],main="posterior for demand (SkyTeam)",type='o',xlab='MCMC sample',ylab='poserior sample value',cex.lab=1.5, cex.main=2, cex.axis=1.5)
plot(posterior_samples[,6],main="preference for cost (SkyTeam)",type='o',xlab='MCMC sample',ylab='poserior sample value',cex.lab=1.5, cex.main=2, cex.axis=1.5)

start = 20000; end = 50000
#plotting all posteriors if dec var is 3
layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
plot(samp_v[start:end],posterior_samples[start:end,1]/100,main="preference for distance (Star)",type='o',xlab='MCMC sample number',ylab='poserior sample value',cex.lab=1.5, cex.main=2, cex.axis=1.5)
plot(samp_v[start:end],posterior_samples[start:end,2],main="posterior for demand (Star)",type='o',xlab='MCMC sample',ylab='poserior sample value',cex.lab=1.5, cex.main=2, cex.axis=1.5)
plot(samp_v[start:end],posterior_samples[start:end,3],main="preference for cost (Star)",type='o',xlab='MCMC sample',ylab='poserior sample value',cex.lab=1.5, cex.main=2, cex.axis=1.5)
plot(samp_v[start:end],posterior_samples[start:end,4]/100,main="preference for distance (SkyTeam)",type='o',xlab='MCMC sample number',ylab='poserior sample value',cex.lab=1.5, cex.main=2, cex.axis=1.5)
plot(samp_v[start:end],posterior_samples[start:end,5],main="posterior for demand (SkyTeam)",type='o',xlab='MCMC sample',ylab='poserior sample value',cex.lab=1.5, cex.main=2, cex.axis=1.5)
plot(samp_v[start:end],posterior_samples[start:end,6],main="preference for cost (SkyTeam)",type='o',xlab='MCMC sample',ylab='poserior sample value',cex.lab=1.5, cex.main=2, cex.axis=1.5)


#plotting after burnin interaction
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
plot(samp_v[40000:50000],posterior_samples[40000:50000,7],main="Interaction parameter (Star)",type='o',xlab='MCMC sample',ylab='poserior sample value',cex.lab=1.5, cex.main=2, cex.axis=1.5)
plot(samp_v[40000:50000],posterior_samples[40000:50000,8],main="Interaction parameter (SkyTeam)",type='o',xlab='MCMC sample',ylab='poserior sample value',cex.lab=1.5, cex.main=2, cex.axis=1.5)

#plotting airport presence
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
plot(posterior_samples[,9],main=expression(paste(alpha[1])),type='o',xlab='MCMC sample',ylab='poserior sample value',cex.lab=1.5, cex.main=2, cex.axis=1.5)
plot(posterior_samples[,10],main=expression(paste(alpha[2])),type='o',xlab='MCMC sample',ylab='poserior sample value',cex.lab=1.5, cex.main=2, cex.axis=1.5)


