C:\Users\joseph\Google Drive\Working folder\discrete games\meetings\Dr D\program\results_20000\normalized_data_compress
burn = 20000; N = 50000; dp = 3 #represent 100point7_all3
burn = 2000;N = 5000;dp=2 #represent 5k20point7_demand_junkcost
burn = 1000;N = 5000;dp=2; MV=10 #compressed 5k20point7_demand_junkcost
cu = 852043; csd = 516957; du = 143610; dsd = 156942
BFv = matrix(data=NA,nrow=3,ncol=1)
validate_eqv = matrix(data=NA,nrow=74,ncol=2*MV)
a1p_v = matrix(data=NA,nrow=74,ncol=MV)
a2p_v = matrix(data=NA,nrow=74,ncol=MV)
for(mv in 1:MV)
{
print(mv)
Pr = matrix(data=NA,nrow=(N-burn+1),ncol=5)
Peq = matrix(data=NA,nrow=(N-burn+1),ncol=4)
pr = matrix(data=NA,nrow=length(a1),ncol=5)
peq = matrix(data=NA,nrow=length(a1),ncol=4)
for (rd in 1:length(a1))
{
#xd_11 = x_11[rd]
xd_11=((x_11[rd]*csd+cu)*(1+0.5*mv/MV)-cu)/csd; xd_12 = x_12[rd]; apd = ap[rd];
xd_21 = x_21[rd]; xd_22 = x_22[rd];
#xd_22 = ((x_22[rd]*dsd+du)*100*mv/MV-du)/dsd
if(dp>2) {xd_13 = x_13[rd]; xd_23 = x_23[rd]}
for (s in burn:N)
{
B11 = posterior_samples[s,1]
B12 = posterior_samples[s,2]
if(dp>2) {B13 = posterior_samples[s,3]}
B21 = posterior_samples[s,(dp+1)]
B22 = posterior_samples[s,(dp+2)]
if(dp>2) {B23 = posterior_samples[s,(2*dp)]}
L1 = posterior_samples[s,(2*dp+1)]
L2 = posterior_samples[s,(2*dp+2)]
A1 = posterior_samples[s,(2*dp+3)]
A2 = posterior_samples[s,(2*dp+4)]
B = c(B11,B12,B21,B22,L1,L2,A1,A2)
if(dp>2) {B = c(B11,B12,B13,B21,B22,B23,L1,L2,A1,A2)}

xB1 = -(xd_11*B11 + xd_12*B12)
xB2 = -(xd_21*B21 + xd_22*B22)
if(dp>2)
{
xB1 = -(xd_11*B11 + xd_12*B12 + xd_13*B13)
xB2 = -(xd_21*B21 + xd_22*B22 + xd_23*B23)
}
xBL1 = xB1-L1
xBL2 = xB2-L2
T = (1+exp(A1+A2*apd))^(-1)
Pr[s-burn+1,1] = pnorm(xB1)*pnorm(xB2)
Pr[s-burn+1,2] = pnorm(-xB1)*pnorm(xB2) +  pnorm(-xBL1)*(pnorm(xBL2)-pnorm(xB2))
Pr[s-burn+1,3] = pnorm(xBL1)*pnorm(-xBL2) +  pnorm(xB1)*(pnorm(xBL2)-pnorm(xB2))
Pr[s-burn+1,4] = pnorm(-xBL1)*pnorm(-xBL2)
Pr[s-burn+1,5] = (pnorm(xBL1)-pnorm(xB1))*(pnorm(xBL2)-pnorm(xB2))
Peq[s-burn+1,1] = Pr[s-burn+1,1]
Peq[s-burn+1,2] = Pr[s-burn+1,2]+T*Pr[s-burn+1,5]
Peq[s-burn+1,3] = Pr[s-burn+1,3]+(1-T)*Pr[s-burn+1,5]
Peq[s-burn+1,4] = Pr[s-burn+1,4]
}
pr[rd,1]=mean(Pr[,1]); pr[rd,2]=mean(Pr[,2]); pr[rd,3]=mean(Pr[,3]); pr[rd,4]=mean(Pr[,4]); pr[rd,5]=mean(Pr[,5])
peq[rd,1]=mean(Peq[,1]); peq[rd,2]=mean(Peq[,2]); peq[rd,3]=mean(Peq[,3]); peq[rd,4]=mean(Peq[,4])
#print(rd)
}
z = cbind(a1,a2)
validate_prob = round(cbind(peq,z),2)
mvalidate_prob = round(cbind(pr,z),2)

eqm = matrix(data=NA,nrow=length(a1),ncol=1)
a1p = matrix(0,nrow=length(a1),ncol=1);a2p = matrix(0,nrow=length(a1),ncol=1)
ex = matrix(0,nrow=length(a1),ncol=1); exp = matrix(0,nrow=length(a1),ncol=1)
err = 0; err1 = 0; err2 = 0; BF = 0;

for (i in 1:nrow(peq))
{
eqm[i,1] = which(peq[i,]==max(peq[i,]))
if(eqm[i,1]==2 || eqm[i,1]==1) {a1p[i,1]=1}; end;
if(eqm[i,1]==3 || eqm[i,1]==4) {a2p[i,1]=1}; end;
if(z_data[i]!=1) {ex[i,1]=1}; end;
if(eqm[i,1]!=1) {exp[i,1]=1}; end;
err1 = err1+abs(a1[i]-a1p[i]); err2 = err2+abs(a2[i]-a2p[i]); 
err = err1+abs(ex[i]-exp[i]); 

for (j in 1:4)
{
if(j==z_data[i]) {BF = BF + (1-peq[i,j])^2}
if(j!=z_data[i]) {BF = BF + peq[i,j]^2}
}
}
validate_eq = cbind(z_data,eqm)
BFv[mv] = BF/length(a1)
a1p_v[,mv]=a1p; a2p_v[,mv]=a2p;
validate_eqv[1:74,(2*mv-1):(2*mv)] = validate_eq
}

#hist(eqm,main="demand +10SDP2",xlab="equilibria",border="blue",col="green")
#write.table(validate_eq,"10SDP2_demand.xlsx")
#write.table(validate_eqv,"d2_eq_validate.xlsx")
#write.table(BFv,"d2_BF.xlsx")
#write.table(a1p_v,"d2_1_routes.xlsx")
#write.table(a2p_v,"d2_2_routes.xlsx")