setwd("W:/current_backup_final/discrete games final/python/Dr D/program/MCMC in R")
rm(list=ls(all=TRUE))

options(scipen=999)

### import and filter all non zero rows
data_2013 = read.csv("final_data_2013.csv", stringsAsFactors = FALSE)
zm234_r <- data_2013[which(data_2013$a1 > "0" | data_2013$a2 > "0"),]
#write.csv(zm234_r,"represent_data.csv")

###find z based on a1 and a2
zm = matrix(data=NA,nrow=nrow(zm234_r),ncol=1)
for (i in 1:nrow(zm234_r))
    {
        
if(zm234_r$a1[i]==0 & zm234_r$a2[i]==0)
{
zm[i,]=1
}
if(zm234_r$a1[i]>0 & zm234_r$a2[i]==0)
{
zm[i,]=2
}
if(zm234_r$a1[i]==0 & zm234_r$a2[i]>0)
{
zm[i,]=3
}
if(zm234_r$a1[i]>0 & zm234_r$a2[i]>0)
{
zm[i,]=4
}
}
zm234_r[,ncol(zm234_r)+1]=zm

#write.csv(zm234_r,"represent_data.csv")

