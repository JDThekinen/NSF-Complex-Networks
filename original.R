# discrete games
#date: "October 12, 2016"

# First defining the prior over the parameters. This joint prior 
#for the independednt variables is coded same as the one given in 
# the paper by Sugaware et al.

setwd("W:/current_backup_final/discrete games final/python/Dr D/program/MCMC in R")
cov_mat = matrix(c(0.00011,0,0,0,0,0,
                   0,0.00016,0,0,0,0,
                   0,0,0.00011,0,0,0,
                   0,0,0,0.00016,0,0,
                   0,0,0,0,72,0,
                   0,0,0,0,0,72),
                 nrow=6,ncol=6,byrow=TRUE)

# Prior over the 6 parameters
prior = function(theta)
{
    library(tmvtnorm)
    #cov_mat = diag(1,6,6)
    value = dtmvnorm(theta, 
                     mean = c(-0.00044,0.00087,-0.00044,0.00087,-100,-100),
                     lower = rep(-Inf,6), 
                     upper = c(Inf,Inf,Inf,Inf,0,0), 
                     sigma = cov_mat, 
                     log=TRUE)
    return(value)
}

# Likelihood using the data for M markets/data points
data = read.csv('joseph_data_old.csv', header=TRUE)
z_data = data$zm
#show(z_data)

# finalizing the values of the decision variables from the data sheet
x_11 = data$distance_1000
x_12 = data$demand_1000
x_21 = data$distance_1000
x_22 = data$demand_1000
x_data = cbind(x_11,x_12,x_21,x_22)
p1 = function(theta,x)
{
    value = pnorm(-t(theta[1:2])%*%x[1:2],log=TRUE) + pnorm(-t(theta[3:4])%*%x[3:4], log=TRUE)
    #if (value==0)
    #  {
    #  value = 0.000001
    #}
    #else{
    #  value = value[1]
    #}
    #print (1)
    #show(t(theta[1:2])%*%x[1:2])
    #show(value)
    #show(log(value))
    return((value[1]))
}

p2  = function(theta,x)
{
    value = pnorm(t(theta[1:2])%*%x[1:2]) - pnorm(t(theta[1:2])%*%x[1:2] + theta[5]) * (pnorm(-t(theta[3:4])%*%x[3:4])) +(pnorm(t(theta[1:2])%*%x[1:2] + theta[5]))*(pnorm(-t(theta[3:4])%*%x[3:4] - theta[6]))
    if (value==0)
    {
        value = 0.000001
    }
    else{
        value = value[1]
    }
    #print (2)
    show(value)
    #show(log(value))
    return((value))
}

p3 = function(theta,x)
{
    
    value = (pnorm(t(theta[1:2])%*%x[1:2]) - pnorm(t(theta[1:2])%*%x[1:2] + theta[5]))*pnorm((t(theta[3:4])%*%x[3:4] + theta[6])) + pnorm(-t(theta[1:2])%*%x[1:2])*(pnorm(t(theta[3:4])%*%x[3:4]))
    if (value==0){
        value = 0.000001
    }
    else{
        value = value[1]
    }
    #print (3)
    #show(t(theta[1:2])%*%x[1:2])
    #show(value)
    #show(log(value))
    return((value))
}

p4 = function(theta,x)
{
    value = pnorm(t(theta[1:2])%*%x[1:2] + theta[5], log=TRUE) + pnorm(t(theta[3:4])%*%x[3:4] + theta[6],log=TRUE)
    #if (value==0){
    #  value = 0.000001
    #}
    #else{
    #  value = value[1]
    #}
    #print (4)
    #show(t(theta[1:2])%*%x[1:2])
    #show(value)
    #show(log(value))
    return((value[1]))
}

p5 = function(theta,x)
{
    value = (pnorm(t(theta[1:2])%*%x[1:2])-pnorm(t(theta[1:2])%*%x[1:2] + theta[5]))*(pnorm(t(theta[3:4])%*%x[3:4]) - pnorm(t(theta[1:2])%*%x[1:2] + theta[6]))
    #print (5)
    if (value==0){
        value = 0.000001
    }
    else{
        value = value[1]
    }
    #show(t(theta[1:2])%*%x[1:2])
    show(value)
    #show(log(value))
    return((value))
}

# Likelihood
likelihood = function(theta,lambda)
{
    # this is the likelihood of z|theta,lambda eqn. A.6 of the paper
    value  <- vector(mode='numeric',length=length(z_data))
    for (i in 1:length(z_data)){
        if (z_data[i]==1)
        {
            likelihood = p1(theta,x_data[i,])
            print (31)
            show(likelihood)
            value[i] = likelihood 
        }
        
        if (z_data[i]==2)
        {
            likelihood = p2(theta,x_data[i,]) + lambda[i]*p5(theta,x_data[i,])
            if (likelihood==0){
                likelihood = 1e-6
            }
            print (32)
            show(log(likelihood))
            value[i] = log(likelihood) 
        }
        
        if (z_data[i]==3)
        {
            likelihood = p3(theta,x_data[i,]) + lambda[i]*p5(theta,x_data[i,])
            if (likelihood==0){
                likelihood = 1e-6
            }
            print (33)
            show(log(likelihood))
            value[i] = log(likelihood) 
        }
        
        if (z_data[i]==4)
        {
            likelihood = p4(theta,x_data[i,])
            print (34)
            show(likelihood)
            value[i] = likelihood 
        }
        #show(likelihood)
        #show(value)
    }
    #show(sum(value))
    #show(dim(value))
    return(sum(value))
}

# Multiplication term
a1 = data$a1
a2 = data$a2

mpl_term = function(lambda,p)
{
    value <- vector(mode='numeric',length=length(z_data))
    for (i in 1:length(z_data))
    {
        mpl_1 = (p[i]^(lambda[i] +a1[i] -1 ))
        #print (11)
        #show(p[i])
        #show(lambda[i])
        mpl_2 = ((1-p[i])^(1-lambda[i]+a2[i]-1))
        if (mpl_1==0){
            mpl_1=1e-6
        }
        if (mpl_2==0){
            mpl_2=1e-6
        }
        mpl = mpl_1*mpl_2
        value[i] = log(mpl)
    }
    show(sum(value))
    return(sum(value))
}

#Coding the joint posterior of the parameters: equation A.7 of the paper.
joint_posterior = function(theta,lambda,p)
{
    value = prior(theta) + likelihood(theta,lambda) 
    #+ mpl_term(lambda,p)
    print (3)
    #show(prior(theta))
    show(likelihood(theta,lambda))
    show(mpl_term(lambda,p))
    return(value)
}


s_likelihood = function(theta,lambda,i)
{
    if (z_data[i]==1)
    {
        likelihood = p1(theta,x_data[i,])
        #show(likelihood)
        value = likelihood 
    }
    if (z_data[i]==2)
    {
        likelihood = p2(theta,x_data[i,]) + lambda*p5(theta,x_data[i,])
        #show(likelihood)
        value = likelihood 
    }
    if (z_data[i]==3)
    {
        likelihood = p3(theta,x_data[i,]) + lambda*p5(theta,x_data[i,])
        #show(likelihood)
        value = likelihood 
    }
    if (z_data[i]==4)
    {
        likelihood = p4(theta,x_data[i,])
        #show(likelihood)
        value = likelihood 
    }
    return(exp(value))
}

# Coding the conditional posterior for lambda for sampling lambda using the gibbs sampler :equation A.8 of the paper.
cond_post_lambda = function(theta,p,z)
{
    lambda <- vector(mode='numeric',length=length(z_data))
    print (1)
    #show(p)
    for (i in 1:length(z_data)){
        q_1 = ((p[i]^(a1[i]))*((1-p[i])^(a2[i]-1))*s_likelihood(theta,lambda=1,i=i))
        q_2 = ((p[i]^(a1[i]))*((1-p[i])^(a2[i]-1))*s_likelihood(theta,lambda=1,i=i)) + ((p[i]^(a1[i]-1))*((1-p[i])^(a2[i]))*s_likelihood(theta,lambda=0,i=i))
        if (q_2==0)
        {
            q_2 = 1e-4
        }
        q_m = q_1/q_2
        print (2)
        show(q_m)
        library(Rlab)
        #show(q_m)
        if (is.na(q_m)){
            q_m=1
        }
        lambda[i] = rbern(1,q_m)
    }
    #print (2)
    #show(lambda)
    return(lambda)
}

# Coding the conditional posterior for p_m for sampling p_m using the gibbs sampler: equation A.9 of the paper
cond_post_p = function(theta,lambda)
{
    p <- vector(mode='numeric',length=length(z_data))
    for (i in 1:length(z_data))
    {
        p[i] = rbeta(1,a1[i] + lambda[i], a2[i]+1-lambda[i])
    }
    return(p)
}

# Defining the jump function sampler for the latent variables
# We pick an initial value for the pm vector and for theta and then sample lambda bsed on that.
jump_function = function(theta_1,theta_2)
{
    mu = theta_2
    sigma = cov_mat
    library(mvtnorm)
    jump_prob = dmvnorm(theta_1,mean=mu,sigma = sigma,log=TRUE)
    return(jump_prob)
}



# M-H & gibbs sampler sampler for the parameters of interest
# Make use of the marginal posterior of the vector theta
# we take single samples of pm and lambda m each step at a time and sample theta wrt that, finally we shall average the samples for the theta vector 
# this shall be the sampled from the marginal posterior of theta
sampler  = function(num_samp,theta_init,p_init)
{
    theta_samp = matrix(data=NA,nrow=n_samp,ncol=6)
    theta_old = theta_init
    acc=0
    p = p_init
    show(theta_old)
    for (i in 1:n_samp)
    {
        show(i);
        lambda = cond_post_lambda(theta_old,p) # sample lambda
        theta_new = MASS::mvrnorm(n=1, mu=theta_old, Sigma=cov_mat/50) # Proposed theta vector 
        r = (joint_posterior(theta_new,lambda,p) - joint_posterior(theta_old,lambda,p) + jump_function(theta_old,theta_new) - jump_function(theta_new,theta_old))
        show(joint_posterior(theta_new,lambda,p))
        #show(r)
        if (r>0)
        {
            theta_samp[i,] = theta_new
            theta_old = theta_new
            acc = acc+1
        }
        else{
            u  = runif(1)
            if(log(u)<r){
                theta_samp[i,] = theta_new
                theta_old = theta_new
                acc = acc+1
            }
            else{
                theta_samp[i,] = theta_old
                theta_old = theta_old
            }
        }
        p = cond_post_p(theta_old,lambda)
    }
    return(theta_samp)
}


# Here, we aim to draw samples from the poterior of the parameter vector theta

theta_init = c(-0.00044,0.00087,-0.00044,0.00087,-100,-100)
p_init = rep(0.326,length.out=length(z_data))

n_samp = 100
burnin = 50
posterior_samples = sampler(num_samp=n_samp,theta_init = theta_init ,p_init =  p_init)
#nett_draws = colSums(posterior_samples[burnin:n_samp,])/(n_samp-burnin)
plot(posterior_samples[,1],main="Scatterplot for first theta component",type='o')
plot(posterior_samples[,2],main="Scatterplot for second theta component",type='o')
plot(posterior_samples[,3],main="Scatterplot for third theta component",type='o')
plot(posterior_samples[,4],main="Scatterplot for fourth theta component",type='o')
plot(posterior_samples[,5],main="Scatterplot for fifth theta component",type='o')
plot(posterior_samples[,6],main="Scatterplot for sixth theta component",type='o')
