require(MCMCpack)

####################################
########### Initialization #########
####################################

K = 6               # category number
rows = 1            # samples
changepoint.num = 0 # changepoint number per row 
#(right now, all the rows has the same changepoint number and position)
paramsnum = K*(changepoint.num+1) 
SampleN = 100       # how many points per section have
rowN = SampleN*(changepoint.num+1)

#####################################
########### Data Simulation #########
#####################################

#simulate the alpha parameter
alpha.true = unlist(matrix(0,rows,paramsnum))
for(i in 1:rows){
    x=c()
    for(j in 1:(changepoint.num+1)){ 
        y=runif(K,0,0.01*10^j)        # control the change rate
        y[which(y<0.003*10^j)]=0
        x = c(x,sort(y*c(1:K)*j))
        #x= rep.int(1:K,changepoint.num+1)
    }
    alpha.true[i,] = x
    
}

#simulate the multinomial data
data = unlist(matrix(0,rows,rowN))
theta.true = unlist(matrix(0,rows,paramsnum))
data.n = unlist(matrix(0,rows,paramsnum))

for(i in 1:rows){
    
    # Dirichlet distribution (from MCMCpack package)
    for(j in 1:(changepoint.num+1)){
        start = K*(j-1)+1 #beginning of jth changepoint section
        end = K*j           #end of the jth changepoint section
        # take the avg of the dirichlet simulation
        # given same alpha,
        # if sample only once, the result of theta has large diff,
        # but if sample many times, the mean of theta is similar
        theta.true[i,start:end] = colMeans(rdirichlet(100, alpha.true[i,start:end]))
    }
    
    
    # Multinomial distribution
    data.n[i,] = t(rmultinom(1, size = rowN, prob = theta.true[i,]))
    
    seeds = rep.int(rep.int(1:K,changepoint.num+1),sapply(data.n[i,],rep))
    
    #sampling to generate the data
    for(j in 1:(changepoint.num+1)){
        start = SampleN*(j-1)+1 # beginning of jth changepoint sample
        end = SampleN*j
        data[i,start:end] = sample(seeds[start:end],SampleN,replace=FALSE) 
    }
}   

# #################################################################
# ################# Plot the Data(1 change points) ################
# #################################################################
# layout(matrix(c(1,1,2,3,4,5,6,7), 4, 2, byrow = TRUE))
# plot(1:rowN,data[i,],main="data",xlab = "", ylab = "")
# hist(data[i,1:SampleN],main="counts",ylim=c(0,max(data.n[i,])),xlim=c(1,K),xlab = "", ylab = "")
# hist(data[i,(SampleN+1):rowN],main="counts",ylim=c(0,max(data.n[i,])),xlim=c(1,K),xlab = "", ylab = "")
# plot(1:K,alpha.true[i,1:K],main="alpha",type="b",ylim=c(0,max(alpha.true[i,])),xlab = "", ylab = "")
# plot(1:K,alpha.true[i,(K+1):paramsnum],main="alpha",type="b",ylim=c(0,max(alpha.true[i,])),xlab = "", ylab = "")
# plot(1:K,theta.true[i,1:K],main="theta",type="b",ylim=c(0,1),xlab = "", ylab = "")
# plot(1:K,theta.true[i,(K+1):paramsnum],main="theta",type="b",ylim=c(0,1),xlab = "", ylab = "")


########################################################
########### optim method to estimate alpha #########
########################################################
# define the function to be optimized
loglikelihood <- function(par,data.n){
    idx = which(data.n>0)
    lgamma(sum(par)) - lgamma(sum(par + data.n)) + 
        sum(lgamma(par[idx] + data.n[idx]) + lgamma(par[idx])) +
        log10(prod(sum(data.n))) + sum(log10(prod(data.n[idx])))
}

gradient <- function(par,data.n){ ##Gradient of loglikelihood
    idx = which(data.n>0)
    gdt = unlist(matrix(0,1,length(par)))
    gdt[idx] = digamma(sum(par)) - digamma(sum(par + data.n)) + 
        digamma(par[idx] + data.n[idx]) - digamma(par[idx])
    gdt
}

# initialize the estimated alpha by random value first
alpha.estimate = unlist(matrix(runif(rows*paramsnum,0.1,1),rows,paramsnum))
theta.estimate = unlist(matrix(0,rows,paramsnum))

Loglikelihood.estimate = unlist(matrix(0,rows,1))
Loglikelihood.true = unlist(matrix(0,rows,1))

for(i in 1:rows){
    Loglikelihood.true[i] =  loglikelihood(alpha.true[i,],data.n[i,])
    
    alpha.estimate[i,which(data.n[i,] == 0)] = 0
    
    res = optim(alpha.estimate[i,],loglikelihood,gradient,data.n = data.n[i,], 
            method="L-BFGS-B", lower=rep(0,paramsnum), upper=rep(Inf,paramsnum),
            control=list(fnscale=-1,ndeps=1e-10,maxit=1e+9))
    alpha.estimate[i,] = res$par
    #smooth
    #alpha.estimate[i,] = alpha.estimate[i,] + 0.01 * min(alpha.estimate[i,which(alpha.estimate[i,]>0)])
    
    Loglikelihood.estimate[i,] = res$value
    theta.estimate[i,] = colMeans(rdirichlet(100, alpha.estimate[i,]))
    
}