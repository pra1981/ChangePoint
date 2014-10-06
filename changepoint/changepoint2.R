require(MCMCpack)

####################################
########### Initialization #########
####################################

N = 100
K = 6
samples = 1
changepoint.num = 0
collen = K*(changepoint.num+1)

#####################################
########### Data Simulation #########
#####################################

#simulate the alpha parameter
alpha.true = unlist(matrix(0,samples,collen))
for(i in 1:samples){
    x=c()
    for(j in 1:(changepoint.num+1)){
        #x = c(x,sort(runif(K,0.5,10)))
        x = rep.int(1,K)
    }
    alpha.true[i,] = x
}

#simulate the multinomial data
data = unlist(matrix(0,samples,N))
theta.true = unlist(matrix(0,samples,collen))
for (i in 1:samples){
    # Dirichlet distribution (from MCMCpack package)
    theta.true[i,] = rdirichlet(1, alpha.true[i,])
    
    # Multinomial distribution
    counts = rmultinom(samples, size = N, prob = theta.true[i,])
    
    #sampling to generate the data
    vector = rep.int(1:K,counts)
    data[i,] = sample(vector,N,replace=FALSE) 
}


################################################
########### Known parameters from Data #########
################################################

# calculate the n,p
data.p = unlist(matrix(0,samples,collen))
data.n = unlist(matrix(0,samples,collen))
true.ctgr = table(factor(rep(1:K,0),levels = 1:K))
for (i in 1:samples){
    curr.data.ctgr = table(data[i,])
    curr.data.Freq = as.data.frame(curr.data.ctgr)$Freq
    # in case some ctgrs are not sampled
    if(nrow(curr.data.ctgr) != K){
        temp = merge(curr.data.ctgr,true.ctgr,by='Var1',all=TRUE)
        temp$Freq.x[is.na(temp$Freq.x)]=0
        temp$sum<-temp$Freq.x + temp$Freq.y
        #sort the level index
        # -- after merge two ctgrs, the levels index would be confused
        # -- need to sort the ctgrs levels as 1:K
        ctgridx = sort(as.numeric(levels(temp$Var1)),index.return=TRUE)$ix
        curr.data.Freq = temp$sum[ctgridx]
    }
    data.n[i,]=curr.data.Freq
    data.p[i,]=curr.data.Freq/N
}

########################################################
########### optim method to estimate alpha #########
########################################################
# define the function to be optimized
loglikelihood <- function(par,data.n){
    log10(gamma(sum(par))) - log10(gamma(sum(par + data.n))) + 
        sum(log10(gamma(par + data.n)) + log10(gamma(par))) +
        log10(prod(sum(data.n))) + sum(log10(prod(data.n)))
}

gradient <- function(par,data.n){ ##Gradient of loglikelihood
    digamma(sum(par)) - digamma(sum(par + data.n)) + 
        digamma(par + data.n) - digamma(par)
}

# initialize the estimated alpha by random value first
alpha.estimate = unlist(matrix(runif(samples*collen,0.1,1),samples,collen))
#constriant: alpha > 0
alpha.estimate[i,which(alpha.estimate[i,] < 0)] = 
    -log10(-alpha.estimate[i,which(alpha.estimate[i,] < 0)])
theta.estimate = unlist(matrix(0,samples,collen))

Loglikelihood.estimate = unlist(matrix(0,samples,1))
Loglikelihood.true = unlist(matrix(0,samples,1))

for(i in 1:samples){
    Loglikelihood.true[i] = log10(gamma(sum(alpha.true[i,]))) - 
                    log10(gamma(sum(alpha.true[i,]) + N)) + 
                    sum(log10(gamma(alpha.true[i,] + data.n[i,])) + 
                    log10(gamma(alpha.true[i,])))
    res = optim(alpha.estimate[i,],loglikelihood,gradient,data.n = data.n[i,], 
            method="L-BFGS-B", lower=rep(0.000000001,length(alpha.estimate[i,])),
            upper=rep(Inf,length(alpha.estimate[i,])),
            control=list(fnscale=-1,ndeps=1e-10,maxit=1e+9))
    alpha.estimate[i,] = res$par
    #smooth
    alpha.estimate[i,] = alpha.estimate[i,] + 0.01 * min(alpha.estimate[i,which(alpha.estimate[i,]>0)])
    
    Loglikelihood.estimate[i,] = res$value
    theta.estimate[i,] = rdirichlet(1,alpha.estimate[i,])
}



