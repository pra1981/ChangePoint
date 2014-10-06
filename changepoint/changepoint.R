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
########### Grediant Descent to estimate alpha #########
########################################################

# initialize the estimated alpha by random value first
alpha.estimate = unlist(matrix(runif(samples*collen,0.1,1),samples,collen))
#constriant: alpha > 0
alpha.estimate[i,which(alpha.estimate[i,] < 0)] = 
    -log10(-alpha.estimate[i,which(alpha.estimate[i,] < 0)])

theta.estimate = unlist(matrix(0,samples,collen))

Loglikelihood.estimate = unlist(matrix(0,samples,1))
Loglikelihood.true = unlist(matrix(0,samples,1))

for(i in 1:samples){
    Loglikelihood.true[i,] = - log10(gamma(sum(alpha.true[i,]))) + log10(gamma(sum(alpha.true[i,]) + N)) - 
        sum(log10(gamma(alpha.true[i,] + data.n[i,])) + log10(gamma(alpha.true[i,]))) -
        log10(prod(sum(data.n))) - sum(log10(prod(data.n)))
    
    prelogL = Inf
    # initialize the gradient value by random value first, it will be updated later
    gradient = runif(K) 
    
    numIterations = 1000000000
    
    for(j in 1:numIterations){

        #current log likelihood
        #currlogL = N * (log10(gamma(sum(alpha.estimate[i,]))) 
         #               - sum(log10(gamma(alpha.estimate[i,]))) 
        #                + sum((alpha.estimate[i,] - 1) * log.theta_t.hat))
        currlogL = - log10(gamma(sum(alpha.estimate[i,]))) + log10(gamma(sum(alpha.estimate[i,]) + N)) - 
                    sum(log10(gamma(alpha.estimate[i,] + data.n[i,])) + log10(gamma(alpha.estimate[i,])))
        
        if(j %% 10000 == 0){
            print(c(currlogL,prelogL,Loglikelihood.true[i,],prelogL - currlogL))
        }
        
        #check whether it's near the optimum
        if(prelogL - currlogL< 0.00000000001 && j > 1)
            break
        
        prelogL = currlogL 
        
        #update the gradient
        #gradient = N * (digamma(sum(alpha.estimate[i,])) 
         #               - digamma(alpha.estimate[i,]) + log.theta_t.hat)
        gradient = - digamma(sum(alpha.estimate[i,])) + digamma(sum(alpha.estimate[i,]) + N) - 
                    digamma(alpha.estimate[i,] + data.n[i,]) + digamma(alpha.estimate[i,])
        
        #update alpha.estimate
        #alpha.estimate[i,] = alpha.estimate[i,]*((digamma(alpha.estimate[i,] + data.n[i,]) - 
        #                    digamma(alpha.estimate[i,]))/(digamma(sum(alpha.estimate[i,]) + 
        #                                            N) - digamma(sum(alpha.estimate[i,]))))
        alpha.estimate[i,] = alpha.estimate[i,] + 0.000000005 * gradient
        #constriant: alpha > 0
        while(length(which(alpha.estimate[i,] < 0)) > 0){
            alpha.estimate[i,which(alpha.estimate[i,] < 0)] = 
                               -log10(-alpha.estimate[i,which(alpha.estimate[i,] < 0)])
        }
        
    }
    Loglikelihood.estimate[i,] = currlogL
    #smooth
    alpha.estimate[i,] = alpha.estimate[i,] + 0.01 * min(alpha.estimate[i,which(alpha.estimate[i,]>0)])
    
    theta.estimate[i,] = rdirichlet(1,alpha.estimate[i,])
}


