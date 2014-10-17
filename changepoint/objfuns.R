# define the function to be optimized

## loglikelihood of DiriMult
loglikelihood <- function(par,data.n){ 
    idx = which(data.n>0)
    lgamma(sum(par)) - lgamma(sum(par + data.n)) + 
        sum(lgamma(par[idx] + data.n[idx]) + lgamma(par[idx])) +
        log10(prod(sum(data.n))) + sum(log10(prod(data.n[idx])))
}

##Gradient of loglikelihood
gradient <- function(par,data.n){
    idx = which(data.n>0)
    gdt = unlist(matrix(0,1,length(par)))
    gdt[idx] = digamma(sum(par)) - digamma(sum(par + data.n)) + 
        digamma(par[idx] + data.n[idx]) - digamma(par[idx])
    gdt
}


## Object function of argmax
argmax <- function(par,data,...){ 
    data.n = getDataFreq(data,length(data))
    N=sum(data.n$alpha)
    
    ### calculate all the logPr for each lambda ###
    logPr = matrix(0,1,N)
    
    ## lamdba = 0 : no change point
    logP0 = log10(1- par[1])
    logP.alpha = loglikelihood(par[2:K+1],data.n$alpha)
    logPr[1] = logP0 + logP.alpha
    
    ## lambda > 0 : change point position
    for(lambda in 1:N-1){
        data.n = getDataFreq(data,lambda)
        logP0 = log10(par[1]/(N-1))
        logP.alpha = loglikelihood(par[2:(K+1)],data.n$alpha)
        logP.beta = loglikelihood(par[(K+2):(2*K+1)],data.n$beta)
        logPr[lambda+1] = logP0 + logP.alpha + logP.beta
    }
    logPmax = max(logPr)
    
    res = logPmax + log10(sum(exp(logPr - logPmax)))
    res
}

getDataFreq <- function(data,lambda,...){
    data.n = vector("list",2)
    names(data.n) = c("alpha","beta")
    
    true.ctgr = table(factor(rep(1:K,0),levels = 1:K)) 
    
    #calculate data.n for alpha
    alpha.ctgr = table(data[1:lambda])
    alpha.Freq = as.data.frame(alpha.ctgr)$Freq
    if(nrow(alpha.ctgr) != K){
        temp = merge(alpha.ctgr,true.ctgr, by='Var1',all=TRUE)
        temp$Freq.x[is.na(temp$Freq.x)]=0
        temp$sum<-temp$Freq.x + temp$Freq.y
        #sort the level index
        # -- after merge two ctgrs, the levels index would be mixed
        # -- need to sort the ctgrs levels as 1:K
        ctgridx = sort(as.numeric(levels(temp$Var1)),index.return=TRUE)$ix
        alpha.Freq = temp$sum[ctgridx]
    }
    
    data.n$alpha = alpha.Freq
    
    #calculate data.n for beta
    if(lambda+1 <= length(data)){
        beta.ctgr = table(data[lambda+1:length(data)])
        beta.Freq = as.data.frame(beta.ctgr)$Freq
        if(nrow(beta.ctgr) != K){
            temp = merge(beta.ctgr,true.ctgr, by='Var1',all=TRUE)
            temp$Freq.x[is.na(temp$Freq.x)]=0
            temp$sum<-temp$Freq.x + temp$Freq.y
            #sort the level index
            # -- after merge two ctgrs, the levels index would be mixed
            # -- need to sort the ctgrs levels as 1:K
            ctgridx = sort(as.numeric(levels(temp$Var1)),index.return=TRUE)$ix
            beta.Freq = temp$sum[ctgridx]
        }
        data.n$beta = beta.Freq
    }
    data.n
}
