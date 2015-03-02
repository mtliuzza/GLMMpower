#http://glmm.wikidot.com/power-analysis
#https://stat.ethz.ch/pipermail/r-sig-mixed-models/2009q1/001790.html

#this script models one random intercept only (VSubj)
#in case you want to model random slopes, add them too, e.g. VbRWA=1, VbGrouo=.5 etc...
library(lme4)

sim1 <- function(J, K, b0= 0, bRWA=0, bGroup1=0, bGroup2=0, bGxR1=0,bGxR2=0,Vsubj=1, Verror=1) {
    Subject <- rep( 1:J)#define subjects
    ID <-as.factor(Subject)#define ID as factor
    #verror 1/44 0.02272727 (ovvero .15 sd)
    #serror 2/44 0.04545455 (0.2132007 sd)
    RWA<-rnorm(J)#create J (N clusters/participants) normally distributed RWAs
    RWA<-rep(RWA,1,each=K)#repeat K (observations) times
    
    #adatto varianza a n soggetti
    Vsubj=Vsubj*J# var=np(1-p)
    Verror=Verror*J
    
    Group <- rep(c('T','G','I'), each=1,J*K/3)#idem for groups
    
    # random effects per subject
    b0 <- rnorm(J, b0, sqrt(Vsubj))#random effect for the intercept
    b0 <-rep(b0,1,each=K)#repeat K (observations) times

    # random effects of the slopes I keep it the same for now, 
    #for future random slopes models it will be tailored on the slope
   
    
    # epsilons -- if we use rbinom later, this is not useful anymore
    eps <- rnorm(J*K, 0, sqrt(Verror))#residual error
    
    
    # put it all together
    pTrust <- 1 / ( 1 + exp(-(b0 + bRWA*RWA+bGroup1*(Group=='I') + #builds probaility of trusting, given the coefficients plus random intercept plus residual error 
    bGroup2*(Group=='G') +bGxR1*(Group=='I')*RWA+bGxR2*(Group=='G')*RWA+eps)))# +S.re)))
    Trust<-rbinom(J*K,1,pTrust)#build a binomial random distribution for J*K trials based on the probability given the intercepts, the random error and the residual
    # put into a data frame
    mydata <- data.frame( ID = ID, 
                          RWA = RWA, Group=Group,
                          Trust = as.factor(Trust))
    
    # analyze looking at interaction term with LR test
    fit1 <- glmer( Trust ~ (RWA*Group) + (1|ID), data=mydata,family=binomial)
    fit2 <- glmer( Trust ~ Group + RWA + (1|ID), data=mydata,family=binomial)
    anova(fit2,fit1)[2,8]
}

#try with .4 as exp(.4)=1.5 odd ratio, that is a small effect,i.e. an increase of 10%

nsims=100
pb <- txtProgressBar(max=nsims,style=3) # or tkProgressBar or txtProgressbar

setTxtProgressBar(pb, 0)
out1 <- replicate( nsims, {setTxtProgressBar(pb, getTxtProgressBar(pb)+1);
                         sim1(40, 180, b0= 0, bRWA=0, bGroup1=-0, bGroup2=0, bGxR1=0,bGxR2=.4,Vsubj=0.02777778, Verror=.25)})
hist(out1)#verror 5.7 because ICC=1/(1+9)=.1, which is similar to the "medium" ICC in Gelman...I made it likely for binary data,where var can be no more than .25
mean( out1 < 0.05 )

## Michele's new TERRIBLE proposal

#library(lme4)
library(MuMIn)
##Nakagawa, S, Schielzeth, H. (2012).
##A general and simple method for obtaining R^2
##from Generalized Linear Mixed-effects Models.
##Methods in Ecology and Evolution: (online)
##doi:10.1111/j.2041-210x.2012.00261.x
## with the r.squaredGLMM we can obtain 2 R^2 indexes:
## the marginal R^2, variance explained by the fixed effects only
## and the conditional R^2 that assesses the variance explained from both random and fixed effects

set.seed(5)
nsim <- 10
## just 10 simulations to try it... it's very slow...

fitsim <- function(i) {
## function that compare the null and alternative hypothesis and gives us the r2c and r2m
  c(anova(fit[[i]]$H1,fit[[i]]$H0)[2,8],r.squaredGLMM(fit[[i]]$H1))
}

## initial values
betas = c(0,.1,.1,.1)## intercept, W, C, V
thetas = c(.3,.7)## ids, subj

start = Sys.time()
fitAll=list()
for(nSubj in c(10,20,30,40,50)){
  numobs = 93*nSubj
  
  pb = txtProgressBar(min=0,max=nsim, style=3)
  fit = list()## within this list there will be all the glmm
  for(j in 1:nsim){
  ## this cycle needs optimization!
  ## because I don't know how to produce a glmm given a R^2 value,
  ## I'm simulating data and fitting glmms until the R^2 is not the desidered value
  
  ## maybe parallelize this part is helpful
    repeat{
      expdata=data.frame(
        W = rnorm(numobs),
        C = rnorm(numobs),
        V = rnorm(numobs),
        subjs = as.factor(rep(1:nSubj,length.out=numobs)),
        ids = as.factor(rep(1:93,length.out=numobs))
      )
      
#       "C"=rnorm(n=1,mean=1/nSubj,sd=1/nSubj),
#       "W"=rnorm(n=1,mean=1/nSubj,sd=1/nSubj),
#       "V"=rnorm(n=1,mean=1/nSubj,sd=1/nSubj))
      
      ## probably we need better betas and thetas value
      ## they need to be scaled, otherwise if they are fixed
      ## with larger samples the R^2 will be too large to fit with the final condition
      
      ## I'm starting with initial values (row 78 and 79), each iteraction
      ## I'm saving the new betas and thetas and each iteraction I'm stricting the SD
      ## however I'm guaranteeing random variations
      beta <- c("(Intercept)"=rnorm(n=1,mean=betas[1],sd=1/(j*10)),
                "C"=rnorm(n=1,mean=betas[2],sd=1/(j*10)),
                "W"=rnorm(n=1,mean=betas[3],sd=1/(j*10)),
                "V"=rnorm(n=1,mean=betas[4],sd=1/(j*10)))
      theta <- c("ids.(Intercept)"=rnorm(n=1,mean=thetas[1],sd=1/(j*10)),
                 "subjs.(Intercept)"=rnorm(n=1,mean=thetas[2],sd=1/(j*10)))
                 
      
      ## I'm using simulate because probably it takes into account all the possible sources of errors
      ## variations and interclass correlations
      expdata$itarum <- unlist(simulate(~C+W+V+ (1 | subjs) + (1|ids),
                     nsim = 1, family = binomial, 
                     #weights = rep(25, nrow(expdat)),
                     newdata = expdata, newparams =
                       list(theta = theta, beta = beta)))
      
      fit[[j]] <- list(
        H1=glmer(itarum~C+W+V+(1 | subjs) + (1|ids), data=expdata,family=binomial),
        H0=glmer(itarum~1+(1 | subjs) + (1|ids),data=expdata,family=binomial)
      )
      r2c=r.squaredGLMM(fit[[j]]$H1)[2]
      
    ## if r^2 is medium let's go...
    ## else let's coninue the simulation
    ## condition too strict?
      if(r2c>=.3^2&&r2c<=.31^2){
      ## saving good betas and thetas
        betas=fit[[1]]$H1@beta
        thetas=fit[[1]]$H1@theta     
        
        setTxtProgressBar(pb, j)
        break
      }
    }
  }
  close(pb)
  
  fitAll[[nSubj]] = lapply(seq(nsim), function(i) fitsim(i)) #<- unlist(clusterApply(cl, seq(nsim), function(i) fitsim(i)))
  
  ## partial ad-interim output
  cat(nSubj,"\n")
  print(c(power=mean(do.call("rbind",fitAll[[nSubj]])[,1]<0.05),
        colMeans(do.call("rbind",fitAll[[nSubj]])[,2:3])))
}
## final table
rbind("10"=c(power=mean(do.call("rbind",fitAll[[10]])[,1]<0.05),colMeans(do.call("rbind",fitAll[[10]])[,2:3])),
      "20"=c(power=mean(do.call("rbind",fitAll[[20]])[,1]<0.05),colMeans(do.call("rbind",fitAll[[20]])[,2:3])),
      "30"=c(power=mean(do.call("rbind",fitAll[[30]])[,1]<0.05),colMeans(do.call("rbind",fitAll[[30]])[,2:3])),
      "40"=c(power=mean(do.call("rbind",fitAll[[40]])[,1]<0.05),colMeans(do.call("rbind",fitAll[[40]])[,2:3])),
      "50"=c(power=mean(do.call("rbind",fitAll[[50]])[,1]<0.05),colMeans(do.call("rbind",fitAll[[50]])[,2:3])))
Sys.Time()-start## blasfemies
