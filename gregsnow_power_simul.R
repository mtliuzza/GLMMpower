#http://glmm.wikidot.com/power-analysis
#https://stat.ethz.ch/pipermail/r-sig-mixed-models/2009q1/001790.html

#this script models one random intercept only (VSubj)
#in case you want to model random slopes, add them too, e.g. VbRWA=1, VbGrouo=.5 etc...
library(lme4)

sim1 <- function(J, K, b0= 0, bRWA=0, bGroup1=0, bGroup2=0, bGxR1=0,bGxR2=0,Vsubj=1, Verror=1) {
    Subject <- rep( 1:J)#define subjects
    ID <-as.factor(Subject)#define ID as factor
    
    RWA<-rnorm(J)#create J (N clusters/participants) normally distributed RWAs
    RWA<-rep(RWA,1,each=K)#repeat K (observations) times
    
    Group <- rep(c('T','G','I'), each=1,J*K/3)#idem for groups
    
    # random effects per subject
    S.re <- rnorm(J, b0, sqrt(Vsubj))#random effect for the intercept
    S.re<-rep(S.re,1,each=K)#repeat K (observations) times
    
    # epsilons -- if we use rbinom later, this is not useful anymore
    #eps <- rnorm(J*K, 0, sqrt(Verror))#residual error
    
    
    # put it all together
    pTrust <- 1 / ( 1 + exp(-(b0 + bRWA*RWA+bGroup1*(Group=='I') + #builds probaility of trusting, given the coefficients plus random intercept plus residual error 
    bGroup2*(Group=='G') +bGxR1*(Group=='I')*RWA+bGxR2*(Group=='G')*RWA+S.re))# + eps)))
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
                         sim1(15, 180, b0= 0, bRWA=0, bGroup1=-0, bGroup2=0, bGxR1=0,bGxR2=-.4,Vsubj=1, Verror=1)})
hist(out1)
mean( out1 < 0.05 )
