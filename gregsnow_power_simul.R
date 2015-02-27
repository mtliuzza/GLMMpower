#http://glmm.wikidot.com/power-analysis
#https://stat.ethz.ch/pipermail/r-sig-mixed-models/2009q1/001790.html
library(lme4)

sim1 <- function(J, K, b0= 0, bRWA=0, bGroup1=0, bGroup2=0, bGxR1=0,bGxR2=0,Vsubj=1, Verror=1) {
    Subject <- rep( 1:J)#define subjects
    ID <-as.factor(Subject)#
    
    RWA<-rnorm(J)
    RWA<-rep(RWA,1,each=K)
    
    Group <- rep(c('T','G','I'), each=1,J*K/3)
    
    # random effects per subject
    S.re <- rnorm(J, b0, sqrt(Vsubj))
    S.re<-rep(S.re,1,each=K)
    
    # epsilons
    eps <- rnorm(J*K, 0, sqrt(Verror))
    
    
    # put it all together
    pTrust <- 1 / ( 1 + exp(-(b0 + bRWA*RWA+bGroup1*(Group=='I') + bGroup2*(Group=='G') +bGxR1*(Group=='I')*RWA+bGxR2*(Group=='G')*RWA+S.re + eps)))
    Trust<-ifelse(pTrust < .5,0,1)
    # put into a data frame
    mydata <- data.frame( ID = ID, 
                          RWA = RWA, Group=Group,
                          Trust = as.factor(Trust))
    
    # analyze looking at interaction term with LR test
    fit1 <- glmer( Trust ~ (RWA*Group) + (1|ID), data=mydata,family=binomial)
    fit2 <- glmer( Trust ~ Group + RWA + (1|ID), data=mydata,family=binomial)
    anova(fit2,fit1)[2,8]
}

#come sono arrivato al .13.
#http://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/effectSize
#partial eta^2 medium = .13
#Cohen's f = Square Root of eta-squared / (1-eta-squared)
#Use the square of a Pearson correlation for effect sizes for partial η2 (R-squared in a multiple regression) 
#giving 0.01 (small), 0.09 (medium) and 0.25 (large) which are intuitively larger values than eta-squared.
#Further to this Cohen, Cohen, West and Aiken (2003) on page 95 of Applied Multiple Regression/Correlation Analysis for the behavioral Sciences third edition for looking at semi-partial effects of single predictors in a regression rather than an overall model R-squared 
#ie looking at sqrt(change in R-squared) from a model with and without the regressor and using the Pearson correlations as a rule of thumb for effect sizes.
#cohen's f =sqrt(.14^2 / (1-.14^2))=.14 ovvero un efetto piccolo

nsims=100
pb <- txtProgressBar(max=nsims,style=3) # or tkProgressBar or txtProgressbar

setTxtProgressBar(pb, 0)
out1 <- replicate( nsims, {setTxtProgressBar(pb, getTxtProgressBar(pb)+1);
                         sim1(15, 180, b0= 0, bRWA=0, bGroup1=-0, bGroup2=0, bGxR1=0,bGxR2=-.4,Vsubj=1, Verror=1)})
hist(out1)
mean( out1 < 0.05 )
