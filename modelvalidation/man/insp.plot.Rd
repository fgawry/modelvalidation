\name{insp.plot}
\alias{insp.plot}
\title{Linear Model Validation Plots}
\usage{
insp.plot(model, data = extract.data(model), bin.size = 10)
}
\arguments{
\item{model}{A fitted model. Accepts lm(), glm(), gls(), lme(), and lmer() models.}

\item{data}{Model data.}

\item{bin.size}{For glm binomial family only. Binned residuals.}
}
\value{
Model specific residuals.
}
\description{
Linear Model Validation Plots for `lm()`, `glm()`, `lmer()`,
`lme()` and `gls()` models.
}
\references{
Zuur et al. 2009
}
\author{
Felipe M. Gawryszewski
}
\examples{
#General Linear Model
n<-100
pred1<-rnorm(n=n, mean=0, sd=1)
pred2<-rnorm(n=n, mean=1, sd=2)
pred3<-gl(n=5,k=n/5)
ERROR<-rnorm(n=n, mean=0, sd=2)
a<-1
b<-2
c<-rep(runif(n=length(levels(pred3)),min = -2, max = 2),
      each=n/length(levels(pred3)))
Y <- a*pred1 + b*pred2 + c + ERROR
Y.true <- a*pred1 + b*pred2 + c
plot(Y~Y.true)

lm.model<-lm(Y ~ pred1 + pred2 + pred3)
insp.plot(model=lm.model, data=extract.data(lm.model))

#Generalised Linear Model - Poisson
#Zuur (2009) page 230; However, it should be noted that we are not
#looking for normality from the Pearson or
#deviance residuals. It is all about lack of
#fit and looking for patterns in the deviance
#or Pearson residuals.
n<-100
pred1<-1:100
a<-0.01
b<-0.03
Y.true <- exp(a + b*pred1)
Y<-Y.true
for (i in 1:length(Y.true)) {
   Y[[i]] <- rpois(n=1,lambda=Y.true[[i]])
}
plot(Y~Y.true)

glm<-glm(Y ~ pred1, family=poisson)
summary(glm)
insp.plot(data=data.frame(Y,pred1), model=glm)


#Generalised Linear Model - Binomial
n<-1000
pred1<-1:1000
a<--3.8
b<-0.03
Y.true <- exp(a + b*pred1)/(1+exp(a + b*pred1))
Y<-Y.true
for (i in 1:length(Y.true)) {
 Y[[i]] <- rbinom(n=1,size=1,prob=Y.true[[i]])
}
plot(Y~Y.true)
glm<-glm(Y ~ pred1, family=binomial)
insp.plot(data=data.frame(Y,pred1), model=glm, bin.size = 100)

\dontrun{
##Generalised least-squares
library(nlme)
n<-100
pred1<-rnorm(n=n, mean=0, sd=1)
pred2<-rnorm(n=n, mean=1, sd=2)
pred3<-gl(n=2,k=n/2)
ERROR1<-rnorm(n=n/2, mean=0, sd=1)
ERROR2<-rnorm(n=n/2, mean=0, sd=2)
ERROR<-c(ERROR1,ERROR2)
a<-1
b<-2
c<-rep(runif(n=length(levels(pred3)),min = 1, max = 5),
      each=n/length(levels(pred3)))
Y <- a*pred1 + b*pred2 + c + ERROR
Y.true <- a*pred1 + b*pred2 + c
plot(Y~Y.true)

gls1<-gls(Y ~ pred1 + pred2 + pred3, weights=varIdent(form=~1|pred3))
insp.plot(data=data.frame(Y,pred1,pred2,pred3), model=gls1)

##General Mixed Model. lme() and lmer().
#random intercept
n<-100
pred1<-seq(-10,10,length.out = n)
pred3<-gl(n=10,k=n/10)
a<-10
b<-0
random<-rnorm(n=length(levels(pred3)),mean=0,sd=1)
random<-rep(random, each=n/length(levels(pred3)))

ERROR<-rnorm(n=n, mean=0, sd=1)

Y <- b + a*pred1 + random + ERROR

#ACF shows a clear correlation of residuals
lm<-lm(Y ~ pred1)
insp.plot(model=lm)

lm<-lm(Y ~ pred1 + pred3)
insp.plot(model=lm)

library(nlme)
lme<-lme(Y ~ pred1, random = ~1|pred3)
insp.plot(data=data.frame(Y,pred1,pred3), model=lme)

library(lme4)
lmer<-lmer(Y ~ pred1 + (1|pred3))
insp.plot(data=data.frame(Y,pred1,pred3), model=lmer)

#random intercept and slope model
n.total<-200
n.groups<-10
size.groups<-n.total/n.groups
pred2<-gl(n=n.groups,k=size.groups)
alpha<-1
beta<-2
#each beach has its own alpha and beta
alpha.beach<-rnorm(n=n.groups,mean=alpha,sd=1)
beta.beach<-rnorm(n=n.groups,mean=beta,sd=1)
#one model for each beach
X<-NULL
Y<-NULL
for (i in 1:n.groups) {
 X.temp<-rnorm(n=size.groups)
 Y.temp <- alpha.beach[[i]]+beta.beach[[i]]*X.temp + rnorm(n=size.groups,mean=0,sd=1)
 X<-c(X.temp,X)
 Y<-c(Y.temp,Y)
}

lm<-lm(Y ~ X)
insp.plot(model=lm)

lm<-lm(Y ~ X*pred2)
insp.plot(model=lm)

library(nlme)
lme1<-lme(Y ~ X, random = ~1|pred2)
insp.plot(data=data.frame(Y,X,pred2), model=lme1)

lme2<-lme(Y ~ X, random = ~1+X|pred2)
insp.plot(data=data.frame(Y,X,pred2), model=lme2)

AIC(lme1,lme2)

library(lme4)
lmer<-lmer(Y ~ X + (1+X|pred2))
insp.plot(data=data.frame(Y,X,pred2), model=lmer)

}
}
