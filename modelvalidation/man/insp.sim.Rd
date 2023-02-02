\name{insp.sim}
\alias{insp.sim}
\title{Model simulation}
\usage{
insp.sim(model, nsim = 19, data = extract.data(model))
}
\arguments{
\item{model}{A fitted model}

\item{nsim}{Number of simulations}

\item{data}{Original model data}
}
\value{
Simulated models
}
\description{
Simulate new models based on `lm()`, `glm()`, `lmer()`, gls() or lme() models. See also `simulate()`, `simulate_gls()` and `simulate_lme()` for further details.
}
\examples{
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
lm.model<-lm(Y ~ pred1 + pred2 + pred3)
insp.sim(model=lm.model,nsim=19,data=extract.data(lm.model))
}
