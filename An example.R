library(dlnm) ; library(mgcv) ; library(splines) ; library(tsModel)
source("my.plot.pred.R")
source("my.crossbasis.R")
source("my.pred.R")
source("helper.R")
source("my.cbPen.R")

x<- chicagoNMMAPS$dptp
x <- (x-min(x))/diff(range(x))*10
u <- chicagoNMMAPS$temp
u <- (u-min(u))/diff(range(u))*10
n<-length(x)

L <- 4
testfun <- fcomplex
cumeff <- gen_cumeff(x,u,L,testfun)
exp(mean(cumeff[(L+1):length(cumeff)]))
alpha <- 5
mu <- exp(cumeff + log(alpha))
suppressWarnings(y <- rpois(n,mu))
cor(mu[(L+1):n],y[(L+1):n])
o3.pred = seq.int(min(x),max(x),length.out=23)
tmp.pred = seq.int(min(u),max(u),length.out=29)

trueeff <- eval_fun_overgrid(x = o3.pred,
                             u = tmp.pred,
                             L,
                             testfun)
plot(trueeff)
cb <- my.crossbasis(x = x,u = u,lag=L,
                    argvarx=list(fun='ps',df=9),
                    argvaru=list(fun='ps',df=9),
                    arglag=list(fun='ps'))

cbgamPen <- my.cbPen(cb)
model <- gam(y~cb,family=poisson,paraPen=list(cb=cbgamPen),method='REML')
cp <- pred(x = o3.pred,
           u = tmp.pred,
           xname = "dptp",
           uname = "temp",
           lag = L,
           cb,
           model,log=FALSE)
plot(cp,truth = trueeff)