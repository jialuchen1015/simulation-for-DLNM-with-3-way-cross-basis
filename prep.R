library(dlnm) ; library(mgcv) ; library(splines) ; library(tsModel)
source("my.plot.pred.R")
source("my.crossbasis.R")
source("my.pred.R")
source("helper.R")
source("my.cbPen.R")

L <- 4
x<- chicagoNMMAPS$dptp
x <- (x-min(x))/diff(range(x))*10
x<-x[1:2000]
u <- chicagoNMMAPS$temp
u <- (u-min(u))/diff(range(u))*10
u<-u[1:2000]
n<-length(x)
qn <- qnorm(0.975)
kmax <- 1
degree <- 1

argvarxlist1 <- c(list(
  list(fun='lin'),list(fun='bs',degree=degree)),
  lapply(seq(kmax), function(k) {
    knots <- equalknots(x,nk=k,fun='bs',degree=degree)
    list(fun='bs',knots=knots,degree=degree)
  })
)
argvarulist1 <- c(list(
  list(fun='lin'),list(fun='bs',degree=degree)),
  lapply(seq(kmax), function(k) {
    knots <- equalknots(u,nk=k,fun='bs',degree=degree)
    list(fun='bs',knots=knots,degree=degree)
  })
)
arglaglist1 <- c(list(
  list(fun='strata',int=T),list(fun='bs',degree=degree,int=T)),
  lapply(seq(kmax), function(k) {
    knots <- equalknots(0:L,nk=k,fun='bs',degree=degree)
    list(fun='bs',knots=knots,degree=degree,int=T)
  })
)


arglaglist1 <- rep(arglaglist1,(kmax+2)^2)
argvarulist1 <- rep(argvarulist1,each=kmax+2,kmax+2)
argvarxlist1 <- rep(argvarxlist1,each=(kmax+2)^2)

# NUMBER OF ITERATIONS (SET TO 1000 TO REPLICATE THE RESULTS)
nsim <- 1

# NUMBER OF SAMPLES USED AS EXAMPLES OF INDIVIDUAL ESTIMATES
nsample <- min(nsim,25)

################################################################################
# COMBINATIONS OF FUNCTIONS USED TO SIMULATE DATA
combsim <- c("fplane","ftemp","fcomplex")
names(combsim) <- c("Plane","Temperature","Complex")

# CREATE THE OBJECT TO STORE RESULTS

# MODELS
models <- c('GAM','GLM-AIC','GAM-AIC','GAM-cr','GAM-ps21','GAM-last','GAM-quad',
            'GAM-exp')
modlab <-list(expression(plain(GAM)),expression(plain(GLM-AIC)),
              expression(plain(GAM-AIC)),expression(plain(GAM-CR)),
              expression(plain(GAM-PS)[plain("2,1")]),expression(plain(GAM-ADD)[plain(LAST)]),
              expression(plain(GAM-ALT)[plain(QUAD)]),expression(plain(GAM-ALT)[plain(EXP)]))

# LISTS FOR STORING THE AVERAGE PREDICTION, BIAS, COVERAGE RMSE AND SAMPLES
# FOR THE WHOLE SURFACE
predmatrix <- lapply(models,function(i) lapply(seq(combsim), function(j) 0))
sample <- lapply(models,function(i) vector('list',length(combsim)))
names(predmatrix) <- names(sample) <- models
for(i in seq(predmatrix)) names(predmatrix[[i]]) <- names(sample[[i]]) <- names(combsim)
bias <- cov <- rmse <- predmatrix

# LISTS OF MATRICES FOR STORING EDF, TIME, AND CONVERGENCE
nonconv <- lapply(models,function(x) 
  matrix(NA,nsim,length(combsim),dimnames=list(seq(nsim),names(predmatrix[[1]]))))
names(nonconv) <- models
edf <- timemodel <- nonconv
base <- c(15,3,5)
o3.pred = seq.int(min(x),max(x),length.out=23)
tmp.pred = seq.int(min(u),max(u),length.out=29)
#################################################################