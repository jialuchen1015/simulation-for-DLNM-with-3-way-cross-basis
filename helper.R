# Functions on the predictor side
flin <- function(x) 0.1*x
fflex <- function(x) {
  coef <- c(1,0.3,0,0,0)
  as.numeric(outer(x,0:4,'^')%*%coef)
}
fdnorm <- function(x) (dnorm(x,1.5,3)+1.5*dnorm(x,5.5,3))


# Functinos on the lag side
wconst <- function(lag) lag-lag+0.20
wdecay <- function(lag) exp(-lag/2)
wpeak1 <- function(lag) 12*dnorm(lag,8,5)
wpeak2 <- function(lag) 15*dnorm(lag,4,10)
wdnorm <- function(lag) 5*(dnorm(lag,4,6)+dnorm(lag,25,4))


# Functions to simulate the exposure-response surface
f1 <- function(x,u,l) 0.1*x*u

fplane <- function(x,u,lag) 1 * (flin(x)-flin(2)) * wconst(lag) *  (flin(u)-flin(2))

ftemp <- function(x,u,lag) 1 * (fflex(x)-fflex(5)) * (fflex(u)-fflex(5)) *
  ifelse(is.na(x),NA,ifelse(x>=5,wdecay(lag),wpeak1(lag)))* 
  ifelse(is.na(u),NA,ifelse(u>=5,wdecay(lag),wpeak1(lag)))

fcomplex <- function(x,u,lag) 1000 * (fdnorm(x)-fdnorm(6)) *  (fdnorm(u)-fdnorm(6))*
  ifelse(is.na(x),NA,ifelse(x>=6,wdnorm(lag),wpeak2(lag)))*
  ifelse(is.na(u),NA,ifelse(u>=6,wdnorm(lag),wpeak2(lag)))


# check dimension compabibility of x and u before applying `fun`
bigF <- function(x,u,l,fun){
  if(length(x) != length(u) ||
     length(x) != length(l)){
    stop("x and u and l should be of same length.")
  }
  return(fun(x,u,l))
}


# Given exposure x and u, and lag, and the function fun,
# generate cumulative effect.
# Note y should be of the same dimension as x and u,
# and that the first L elements of y is NA
gen_cumeff <- function(x,u,L,fun){
  if(length(L)!=1){
    stop("L in gen_y should be integer.")
  }

  nobs <- length(x)
  if(TRUE){
      x <- Lag(x,0:L)
      u <- Lag(u,0:L)
  }

  l.hist <- as.vector(0:L)
  y <- vector(mode='numeric',length=nobs)
  for(i in seq(y)){
    y[i] <- sum(bigF(x[i,],
                     u[i,],
                     l.hist,
                     fun))
  }
  return(y)
}


# Generate a 3-D array, containing the values
# of function `fun` evaluated on a 3D grid
eval_fun_overgrid <- function(x,u,L,fun,xname="x",uname="u"){
  if(length(L)!=1){
    stop("L in gen_y should be integer.")
  }
  nx <- length(x)
  nu <- length(u)
  nl <- L + 1

  oldx <- x
  oldu <- u

  l <- rep(0:L,nu*nx)
  u <- rep(u,each = nl,nx)
  x <- rep(x,each = nl*nu)
  res <- fun(x,u,l)

  ar <- array(res,dim=c(nl,nu,nx),
              dimnames = list(
                paste0("lag",seq(nl)),
                paste0(uname,seq(nu)),
                paste0(xname,seq(nx))
              ))

  se <- array(0,dim=c(nl,nu,nx))
  res <- list(ar = ar,
              se = se,
              h = ar,
              l = ar,
              x = oldx,
              u = oldu,
              xname = xname,
              uname = uname,
              lag = L
  )

  class(res) <- "mypred"
  return(res)
}

strata <- function (x, df = 1, breaks = NULL, ref = 1, intercept = FALSE) 
{
  nx <- names(x)
  x <- as.vector(x)
  range <- range(x, na.rm = TRUE)
  if (!is.null(breaks)) {
    breaks <- sort(unique(breaks))
  }
  else if (df - intercept > 0) 
    breaks <- quantile(x, 1/(df - intercept + 1) * 1:((df - 
                                                         intercept)), na.rm = TRUE)
  df <- length(breaks) + intercept
  xcat <- cut(x, c(range[1] - 1e-04, breaks, range[2] + 1e-04), 
              right = FALSE)
  basis <- matrix(outer(xcat, levels(xcat), "==") + 0, ncol = length(levels(xcat)))
  if (!ref %in% seq(0, ncol(basis))) 
    stop("wrong value in 'ref' argument. See help('strata')")
  if (!intercept && ref == 0) 
    ref <- 1
  if (!is.null(breaks)) {
    if (ref != 0) 
      basis <- basis[, -ref, drop = FALSE]
    if (intercept && ref != 0) 
      basis <- cbind(1, basis)
  }
  dimnames(basis) <- list(nx, seq(ncol(basis)))
  attributes(basis) <- c(attributes(basis), list(df = df, breaks = breaks, 
                                                 ref = ref, intercept = intercept))
  class(basis) <- c("strata", "matrix")
  return(basis)
}

lin <- function (x, intercept = FALSE) 
{
  nx <- names(x)
  x <- as.vector(x)
  basis <- as.matrix(x)
  if (intercept) 
    basis <- cbind(1, basis)
  dimnames(basis) <- list(nx, seq(ncol(basis)))
  attributes(basis) <- c(attributes(basis), list(intercept = intercept))
  class(basis) <- c("lin", "matrix")
  return(basis)
}
