pred <- function(x,u,
                 xname = "x",
                 uname = "u",
                 lag,basis,coef=NULL,
                 ci.level = 0.95,
                 log = TRUE){

    argvarx <- attr(basis,"argvarx")
    argvaru <- attr(basis,"argvaru")
    l <- as.vector(0:lag)
    basisvarx <- do.call("onebasis",c(list(x=x),attr(basis,"argvarx")))
    basisvaru <- do.call("onebasis",c(list(x=u),attr(basis,"argvaru")))
    basislag <- do.call("onebasis",c(list(x=l),attr(basis,"arglag")))
    mat1 <- basisvarx %x% basisvaru
    mat2 <- mat1 %x% basislag

    nx <- length(x)
    nu <- length(u)
    nl <- lag + 1
    vx <- ncol(basisvarx)
    vu <- ncol(basisvaru)
    vl <- ncol(basislag)

    coef <- as.matrix(model$coefficients[2:(vx*vu*vl+1)])
    vcov <- as.matrix(vcov(model)[2:(vx*vu*vl+1),2:(vx*vu*vl+1)])
    z <- qnorm(1-(1-ci.level)/2)

    ff <- mat2 %*% coef
    ar <- array(ff,
                dim=c(nl,nu,nx),
                dimnames = list(
                    paste0("lag=",round(l,2)),
                    paste0(uname,"=",round(u,2)),
                    paste0(xname,"=",round(x,2))
                    )
                )

    ff.se <- sqrt(pmax(0,rowSums((mat2 %*% vcov) * mat2)))
    ff.se.mat <- array(ff.se,dim=c(nl,nu,nx),
                       dimnames = list(
                         paste0("lag=",round(l,2)),
                         paste0(uname,"=",round(u,2)),
                         paste0(xname,"=",round(x,2))
                       ))

    high <- ar + z * ff.se.mat
    low <- ar - z * ff.se.mat

    if(log){
        ar <- exp(ar)
        ff.se.mat <- exp(ff.se.mat)
        high <- exp(high)
        low <- exp(low)
    }
    res <- list(ar = ar,
                coef = coef,
                se = ff.se.mat,
                h = high,
                l = low,
                x = x,
                u = u,
                xname = xname,
                uname = uname,
                lag = lag
                )

    class(res) <- "mypred"
    return(res)
}
