my.crossbasis <- function (x, u, lag,
                           argvarx = list(),
                           argvaru = list(),
                           arglag = list(),
                           group = NULL,
                           ...)
{

    lag <- mklag(lag)
    require(dlnm)

    # variables
    x <- as.matrix(x)
    u <- as.matrix(u)
    dimx <- dim(x)
    dimu <- dim(u)
    if (!dimx[2] %in% c(1L, diff(lag) + 1L)){
        stop("NCOL(x) must be equal to 1 (if x is a time series vector), ",
             "otherwise to the lag period (for x as a matrix of lagged occurrences)")}
    basisvarx <- do.call("onebasis", modifyList(argvarx, list(x = as.numeric(x))))
    basisvaru <- do.call("onebasis", modifyList(argvaru, list(x = as.numeric(u))))
    basisvarxu <- cross.row.prod(basisvarx,basisvaru)


    # lag
    if (length(arglag) == 0L || diff(lag) == 0L){
        arglag <- list(fun = "strata", df = 1, intercept = TRUE)
    }
    if ((is.null(arglag$fun) ||
         "intercept" %in% names(formals(arglag$fun))) &&
        sum(pmatch(names(arglag), "intercept", nomatch = 0)) == 0){
            arglag$intercept <- TRUE
    }
    arglag$cen <- NULL
    basislag <- do.call("onebasis", modifyList(arglag, list(x = seqlag(lag))))


    # cross basis
    crossbasis <- matrix(0, nrow = dimx[1],
                         ncol = ncol(basisvarxu) * ncol(basislag))
    for (v in seq(length = ncol(basisvarxu))) {

        mat <- as.matrix(Lag(basisvarxu[, v], seqlag(lag)))

        for (l in seq(length = ncol(basislag))) {
            crossbasis[, ncol(basislag) * (v - 1) + l] <- mat %*%
                (basislag[, l])
        }
    }


    # names of the cross basis matrix
    cxcu <- paste0("bx", rep(seq(ncol(basisvarx)), each = ncol(basisvaru)),
                   ".bu", rep(seq(ncol(basisvaru)), ncol(basisvarx)))
    cn <- paste0(rep(cxcu, each = ncol(basislag)),
                 ".bl", rep(seq(ncol(basislag)), ncol(basisvarx)))
    dimnames(crossbasis) <- list(rownames(x), cn)


    # info of the onebasis
    ind <- match(names(formals(attributes(basisvarx)$fun)), names(attributes(basisvarx)),
                 nomatch = 0)
    argvarx <- c(attributes(basisvarx)["fun"], attributes(basisvarx)[ind])

    ind <- match(names(formals(attributes(basisvaru)$fun)), names(attributes(basisvaru)),
                 nomatch = 0)
    argvaru <- c(attributes(basisvaru)["fun"], attributes(basisvaru)[ind])

    ind <- match(names(formals(attributes(basislag)$fun)), names(attributes(basislag)),
                 nomatch = 0)
    arglag <- c(attributes(basislag)["fun"], attributes(basislag)[ind])

    argvarx$cen <- attributes(basisvarx)$cen
    argvaru$cen <- attributes(basisvaru)$cen


    attributes(crossbasis) <- c(attributes(crossbasis),
                                list(df = c(ncol(basisvarx),ncol(basisvaru),ncol(basislag)),
                                     xrange = range(x, na.rm = T),
                                     urange = range(u, na.rm = T),
                                     lag = lag,
                                     argvarx = argvarx,
                                     argvaru = argvaru,
                                     arglag = arglag)
                                )
    class(crossbasis) <- c("crossbasis", "matrix")
    return(crossbasis)
}

mklag <- function (lag)
{
    if (any(!is.numeric(lag)) || length(lag) > 2)
        stop("'lag' must a integer vector or length 2 or 1")
    if (length(lag) == 1L)
        lag <- if (lag < 0L)
            c(lag, 0L)
    else c(0L, lag)
    if (diff(lag) < 0L)
        stop("lag[1] must be <= lag[2]")
    return(round(lag[1L:2L]))
}

seqlag <- function (lag, by = 1){
    seq(from = lag[1], to = lag[2], by = by)
}

cross.row.prod <- function(x,u){
    if(nrow(x) != nrow(u)){
        stop("x and u should have the same number of rows.")
    }
    res <- matrix(0,nrow(u),ncol(u)*ncol(x))
    for(i in 1:ncol(x)){
        for(j in 1:ncol(u)){
            res[,(i-1)*ncol(u)+j] = x[,i] * u[,j]
        }
    }
    return(res)
}

Lag <- function (v, k, group = NULL) 
{
  stopifnot(length(k) > 0)
  v <- as.numeric(v)
  if (max(abs(k)) >= length(v)) 
    stop("largest lag in 'k' must be less than 'length(v)'")
  lag.f <- function(x) {
    lagmat <- matrix(nrow = length(x), ncol = length(k))
    n <- length(x)
    for (i in seq(along = k)) {
      lag <- k[i]
      if (lag > 0) 
        lagmat[, i] <- c(rep(NA, lag), x[1:(n - lag)])
      else if (lag < 0) 
        lagmat[, i] <- c(x[(-lag + 1):n], rep(NA, -lag))
      else lagmat[, i] <- x
    }
    lagmat
  }
  lagmat <- if (!is.null(group)) {
    groupLag <- tapply(v, group, lag.f)
    lagmat <- matrix(nrow = length(v), ncol = length(k))
    split(lagmat, group) <- groupLag
    lagmat
  }
  else lag.f(v)
  colnames(lagmat) <- as.character(k)
  drop(lagmat)
}
