my.cbPen <- function (cb, sp = -1, addSlag = NULL)
{
  if (all(class(cb) != "crossbasis") & all(class(cb) != "onebasis")) 
    stop("first argument must be object of class 'crossbasis' or 'onebasis")
  attr <- attributes(cb)
  fff <- c(attr$argvarx$fun,attr$argvaru$fun,attr$arglag$fun)
  fx <- c(!fff[1] %in% c("ps", "cr") || attr$argvarx$fx,
          !fff[2] %in% c("ps", "cr") || attr$argvaru$fx,
          !fff[3] %in% c("ps", "cr") || attr$arglag$fx)
  Slist <- list()
  if (!fx[1]) 
    Slist <- c(Slist, list(Svarx = attr$argvarx$S %x% diag(attr$df[3]*attr$df[2]))) #diag(n) generate n-dimensional identity matrix 
  if (!fx[2]) 
    Slist <- c(Slist, list(Svaru = diag(attr$df[1])%x% attr$argvaru$S %x% diag(attr$df[3])))
  if (!fx[3]) 
    Slist <- c(Slist, list(Slag = diag(attr$df[1]*attr$df[2]) %x% attr$arglag$S))
  Slist <- lapply(Slist, function(X) X/eigen(X, symmetric = TRUE, 
                                             only.values = TRUE)$values[1])
  one <- any(class(cb) == "onebasis")
  if (one & !is.null(addSlag)) 
    stop("penalties on lag not allowed for class 'onebasis")
  if (!is.null(addSlag)) 
    Slist <- c(Slist, mkaddSlag(addSlag, attr$df))
  rank <- sapply(Slist, findrank)
  npen <- length(Slist)
  if (npen == 0L) 
    stop("no penalization defined")
  if (length(sp) == 1L) 
    sp <- rep(sp, npen)
  if (!is.numeric(sp) || length(sp) != npen) 
    stop("'sp' must be numeric and consistent with number of penalty terms")
  names(sp) <- names(Slist)
  res <- c(Slist, list(rank = rank, sp = sp))
  return(res)
}
mkaddSlag <-function (addSlag, d) 
{
  Slist <- if (is.list(addSlag)) 
    addSlag
  else list(addSlag)
  if (!all(sapply(Slist, is.numeric))) 
    stop("non-numeric values supplied in 'addSlag'")
  Slist <- lapply(Slist, function(x) if (is.matrix(x)) 
    x
    else diag(x))
  for (i in seq(Slist)) if (any(dim(Slist[[i]]) != d[3])) 
    stop("terms in addSlag with dimensions not consistent with basis for lag")
  Slist <- lapply(Slist, function(X) {
    X <- X/eigen(X, symmetric = TRUE, only.values = TRUE)$values[1]
    diag(d[1]*d[2]) %x% X
  })
  names(Slist) <- paste0("Slag", seq(Slist) + 1)
  return(Slist)
}

findrank <- function (X) 
{
  ev <- eigen(X, symmetric = TRUE, only.values = TRUE)$values
  rank <- sum(ev > max(ev) * .Machine$double.eps * 10)
  return(rank)
}