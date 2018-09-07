time <- proc.time()
# LOOP ACROSS SIMULATED EXPOSURE-LAG-RESPONSES
#for(j in seq(nrow(combsim))) {
for(j in seq(combsim)) {
  
  # PRINT
  cat("\n\n ",names(combsim)[j],"\n")
  trueeff <- eval_fun_overgrid(x = o3.pred,
                               u = tmp.pred,
                               L,
                               eval(parse(text=combsim[j])))
  #plot(trueeff)
  ################################################################################
  # LOOP ACROSS RANDOMLY SIMULATED DATA
  
  for(i in seq(nsim)) {
    
    # PRINT
    cat(i,"")
    
    # SET THE SEED
    seed <- 13041975 + i
    set.seed(seed)
    
    # SIMULATE THE DATA
    cumeff <- gen_cumeff(x,u,L,eval(parse(text=combsim[j])))
    suppressWarnings(y <- rpois(length(x),exp((log(base[j])+cumeff))))
    
    ################################################################################
    # FIT GAM-REML
    
    p <- 1
    
    # RUN THE MODEL
    mtime <- proc.time()
    cb <- my.crossbasis(x = x,u = u,lag=L,
                        argvarx=list(fun='ps',df=9),
                        argvaru=list(fun='ps',df=9),
                        arglag=list(fun='ps'))
    
    cbgamPen <- my.cbPen(cb)
    model <- gam(y~cb,family=poisson,paraPen=list(cb=cbgamPen),method='REML')
    cp <- pred(x = o3.pred,
              u = tmp.pred,
              xname = "o3",
              uname = "temp",
              lag = L,
              cb,
              model,log=FALSE)
    
    
    #plot(cp,truth = trueeff)
    # STORE THE TIME
    timemodel[[p]][i,j] <- (proc.time()-mtime)[3]
    
    # MEASURE SIGNAL-TO-NOISE
    #cor(na.omit(log(y+1)),log(predict(model,type="response")+1))
 
    
    # STORE THE RESULTS
    predmatrix[[p]][[j]] <- predmatrix[[p]][[j]] + cp$ar
    bias[[p]][[j]] <- bias[[p]][[j]] + (cp$ar - trueeff[[j]])
    cov[[p]][[j]] <- cov[[p]][[j]] + (trueeff[[j]] >= cp$ar-qn*cp$se & 
                                        trueeff[[j]] <= cp$ar+qn*cp$se)
    rmse[[p]][[j]] <- rmse[[p]][[j]] + (cp$ar - trueeff[[j]])^2
    if(i<=nsample) sample[[p]][[j]] <- c(sample[[p]][[j]],list(cp$ar))
    edf[[p]][i,j] <- sum(model$edf[-1])
    nonconv[[p]][i,j] <- !model$converged
    
    ################################################################################
    # FIT GLM-AIC
    
    p <- 2
    
    # LOOP ACROSS GLM WITH DIFFERENT NUMBER OF KNOTS
    modellist <- list()
    mtime <- proc.time()
    for(k in seq(argvarxlist1)) {
      
      # DEFINE THE CROSS-BASIS
      # NB: USE FIRST COLUMN OF Q AS TIME SERIES DATA -> FASTER
      cb <- my.crossbasis(x = x,
                          u = u,
                          lag = L,argvarx=argvarxlist1[[k]],argvaru=argvarulist1[[k]],
                       arglag=arglaglist1[[k]])
      
      # RUN THE MODEL, SAVING IT IN THE LIST WITH MINIMAL INFO (SAVE MEMORY)
      modellist[[k]] <- glm(y~cb,family=poisson,model=F)
    }
    
    # DETERMINE THE BEST GLM FOR AIC
    best <- which.min(sapply(modellist,AIC))
    cb <- my.crossbasis(x = x,
                        u = u,
                        lag = L,argvarx=argvarxlist1[[best]],argvaru=argvarulist1[[best]],
                        arglag=arglaglist1[[best]])
    timemodel[[p]][i,j] <- (proc.time()-mtime)[3]
    
    # PREDICT
    cp <- pred(x = o3.pred,
               u = tmp.pred,
               xname = "o3",
               uname = "temp",
               lag = L,
               cb,
               model,log=FALSE)
    
    # STORE THE RESULTS
    
    predmatrix[[p]][[j]] <- predmatrix[[p]][[j]] + cp$ar
    bias[[p]][[j]] <- bias[[p]][[j]] + (cp$ar - trueeff[[j]])
    cov[[p]][[j]] <- cov[[p]][[j]] + (trueeff[[j]] >= cp$ar-qn*cp$se & 
                                        trueeff[[j]] <= cp$ar+qn*cp$se)
    rmse[[p]][[j]] <- rmse[[p]][[j]] + (cp$ar - trueeff[[j]])^2
    if(i<=nsample) sample[[p]][[j]] <- c(sample[[p]][[j]],list(cp$ar))
    edf[[p]][i,j] <- length(cp$coef) 
    nonconv[[p]][i,j] <- !modellist[[best]]$converged  #k change to best
    
    ################################################################################
    # FIT GAM-AIC
    
    p <- 3
    
    # RUN THE MODEL
    mtime <- proc.time()
    cb <- my.crossbasis(x = x,
                  u = u,lag=L,argvarx=list(fun='ps',df=9),argvaru=list(fun='ps',df=9),
                  arglag=list(fun='ps'))
    
    cbgamPen <- my.cbPen(cb)
    model <- gam(y~cb,family=poisson,paraPen=list(cb=cbgamPen),method="GCV.Cp")
    
    # STORE THE TIME
    timemodel[[p]][i,j] <- (proc.time()-mtime)[3]
    
    # PREDICT
    cp <- pred(x = o3.pred,
               u = tmp.pred,
               xname = "o3",
               uname = "temp",
               lag = L,
               cb,
               model,log=FALSE)
    
    # STORE THE RESULTS
    predmatrix[[p]][[j]] <- predmatrix[[p]][[j]] + cp$ar
    bias[[p]][[j]] <- bias[[p]][[j]] + (cp$ar - trueeff[[j]])
    cov[[p]][[j]] <- cov[[p]][[j]] + (trueeff[[j]] >= cp$ar-qn*cp$se & 
                                        trueeff[[j]] <= cp$ar+qn*cp$se)
    rmse[[p]][[j]] <- rmse[[p]][[j]] + (cp$ar - trueeff[[j]])^2
    if(i<=nsample) sample[[p]][[j]] <- c(sample[[p]][[j]],list(cp$ar))
    edf[[p]][i,j] <- sum(model$edf[-1])
    nonconv[[p]][i,j] <- !model$converged
    
    
    ################################################################################
    # FIT GAMcr
    
    p <- 4
    
    # RUN THE MODEL
    mtime <- proc.time()
    cb <- my.crossbasis(x = x,
                        u = u,lag=L,argvarx=list(fun='cr',df=9),argvaru=list(fun='cr',df=9),
                        arglag=list(fun='cr'))
    cbgamPen <- my.cbPen(cb)
    model <- gam(y~cb,family=poisson,paraPen=list(cb=cbgamPen),method='REML')
    # STORE THE TIME
    timemodel[[p]][i,j] <- (proc.time()-mtime)[3]
    
    # PREDICT
    cp <- pred(x = o3.pred,
               u = tmp.pred,
               xname = "o3",
               uname = "temp",
               lag = L,
               cb,
               model,log=FALSE)
    
    # STORE THE RESULTS
    predmatrix[[p]][[j]] <- predmatrix[[p]][[j]] + cp$ar
    bias[[p]][[j]] <- bias[[p]][[j]] + (cp$ar - trueeff[[j]])
    cov[[p]][[j]] <- cov[[p]][[j]] + (trueeff[[j]] >= cp$ar-qn*cp$se & 
                                        trueeff[[j]] <= cp$ar+qn*cp$se)
    rmse[[p]][[j]] <- rmse[[p]][[j]] + (cp$ar - trueeff[[j]])^2
    if(i<=nsample) sample[[p]][[j]] <- c(sample[[p]][[j]],list(cp$ar))
    edf[[p]][i,j] <- sum(model$edf[-1])
    nonconv[[p]][i,j] <- !model$converged
    
    ################################################################################
    # FIT GAMps21
    
    p <- 5
    
    # RUN THE MODEL
    mtime <- proc.time()
    cb <- my.crossbasis(x = x,u = u,lag=L,
                        argvarx=list(fun='ps',df=9),
                        argvaru=list(fun='ps',df=9),
                        arglag=list(fun='ps',diff=1))
    cbgamPen <- my.cbPen(cb)
    model <- gam(y~cb,family=poisson,paraPen=list(cb=cbgamPen),method='REML')
    # STORE THE TIME
    timemodel[[p]][i,j] <- (proc.time()-mtime)[3]
    
    # PREDICT
    cp <- pred(x = o3.pred,
               u = tmp.pred,
               xname = "o3",
               uname = "temp",
               lag = L,
               cb,
               model,log=FALSE)
    
    # STORE THE RESULTS
    predmatrix[[p]][[j]] <- predmatrix[[p]][[j]] + cp$ar
    bias[[p]][[j]] <- bias[[p]][[j]] + (cp$ar - trueeff[[j]])
    cov[[p]][[j]] <- cov[[p]][[j]] + (trueeff[[j]] >= cp$ar-qn*cp$se & 
                                        trueeff[[j]] <= cp$ar+qn*cp$se)
    rmse[[p]][[j]] <- rmse[[p]][[j]] + (cp$ar - trueeff[[j]])^2
    if(i<=nsample) sample[[p]][[j]] <- c(sample[[p]][[j]],list(cp$ar))
    edf[[p]][i,j] <- sum(model$edf[-1])
    nonconv[[p]][i,j] <- !model$converged
    
    ################################################################################
    # FIT GAM-last
    
    p <- 6
    
    # RUN THE MODEL
    mtime <- proc.time()
    cb <- my.crossbasis(x = x,u = u,lag=L,
                        argvarx=list(fun='ps',df=9),
                        argvaru=list(fun='ps',df=9),
                        arglag=list(fun='ps'))
    cbgamPen <- my.cbPen(cb,addSlag=rep(0:1,c(6,4)))
    model <- gam(y~cb,family=poisson,paraPen=list(cb=cbgamPen),method='REML')
    # STORE THE TIME
    timemodel[[p]][i,j] <- (proc.time()-mtime)[3]
    
    # PREDICT
    cp <- pred(x = o3.pred,
               u = tmp.pred,
               xname = "o3",
               uname = "temp",
               lag = L,
               cb,
               model,log=FALSE)
    
    # STORE THE RESULTS
    predmatrix[[p]][[j]] <- predmatrix[[p]][[j]] + cp$ar
    bias[[p]][[j]] <- bias[[p]][[j]] + (cp$ar - trueeff[[j]])
    cov[[p]][[j]] <- cov[[p]][[j]] + (trueeff[[j]] >= cp$ar-qn*cp$se & 
                                        trueeff[[j]] <= cp$ar+qn*cp$se)
    rmse[[p]][[j]] <- rmse[[p]][[j]] + (cp$ar - trueeff[[j]])^2
    if(i<=nsample) sample[[p]][[j]] <- c(sample[[p]][[j]],list(cp$ar))
    edf[[p]][i,j] <- sum(model$edf[-1])
    nonconv[[p]][i,j] <- !model$converged
    
    ################################################################################
    # FIT GAM-quad
    
    p <- 7
    
    # RUN THE MODEL
    mtime <- proc.time()
    cb <- my.crossbasis(x = x,u = u,lag=L,
                        argvarx=list(fun='ps',df=9),
                        argvaru=list(fun='ps',df=9),
                        arglag=list(fun='ps',fx=T))
    
    C <- do.call(onebasis,c(list(x=0:L),attr(cb,"arglag")))
    D <- diff(diag(L+1),diff=2)
    P <- diag((0:(L-2))^2)# has changed
    Slag2 <- t(C)%*%t(D)%*%P%*%D%*%C
    cbgamPen <- my.cbPen(cb,addSlag=Slag2)
    model <- gam(y~cb,family=poisson,paraPen=list(cb=cbgamPen),method='REML')
    # STORE THE TIME
    timemodel[[p]][i,j] <- (proc.time()-mtime)[3]
    
    # PREDICT
    cp <- pred(x = o3.pred,
               u = tmp.pred,
               xname = "o3",
               uname = "temp",
               lag = L,
               cb,
               model,log=FALSE)
    
    # STORE THE RESULTS
    predmatrix[[p]][[j]] <- predmatrix[[p]][[j]] + cp$ar
    bias[[p]][[j]] <- bias[[p]][[j]] + (cp$ar - trueeff[[j]])
    cov[[p]][[j]] <- cov[[p]][[j]] + (trueeff[[j]] >= cp$ar-qn*cp$se & 
                                        trueeff[[j]] <= cp$ar+qn*cp$se)
    rmse[[p]][[j]] <- rmse[[p]][[j]] + (cp$ar - trueeff[[j]])^2
    if(i<=nsample) sample[[p]][[j]] <- c(sample[[p]][[j]],list(cp$ar))
    edf[[p]][i,j] <- sum(model$edf[-1])
    nonconv[[p]][i,j] <- !model$converged
    
    ################################################################################
    # FIT GAM-exp
    
    p <- 8
    
    # RUN THE MODEL
    mtime <- proc.time()
    cb <- my.crossbasis(x = x,u = u,lag=L,
                        argvarx=list(fun='ps',df=9),
                        argvaru=list(fun='ps',df=9),
                        arglag=list(fun='ps',fx=T))
    cbgamPen <- my.cbPen(cb,addSlag=exp(0:9))
    model <- gam(y~cb,family=poisson,paraPen=list(cb=cbgamPen),method='REML')
    # STORE THE TIME
    timemodel[[p]][i,j] <- (proc.time()-mtime)[3]
    
    # PREDICT
    cp <- pred(x = o3.pred,
               u = tmp.pred,
               xname = "o3",
               uname = "temp",
               lag = L,
               cb,
               model,log=FALSE)
    
    # STORE THE RESULTS
    predmatrix[[p]][[j]] <- predmatrix[[p]][[j]] + cp$ar
    bias[[p]][[j]] <- bias[[p]][[j]] + (cp$ar - trueeff[[j]])
    cov[[p]][[j]] <- cov[[p]][[j]] + (trueeff[[j]] >= cp$ar-qn*cp$se & 
                                        trueeff[[j]] <= cp$ar+qn*cp$se)
    rmse[[p]][[j]] <- rmse[[p]][[j]] + (cp$ar - trueeff[[j]])^2
    if(i<=nsample) sample[[p]][[j]] <- c(sample[[p]][[j]],list(cp$ar))
    edf[[p]][i,j] <- sum(model$edf[-1])
    nonconv[[p]][i,j] <- !model$converged
  }
  
  ################################################################################
  # COMPUTE THE AVERAGE FOR THE STATISTICS OF THE WHOLE SURFACE
  
  for(i in seq(predmatrix)) {
    predmatrix[[i]][[j]] <- predmatrix[[i]][[j]]/nsim
    bias[[i]][[j]] <- bias[[i]][[j]]/nsim
    cov[[i]][[j]] <- cov[[i]][[j]]/nsim
    rmse[[i]][[j]] <- sqrt(rmse[[i]][[j]]/nsim)
  }
  
}
(tottime <- proc.time()-time)

################################################################################
# SAVE

save.image("simext.RData")

#
