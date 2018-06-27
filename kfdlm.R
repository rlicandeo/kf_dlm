runFKFdlm = function(rt,st,init_sd_obs_kf=0.2,init_sd_alpha_kf=1.0,init_sd_beta_kf=1.0,init_alpha_lm=2,init_beta_lm=-1,
                       true_alpha=NA,true_beta=NA,fkf_estimator,scaler4beta=1,true_sd_alpha=NA,true_sd_beta=NA,true_sd_obs=NA,
                       getplot=T,plot_path=NA,txtlab="",ylim_alpha=10,ylim_beta=3,hess=F,getinit=F) {
    
    require(PerformanceAnalytics, quietly = TRUE,  warn.conflicts = FALSE)
    require(ggplot2)
    library(gridExtra)
    # require("ggpubr")
    # require("magrittr")
    require(dlm)
    
    y = log( (rt/scaler4beta)/(st/scaler4beta) )
    x = st/scaler4beta
     

    if(getinit==T)
    {
      
      ## lm( log(rt/st) ~ st )
      fit = nls(y ~ a + b*x, start=list(a=init_alpha_lm,b=init_beta_lm),
                model=T, control=list(maxiter = 1000, warnOnly = F, tol = 1e-03))
      
      init_sd_obs = as.numeric(summary(fit)[3])
      init_alpha0 = coef(fit)[[1]]
      init_beta0 = coef(fit)[[2]]
      conv_lm = fit$convInfo$finTol
      conv_inf_lm = fit$convInfo$stopMessage
      cat(paste("lm_estimates:","alpha=",init_alpha0,"beta=",init_beta0,"sd_obs=",init_sd_obs,"conv=",conv_inf_lm),"\n")
      init_sd_alpha = init_sd_alpha_kf    # Variance of the alpha regression parameter
      init_sd_beta = init_sd_beta_kf      # Variance of the beta regression parameter  
      
    } else {
      
      ## Specifying initial Variance for KF 
      init_sd_obs = init_sd_obs_kf        # Variance of observation equation
      init_sd_alpha = init_sd_alpha_kf    # Variance of the alpha regression parameter
      init_sd_beta = init_sd_beta_kf      # Variance of the beta regression parameter  
      init_alpha0 = init_alpha_lm
      init_beta0 = init_beta_lm
    }
    
    
    
    
    
    switch(fkf_estimator,
           
           "fitVobs" = {
             
             start.vals = log(init_sd_obs^2)

             buildTVP <- function(parm, x.mat){
               parm <- exp(parm)
               return( dlmModReg(X=x.mat, dV=parm[1], dW=c(0.,0.)) )
             }
             
             TVP.mle = dlmMLE(y=y, parm=start.vals, x.mat=x, build=buildTVP,
                              hessian=hess,
                      method = "Nelder-Mead", control = list(maxit = 10000))
             sd_estimates <- sqrt(exp(TVP.mle$par))
             
             # Build fitted ss model, passing to it x as the matrix X in the model
             TVP.dlm <- buildTVP(TVP.mle$par, x)
             TVP.f <- dlmFilter(y = y, mod = TVP.dlm)
             # Optimal estimates of theta_t given information available at time T.
             TVP.s <- dlmSmooth(TVP.f)
             
             # extract smoothed states - intercept and slope coefs
             alpha.s = TVP.s$s[-1,1,drop=FALSE]
             beta.s  = TVP.s$s[-1,2,drop=FALSE]
             colnames(alpha.s) = "alpha"
             colnames(beta.s)  = "beta"
             nT = length(alpha.s)
             
             if(hess==T) {
               # extract std errors - dlmSvd2var gives list of MSE matrices
               mse.list = dlmSvd2var(TVP.s$U.S, TVP.s$D.S)
               se.mat = t(sapply(mse.list, FUN=function(x) sqrt(diag(x))))
               se.xts = se.mat[-1, ]
               colnames(se.xts) = c("alpha", "beta")
               a.u = alpha.s + 1.96*se.xts[, "alpha"]
               a.l = alpha.s - 1.96*se.xts[, "alpha"]
               b.u = beta.s  + 1.96*se.xts[, "beta"]
               b.l = beta.s  - 1.96*se.xts[, "beta"]
             } else {
               se.xts = matrix(0,nT,2)
               colnames(se.xts) = c("alpha", "beta")
               a.u = alpha.s
               a.l = alpha.s
               b.u = beta.s
               b.l = beta.s
             }
             
             sd_alpha = 0.
             sd_beta = 0.
             sd_obs = sd_estimates
             conv_kf = TVP.mle$convergence
             nloglik = TVP.mle$value
             npar = length(sd_estimates)
             aic = npar + 2*nloglik
             
           },
           
           "fitVobsVa" = {

             start.vals = c(log(init_sd_obs^2),log(init_sd_alpha^2))

             buildTVP <- function(parm, x.mat){
              parm <- exp(parm)
              return( dlmModReg(X=x.mat, dV=parm[1], dW=c(parm[2],0.)) )
             }
             
             TVP.mle = dlmMLE(y=y, parm=start.vals, x.mat=x, build=buildTVP,
                            hessian = hess,
                      method = "Nelder-Mead", control = list(maxit = 10000))
             sd_estimates <- sqrt(exp(TVP.mle$par))
             
             # Build fitted ss model, passing to it x as the matrix X in the model
             TVP.dlm <- buildTVP(TVP.mle$par, x)
             TVP.f <- dlmFilter(y = y, mod = TVP.dlm)
             TVP.s <- dlmSmooth(TVP.f)            
             
             # extract smoothed states - intercept and slope coefs
             alpha.s = TVP.s$s[-1,1,drop=FALSE]
             beta.s  = TVP.s$s[-1,2,drop=FALSE]
             colnames(alpha.s) = "alpha"
             colnames(beta.s)  = "beta"
             nT = length(alpha.s)
             
             if(hess==T) {
               # extract std errors - dlmSvd2var gives list of MSE matrices
               mse.list = dlmSvd2var(TVP.s$U.S, TVP.s$D.S)
               se.mat = t(sapply(mse.list, FUN=function(x) sqrt(diag(x))))
               se.xts = se.mat[-1, ]
               colnames(se.xts) = c("alpha", "beta")
               a.u = alpha.s + 1.96*se.xts[, "alpha"]
               a.l = alpha.s - 1.96*se.xts[, "alpha"]
               b.u = beta.s  + 1.96*se.xts[, "beta"]
               b.l = beta.s  - 1.96*se.xts[, "beta"]
             } else {
               se.xts = matrix(0.,nT,2.)
               colnames(se.xts) = c("alpha", "beta")
               a.u = alpha.s
               a.l = alpha.s
               b.u = beta.s
               b.l = beta.s
             }
             
             sd_alpha = sd_estimates[2]
             sd_beta = 0.
             sd_obs = sd_estimates[1]
             conv_kf = TVP.mle$convergence
             nloglik = TVP.mle$value
             npar = length(sd_estimates)
             aic = npar + 2*nloglik
             
             
           },
           
           "fitVobsVb" = {

             start.vals = c(log(init_sd_obs^2),log(init_sd_beta^2))

             buildTVP <- function(parm, x.mat){
              parm <- exp(parm)
              return( dlmModReg(X=x.mat, dV=parm[1], dW=c(0.,parm[2])) )
             }
             
             TVP.mle = dlmMLE(y=y, parm=start.vals, x.mat=x, build=buildTVP,
                            hessian = hess,
                    # method = "L-BFGS-B", lower=c(log(1e-9),log(1e-9)), upper=c(log(0.6),log(0.6)) )
                    # method = "BFGS") ## dont work very well with beta_RS (dWbeta goes to zero and dW does the opposite)
                    method = "Nelder-Mead", control = list(maxit = 10000))

             sd_estimates <- sqrt(exp((TVP.mle$par)))
             
             ## Build fitted ss model, passing to it x as the matrix X in the model
             TVP.dlm <- buildTVP(TVP.mle$par, x)
             TVP.f <- dlmFilter(y = y, mod = TVP.dlm)
             TVP.s <- dlmSmooth(TVP.f)       
             
             ## extract smoothed states - intercept and slope coefs
             alpha.s = TVP.s$s[-1,1,drop=FALSE]
             beta.s  = TVP.s$s[-1,2,drop=FALSE]
             colnames(alpha.s) = "alpha"
             colnames(beta.s)  = "beta"
             nT = length(alpha.s)
             
             if(hess==T) {
               # extract std errors - dlmSvd2var gives list of MSE matrices
               mse.list = dlmSvd2var(TVP.s$U.S, TVP.s$D.S)
               se.mat = t(sapply(mse.list, FUN=function(x) sqrt(diag(x))))
               se.xts = se.mat[-1, ]
               colnames(se.xts) = c("alpha", "beta")
               a.u = alpha.s + 1.96*se.xts[, "alpha"]
               a.l = alpha.s - 1.96*se.xts[, "alpha"]
               b.u = beta.s  + 1.96*se.xts[, "beta"]
               b.l = beta.s  - 1.96*se.xts[, "beta"]
             } else {
               se.xts = matrix(0,nT,2)
               colnames(se.xts) = c("alpha", "beta")
               a.u = alpha.s
               a.l = alpha.s
               b.u = beta.s
               b.l = beta.s
             }
             
             sd_alpha = 0.
             sd_beta = sd_estimates[2]
             sd_obs = sd_estimates[1]
             conv_kf = TVP.mle$convergence
             nloglik = TVP.mle$value
             npar = length(sd_estimates)
             aic = npar + 2*nloglik
             
             
           },
           
           "fitAll" = {

             start.vals = c(log(init_sd_obs^2),log(init_sd_alpha^2),log(init_sd_beta^2))

             buildTVP <- function(parm, x.mat){
                parm <- exp(parm)
                return( dlmModReg(X=x.mat, dV=parm[1], dW=c(parm[2],parm[3])) )
             }
             
             TVP.mle = dlmMLE(y=y, parm=start.vals, x.mat=x, build=buildTVP,
                              hessian = hess,
                      method = "Nelder-Mead", control = list(maxit = 10000))
             
             sd_estimates <- sqrt(exp(TVP.mle$par))
             
             # Build fitted ss model, passing to it x as the matrix X in the model
             TVP.dlm <- buildTVP(TVP.mle$par, x)
             TVP.f <- dlmFilter(y = y, mod = TVP.dlm)
             
             # Optimal estimates of ??_t given information available at time T.
             TVP.s <- dlmSmooth(TVP.f)
             
             # extract smoothed states - intercept and slope coefs
             alpha.s = TVP.s$s[-1,1,drop=FALSE]
             beta.s  = TVP.s$s[-1,2,drop=FALSE]
             colnames(alpha.s) = "alpha"
             colnames(beta.s)  = "beta"
             nT = length(alpha.s)
             
             if(hess==T) {
               # extract std errors - dlmSvd2var gives list of MSE matrices
               mse.list = dlmSvd2var(TVP.s$U.S, TVP.s$D.S)
               se.mat = t(sapply(mse.list, FUN=function(x) sqrt(diag(x))))
               se.xts = se.mat[-1, ]
               colnames(se.xts) = c("alpha", "beta")
               a.u = alpha.s + 1.96*se.xts[, "alpha"]
               a.l = alpha.s - 1.96*se.xts[, "alpha"]
               b.u = beta.s  + 1.96*se.xts[, "beta"]
               b.l = beta.s  - 1.96*se.xts[, "beta"]
             } else {
               se.xts = matrix(0,nT,2)
               colnames(se.xts) = c("alpha", "beta")
               a.u = alpha.s
               a.l = alpha.s
               b.u = beta.s
               b.l = beta.s
             }
             
             sd_alpha = sd_estimates[2]
             sd_beta = sd_estimates[3]
             sd_obs = sd_estimates[1]
             conv_kf = TVP.mle$convergence
             nloglik = TVP.mle$value
             npar = length(sd_estimates)
             aic = npar + 2*nloglik
             
             
           }, 
           

           "dlmkf" = {

             mod <- dlmModReg(x, dV=init_sd_obs^2, m0=c(init_alpha0,init_beta0), dW=c(init_sd_alpha^2,init_sd_beta^2), C0=diag(c(1.^2,1.^2)))
             outF <- dlmFilter(y, mod)
             
             alpha.s = as.matrix(data.frame("alpha"=outF$m[-1,1]))
             beta.s = as.matrix(data.frame("beta.s"=outF$m[-1,2]))
             nT = dim(alpha.s)[1]
             se.xts = matrix(0,nT,2)
             colnames(se.xts) = c("alpha", "beta")
             a.u = alpha.s
             a.l = alpha.s
             b.u = beta.s
             b.l = beta.s
             sd_alpha = NA
             sd_beta = NA
             sd_obs = NA
             conv_kf = NA
             nloglik = NA
             aic = NA
             
           } 
           
    )
    
   
      
      alpha.df <- data.frame(dateTime = index(se.xts), alpha = exp(alpha.s), upr = exp(a.u), lwr = exp(a.l))
      names(alpha.df) <- c("Time", "alpha", "upr", "lwr")
      alpha.df$true_alpha = true_alpha[1:nT]
      
      beta.df  <- data.frame(dateTime = index(se.xts), beta = -1*beta.s, upr = -1*b.u, lwr = -1*b.l)
      names(beta.df) <- c("Time", "beta", "upr", "lwr")
      beta.df$true_beta = 1*true_beta[1:nT]
      
      if(getplot==T) {
      
        
      ## formatting beta parameter
      if(fkf_estimator=="dlmkf") {sd_beta_label=NA} else {sd_beta_label=formatC(sd_beta, format="e", digits=2)}
      # if(fkf_estimator=="dlmkf") {sd_alpha_label=NA} else {sd_alpha_label=formatC(sd_alpha, format="e", digits=2)}
 
      estimates = substitute(
        paste( hat(sigma)[alpha],"=", sd_alpha,"  ", hat(sigma)[beta],"=", sd_beta, "  ", hat(sigma)[obs],"=",sd_obs  ),
        list(sd_alpha = round(sd_alpha,6), sd_beta = sd_beta_label, sd_obs = round(sd_obs,3) )
      )
      
      true = substitute(
        paste( sigma[alpha],"=", true_sd_alpha,"  ", sigma[beta],"=", true_sd_beta, "  ", sigma[obs],"=",true_sd_obs  ),
        list(true_sd_alpha = round(true_sd_alpha,3), true_sd_beta = round(true_sd_beta,3) , true_sd_obs = round(true_sd_obs,3))
      )
      
      ## Plotting alpha
      # par(mfrow=c(2,1), mar=c(4,5,0,1), oma=c(0,0,0,0))
      pa = ggplot(data = alpha.df, aes(Time, alpha) ) + geom_point(colour="red", alpha=0.3) + geom_line(colour = "red", alpha=0.3) + geom_ribbon(data=alpha.df, aes(ymin=lwr,ymax=upr), alpha=0.1, fill="blue")
      pa = pa + labs(x = "year", y = expression(hat(alpha)), title = estimates, subtitle= paste0("nll=",round(nloglik,1),"  ", "aic=",round(aic,2)) )
      if( !is.na(true_alpha[1]) ) pa = pa + geom_line(data=alpha.df,aes(Time, true_alpha, colour = "true")) + scale_colour_manual("",values="blue") else pa = pa
      pa <- pa + ylim(0, ylim_alpha)
      pa <- pa + theme(axis.text = element_text(face="bold", size=15),
                       axis.text.y= element_text(angle=0, hjust = 0.5, size=15),
                       axis.text.x= element_text(angle=0, hjust = 0.5, size=15),
                       axis.title = element_text(face="bold", size=20),
                       strip.text = element_text(face="bold", size=15)
                       )
      pa = pa + theme_light()
      pa = pa + theme(legend.justification = c("top"), legend.position = c(.85,.95) )
      # print(pa)
      
      ## Plotting beta
      pb = ggplot(data = beta.df, aes(Time, beta) ) + geom_point(colour="red", alpha=0.3 ) + geom_line(colour="red", alpha=0.3) + geom_ribbon(data=beta.df , aes(ymin=lwr,ymax=upr), alpha=0.1, fill="blue")
      pb = pb + labs(x = "year", y = expression(hat(beta)), title = true, subtitle="")
      if( !is.na(true_beta[1]) )  pb = pb + geom_line(data=beta.df, aes(Time, true_beta, colour = "true")) + scale_colour_manual("",values="blue") else pb = pb
      pb <- pb + ylim(0, ylim_beta)
      pb <- pb + theme(axis.text = element_text(face="bold", size=15),
                       axis.text.y= element_text(angle=0, hjust = 0.5, size=15),
                       axis.text.x= element_text(angle=0, hjust = 0.5, size=15),
                       axis.title = element_text(face="bold", size=20),
                       strip.text = element_text(face="bold", size=15)
                       )
      pb = pb + theme_light()
      pb = pb + theme(legend.justification = c("top"), legend.position = c(.85,.95) )
      # print(pb)
      
      # Q-Q plot
      # qqnorm(TVP.res)
      # qqline(TVP.res)
      # # Plotting Diagnostics for Time Series fits
      # tsdiag(TVP.f)
      
      label = paste0(txtlab,"_",Sys.Date(),".png",sep="")
      ggsave(grid.arrange(pa, pb, nrow = 1), filename=label,
             path = plot_path,
             width = 20, height = 10,
             units = c("cm"),
             dpi = 400)
      
      ## label = paste0(txtlab,"_",Sys.Date(),".pdf",sep="")
      ## pab  = grid.arrange(pa, pb, nrow = 1,  newpage = TRUE)
      ## ggsave(pab, filename=label)
      
      
    }
    
    return( list( arecnT=exp(alpha.s[nT]), brecnT=-1*beta.s[nT], kfconv=conv_kf,
                  sd_obs=sd_obs, sd_alpha=sd_alpha, sd_beta=sd_beta,
                  arectv=exp(alpha.s), brectv=-1*beta.s, alpha.df=alpha.df, beta.df=beta.df)  )
    
  }
  