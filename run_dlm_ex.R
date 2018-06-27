rm(list=ls())
options(max.print=999999)
library(dlm)
library(ggplot2, warn.conflicts = FALSE)

dir = "P:/pCloud Sync/MSE_RS_test"
setwd(dir)
data = read.csv("CarlAnchovetaData.csv")
source('P:/pCloud Sync/MSE_RS_test/kfdlm.R')

## run standard ricker 
# ricker =  runFKFdlm(rt=data$rt,st=data$sbt,
#                     init_sd_obs_kf=0.6,init_sd_alpha_kf=0.01,init_sd_beta_kf=0.2,init_alpha_lm=1.0,init_beta_lm=-0.23,
#                     true_alpha=NA,true_beta=NA,fkf_estimator="fitVobs",scaler4beta=1,true_sd_alpha=NA,true_sd_beta=NA,true_sd_obs=NA,
#                     getplot=T,plot_path=dir,txtlab="anchovy",ylim_alpha=2,ylim_beta=1,hess=T,getinit=T )
                    




    for( j in c("fitVobs","fitVobsVa","fitVobsVb","fitAll",
                "dlmkf_nopr","dlmkf_pr_ars","dlmkf_pr_brs",
                "dlmkf_pr_Vb0","dlmkf_pr_Va0","dlmkf_pr_Vab0"
               )) {
    
      switch(j,
             "fitVobs"      = { s2_obs=0.6; s2_alpha=0.2; s2_beta=0.2; a0=2; b0=-1; hess=T },
             "fitVobsVa"    = { s2_obs=0.6; s2_alpha=0.2; s2_beta=0.2; a0=2; b0=-1; hess=T },
             "fitVobsVb"    = { s2_obs=0.6; s2_alpha=0.2; s2_beta=0.2; a0=2; b0=-1; hess=T },
             "fitAll"       = { s2_obs=0.6; s2_alpha=0.2; s2_beta=0.2; a0=2; b0=-1; hess=T },
             "dlmkf_nopr"   = { s2_obs=0.6; s2_alpha=0.6; s2_beta=0.6; a0=2; b0=-1; hess=F  },
             "dlmkf_pr_ars" = { s2_obs=0.6; s2_alpha=0.2; s2_beta=0.2/4; a0=2; b0=-1; hess=F },
             "dlmkf_pr_brs" = { s2_obs=0.6; s2_alpha=0.2/4; s2_beta=0.2; a0=2; b0=-1; hess=F },
             "dlmkf_pr_Vb0" = { s2_obs=0.6; s2_alpha=0.2; s2_beta=0.; a0=2; b0=-1; hess=F },
             "dlmkf_pr_Va0" = { s2_obs=0.6; s2_alpha=0.; s2_beta=0.2; a0=2; b0=-1; hess=F },
             "dlmkf_pr_Vab0" = { s2_obs=0.6; s2_alpha=0.; s2_beta=0.; a0=2; b0=-1;hess=F }
             )

      if( gsub("_.*","",j) =="dlmkf") fkf_estimator = "dlmkf" else fkf_estimator = j
      
    runFKFdlm(rt=data$rt,st=data$sbt,
              init_sd_obs_kf=s2_obs,
              init_sd_alpha_kf=s2_alpha,
              init_sd_beta_kf=s2_beta,
              init_alpha_lm=1.0,init_beta_lm=-0.23,
              true_alpha=NA,true_beta=NA,fkf_estimator=fkf_estimator,scaler4beta=1,true_sd_alpha=NA,true_sd_beta=NA,true_sd_obs=NA,
              getplot=T,plot_path=dir,txtlab=paste0("anchovy","_",j),ylim_alpha=10,ylim_beta=5,hess=hess,getinit=T )
  }


