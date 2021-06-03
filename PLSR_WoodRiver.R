#####################################################################################
### PLSR models AISA ################################################################
### Anna K Schweiger July 31 2020 ###################################################
### with lines of code and the VIPjh function from Shawn P. Serbin's awesome GitHub repo # 
### https://github.com/serbinsh/Spectra_PLSR_models #################################

library(pls)
library(reshape)
library(agricolae)
library(caret)

VIPjh=function(object, j, h) { 
  if (object$method != "oscorespls") stop("Only implemented for orthogonal scores algorithm. Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1) stop("Only implemented for single-response models")
  b=c(object$Yloadings)[1:h]
  T=object$scores[,1:h, drop = FALSE]
  SS=b^2 * colSums(T^2)
  W=object$loading.weights[,1:h, drop = FALSE]
  Wnorm2=colSums(W^2)
  VIP=sqrt(nrow(W) * sum(SS * W[j,]^2 / Wnorm2) / sum(SS)) 
  return(VIP)
}

### Read data 
# allData <- read.csv("./plot_chem_spec.csv")
allData <- read.csv("./subplot_chem_spec.csv")
# allData <- read.csv("./transect_chem_spec.csv")

inBands <- names(allData)[10:187] # bands to work with
inVars <- names(allData[c(3:9)]) # traits to work with


###############################################
### Find number of components
for (i in 1:length(inVars)){
  inVar <- inVars[i]
  iForm <- paste(inBands,collapse="+")
  iForm <- as.formula(paste(inVar,"~",iForm))
  
  ### Separate data into model CAL/VAL and IVAL (independent or external validation) data
  set.seed(19841002+i)
  rand <- createDataPartition(allData[, inVars[i]], p=0.75, list=F) # here: 75% data for model CAL/VAL
  subData <- allData[rand,] ## CAL
  tstData <- allData[-rand,] ## IVAL
  
  ### Calibrate the model
  nComps <- 15 # max no of components to be evaluated 
  nsims <- 50 # no of simulations
  outMat <- matrix(data=NA,nrow=nsims,ncol=nComps)
  outMatRM <- matrix(data=NA,nrow=nsims,ncol=nComps)
  
  for (nsim in seq(nsims)){
    print(nsim)
    flush.console()
    segs <- cvsegments(N = nrow(subData), k = 5, type="random") # k fold CV, 5 folds = 80:20 split
    resNCOMP <- plsr(iForm,data=subData,ncomp=nComps, 
                     validation="CV", segments=segs,
                     method="oscorespls")
    resPRESS <- as.vector(resNCOMP$validation$PRESS)
    outMat[nsim,seq(resNCOMP$validation$ncomp)] <-resPRESS
    resRMSEP <- as.numeric(RMSEP(resNCOMP,estimate="CV",intercept=F)$val)
    outMatRM[nsim,seq(resNCOMP$validation$ncomp)] <-resRMSEP
  }
  
  ### PRESS statistic: Tukey test - sign diff among no of components?
  pressDF <- as.data.frame(outMat)
  names(pressDF) <- as.character(seq(nComps))
  pressDFres <- melt(pressDF)
  
  modi <- lm (value~variable, pressDFres) 
  tuk <- HSD.test (modi,"variable")
  tuk_dat <- as.data.frame(tuk$groups)
  tuk_dat$var <- as.numeric(row.names(tuk_dat))
  tuk_dat <- tuk_dat[order(tuk_dat$var,decreasing = F),]
  letters <- as.character(tuk_dat$groups)
  
  jpeg(paste("./", inVar, "_PRESS.jpg",sep=""),width=6,height=4.5,units="in",res=200)
  par(bty="l")
  boxplot(pressDFres$value~pressDFres$variable, xlab="n Components",ylab="PRESS",main=inVar)
  text(x=1:max(as.numeric(pressDFres$variable)),  y=rep(max(pressDFres$value),15),letters)
  dev.off()
  
  ### RMSEP: Tukey test 
  RMDF <- as.data.frame(outMatRM)
  names(RMDF) <- as.character(seq(nComps))
  RMDFres <- melt(RMDF)
  
  modi <- lm (value~variable, RMDFres) 
  tuk <- HSD.test (modi,"variable")
  tuk_dat <- as.data.frame(tuk$groups)
  tuk_dat$var <- as.numeric(row.names(tuk_dat))
  tuk_dat <- tuk_dat[order(tuk_dat$var,decreasing = F),]
  letters <- as.character(tuk_dat$groups)
  
  jpeg(paste("./", inVar, "_RMSEP.jpg",sep=""),width=6,height=4.5,units="in",res=200)
  par(bty="l")
  boxplot(RMDFres$value~RMDFres$variable, xlab="n Components",ylab="RMSEP",main=inVar)
  text(x=1:max(as.numeric(RMDFres$variable)), y=rep(max(RMDFres$value),15),letters)
  dev.off()
}


#######################################
###### Final models ###################

### Make dataframe containing selected no of components (=data)  per trait
compis <- as.data.frame(matrix(nrow=1, ncol=length(inVars),dimnames = list(NULL, c(inVars)), 
                               data=c(3,4,2,4,2,3,2))) 

### One trait at a time ...
k <- 1 
inVar <- inVars[k]  
nCompss <- compis[, inVar] # nCompss could also be a vector of a number of options for ncomps

### Separate data into CAL and IVAL data, use same seed same as above
set.seed(19841002+k)
rand <- createDataPartition(allData[, inVars[k]], p=0.75, list=F) ## 75% data for calibration
subData <- allData[rand,] ## CAL
tstData <- allData[-rand,] ## IVAL

iForm <- paste(inBands,collapse="+")
iForm <- as.formula(paste(inVar,"~",iForm))

set.seed(1840)
for (i in 1:length(nCompss)){
  nsims <- 500
  nComps <- nCompss[i]
  nCutmod <- floor(0.8*nrow(subData)) ## 80% data for calibration, 20% for CV
  
  coefMat <- matrix(data=NA,nrow=nsims,ncol=length(inBands)+1)
  coefStd <- matrix(data=NA,nrow=nsims,ncol=length(inBands)+1)
  vipMat <- matrix(data=NA,ncol=nsims,nrow=length(inBands))
  statMat <- matrix(data=NA,nrow=nsims,ncol=10)
  
  for (nsim in seq(nsims)){
    print(nsim)
    flush.console()
    
    set.seed (19818411+nsim)
    randx <-createDataPartition(subData[, inVars[k]], p=0.8, list=F)
    calData <- allData[randx,] ## CAL
    valData <- allData[-randx,] ## VAL (internal)
    
    resX <- plsr(iForm,data=calData,ncomp=nComps,method="oscorespls") ###
    resS <- plsr(iForm,data=calData,ncomp=nComps, method="oscorespls",scale=T) ### scaled by SD
    
    ### Coefficients (raw and standardized)
    coefs <- as.vector(coef(resX,ncomp=nComps,intercept=T))
    zcoef <- as.vector(coef(resS,ncomp=nComps,intercept=T))
    coefMat[nsim,] <- coefs
    coefStd[nsim,] <- zcoef ### standardized coeffis for importance of wvls
    
    ### VIP
    vip <- c()
    for (j in seq(length(inBands))){
      vip <- c(vip,VIPjh(resS,j,nComps))
    }
    vipMat[,nsim] <- vip
    
    ### Model stats
    fitX <- as.vector(unlist(resX$fitted.values[,1,nComps])) ### fitted values
    preX <- as.vector(unlist(predict(resX,ncomp=nComps,newdata=valData))) ### internal val
    fitBias <- mean(calData[,inVar]-fitX)
    valBias <- mean(valData[,inVar]-preX)
    fitR2 <- summary(lm(calData[,inVar]~fitX))$r.squared
    valR2 <- summary(lm(valData[,inVar]~preX))$r.squared
    fitP <- summary(lm(calData[,inVar]~fitX))[[4]][[8]]
    valP <- summary(lm(valData[,inVar]~preX))[[4]][[8]]
    fitRMSE <- sqrt(mean((calData[,inVar]-fitX)^2))
    valRMSE <- sqrt(mean((valData[,inVar]-preX)^2))
    fitRMSEperc <- fitRMSE/(max(calData[,inVar])-min(subData[,inVar]))*100
    valRMSEperc <- valRMSE/(max(valData[,inVar])-min(valData[,inVar]))*100
    outVec <- c(fitR2,fitP, fitRMSE,fitBias,fitRMSEperc,
                valR2,valP, valRMSE,valBias,valRMSEperc)
    statMat[nsim,] <- outVec
  }
  
  statMat <- as.data.frame(statMat)
  names(statMat) <- c("fitR2", "fitP","fitRMSE","fitBias","fitRMSEperc",
                      "valR2","valP","valRMSE","valBias","valRMSEperc")
  
  meanStat <- as.data.frame(t(colMeans(statMat)))
  
  write.csv(statMat,paste("./",inVar,"_", nComps,"comps_stats.csv", sep=""),row.names=FALSE)
  write.csv(meanStat, paste("./",inVar,"_", nComps, "_mean_stats.csv", sep=""))
  
  
  ### Coefficients
  coeffis <- data.frame(matrix(nrow = length(inBands)+1, ncol = 3))
  names(coeffis) <- c("bands", "mean","stdv")
  coeffis$bands <- c("Intercept",inBands)
  coeffis$mean <- apply(coefMat,MAR=2,FUN=mean)
  coeffis$stdv <- apply(coefMat,MAR=2,FUN=sd)
  write.table(coeffis, paste("./", inVar, "_",nComps,"comps_coeffMEAN.csv", sep=""),
              sep=",",col.names=T,row.names=F)
  
  ### Predictions
  specMat <- subData[,inBands]
  specMat <- cbind(rep(1,nrow(specMat)),specMat)
  specMat <- as.matrix(specMat)
  
  predMat <- specMat%*%t(coefMat)
  predMean <- apply(predMat,FUN=mean,MAR=1)
  predStdv <- apply(predMat,FUN=sd,MAR=1)
  
  preds <- subData[, !(names(subData) %in% inBands|names(subData)=="RAND")]
  preds[,paste("predMean_",inVar,sep="")] <- predMean
  preds[,paste("predStdv_",inVar,sep="")] <- predStdv
  write.csv(preds,paste("./", inVar, "_", nComps,"comps_preds_mod.csv", sep=""), row.names=FALSE)
  
  ### Plot predictions 
  modCI <- quantile(statMat$fitR2, probs=c(0.05,0.95))
  modCIval <- quantile(statMat$valR2, probs=c(0.05,0.95))
  RMCI <- quantile(statMat$fitRMSE, probs=c(0.05,0.95))
  RMCIperc <- quantile(statMat$fitRMSEperc, probs=c(0.05,0.95))
  RMCIval <- quantile(statMat$valRMSE, probs=c(0.05,0.95))
  RMCIvalperc <- quantile(statMat$valRMSEperc, probs=c(0.05,0.95))
  
  formi <- as.formula(paste(paste(inVar," ~ predMean_",inVar, sep="")))
  lmformi <-   as.formula(paste(paste("predMean_",inVar," ~ ", inVar, sep="")))
  
  jpeg(paste("./",inVar,"_", nComps, "comp_predplot.jpg",sep=""),
       width=6,height=5,units="in",res=300)
  plot(formi, data= preds, pch=16,cex=0.8,ylab="measured",xlab="predicted",
       main=inVar,xlim=c(min(preds[,grepl(paste0("predMean_", inVar),names(preds))] -
                               preds[,grepl(paste0("predStdv_", inVar),names(preds))]),
                         max(preds[,grepl(paste0("predMean_", inVar),names(preds))] +
                               preds[,grepl(paste0("predStdv_", inVar),names(preds))])))
  
  abline(lm(lmformi, data= preds))
  abline(a = 0,b = 1, lty=2)
  arrows(preds[,grepl(paste0("predMean_", inVar),names(preds))],
         preds[,match(inVar,names(preds))],
         preds[,grepl(paste0("predMean_", inVar),names(preds))]+preds[,grepl(paste0("predStdv_", inVar),names(preds))],
         preds[,match(inVar,names(preds))],angle=90,length=0.05, lwd=0.8)
  arrows(preds[,grepl(paste0("predMean_", inVar),names(preds))],
         preds[,match(inVar,names(preds))],
         preds[,grepl(paste0("predMean_", inVar),names(preds))]-preds[,grepl(paste0("predStdv_", inVar),names(preds))],
         preds[,match(inVar,names(preds))],angle=90,length=0.05, lwd=0.8)
   legend("topleft", bty="n", cex=0.8,
         c(paste("R² cal= ", sprintf("%.2f",signif (mean(statMat$fitR2),3))," [",signif(modCI[1],2),",",signif(modCI[2],2),"]", sep = ""),
           paste("R² val= ", sprintf("%.2f",signif (mean(statMat$valR2),3))," [",signif(modCIval[1],2),",",signif(modCIval[2],2),"]", sep = ""),
           paste("RMSEP % cal= ", sprintf("%.2f",signif (mean(statMat$fitRMSEperc),3)), " [",signif(RMCIperc[1],3),",",signif(RMCIperc[2],3),"]", sep=""),
           paste("RMSEP % val= ", sprintf("%.2f",signif (mean(statMat$valRMSEperc),3)), " [",signif(RMCIvalperc[1],3),",",signif(RMCIvalperc[2],3),"]", sep=""),
           paste("ncomps =", nComps, sep=" ")))
  dev.off()
  
  ### VIP for plotting
  vipAggr <- as.data.frame(t(apply(vipMat,MAR=1,FUN=quantile,probs=c(0.05,0.5,0.95))))
  vipAggr$mean_VIP <- apply(vipMat,MAR=1,FUN=mean)
  vipAggr$stdv <- apply(vipMat,MAR=1,FUN=sd)
  serr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))  ## standard error of mean 
  vipAggr$se <- apply(vipMat,MAR=1,FUN=serr)
  vipAggr$band <- inBands
  
  ### Standardized coefficients for plotting
  coeff_std <- data.frame(matrix(nrow = length(inBands)+1, ncol = 3))
  names(coeff_std) <- c("bands", "mean","stdv")
  coeff_std$bands <- c("Intercept",inBands)
  coeff_std$mean <- apply(coefStd,MAR=2,FUN=mean)
  coeff_std$stdv <- apply(coefStd,MAR=2,FUN=sd)
  
  ### Plot VIP and standardized coefficients
  jpeg(paste("./", inVar, "_", nComps, "comp_varimp.jpg",sep=""),
       width=6,height=7,units="in",res=300)
  par(mfrow=c(2,1), mar=c(1.5,4,2.5,1.5), oma=c(3,0,0,0))
  plot(coeff_std$mean[-1]~as.numeric(substr(coeff_std$bands[-1],2,nchar(coeff_std$bands[-1]))),
       type="p",pch=19, xlab="",ylab="coeff_stdmean",main=paste(inVar,nComps,"comps",sep = "_"),
       ylim=c(-max(abs(coeff_std$mean[-1])),max(abs(coeff_std$mean[-1]))), bty="l")
  abline(h=0)
  points(abs(coeff_std$mean)[-1]~as.numeric(substr(coeff_std$bands[-1],2,nchar(coeff_std$bands[-1]))),
         xlab="wvl",ylab="coeff_stdmean", col=2, pch=16, cex=0.8)
  lines(abs(coeff_std$mean)[-1]~as.numeric(substr(coeff_std$band[-1],2,nchar(coeff_std$band[-1]))), col=2)
  plot(as.numeric(substr(vipAggr$band,2,nchar(vipAggr$band))),vipAggr$mean_VIP, type="l",
       xlab = "wvl",ylab = "VIP", bty="l")
  polygon(x=c(as.numeric(substr(vipAggr$band,2,nchar(vipAggr$band))),
              rev(as.numeric(substr(vipAggr$band,2,nchar(vipAggr$band))))),
          y=c(vipAggr$mean_VIP+vipAggr$stdv*1.96, rev(vipAggr$mean_VIP-vipAggr$stdv*1.96)),  
          col =  adjustcolor("red", alpha.f = 0.2), border = NA)
  mtext("wavelength(nm)",1,outer = T,line = 1)
  dev.off()
  
  ######################################
  ### Independent (external) validation
  specMatIVAL <- tstData[,inBands]
  specMatIVAL <- cbind(rep(1,nrow(specMatIVAL)),specMatIVAL)
  specMatIVAL <- as.matrix(specMatIVAL)
  
  predMatIVAL <- specMatIVAL%*%t(coefMat)
  predMeanIVAL <- apply(predMatIVAL,FUN=mean,MAR=1)
  predStdvIVAL <- apply(predMatIVAL,FUN=sd,MAR=1)
  
  predsIVAL <- tstData[, !(names(tstData) %in% inBands|names(tstData)=="RAND")]
  predsIVAL[,paste("predMeanIVAL_",inVar,sep="")] <- predMeanIVAL
  predsIVAL[,paste("predStdvIVAL_",inVar,sep="")] <- predStdvIVAL
  write.csv(predsIVAL,paste("./", inVar, "_", nComps, "comps_preds_IVAL.csv", sep=""), 
            row.names=FALSE)
  
  ### Model stats
  statMatIVAL <- matrix(data=NA,nrow=nsims,ncol=5)
  for (nsim in seq(nsims)){
    valBias <- mean(tstData[,inVar]-predMatIVAL[,nsim])
    valR2 <- summary(lm(tstData[,inVar]~predMatIVAL[,nsim]))$r.squared
    valRMSE <- sqrt(mean((tstData[,inVar]-predMatIVAL[,nsim])^2))
    valRMSEperc <- valRMSE/(max(tstData[,inVar])-min(tstData[,inVar]))*100
    valP <- summary(lm(tstData[,inVar]~predMatIVAL[,nsim]))[[4]][[8]]
    outVec <- c(valR2,valP,valRMSE,valBias,valRMSEperc)
    statMatIVAL[nsim,] <- outVec
  }
  
  statMatIVAL <- as.data.frame(statMatIVAL)
  names(statMatIVAL) <- c("IvalR2","IvalP","IvalRMSE","IvalBias","IvalRMSEperc")
  meanIVAL <- as.data.frame(t(colMeans(statMatIVAL)))
  
  write.csv(statMatIVAL,paste("./",inVar,"_", nComps, "comps_stats_IVAL.csv", sep=""),
            row.names=FALSE)
  write.csv(meanIVAL, paste("./",inVar,"_", nComps, "_mean_stats_IVAL.csv", sep=""))
  
  ### Plot
  modCIval <- quantile(statMatIVAL$IvalR2, probs=c(0.05,0.95))
  RMCIval <- quantile(statMatIVAL$IvalRMSE, probs=c(0.05,0.95))
  RMCIvalperc <- quantile(statMatIVAL$IvalRMSEperc, probs=c(0.05,0.95))
  
  formi <- as.formula(paste(paste(inVar," ~ predMeanIVAL_",inVar, sep="")))
  lmformi <-   as.formula(paste(paste("predMeanIVAL_",inVar," ~ ", inVar, sep="")))
  
  jpeg(paste("./", inVar,"_", nComps, "comp_predplot_IVAL.jpg",sep=""),
       width=6,height=5,units="in",res=300)
  plot(formi, data= predsIVAL, pch=16,cex=0.8,ylab="measured",xlab="predicted",
       main=inVar,xlim=c(min(predsIVAL[,grepl(paste0("predMeanIVAL_", inVar),names(predsIVAL))] -
                               predsIVAL[,grepl(paste0("predStdvIVAL_", inVar),names(predsIVAL))]),
                         max(predsIVAL[,grepl(paste0("predMeanIVAL_", inVar),names(predsIVAL))] +
                               predsIVAL[,grepl(paste0("predStdvIVAL_", inVar),names(predsIVAL))])))
  
  abline(lm(lmformi, data= predsIVAL))
  abline(a = 0,b = 1, lty=2)
  arrows(predsIVAL[,grepl(paste0("predMeanIVAL_", inVar),names(predsIVAL))],
         predsIVAL[,match(inVar,names(predsIVAL))],
         predsIVAL[,grepl(paste0("predMeanIVAL_", inVar),names(predsIVAL))]+predsIVAL[,grepl(paste0("predStdvIVAL_", inVar),names(predsIVAL))],
         predsIVAL[,match(inVar,names(predsIVAL))],angle=90,length=0.05, lwd=0.8)
  arrows(predsIVAL[,grepl(paste0("predMeanIVAL_", inVar),names(predsIVAL))],
         predsIVAL[,match(inVar,names(predsIVAL))],
         predsIVAL[,grepl(paste0("predMeanIVAL_", inVar),names(predsIVAL))]-predsIVAL[,grepl(paste0("predStdvIVAL_", inVar),names(predsIVAL))],
         predsIVAL[,match(inVar,names(predsIVAL))],angle=90,length=0.05, lwd=0.8)
  legend("topleft", bty="n", cex=0.8,
        c(paste("R² = ", sprintf("%.2f",signif (mean(statMatIVAL$IvalR2),3))," [",signif(modCIval[1],2),",",signif(modCIval[2],2),"]", sep = ""),
          paste("RMSEP % = ", sprintf("%.2f",signif (mean(statMatIVAL$IvalRMSEperc),3)), " [",signif(RMCIvalperc[1],2),",",signif(RMCIvalperc[2],2),"]", sep=""),
          paste("ncomps =", nComps, sep=" ")))
  dev.off()
}

###### END #########

  