#################################################################
## Compare with splines on prediction with n=10 ###################
library(splines)
library(tidymodels)
library(pROC)
library(plgp)
library(mvtnorm)
library(dplyr)
library(MASS)
library(pracma)
library(Matrix)
library(condmixt)
library(ggplot2)
library("wesanderson")
library(latex2exp)
source("script/TMGGP_functions.R")
n=20
ntest = 50
res.cor10 <- pbsapply(1:10, function(i){
  set.seed(665+i)
  Y.sep.ts20 <- simulateData2(b,a,s=c(0,2),tau2,n,kernel = "seperate")
  train.index = sample(1:n,n/2)
  Y.sep.ts.train10 <- Y.sep.ts20[c(train.index,train.index+n),]
  # plotSimulate(Y.sep.ts.train10,i=1,s= c(0,2))
  Y_exp20 <- simulateData2(b,a,s=c(0,2),tau2,n,kernel = "exp")
  Y_exp.train10 <- Y_exp20[c(train.index,train.index+n),]
  # plotSimulate(Y_exp.train10,i=1,s= c(0,2))
  sep.mtsgp.train10 <- deriveEstimation(Y.sep.ts.train10, "exp", sl = -1, su = 5)
  exp.mtsgp.train10 <- deriveEstimation(Y_exp.train10, "exp", sl = -1, su = 5)

  sep.sep.train10 <- deriveEstimation(Y.sep.ts.train10, "sep.exp")
  exp.sep.train10 <- deriveEstimation(Y_exp.train10, "sep.exp")

  # ttest = rownames(Y.sep.ts.train10) %>% as.numeric()
  t <- (rownames(Y.sep.ts.train10) %>% as.numeric())
  corr.sep10 <- sapply(1:100, function(i){
    Y = Y.sep.ts.train10[,i]
    ttest = rep(seq(min(Y), max(Y), length.out = ntest),2)
    params = sep.mtsgp.train10[i,]
    yhat <- interpolation_sim(t,Y,params,ttest,prediction = FALSE)
    res.mtsgp <- cor(yhat[1:ntest],yhat[(ntest+1):(2*ntest)])
    # res.mtsgp.a = yhat[2*ntest+1]
    # fit splines
    ntrain = n/2
    t1 <- t[1:ntrain]
    t2 <- t[(ntrain+1):(2*ntrain)]
    fit = lm(Y[1:ntrain] ~ ns(t1, knots=3))
    y1hat = predict(fit,data.frame(t1=ttest[1:ntest]))
    fit2 = lm(Y[(ntrain+1):(2*ntrain)] ~ ns(t2, knots=3))
    y2hat = predict(fit2,data.frame(t2=ttest[(ntest+1):(2*ntest)]))
    res.splines = cor(y1hat,y2hat)

    # plot(t1,Y[1:n],col="grey",xlab="Age",ylab="Wages")
    # points(ttest[1:n],yhat[1:n],col="darkgreen")
    # points(ttest[1:n],Y.sep.ts.test[1:n,i])

    # fit sep GP
    params = sep.sep.train10[i,]
    yhat <- interpolation_sim(t,Y,params,ttest)
    res.sep <- cor(yhat[1:ntest],yhat[(ntest+1):(2*ntest)])
    return(c(res.mtsgp, res.splines,res.sep))
  })

  t <- (rownames(Y_exp.train10) %>% as.numeric())
  corr.exp10 <- sapply(1:100, function(i){
    Y = Y_exp.train10[,i]
    ttest = rep(seq(min(Y), max(Y), length.out = ntest),2)
    params = exp.mtsgp.train10[i,]
    yhat <- interpolation_sim(t,Y,params,ttest,prediction = FALSE)
    res.mtsgp <- cor(yhat[1:ntest],yhat[(ntest+1):(2*ntest)])


    # fit splines
    ntrain = n/2
    t1 <- t[1:ntrain]
    t2 <- t[(ntrain+1):(2*ntrain)]
    fit = lm(Y[1:ntrain] ~ ns(t1, knots=3))
    y1hat = predict(fit,data.frame(t1=ttest[1:ntest]))
    fit2 = lm(Y[(ntrain+1):(2*ntrain)] ~ ns(t2, knots=3))
    y2hat = predict(fit2,data.frame(t2=ttest[(ntest+1):(2*ntest)]))
    res.splines = cor(y1hat,y2hat)

    # plot(t1,Y[1:n],col="grey",xlab="Age",ylab="Wages")
    # points(ttest[1:n],y1hat,col="darkgreen")
    # points(ttest[1:n],Y_exp.test[1:n,i])

    # fit sep GP
    params = exp.sep.train10[i,]
    yhat <- interpolation_sim(t,Y,params,ttest)
    res.sep <- cor(yhat[1:ntest],yhat[(ntest+1):(2*ntest)])
    return(c(res.mtsgp,res.splines,res.sep))
  })

  dat3.1 <- t(corr.exp10) %>% as.data.frame() %>% `colnames<-`(c("MTSGP","splines","sep")) %>%
    tidyr::pivot_longer(everything(),names_to = "model")

  dat3.2 <- t(corr.sep10) %>% as.data.frame() %>% `colnames<-`(c("MTSGP","splines","sep")) %>%
    tidyr::pivot_longer(everything(),names_to = "model")
  dat <- rbind(dat3.1,dat3.2)%>% mutate(generation=c(rep("mtsgp",300),rep("gp",300)))

  # dat %>%
  #   ggplot(aes(x=generation,y=value, fill = model))+
  #   geom_boxplot() +
  #   theme(panel.background = element_rect(fill = "white", colour = "grey20")) +
  #   xlab("Model from which data is generated") +
  #   ylab("AUC") +
  #   scale_x_discrete(breaks=c("gp","mtsgp"),
  #                    labels=c("seperate GPs \n s=2","MTSGP \n s = 2 "))

  dat$generation<-as.factor(dat$generation)
  res.cor <- dat %>%
    group_by(model) %>%
    roc_auc(truth = generation, value, event_level = "second")

  # compare with cross-correlation
  m = nrow(Y.sep.ts.train10)/2
  sep.cc <- sapply(1:ncol(Y.sep.ts.train10), function(i){
    cc <- ccf(Y.sep.ts.train10[1:m,i],Y.sep.ts.train10[(m+1):(2*m),i],plot=FALSE)
    max(cc$acf)
  })

  mtsgp.cc <- sapply(1:ncol(Y_exp.train10), function(i){
    cc <- ccf(Y_exp.train10[1:m,i],Y_exp.train10[(m+1):(2*m),i],plot=FALSE)
    max(cc$acf)
  })

  dat2 <- data.frame(mtsgp = mtsgp.cc, gp = sep.cc) %>%
    tidyr::pivot_longer(everything(),names_to = "generation")

  cccor <- roc(dat2$generation ~ dat2$value, plot = TRUE) # 0.5579

  dat1 <- data.frame(mtsgp = exp.mtsgp.train10[,2], gp = sep.mtsgp.train10[,2]) %>%
    tidyr::pivot_longer(everything(),names_to = "generation")
  mtsgp.a <- roc(dat1$generation ~ log2(dat1$value), plot = TRUE) # 0.673

  return(c(res.cor$.estimate,cccor$auc, mtsgp.a$auc))
}, cl= 4)

res.cor10
rowMeans(res.cor10[,-9])
rowSds(res.cor10[,-9])

#################################################################
## Compare with splines on prediction with n=25 ###################
n=50
ntest = 50
res.cor25 <- pbsapply(1:10, function(i){
  set.seed(665+i)
  Y.sep.ts50 <- simulateData2(b,a,s=c(0,2),tau2,n,kernel = "seperate")
  train.index = sample(1:n,n/2)
  Y.sep.ts.train25 <- Y.sep.ts50[c(train.index,train.index+n),]
  # plotSimulate(Y.sep.ts.train25,i=1,s= c(0,2))
  Y_exp50 <- simulateData2(b,a,s=c(0,2),tau2,n,kernel = "exp")
  Y_exp.train25 <- Y_exp50[c(train.index,train.index+n),]
  # plotSimulate(Y_exp.train25,i=1,s= c(0,2))
  sep.mtsgp.train25 <- deriveEstimation(Y.sep.ts.train25, "exp", sl = -1, su = 5)
  exp.mtsgp.train25 <- deriveEstimation(Y_exp.train25, "exp", sl = -1, su = 5)

  sep.sep.train25 <- deriveEstimation(Y.sep.ts.train25, "sep.exp")
  exp.sep.train25 <- deriveEstimation(Y_exp.train25, "sep.exp")

  # ttest = rownames(Y.sep.ts.train25) %>% as.numeric()
  t <- (rownames(Y.sep.ts.train25) %>% as.numeric())
  corr.sep25 <- sapply(1:100, function(i){
    Y = Y.sep.ts.train25[,i]
    ttest = rep(seq(min(Y), max(Y), length.out = ntest),2)
    params = sep.mtsgp.train25[i,]
    yhat <- interpolation_sim(t,Y,params,ttest,prediction = FALSE)
    res.mtsgp <- cor(yhat[1:ntest],yhat[(ntest+1):(2*ntest)])
    # res.mtsgp.a = yhat[2*ntest+1]
    # fit splines
    ntrain = n/2
    t1 <- t[1:ntrain]
    t2 <- t[(ntrain+1):(2*ntrain)]
    fit = lm(Y[1:ntrain] ~ ns(t1, knots=3))
    y1hat = predict(fit,data.frame(t1=ttest[1:ntest]))
    fit2 = lm(Y[(ntrain+1):(2*ntrain)] ~ ns(t2, knots=3))
    y2hat = predict(fit2,data.frame(t2=ttest[(ntest+1):(2*ntest)]))
    res.splines = cor(y1hat,y2hat)

    # plot(t1,Y[1:n],col="grey",xlab="Age",ylab="Wages")
    # points(ttest[1:n],yhat[1:n],col="darkgreen")
    # points(ttest[1:n],Y.sep.ts.test[1:n,i])

    # fit sep GP
    params = sep.sep.train25[i,]
    yhat <- interpolation_sim(t,Y,params,ttest)
    res.sep <- cor(yhat[1:ntest],yhat[(ntest+1):(2*ntest)])
    return(c(res.mtsgp, res.splines,res.sep))
  })

  t <- (rownames(Y_exp.train25) %>% as.numeric())
  corr.exp25 <- sapply(1:100, function(i){
    Y = Y_exp.train25[,i]
    ttest = rep(seq(min(Y), max(Y), length.out = ntest),2)
    params = exp.mtsgp.train25[i,]
    yhat <- interpolation_sim(t,Y,params,ttest,prediction = FALSE)
    res.mtsgp <- cor(yhat[1:ntest],yhat[(ntest+1):(2*ntest)])


    # fit splines
    ntrain = n/2
    t1 <- t[1:ntrain]
    t2 <- t[(ntrain+1):(2*ntrain)]
    fit = lm(Y[1:ntrain] ~ ns(t1, knots=3))
    y1hat = predict(fit,data.frame(t1=ttest[1:ntest]))
    fit2 = lm(Y[(ntrain+1):(2*ntrain)] ~ ns(t2, knots=3))
    y2hat = predict(fit2,data.frame(t2=ttest[(ntest+1):(2*ntest)]))
    res.splines = cor(y1hat,y2hat)

    # plot(t1,Y[1:n],col="grey",xlab="Age",ylab="Wages")
    # points(ttest[1:n],y1hat,col="darkgreen")
    # points(ttest[1:n],Y_exp.test[1:n,i])

    # fit sep GP
    params = exp.sep.train25[i,]
    yhat <- interpolation_sim(t,Y,params,ttest)
    res.sep <- cor(yhat[1:ntest],yhat[(ntest+1):(2*ntest)])
    return(c(res.mtsgp,res.splines,res.sep))
  })

  dat3.1 <- t(corr.exp25) %>% as.data.frame() %>% `colnames<-`(c("MTSGP","splines","sep")) %>%
    tidyr::pivot_longer(everything(),names_to = "model")

  dat3.2 <- t(corr.sep25) %>% as.data.frame() %>% `colnames<-`(c("MTSGP","splines","sep")) %>%
    tidyr::pivot_longer(everything(),names_to = "model")
  dat <- rbind(dat3.1,dat3.2)%>% mutate(generation=c(rep("mtsgp",300),rep("gp",300)))

  dat %>%
    ggplot(aes(x=generation,y=value, fill = model))+
    geom_boxplot() +
    theme(panel.background = element_rect(fill = "white", colour = "grey20")) +
    xlab("Model from which data is generated") +
    ylab("AUC") +
    scale_x_discrete(breaks=c("gp","mtsgp"),
                     labels=c("seperate GPs \n s=2","MTSGP \n s = 2 "))

  dat$generation<-as.factor(dat$generation)
  res.cor <- dat %>%
    group_by(model) %>%
    roc_auc(truth = generation, value, event_level = "second")

  # compare with cross-correlation
  m = nrow(Y.sep.ts.train25)/2
  sep.cc <- sapply(1:ncol(Y.sep.ts.train25), function(i){
    cc <- ccf(Y.sep.ts.train25[1:m,i],Y.sep.ts.train25[(m+1):(2*m),i],plot=FALSE)
    max(cc$acf)
  })

  mtsgp.cc <- sapply(1:ncol(Y_exp.train25), function(i){
    cc <- ccf(Y_exp.train25[1:m,i],Y_exp.train25[(m+1):(2*m),i],plot=FALSE)
    max(cc$acf)
  })

  dat2 <- data.frame(mtsgp = mtsgp.cc, gp = sep.cc) %>%
    tidyr::pivot_longer(everything(),names_to = "generation")

  cccor <- roc(dat2$generation ~ dat2$value, plot = TRUE) # 0.5579

  dat1 <- data.frame(mtsgp = exp.mtsgp.train25[,2], gp = sep.mtsgp.train25[,2]) %>%
    tidyr::pivot_longer(everything(),names_to = "generation")
  mtsgp.a <- roc(dat1$generation ~ log2(dat1$value), plot = TRUE) # 0.673

  return(c(res.cor$.estimate,cccor$auc, mtsgp.a$auc))
}, cl= 8)

write_csv(res.cor25 %>% as.data.frame(),file="data/res.cor25.csv")
res.cor25.2 = res.cor25
res.cor25.2[which(res.cor25<0.5)]<-1-res.cor25[which(res.cor25<0.5)]
rowMeans(res.cor25.2[,])
rowSds(res.cor25.2[,])

t1 = rownames(Y_exp.train25) %>% as.numeric()
dat.exp<-Y_exp.train25[,1:100] %>% as.data.frame() %>% `colnames<-`(paste0("ts",1:100)) %>%
  mutate(feature=rep(1:2,each=25), t= c(t1[1:25],t1[1:25]+2)) %>%
  pivot_longer(ts1:ts100,names_to = "sample")

ggplot(aes(x=t,y=value, col= sample))+
  geom_smooth(se=FALSE)+
  facet_wrap(~feature)+
  theme(panel.background = element_rect(fill = "white", colour = "grey20"),
        legend.position ="none")
t2 = rownames(Y.sep.ts.train25) %>% as.numeric()
dat.sep<-Y.sep.ts.train25[,1:100] %>% as.data.frame() %>% `colnames<-`(paste0("ts",1:100)) %>%
  mutate(feature=rep(1:2,each=25), t= c(t2[1:25],t2[1:25]+2)) %>%
  pivot_longer(ts1:ts100,names_to = "sample")

dat <- rbind(dat.sep,dat.exp)%>% mutate(generation=c(rep("SGP",5000),rep("MTSGP",5000))) %>%
  filter(sample == "ts10")
dat$feature <- as.factor(dat$feature)
ggplot(dat,aes(x=t,y=value, col= feature))+
  geom_smooth()+
  facet_grid(~generation)+
  theme(panel.background = element_rect(fill = "white", colour = "grey20"),
        legend.position ="none")
index = order(t1[1:25])
y=Y_exp.train25[1:25,1:3]
y2=Y_exp.train25[26:50,1:3]
library("wesanderson")
matplot(t1[1:25][index],y[index,],pch = 16, col=wes_palette("Darjeeling1"), lty=1, xlab="",type = "l", ylab="",ylim = c(-8,8),lwd=3)
lines(t1[1:25][index],y2[index,1],pch = 16, col=wes_palette("Darjeeling1")[1],lwd=3)

points(t1,yintro[(cumsum.index[1]+1):(cumsum.index[2])],pch=18, col = wes_palette("Darjeeling1")[2],cex=2)
lines(t1,yintro[(cumsum.index[1]+1):(cumsum.index[2])],lty=2,col = wes_palette("Darjeeling1")[2],lwd=3)
points(t1,yintro[(cumsum.index[2]+1):(cumsum.index[3])],pch = 8, col= wes_palette("Darjeeling1")[3],cex=2)
lines(t1,yintro[(cumsum.index[2]+1):(cumsum.index[3])],lty=3, col = wes_palette("Darjeeling1")[3],lwd=3)
