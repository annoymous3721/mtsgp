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
## Simulation 2 - evaluate a is close to simularity
n = 50
## set truth b=0.3, a=1,
tau2 <- 4 # scale =2^2, so the 95% amplitade within -4-4
b = 1
a = 1

## data generated from seperate kernel with no time shift
set.seed(123)
Y_sep = simulateData(b =1,a,s = c(0,0), tau2, n, kernel = "seperate",timemax = n)

ll.sep.sep.exp = deriveEstimation(Y_sep, "sep.exp")
summary(ll.sep.sep.exp[ll.sep.sep.exp<0]) # -108.55

# perform MTSE
sep.mtsgp <- deriveEstimation(Y_sep, "exp") #
sep.mtsgp.filter <- sep.mtsgp[sep.mtsgp[,1] > 0.1 ,] %>% as.data.frame() %>% `colnames<-`(c("b","a","s","tau2","ll"))
summary(sep.mtsgp.filter[,2])#  check estimation of a
boxplot(log10(sep.mtsgp.filter[,2]))

boxplot(sep.mtsgp.filter[,4])
summary(sep.mtsgp.filter[,5]) # check estimation of ll
boxplot(sep.mtsgp.filter[,5])
ggplot(sep.mtsgp.filter,aes(ll,log10(a)))+
  geom_point()

sep.mggp <- deriveEstimation(Y_sep, "mggp.exp") # estimate t
sep.mggp.filter <- sep.mggp[sep.mggp[,1] > 0.1 ,]
boxplot(sep.mggp.filter[,4])
summary(sep.mggp.filter[,4])

## data generated from seperate kernel with time shift=2
Y.sep.ts = simulateData(b,a,s = c(0,2), tau2, n, kernel = "seperate", timemax = n)
plotSimulate(Y.sep.ts,i = 1, s = c(0,2)) # visualize one example
ll.sep.sep.ts = deriveEstimation(Y.sep.ts, "sep.exp")
summary(ll.sep.sep.ts[ll.sep.sep.ts<0])

sep.mtsgp.ts <- deriveEstimation(Y.sep.ts, "exp") # estimate t
sep.mtsgp.ts.filter <- sep.mtsgp.ts[sep.mtsgp.ts[,1] >0.1,]
summary(sep.mtsgp.ts.filter[,2]) # check estimation of a
summary(sep.mtsgp.ts.filter[,5]) # check estimation of ll # -61
boxplot(summary(sep.mtsgp.ts.filter[,5]) )

sep.mggp.ts <- deriveEstimation(Y.sep.ts, "mggp.exp") # estimate t
sep.mggp.ts.filter <- sep.mggp.ts[sep.mggp.ts[,1] >0.1,]
summary(sep.mggp.ts.filter[,2]) # check estimation of a
# summary(sep.mggp.ts.filter[,5]) # check estimation of ll # -61
# boxplot(summary(sep.mggp.ts.filter[,5]) )

## data generated from mggp kernel
Y_mggp = simulateData(b,a = 0.3,s = c(0,0), tau2, n, kernel = "mggp",timemax = n)
# return the likelihood of seperated GP
ll.mggp.sep = deriveEstimation(Y_mggp, "sep.exp")
boxplot(ll.exp.sep)
summary(ll.exp.sep[ll.exp.sep<0]) # -109

start <- Sys.time()
# Y_mggp2 = simulateData(b,a = 0.3,s = c(0,0), tau2, n, kernel = "mggp.rbf")
mggp.mggp = deriveEstimation(Y_mggp, "mggp.exp")
end <- Sys.time()
spend <- end - start
mggp.mggp.filter <- mggp.mggp[mggp.mggp[,1] >0.1,]
boxplot(mggp.mggp.filter[,4])
summary(mggp.mggp.filter[,4])

mggp.mtsgp <- deriveEstimation(Y_mggp, "exp", sl = -1, su = 2) # estimate t
mggp.mtsgp.filter <- mggp.mtsgp[mggp.mtsgp[,1] <3,]
summary(mggp.mtsgp.filter[,2])
boxplot(mggp.mtsgp.filter[,3])
summary(mggp.mtsgp.filter[,3]) # median 0.09
summary(mggp.mtsgp.filter[,5])
boxplot(mggp.mtsgp.filter[,5])

## data generated from mtsgp kernel
Y_exp = simulateData(b = 1,a = 0.3,s = c(0,2), tau2, n, kernel = "exp", timemax = n)
plotSimulate(Y_exp,i=1,s= c(0,2))
# Y_rbf = simulateData(b,a,s = c(0,8), tau2, n, kernel = "rbf")

# return the likelihood of seperated GP
ll.exp.sep = deriveEstimation(Y_exp, "sep.exp")
boxplot(ll.exp.sep)
summary(ll.exp.sep[ll.exp.sep<0]) # -109

# perform mggp - return b,a,tau2, ll in each column
exp.mggp = deriveEstimation(Y_exp, "mggp.exp")
exp.mggp.filter <- exp.mggp[exp.mggp[,1] > 0.1,]
boxplot(exp.mggp.filter[,4])
summary(exp.mggp.filter[,4])

# perform MTSE - return b,a,s, tau2, ll in each column
# sl represent lower bound, su represent upper bound
exp.mtsgp <- deriveEstimation(Y_exp, "exp", sl = -1, su = 2)
exp.mtsgp.filter <- exp.mtsgp[exp.mtsgp[,1] > 0.1,]
summary(exp.mtsgp.filter[,2]) # check estimation of a
boxplot(log10(exp.mtsgp.filter[,2]))
boxplot(exp.mtsgp.filter[,3]) # check estimation of s
summary(exp.mtsgp.filter[,3])
summary(exp.mtsgp.filter[,5]) # check estimation of ll


generation = c("gp.ts0","gp.ts2","mggp","mtsgp.ts2")


wes_palette("Royal1")

df5.1 <- data.frame(SGP = c(ll.sep.sep.exp[,3], ll.sep.sep.ts[,3], ll.mggp.sep[,3], ll.exp.sep[,3]),
                    MGGP = c(sep.mggp[,4], sep.mggp.ts[,4], mggp.mggp[,4],exp.mggp[,4]),
                    MTSGP = c(sep.mtsgp[,5], sep.mtsgp.ts[,5], mggp.mtsgp[,5], exp.mtsgp[,5]),
                    generation = c(rep("gp.ts0",100),rep("gp.ts2",100),rep("mggp",100), rep("mtsgp.ts2", 100))) %>%
  tidyr::pivot_longer(SGP:MTSGP,names_to = "model") %>%
  group_by(generation,model) %>%
  mutate(outlier.p = is.outlier(value)) %>%
  ungroup()
supp.labs <- c("seperate GPs \n s=0", "seperate GPs \n s=2", "MGGP \n s = 0", "MTSGP \n s = 2 ")
names(supp.labs) <- c("gp.ts0", "gp.ts2", "mggp", "mtsgp.ts2")
p5 <- ggplot(df5.1,aes(x = model, y = value, fill = model))+
  geom_boxplot()+
  geom_point(data = df5.1[df5.1$outlier.p,],aes( col = model))+
  facet_wrap(~generation, labeller = labeller(generation = supp.labs))+
  scale_alpha(0.3)+
  # scale_fill_npg()+
  scale_fill_manual(name = "Models",values = wes_palette("Royal1"), labels = c("Seperate GPs", "MGGP", "MTSGP"))+
  scale_colour_manual(values = wes_palette("Royal1"))+
  # scale_fill_manual(name = "Models", values = mypal[c(1,3,4)], labels = c("Seperate GPs", "MGGP", "MTSGP"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey20"),
        legend.position = "none") +
  ylab(TeX("$\\log~p(Y)$"))
# xlab("Model from which data is generated") +
#
#   scale_x_discrete(breaks=model,
#                    labels=c("seperate GPs", "MGGP", "MTSGP"))

jpeg(file="plots/simulation2.jpg",width = 10, height = 7,units = "in",res=450)
p5
dev.off()

## Figure 3B
## Prediction tasks on irregualr data
set.seed(666)
train.index = sample(1:50,25)
Y.sep.ts.train <- Y.sep.ts[c(train.index,train.index+50),]
plotSimulate(Y.sep.ts.train,i=1,s= c(0,2))
Y.mtsgp.train <- Y_exp[c(train.index,train.index+50),]

# perform MTSE - return b,a,s, tau2, ll in each column
sgp.mtsgp.train <- deriveEstimation(Y.sep.ts.train, "exp", sl = -1, su = 5)
mtsgp.mtsgp.train <- deriveEstimation(Y.mtsgp.train, "exp", sl = -1, su = 5)

# perform SGP
sgp.sgp.train <- deriveEstimation(Y.sep.ts.train, "sep.exp")
mtsgp.sgp.train <- deriveEstimation(Y.mtsgp.train, "sep.exp")

# derive prediction
Y.mtsgp.test <- Y_exp[-c(train.index,train.index+50),]
t <- (rownames(Y.mtsgp.train) %>% as.numeric())
ttest = (rownames(Y.mtsgp.test) %>% as.numeric())
ntest = length(ttest)/2

mse.mtsgp <- sapply(1:100, function(i){
  Y = Y.mtsgp.train[,i]
  params = mtsgp.mtsgp.train[i,]
  yhat <- interpolation_sim(t,Y,params,ttest)
  res.mtsgp <- mean((Y.mtsgp.test[,i] - yhat)^2)

  # fit splines
  n = length(Y)/2
  t1 <- t[1:n]
  t2 <- t[(n+1):(2*n)]
  fit = lm(Y[1:n] ~ ns(t1, knots=NULL))
  ntest = length(ttest)/2
  y1hat = predict(fit,data.frame(t1=ttest[1:ntest]))
  fit2 = lm(Y[(n+1):(2*n)] ~ ns(t2, knots=NULL))
  y2hat = predict(fit2,data.frame(t2=ttest[(ntest+1):(2*ntest)]))
  res.splines = mean((Y.mtsgp.test[,i] - c(y1hat,y2hat))^2)

  # fit sep GP
  params = mtsgp.sgp.train[i,]
  yhat <- interpolation_sim(t,Y,params,ttest)
  res.sep <- mean((Y.mtsgp.test[,i] - yhat)^2)

  # fit MAGMAclust
  magma_train <- tibble(ID = rep(c("1","2"), each = n), Input = t, Output = Y)
  magma_pred1 <- tibble(ID = rep("1", n), Input = ttest[1:n], Output = unname(Y.mtsgp.test[1:n,i]))
  magma_pred2 <- tibble(ID = rep("2", n), Input = ttest[(n+1):(2*n)], Output = unname(Y.mtsgp.test[(n+1):(2*n),i]))
  model <- train_magma(data = magma_train, common_hp = T)
  pred1  <- pred_magma(data = magma_pred1,
                       trained_model = model,
                       grid_inputs = ttest[1:n])
  pred1 <- pred1[order(pred1$Input),]
  y1hat = pred1$Mean
  pred2  <- pred_magma(data = magma_pred1,
                       trained_model = model,
                       grid_inputs = ttest[(n+1):(2*n)])
  pred2 <- pred2[order(pred2$Input),]
  y2hat = pred2$Mean
  res.magma = mean((Y.mtsgp.test[,i] - c(y1hat,y2hat))^2)

  return(c(res.mtsgp,res.splines,res.sep, res.magma))
})

Y.sep.ts.test <- Y.sep.ts[-c(train.index,train.index+50),]
t <- (rownames(Y.sep.ts.train) %>% as.numeric())
ttest = (rownames(Y.sep.ts) %>% as.numeric())[-c(train.index,train.index+50)]
mse.sgp <- sapply(1:100, function(i){
  Y = Y.sep.ts.train[,i]
  params = sgp.mtsgp.train[i,]
  # params = c(1,0.3,2,4,-49)
  yhat <- interpolation_sim(t,Y,params,ttest)
  res.mtsgp <- mean((Y.sep.ts.test[,i] - yhat)^2)

  # fit splines
  n = length(Y)/2
  t1 <- t[1:n]
  t2 <- t[(n+1):(2*n)]
  fit = lm(Y[1:n] ~ ns(t1, knots=NULL))
  y1hat = predict(fit,data.frame(t1=ttest[1:ntest]))
  fit2 = lm(Y[(n+1):(2*n)] ~ ns(t2, knots=NULL))
  y2hat = predict(fit2,data.frame(t2=ttest[(ntest+1):(2*ntest)]))
  res.splines = mean((Y.sep.ts.test[,i] - c(y1hat,y2hat))^2)

  # plot(t1,Y[1:n],col="grey",xlab="Age",ylab="Wages")
  # points(ttest[1:n],yhat[1:n],col="darkgreen")
  # points(ttest[1:n],Y.sep.ts.test[1:n,i])

  # fit sep GP
  params = sgp.sgp.train[i,]
  params = c(1,1, -46)
  yhat <- interpolation_sim(t,Y,params,ttest)
  res.sep <- mean((Y.sep.ts.test[,i] - yhat)^2)

  # fit MAGMAclust
  magma_train <- tibble(ID = rep(c("1","2"), each = n), Input = t, Output = Y)
  magma_pred1 <- tibble(ID = rep("1", n), Input = ttest[1:n], Output = unname(Y.sep.ts.test[1:n,i]))
  magma_pred2 <- tibble(ID = rep("2", n), Input = ttest[(n+1):(2*n)], Output = unname(Y.sep.ts.test[(n+1):(2*n),i]))
  model <- train_magma(data = magma_train, common_hp = T)
  pred1  <- pred_magma(data = magma_pred1,
                       trained_model = model,
                       grid_inputs = ttest[1:n])
  pred1 <- pred1[order(pred1$Input),]
  y1hat = pred1$Mean
  pred2  <- pred_magma(data = magma_pred2,
                       trained_model = model,
                       grid_inputs = ttest[(n+1):(2*n)])
  pred2 <- pred2[order(pred2$Input),]
  y2hat = pred2$Mean
  res.magma = mean((Y.sep.ts.test[,i] - c(y1hat,y2hat))^2)

  return(c(res.mtsgp,res.splines,res.sep, res.magma))
})

dat3.1 <- t(mse.mtsgp) %>% as.data.frame() %>% `colnames<-`(c("MTSGP","splines","SGP", "MAGMA")) %>%
  tidyr::pivot_longer(everything(),names_to = "model")

dat3.2 <- t(mse.sgp) %>% as.data.frame() %>% `colnames<-`(c("MTSGP","splines","SGP", "MAGMA")) %>%
  tidyr::pivot_longer(everything(),names_to = "model")

df6.1 <- rbind(dat3.1,dat3.2)%>% mutate(generation=c(rep("mtsgp",400),rep("gp",400))) %>%  group_by(generation,model) %>%
  mutate(outlier.p = is.outlier(value)) %>%
  ungroup()
supp.labs <- c("seperate GPs \n s=2", "MTSGP \n s = 2 ")
names(supp.labs) <- c("gp", "mtsgp")
p6 <- ggplot(df6.1,aes(x=model,y=value, fill = model))+
  geom_boxplot() +
  geom_point(data = df6.1[df6.1$outlier.p,],aes( col = model))+
  facet_wrap(~generation, labeller = labeller(generation = supp.labs))+
  scale_alpha(0.3)+
  # scale_fill_npg()+
  scale_fill_manual(name = "Models",values = wes_palette("Royal1")[c(2,3,4,5)], labels = c("MTSGP","Seperate GPs", "Splines"))+
  scale_colour_manual(values = wes_palette("Royal1")[c(2,3,4,5)])+
  theme(panel.background = element_rect(fill = "white", colour = "grey20"),
        legend.position = "none") +
  # scale_fill_manual(name = "Models", values = mypal[c(4,1,2)], labels = c("MTSGP","Seperate GPs", "Splines"))+
  # xlab("Model from which data is generated") +
  ylab("Prediction MSE")
p6
# scale_x_discrete(breaks=c("gp","mtsgp"),
#                  labels=c("seperate GPs \n s=2","MTSGP \n s = 2 "))
df6.1 %>% group_by(generation, model) %>%
  summarise(mean = mean(value), sd = sd(value))
