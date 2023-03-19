library(plgp)
library(mvtnorm)
library(dplyr)
library(MASS)
library(pracma)
library(Matrix)
library(condmixt)
library(ggplot2)
source("script/TMGGP_functions.R")
## Simulation 3
n <- 50
X <- matrix(seq(0, 2*pi, length=n), ncol=1)
X <- matrix(seq(-2, 2, length=n), ncol=1)
k=0.01
y1 <- atan(k*X)/atan(k)
plot(X,y1)
k=1
y2 <- atan(k*sin(X))/atan(k)
k=0.01
y2 <- atan(k*sin(X))/atan(k)

dat <- data.frame(x= X, y1= y1, y2 = y2)
ggplot(dat, aes(x = X))+
  geom_point(aes(y=y1,col = "blue"))+
  geom_point(aes(y=y2,col = "red"))+
  geom_line(aes(y=y1,col = "blue"))+
  geom_line(aes(y=y2,col = "red"))
# geom_line(aes(y=y3,col = "black"))
delta = Matrix(1-bdiag(replicate(2,matrix(1,n,n),simplify=FALSE)),sparse = TRUE)
group.index <- rep(c(1,2),each=n)
group.size <- rle(group.index)
out <- optim(c(1,1), nl_mggp, method="L-BFGS-B",lower=c(-10,-10),
             upper=c(20,Inf),
             t= c(X,X), Y=c(y1,y2), delta=delta,group.size = group.size)
bhat = softplus(out$par[1])
ahat = softplus(out$par[2])

k_across = rep(c(0.01,1,10),3)
s_across = rep(c(0,0.5,1),each=3)

yother <- sapply(1:length(s_across), function(i){
  k = k_across[i]
  s = s_across[i]
  # y2 <- atan(k*sin(X+s))/atan(k)
  atan(k*(X+s))/atan(k)
})


dat1 <- cbind(X, yother) %>% as.data.frame() %>% `colnames<-`(c("t",paste0("ts",1:9))) %>%
  tidyr::pivot_longer(!t, names_to = "groups") %>%
  mutate(k = rep(k_across %>% as.factor(), n),
         s = rep(s_across %>% as.factor(), n))
dat1$label <- sprintf(
  "k = %s, s = %s",
  dat1$k,
  dat1$s
)
mypal = pal_npg("nrc", alpha = 0.8)(9)
mypal = pal_d3("category10")(9)
p1 <- ggplot() +
  geom_line(data = dat1,aes(x = t, y = value,col = groups),linewidth = 2)+
  geom_line(data = data.frame(t = X, target = y1), aes(x = t,y = target), linewidth = 2)+
  facet_wrap(~groups) +
  scale_color_manual(name = "Groups", values = mypal) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.position="none")+
  geom_text(
    size    = 4,
    data    = dat1,
    mapping = aes(x = Inf, y = Inf, label = label),
    hjust   = 1.0,
    vjust   = 1.5
  )
# xlab("Kernel from which data is generated")+
# ylab("Fitted value of a")

Xtest = matrix(seq(-2, 2, length=200), ncol=1)
a_across<- sapply(1:length(s_across), function(i){
  y2 <- yother[,i]
  out <- optim(c(1,1,0), nl_exp, method="L-BFGS-B",lower=c(-10,-10,-1),
               upper=c(20,Inf,4),
               t= c(X,X), Y=c(y1,y2), delta=delta,group.size = group.size)
  # params = c(softplus(out$par[1:2]),out$par[3])
  # yhat <- interpolation_sim(t = c(X,X),Y=c(y1,y2), params,ttest = c(Xtest, Xtest))
  # cor = cor(yhat[1:200,1],yhat[201:400,1])
  return(c(softplus(out$par[2]),out$par[3]))
})

dat <- data.frame(k= k_across %>% as.factor(), s =s_across %>% as.factor(), a = a_across[1,], shat = a_across[2,],
                  groups = paste0("ts",1:9))

p2 <- ggplot(dat, aes(x= k,y=a, col=groups, shape=s))+
  geom_point(size = 5) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.position = c(.05, .99),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))+
  scale_color_manual(name = "Groups", values = mypal, guide = "none")
p1+p2

jpeg(file="plots/simulation3.jpg",width = 10, height = 6,units = "in",res=450)
p1 + p2 + plot_annotation(tag_levels = 'A') &
  theme(text = element_text(size = 18),
        plot.tag = element_text(size = 14)
  )
dev.off()

write.csv(yother, file = "cluster.csv",row.names = F)
write.csv(y1, file = "target.csv",row.names = F)


