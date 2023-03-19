library(plgp)
library(mvtnorm)
library(dplyr)
library(MASS)
library(pracma)
library(Matrix)
library(condmixt)
library(ggplot2)
source("script/TMGGP_functions.R")
## simulate t with shift
n <- 20
## set truth b=0.3, a=1,
tau2 <- 4 # scale =2^2, so the 95% amplitade within -4-4
b = 0.3
a = 1
s <- c(0,2)
## Simulate data for different kernels
Y_exp = simulateData(b,a,s,tau2, n, kernel = "exp")
Y_rbf = simulateData(b,a,s,tau2, n, kernel = "rbf")
Y_matern = simulateData(b,a,s,tau2, n, kernel = "matern")

### 1. same group covariance
eps <- sqrt(.Machine$double.eps)
# microbenchmark("plgp"={sqrt(plgp::distance(tinput))}, "dist" = {dist(tinput, diag = T,upper = T)})
i=1
e <- eigen(Vk)
e <- eigen(sigma)
e$values

res <- list()
res[["rbf_rbf"]] <- deriveEstimation(Y_rbf, nl_rbf, "rbf")
res[["exp_exp"]] <- deriveEstimation(Y_exp, nl_exp, "exp")
res[["matern_matern"]] <- deriveEstimation(Y_matern, nl_matern, "matern")
res_filter <- lapply(1:3,function(i){
  res[[i]][res[[i]][,1] != softplus(-10) & res[[i]][,2] < 10,]
})
kernel = c("rbf","exp","matern3/2")
df <- do.call("rbind", res_filter) %>% as.data.frame() %>%
  `colnames<-`(c("bhat","ahat","shat","tau2", "ll")) %>%
  mutate(n = 20,
         kernel = rep(kernel, unname(unlist(lapply(res_filter,nrow)))))
# generation = rep(c("rbf","exp","matern3/2"),each = 3)


## larger sample size
n = 10
Y_rbf10 = simulateData(b,a,s,tau2, n, kernel = "rbf")
Y_exp10 = simulateData(b,a,s,tau2, n, kernel = "exp")
Y_matern10 = simulateData(b,a,s,tau2, n, kernel = "matern")
res10 <- list()
res10[["rbf_rbf"]] <- deriveEstimation(Y_rbf10, nl_rbf, "rbf")
res10[["exp_exp"]] <- deriveEstimation(Y_exp10, nl_exp, "exp")
res10[["matern_matern"]] <- deriveEstimation(Y_matern10, nl_matern, "matern")
res_filter10 <- lapply(1:3,function(i){
  res10[[i]][res10[[i]][,1] != softplus(-10) & res10[[i]][,2] < 10 & res10[[i]][,2] != softplus(-10),]
})

n = 50
Y_rbf50 = simulateData(b,a,s,tau2, n, kernel = "rbf")
Y_exp50 = simulateData(b,a,s,tau2, n, kernel = "exp")
Y_matern50 = simulateData(b,a,s,tau2, n, kernel = "matern")
res50 <- list()
res50[["rbf_rbf"]] <- deriveEstimation(Y_rbf50, nl_rbf, "rbf")
res50[["exp_exp"]] <- deriveEstimation(Y_exp50, nl_exp, "exp")
res50[["matern_matern"]] <- deriveEstimation(Y_matern50, nl_matern, "matern")
res_filter50 <- lapply(1:3,function(i){
  res50[[i]][res50[[i]][,1] != softplus(-10) & res50[[i]][,2] < 10,]
})

n = 100
Y_exp100 = simulateData(b,a,s,tau2, n, kernel = "exp")
Y_matern100 = simulateData(b,a,s,tau2, n, kernel = "matern")
res100 <- list()
res100[["exp_exp"]] <- deriveEstimation(Y_exp100, nl_exp, "exp")
res100[["matern_matern"]] <- deriveEstimation(Y_matern100, nl_matern, "matern")
res_filter100 <- lapply(1:3,function(i){
  res100[[i]][res100[[i]][,1] != softplus(-10) & res100[[i]][,2] < 10,]
})

n = 200
Y_exp200 = simulateData(b,a,s,tau2, n, kernel = "exp")
res200 <- list()
res200[["exp_exp"]] <- deriveEstimation(Y_exp200, nl_exp, "exp")
res_filter200 <- lapply(1:1,function(i){
  res200[[i]][res200[[i]][,1] != softplus(-10) & res200[[i]][,2] < 10,]
})

save(res10, res,res50,res100,res200, file = "data/simulation.rda")
df10 <- do.call("rbind", res_filter10) %>% as.data.frame() %>%
  `colnames<-`(c("bhat","ahat","shat","tau2", "ll")) %>%
  mutate(n = 10,
         kernel = rep(kernel, unname(unlist(lapply(res_filter10,nrow)))))

df50 <- do.call("rbind", res_filter50) %>% as.data.frame() %>%
  `colnames<-`(c("bhat","ahat","shat","tau2", "ll")) %>%
  mutate(n = 50,
         kernel = rep(kernel[c(2,3,1)], unname(unlist(lapply(res_filter50,nrow)))))

df100 <- do.call("rbind", res_filter100) %>% as.data.frame() %>%
  `colnames<-`(c("bhat","ahat","shat","tau2", "ll")) %>%
  mutate(n = 100,
         kernel = rep(kernel, unname(unlist(lapply(res_filter100,nrow)))))

df200 <- do.call("rbind", res_filter200) %>% as.data.frame() %>%
  `colnames<-`(c("bhat","ahat","shat","tau2", "ll")) %>%
  mutate(n = 200,
         kernel = rep(kernel[2], unname(unlist(lapply(res_filter200,nrow)))))
dat <- rbind(df,df50, df100,df200)
dat$n <- as.factor(dat$n)

library(ggsci)
library(patchwork)
library(wesanderson)
is.outlier <- function (x) {
  x < quantile(x, .25) - 1.5 * IQR(x) |
    x > quantile(x, .75) + 1.5 * IQR(x)
}
df100 %>% group_by(kernel) %>%
  mutate(outlier.p = is.outlier(ahat),outlier.s = is.outlier(shat)) %>%
  ungroup() -> df100.out
p1 = ggplot(df100.out,aes(y = ahat, x= kernel, fill = kernel))+
  geom_boxplot(outlier.size = 0.5)+
  geom_point(data = df100.out[df100.out$outlier.p,], aes(col = kernel))+
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1, col = wes_palette("Rushmore1")[5]) +
  scale_fill_manual(values = wes_palette("Rushmore1"))+
  scale_colour_manual(values = wes_palette("Rushmore1"))+
  # scale_fill_startrek() +
  theme(panel.background = element_rect(fill = "white", colour = "grey20"),
        legend.position="none") +
  xlab("Kernel from which data is generated")+
  ylab("Fitted value of a") +
  ylim(-1,3)
# scale_fill_simpsons()
p2 = ggplot(df100.out,aes(y = shat, x= kernel, fill = kernel))+
  geom_boxplot(outlier.size = 0.5)+
  geom_point(data = df100.out[df100.out$outlier.s,], aes(col = kernel))+
  geom_hline(yintercept = 2, linetype = "dashed", linewidth = 1, col = wes_palette("Rushmore1")[5]) +
  scale_fill_manual(values = wes_palette("Rushmore1"))+
  scale_colour_manual(values = wes_palette("Rushmore1"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey20"),
        legend.position="none") +
  xlab("Kernel from which data is generated")+
  ylab("Fitted value of s ")

df_matern <- df100 %>% filter(kernel != "rbf") %>%
  mutate(multipllicity = tau2* (bhat^ifelse(kernel=="exp",1,3)))%>% group_by(kernel) %>%
  mutate(outlier.p = is.outlier(multipllicity)) %>%
  ungroup()


p3 <- ggplot(df_matern,aes(y = multipllicity, x= kernel, fill = kernel))+
  geom_boxplot(outlier.size = 0.5)+
  geom_point(data = df_matern[df_matern$outlier.p,], aes(col = kernel))+
  # geom_segment(aes(x="exp",xend="matern3/2",y=4*(0.3^3),yend=4*(0.3^3)),linetype = "dashed", linewidth = 1,col = wes_palette("Rushmore1")[5])
  geom_hline(yintercept = 4*(0.3^3), linetype = "dashed", linewidth = 1, col = wes_palette("Rushmore1")[5]) +
  geom_hline(yintercept = 4*(0.3^1), linetype = "dashed", linewidth = 1, col = wes_palette("Rushmore1")[5]) +
  scale_fill_manual(values = wes_palette("Rushmore1"))+
  scale_colour_manual(values = wes_palette("Rushmore1"))+
  ylim(0,2.5)+
  theme(panel.background = element_rect(fill = "white", colour = "grey20"),
        legend.position="none") +
  xlab("Kernel from which data is generated") +
  ylab(TeX("$\\sigma^2b^{2\\nu}$"))

library(ggsci)
library(patchwork)
library(latex2exp)
dat_exp <- dat %>% filter(kernel == "exp") %>%
  dplyr::select(ahat,shat,n) %>%
  tidyr::pivot_longer(ahat:shat,names_to ="paramater")%>% group_by(paramater,n) %>%
  mutate(outlier.p = is.outlier(value)) %>%
  ungroup()
p4 <- ggplot(dat_exp,aes(y = value, x= paramater, fill = n))+
  geom_boxplot(outlier.size = 0.5)+
  # geom_point(data = dat_exp[dat_exp$outlier.p,], aes(fill = n, col = n))+
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1, col = wes_palette("Rushmore1")[5]) +
  geom_hline(yintercept = 2, linetype = "dashed", linewidth = 1, col =  wes_palette("Rushmore1")[5]) +
  scale_fill_startrek() +
  scale_color_startrek() +
  theme(panel.background = element_rect(fill = "white", colour = "grey10"),
        legend.position = c(.95, 0.05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) +
  xlab("Paramater")+
  ylab("Fitted value")+
  ylim(-1,4)
# scale_fill_simpsons()

jpeg(file="plots/simulation_across_kernels.jpg",width = 12.5, height = 10,units = "in",res=450)
(p1 + p2) / (p3 +p4) + plot_annotation(tag_levels = 'A') &
  theme(text = element_text(size = 24),
        plot.tag = element_text(size = 16)
  )
dev.off()
