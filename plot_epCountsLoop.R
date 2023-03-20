res2<-read.csv(file="LIMA_binnedEPCounts_nonstatic.rbf.estimated.param.csv")
ci <- read.csv(file="LIMA_binnedEPCounts_nonstatic.rbf.ci.csv")
load(file = "LIMA_binnedEPCounts_nonstatic.rda")

negative_index = floor(10*log10(n/2))+1
shiftcorr <- pbsapply(1:ncol(Y.shift), function(i){
  cc <- ccf(Y.shift[1:n,i],Y.shift[(n+1):(2*n),i],plot=FALSE)
  max(cc$acf[negative_index:(2*negative_index-2)])
},cl=18)

pos <- order(shiftcorr,decreasing = TRUE)
keep = shiftprior!=0
pos2 <- pos[keep[pos]]
pdf("plots/Top-correlation-pairs.pdf",         # File name
    width = 8, height = 5, # Width and height in inches
    colormodel = "cmyk")
par(mfrow = c(2, 3))
for(j in 1: 5){
  i <- pos2[j]
  bound <-quantile(Y.shift[,i],probs = c(0,1))
  matplot(t,Y.shift[1:n,i],type="l", col="gray", lty=1, xlab="t(h)", ylab="Expression",ylim = bound)
  title(main=colnames(Y.shift)[i],
        sub = paste0("cor = ", round(ci[,i],3), ";  ", "a =", round(res2[2,i],3)))

  points(t,Y.shift[1:n,i], pch=20, cex=1)

  lines(t,Y.shift[(n+1):(2*n),i],lty=2,col="green")
  points(t,Y.shift[(n+1):(2*n),i], pch=20, cex=1,col="green")
  # lines(t,c[,i]-mean(c[,i]),lty=2,col="blue")
  # points(t,c[,i]-mean(c[,i]), pch=20, cex=1,col="blue")
  # lines(t,Y.shift[(n+1):(2*n),i]-res2[4,i],lty=2,col="blue")
  # points(t,Y.shift[(n+1):(2*n),i]-res2[4,i], pch=20, cex=1,col="blue")

  lines(t+res2[3,i],Y.shift[(n+1):(2*n),i]-res2[4,i],lty=2,col="red")
  points(t+res2[3,i],Y.shift[(n+1):(2*n),i]-res2[4,i], pch=20, cex=1,col="red")

  lines(t+shiftprior[i],Y.shift[(n+1):(2*n),i],lty=2,col="orange3")
  points(t+shiftprior[i],Y.shift[(n+1):(2*n),i], pch=20, cex=1,col="orange3")

  legend("bottomright", legend=c("gene", "peak","peak w/ time & y shift"),col=c("black", "green","red","blue"), lty=c(1,2,2,2), cex=0.5)

}
dev.off()


library(nullranges)

## validate with HiC data
epCounts_filter$obs_exp_CMB_omegaMap_inter_withNorms.hic <- 1 + (epCounts_filter$LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic - epCounts_filter$expectedCounts) / (epCounts_filter$expectedCounts + 1)
hist(epCounts_filter$obs_exp_CMB_omegaMap_inter_withNorms.hic)

# dat_plot <- data.frame(ranges = start(anchors(epCounts_filter)$first),
#                        correlation = ci[1,],
#                        cor_0.9 = ifelse(ci[1,]>0.9,1,0) %>% as.factor(),
#                        lower = ci[2,],
#                        upper = ci[3,],
#                        similarity = res2[2,] %>% as.numeric())
# p<-ggplot(dat_plot, aes(x=ranges, y=correlation,color = cor_0.9)) +
#   geom_pointrange(aes(ymin=lower, ymax=upper),position = "jitter",size=0.1)+
#   xlim(36800000,37500000)
# p2 <- ggplot(dat_plot,aes(x=ranges,y=similarity,color = cor_0.9))+
#   geom_point(size=0.6)+
#   scale_y_reverse()
# library(patchwork)
# p/p2

## Use observed / expected
hic = mcols(epCounts_filter)[,183:189] %>% as.matrix()

## filter not converged pairs; when a is small, all cor variance come from noise
# when a is too small or noise is large, cor will either =1 or cor not fall in corrlow and corrhigh
index1 <- (ci[2,]<= ci[1,] & ci[3,] >= ci[1,] & res2[5,] <1e-2) %>% as.vector()
# index1 <-  (res2[5,] <1e-4 & res2[2,] > 1e-4)  %>% as.vector()
table(index1)
plotPair(Y.shift,1,t, s=c(0,res2[3,i]), mu=res2[4,i])

## second filter to high dynamic change
peak_max = apply(peak,2,max)
peak_min = apply(peak,2,min)
peak_diff = peak_max-peak_min
peak_index = peak_diff > 1
gene_max = apply(gene,2,max)
gene_min = apply(gene,2,min)
gene_diff = gene_max-gene_min
gene_index = gene_diff > 1

hic_index = epCounts_filter$LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic>100

filter2 = peak_index & gene_index & index1 & hic_index
# filter2 = peak_index & gene_index
epCounts_filter2 = epCounts_filter[filter2] # 3078

Y.shift2 <- Y.shift[,filter2]

data = data.frame(
  corr = ci[1,] %>% unlist(),
  corrlow = ci[2,] %>% unlist(),
  corrhigh = ci[3,] %>% unlist(),
  hic = epCounts_filter$obs_exp_CMB_omegaMap_inter_withNorms.hic,
  hic_max = epCounts_filter$hic_max,
  obs_hic = epCounts_filter$LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic,
  obs_hic_max = rowMaxs(hic),
  pvalue_adj = epCounts_filter$HIC_adjusted_pval,
  loop = !is.na(epCounts_filter$HIC_value),
  similarity = res2[2,] %>% unlist(),
  noise = res2[5,] %>% unlist(),
  bhat = res2[1,] %>% unlist(),
  binit = theta,
  time_shift = res2[3,] %>% unlist(),
  spatial_variance = res2[6,] %>% unlist(),
  y_shift = res2[4,] %>% unlist()
)

data2 <- data[filter2,]
ggplot(data2, aes(hic,corr)) + geom_point()
ggplot(data2, aes(hic_max,corr)) + geom_point()
ggplot(data2, aes(obs_hic,corr)) + geom_point()
ggplot(data2, aes(obs_hic_max,corr)) + geom_point()
ggplot(data2, aes(hic,similarity)) + geom_point()
data2 %>% mutate(log10a = -log10(similarity)) %>%
  ggplot(aes(obs_hic,log10a)) + geom_point()
data2 %>% mutate(log10a = -log10(similarity)) %>%
  ggplot(aes(hic,log10a)) + geom_point()

epCounts_filter2$candidate1 = ifelse(ci[1,filter2]>0.95 & (ci[2,filter2]<ci[1,filter2] & ci[3,filter2] > ci[1,filter2] & res2[5,filter2] <1e-2) %>% as.vector(), TRUE, FALSE)
epCounts_filter2$candidate1 = ifelse(ci[2,filter2]>0.8, TRUE, FALSE)

epCounts_filter2$candidate2 = ifelse(res2[2,filter2] %>% as.numeric()<0.2 & (ci[2,filter2]<ci[1,filter2] & ci[3,filter2] > ci[1,filter2] & res2[5,filter2] <1e-2) %>% as.vector(),TRUE, FALSE)
epCounts_filter2$candidate1 = ifelse(ci[1,filter2]>0.65, TRUE, FALSE)
epCounts_filter2$candidate2 = ifelse(res2[2,filter2] %>% as.numeric()<0.2,TRUE, FALSE)
table(epCounts_filter2$candidate1) # 912
table(epCounts_filter2$candidate2) # 714

set.seed(2023)
mgi <- matchRanges(focal = epCounts_filter2[epCounts_filter2$candidate1],
                   pool = epCounts_filter2[!epCounts_filter2$candidate1],
                   covar = ~ distance,
                   method = 'stratified')

mgi <- matchRanges(focal = epCounts_filter2[epCounts_filter2$candidate2],
                   pool = epCounts_filter2[!epCounts_filter2$candidate2],
                   covar = ~ distance,
                   method = 'stratified')
ov <- overview(mgi)
ov
plotPropensity(mgi, sets = c('f', 'p', 'm'))
library(patchwork)
plots <- lapply(covariates(mgi), plotCovariate, x=mgi, sets = c('f', 'm', 'p'))
Reduce('/', plots)

## Generate a randomly selected set from all binPairs
all <- c(focal(mgi), pool(mgi))
set.seed(123)
random <- all[sample(1:length(all), length(mgi), replace = FALSE)]
## Calculate the percent of convergent CTCF sites for each group
g1 <- mean(focal(mgi)$LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic)
g2 <- mean(pool(mgi)$LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic)
g3 <- mean(random$LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic)
g4 <- mean(mgi$LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic)
g1-g4


wilcox.test(focal(mgi)$LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic, mgi$LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic, paired=TRUE)
# Perm.test.stat <- pbsapply(1:1000,function(i){
#   set.seed(i+2000)
#   mgi <- matchRanges(focal = epCounts_filter2[epCounts_filter2$candidate1],
#                      pool = epCounts_filter2[!epCounts_filter2$candidate1],
#                      covar = ~ distance,
#                      method = 'stratified')
#   # all.hic = c(focal(mgi)$LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic, mgi$LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic)
#   # random = sample(all.hic, length(all.hic), replace = FALSE)
#   # Perm.res <- mean(random[1:length(mgi)]) - mean(random[(length(mgi)+1):(2*length(mgi))])
#   Perm.res <- mean(focal(mgi)$LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic) - mean(mgi$LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic)
#   return(Perm.res)
# },cl=8)
# quantile(Perm.test.stat,c(0.05)) # if <5% difference<0
# sum(Perm.test.stat<=0)/1000 # 0.8 = 4.9%; 0.7=
# sum((g1-g4)>Perm.test.stat)/10

## Visualize
barplot(height = c(g1, g2, g3, g4),
        names.arg = c('candidate\n(focal)', 'unlooped\n(pool)',
                      'all pairs\n(random)', 'uncorrelated\n(matched)'),
        col = c('#1F78B4', '#33A02C', 'orange2', '#A6CEE3'),
        ylab = "Normalized hic score",
        main = "Validation with Hi-C",
        border = NA,
        las = 1)



### use obs / exp but control for distance
g1 <- mean(focal(mgi)$obs_exp_CMB_omegaMap_inter_withNorms.hic)
g2 <- mean(pool(mgi)$obs_exp_CMB_omegaMap_inter_withNorms.hic)
g3 <- mean(random$obs_exp_CMB_omegaMap_inter_withNorms.hic)
g4 <- mean(mgi$obs_exp_CMB_omegaMap_inter_withNorms.hic)

g1-g4


p1 <- ggplot(data.frame(hic = c(g1, g2, g3, g4), names = c('candidate\n(focal)','uncorrelated\n(pool)','all pairs\n(random)','uncorrelated\n(matched)')), aes(x = names, y = hic, fill = names)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(hic,3)), vjust = 0) +
  ylab("expted /observed counts") +
  theme(panel.background = element_rect(fill = "white", colour = "grey20"),
        legend.position="none") +
  xlab("Model from which data is generated")

# p-value = 0.006441 when no noise filter cor = 0.8; p=0.04 when noise filter a<0.7
wilcox.test(focal(mgi)$obs_exp_CMB_omegaMap_inter_withNorms.hic,mgi$obs_exp_CMB_omegaMap_inter_withNorms.hic)

# removing non-converge pairs : p-value 0.0047 when noise filter cor>0.5; p=0.0096 with cor>0.65
# don't remove non-converge pairs : p-value 0.01365 when noise filter cor>0.8/0.75; p = 0.057 cor>0.9
library(ggsci)
p1 <- ggplot(data = data.frame(hic = c(epCounts_filter2[epCounts_filter2$candidate1]$obs_exp_CMB_omegaMap_inter_withNorms.hic,
                          epCounts_filter2[!epCounts_filter2$candidate1]$obs_exp_CMB_omegaMap_inter_withNorms.hic),
                  groups = c(rep('High correlation pairs',sum(epCounts_filter2$candidate1)),rep('Others',length(epCounts_filter2) - sum(epCounts_filter2$candidate1)))),
             aes(x = groups, y = hic, fill = groups)) +
  geom_violin() +
  # geom_boxplot()+
  theme_bw()+
  scale_fill_d3()+
  # geom_text(aes(label = round(hic,3)), vjust = 0) +
  ylab("expted /observed Hic")

## Visualize
barplot(height = c(g1, g2, g3, g4),
        names.arg = c('candidate\n(focal)', 'unlooped\n(pool)',
                      'all pairs\n(random)', 'uncorrelated\n(matched)'),
        col = c('#1F78B4', '#33A02C', 'orange2', '#A6CEE3'),
        ylab = "Normalized hic score",
        main = "Testing the Convergence Rule",
        border = NA,
        las = 1)

dat_plot <- data.frame(ranges = start(anchors(epCounts_filter2)$first),
                       correlation = ci[1,],
                       hic = log10(hic),
                       cor_0.99 = ifelse(ci[1,]>0.99,1,0) %>% as.factor(),
                       lower = ci[2,],
                       upper = ci[3,],
                       similarity = res2[2,] %>% as.numeric())
ggplot(dat_plot,aes(hic,y=correlation))+
  geom_point()+
  geom_smooth()
ggplot(dat_plot,aes(hic,y=similarity))+
  geom_point()+
  geom_smooth()


## Enrichment with loop # no identify
epCounts_filter2$candidate1 = ifelse(ci[2,filter2]>0.8, TRUE, FALSE)
epCounts_filter2$loop = data2$loop
epCounts_filter2$corr = data2$corrlow
epCounts_filter2$similarity = data2$similarity
epCounts_filter2$candidate2 = ifelse(res2[2,filter2] %>% as.numeric()<1,TRUE, FALSE)

table(epCounts_filter2$candidate2)
table(epCounts_filter2$candidate3)
plotPair(Y.shift2,2,t, s=c(0,(res2[3,filter2])[2] %>% as.numeric()), mu=(res2[4,filter2])[2] %>% as.numeric())

set.seed(2023)
mgi <- matchRanges(focal = epCounts_filter2[epCounts_filter2$candidate1],
                   pool = epCounts_filter2[!epCounts_filter2$candidate1],
                   covar = ~ distance,
                   method = 'stratified')
mgi <- matchRanges(focal = epCounts_filter2[epCounts_filter2$candidate3],
                   pool = epCounts_filter2[!epCounts_filter2$candidate3],
                   covar = ~ distance,
                   method = 'stratified')
mgi <- matchRanges(focal = epCounts_filter2[epCounts_filter2$candidate2],
                   pool = epCounts_filter2[!epCounts_filter2$candidate2],
                   covar = ~ obs_exp_CMB_omegaMap_inter_withNorms.hic,
                   method = 'stratified')
mgi <- matchRanges(focal = epCounts_filter2[epCounts_filter2$loop],
                   pool = epCounts_filter2[!epCounts_filter2$loop],
                   covar = ~ distance,
                   method = 'stratified')
all <- c(focal(mgi), mgi)
table <- table(all$candidate1, all$loop) #
table <- table(all$candidate2, all$loop)
table <- table(all$candidate3, all$loop)
## p-value = 0.003567 when cor > 0.9 when filter non-converged pairs;
chisq.test(table)
fisher.test(table)
wilcox.test(all[which(all$loop==TRUE)]$similarity,all[which(all$loop==FALSE)]$similarity) # p = 3.984e-5
wilcox.test(all[which(all$loop==TRUE)]$corr,all[which(all$loop==FALSE)]$corr)
p2 <- ggplot(data.frame(corr = all$corr, loop = all$loop), aes(x = loop, y = corr, fill = loop))+
  geom_violin()+
  geom_boxplot()
  ylab("correlation") +
  theme_bw()+
  scale_fill_d3()

index1 <- (ci[2,]<= ci[1,] & ci[3,] >= ci[1,] & res2[5,] <1e-2 & res2[2,] >1e-3) %>% as.vector()
filter2 = peak_index & gene_index & index1 & hic_index
data2 <- data[filter2,]
wilcox.test(data2[which(data2$loop==TRUE),"corr"],data2[which(data2$loop==FALSE),"corr"]) # 0.0001366
wilcox.test(log2(data2[which(data2$loop==TRUE),"similarity"]),log2(data2[which(data2$loop==FALSE),"similarity"])) # 0.0001797

p3 <- ggplot(data2,aes(x = loop, y = corr, fill = loop))+
  geom_violin()+
  ylab("correlation") +
  theme_bw()+
  scale_fill_d3()

p3 <- ggplot(data2,aes(x = loop, y = -log2(similarity), fill = loop))+
  geom_violin()+
  ylab("correlation") +
  theme_bw()+
  scale_fill_d3()

library(patchwork)
p1 / (p2+p3)
pdf("plots/Top-correlation-pairs.pdf",         # File name
    width = 8, height = 5, # Width and height in inches
    colormodel = "cmyk")
par(mfrow = c(2, 3))
data3 = data2[order(data2$corrlow,decreasing = TRUE),]
pair.names <- rownames(data3[which(data3$corrlow>0.99),])
pair.names2 <- str_replace(pair.names,"\\.","-")
for(j in 1: length(pair.names)){
  i <- which(colnames(dat_filter) == pair.names2[j])
  bound <-quantile(Y.shift[,i],probs = c(0,1))
  matplot(t,Y.shift[1:n,i],type="l", col="gray", lty=1, xlab="t(h)", ylab="Expression",ylim = bound)
  title(main=colnames(Y.shift)[i],
        sub = paste0("cor = ", round(data3[pair.names[j],"corrlow"],3), ";  ", " TLCC = ", round(shiftcorr[i],3), "; ", "a =", round(data3[pair.names[j],"similarity"],3), ";  "
                     # "hic = ", round(data3[pair.names[j],"hic"],3), "obs hic = ", round(data3[pair.names[j],"obs_hic"],3)
                     ))

  points(t,Y.shift[1:n,i], pch=20, cex=1)

  lines(t,Y.shift[(n+1):(2*n),i],lty=2,col="green")
  points(t,Y.shift[(n+1):(2*n),i], pch=20, cex=1,col="green")
  # lines(t,c[,i]-mean(c[,i]),lty=2,col="blue")
  # points(t,c[,i]-mean(c[,i]), pch=20, cex=1,col="blue")
  # lines(t,Y.shift[(n+1):(2*n),i]-res2[4,i],lty=2,col="blue")
  # points(t,Y.shift[(n+1):(2*n),i]-res2[4,i], pch=20, cex=1,col="blue")

  lines(t+res2[3,i],Y.shift[(n+1):(2*n),i]-res2[4,i],lty=2,col="red")
  points(t+res2[3,i],Y.shift[(n+1):(2*n),i]-res2[4,i], pch=20, cex=1,col="red")
  lines(t+shiftprior[i],Y.shift[(n+1):(2*n),i],lty=2,col="orange3")
  points(t+shiftprior[i],Y.shift[(n+1):(2*n),i], pch=20, cex=1,col="orange3")

  legend("bottomright", legend=c("gene", "peak","MTSGP", "TLCC"),col=c("black","green","red","orange3"), lty=c(1,2,2,2), cex=0.5)

}
dev.off()

pdf("plots/Top-similarity-pairs.pdf",         # File name
    width = 8, height = 5, # Width and height in inches
    colormodel = "cmyk")
par(mfrow = c(2, 3))
data3 = data2[order(data2$similarity),]
pair.names <- rownames(data3[which(data3$similarity < 0.1 & data3$time_shift != 2 & data3$y_shift != -7),])
pair.names2 <- str_replace_all(pair.names,"\\.","-")
for(j in 1: length(pair.names)){
  i <- which(colnames(dat_filter) == pair.names2[j])
  bound <-quantile(Y.shift[,i],probs = c(0,1))
  matplot(t,Y.shift[1:n,i],type="l", col="gray", lty=1,xlab="t(h)", ylab="Expression",ylim = bound)
  title(main=colnames(Y.shift)[i],
        sub = paste0("correlation = ", round(data3[pair.names[j],"corr"],3), ";  ", "a =", round(data3[pair.names[j],"similarity"],3), ";  ", "obs/exp hic = ", round(data3[pair.names[j],"hic"],3)))

  points(t,Y.shift[1:n,i], pch=20, cex=1)

  lines(t,Y.shift[(n+1):(2*n),i],lty=2,col="green")
  points(t,Y.shift[(n+1):(2*n),i], pch=20, cex=1,col="green")

  # lines(t,Y.shift[(n+1):(2*n),i]-res2[4,i],lty=2,col="blue")
  # points(t,Y.shift[(n+1):(2*n),i]-res2[4,i], pch=20, cex=1,col="blue")

  lines(t+res2[3,i],Y.shift[(n+1):(2*n),i]-res2[4,i],lty=2,col="red")
  points(t+res2[3,i],Y.shift[(n+1):(2*n),i]-res2[4,i], pch=20, cex=1,col="red")

  legend("bottomright", legend=c("gene", "peak","peak w/ time & y shift"),col=c("black", "green","red"), lty=c(1,2,2), cex=0.5)

}
dev.off()


pdf("plots/Top-hic-pairs.pdf",         # File name
    width = 8, height = 5, # Width and height in inches
    colormodel = "cmyk")
par(mfrow = c(2, 3))
data3 = data2[order(data2$hic,decreasing = TRUE),]
pair.names <- rownames(data3[which(data3$hic>5),])
pair.names2 <- str_replace(pair.names,"\\.","-")
for(j in 1: length(pair.names)){
  i <- which(colnames(dat_filter) == pair.names2[j])
  bound <-quantile(Y.shift[,i],probs = c(0,1))
  matplot(t,Y.shift[1:n,i],type="l", col="gray", lty=1,xlab="t(h)", ylab="Expression", ylim = bound)

  title(main=colnames(Y.shift)[i],
        sub = paste0("correlation = ", round(data3[pair.names[j],"corr"],3), ";  ", "a =", round(data3[pair.names[j],"similarity"],3), ";  ", "obs/exp hic = ", round(data3[pair.names[j],"hic"],3)))

  points(t,Y.shift[1:n,i], pch=20, cex=1)

  lines(t,Y.shift[(n+1):(2*n),i],lty=2,col="green")
  points(t,Y.shift[(n+1):(2*n),i], pch=20, cex=1,col="green")

  # lines(t,Y.shift[(n+1):(2*n),i]-res2[4,i],lty=2,col="blue")
  # points(t,Y.shift[(n+1):(2*n),i]-res2[4,i], pch=20, cex=1,col="blue")

  lines(t+res2[3,i],Y.shift[(n+1):(2*n),i]-res2[4,i],lty=2,col="red")
  points(t+res2[3,i],Y.shift[(n+1):(2*n),i]-res2[4,i], pch=20, cex=1,col="red")

  legend("bottomright", legend=c("gene", "peak","peak w/ time & y shift"),col=c("black", "green","red"), lty=c(1,2,2), cex=0.5)

}
dev.off()

hist(epCounts_filter2$obs_exp_CMB_omegaMap_inter_withNorms.hic,breaks=50)
hist(epCounts_filter2$LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic,breaks=50)


## enrichment with eQTL
library(vroom)
eQTL <- vroom("Alasoo_immune_response_eQTL_caQTL_coloc.tsv")
epCounts_filter3 <- epCounts_filter2[epCounts_filter2$anchor2.RNA_ensembl %in% eQTL$gene_id]

load("data/eQTL/naive_qtls.rda")
load("data/eQTL/ifng_qtls.rda")
load("data/eQTL/ciff.rda")

## Filter for very significant eQTLS (not necessary for ciff)
naive <- naive[naive$anchor1.pvalue <= 1e-08,]
ifng <- ifng[ifng$anchor1.pvalue <= 1e-08,]

## Filter out interchromosomal & less than 1e06 distances
naive <- naive[!is.na(pairdist(naive)) & pairdist(naive) <= 1e06,]
ifng <- ifng[!is.na(pairdist(ifng)) & pairdist(ifng) <= 1e06,]
ciff <- ciff[!is.na(pairdist(ciff)) & pairdist(ciff) <= 1e06,]

index1 <- (ci[2,]<= ci[1,] & ci[3,] >= ci[1,] & res2[5,] <1e-2 & res2[2,] >1e-4) %>% as.vector()
filter2 = peak_index & gene_index & index1 & hic_index
epCounts_filter2 <- epCounts_filter[filter2]
data2 <- data[filter2,]
epCounts_filter2$cd = (epCounts_filter2$obs_exp_CMB_omegaMap_inter_withNorms.hic * ci[1,filter2] )%>% unlist()
epCounts_filter2$cd2 = (epCounts_filter2$hic_max * ci[1,filter2] )%>% unlist()
peakmax = colMaxs(peak)
epCounts_filter2$abcd = (epCounts_filter2$obs_exp_CMB_omegaMap_inter_withNorms.hic * ci[1,filter2] * peakmax[filter2])%>% unlist()

epCounts_filter2$abcd2 = (epCounts_filter2$hic_max * ci[1,filter2] * peakmax[filter2])%>% unlist()
epCounts_filter2$abc = (epCounts_filter2$obs_exp_CMB_omegaMap_inter_withNorms.hic * peakmax[filter2])%>% unlist()

epCounts_filter2$candidate1 = ifelse(ci[1,filter2]>0.98, TRUE, FALSE)
epCounts_filter2$candidate2 = ifelse(epCounts_filter2$cd > quantile(epCounts_filter2$cd,0.8), TRUE, FALSE)
epCounts_filter2$candidate2.2 = ifelse(epCounts_filter2$cd2 > quantile(epCounts_filter2$cd2,0.8), TRUE, FALSE)
epCounts_filter2$candidate3 = ifelse(epCounts_filter2$abcd > quantile(epCounts_filter2$abcd,0.9), TRUE, FALSE)
epCounts_filter2$candidate3.2 = ifelse(epCounts_filter2$abcd2 > quantile(epCounts_filter2$abcd2,0.75), TRUE, FALSE)
table(epCounts_filter2$candidate1)
table(epCounts_filter2$candidate2)
table(epCounts_filter2$candidate3)
epCounts_filter2$naive = rep(0,length(epCounts_filter2))
epCounts_filter2$ifng = rep(0,length(epCounts_filter2))
epCounts_filter2$ciff = rep(0,length(epCounts_filter2))

## look at correlation in eQTL
corValidNaive <- data2$corr[countOverlaps(epCounts_filter2, naive) > 0]
corValidIfng <- data2$corr[countOverlaps(epCounts_filter2, ifng) > 0]
corValidCiff <- data2$corr[countOverlaps(epCounts_filter2, ciff) > 0]
## Visualize ABC scores among eQTL-validated subsets
plot(density(data2$corr), main = "Correlation Distributions")
plot(density(data2$corr[countOverlaps(epCounts_filter2, ifng) == 0]), main = "Correlation Distributions")
lines(density(corValidNaive), col="blue")
lines(density(corValidIfng), col="darkgreen")
lines(density(corValidCiff), col="orange3")
legend(x = "topleft",
       legend = c("All ABC Scores",
                  # "Naive eQTL-Validated ABC Scores",
                  "Ifng eQTL-Validated ABC Scores"
                  ),
       text.col = c("black", "darkgreen"), bty = 'n')
legend(x = "topright",
       legend = c("All ABC Scores",
                  # "Naive eQTL-Validated ABC Scores",
                  "Ifng eQTL-Validated ABC Scores",
                  "CRISPRi-FlowFISH Validated"),
       text.col = c("black", "blue", "darkgreen", "orange3"), bty = 'n')
wilcox.test(corValidCiff,data2$corr[countOverlaps(epCounts_filter2, ciff) == 0]) # 0.016
wilcox.test(corValidNaive,data2$corr[countOverlaps(epCounts_filter2, naive) == 0]) # 0.016
wilcox.test(corValidIfng,data2$corr[countOverlaps(epCounts_filter2, ifng) == 0]) # 0.016
plot(density(corValidIfng), main = "Correlation Distributions")
lines(density(data2$corr[countOverlaps(epCounts_filter2, ifng) == 0]), col="blue")

overpos <- findOverlaps(epCounts_filter2, ciff)
epCounts_filter2$ciff[overpos@from] = 1
fisher.test(table(epCounts_filter2$candidate1,epCounts_filter2$ciff))

overpos <- findOverlaps(epCounts_filter2, naive)
epCounts_filter2$naive[overpos@from] = 1
fisher.test(table(epCounts_filter2$candidate1,epCounts_filter2$naive))

overpos <- findOverlaps(epCounts_filter2, ifng)
epCounts_filter2$ifng[overpos@from] = 1
tab<-table(epCounts_filter2$candidate1,epCounts_filter2$ifng)
pos <- which(epCounts_filter2$candidate1==1 & epCounts_filter2$ifng==1)
order(data2$corr[pos],decreasing = TRUE)
fisher.test(tab) # 0.0007328

## look at correlation * hic in eQTLs
cdValidNaive <- epCounts_filter2$cd[countOverlaps(epCounts_filter2, naive) > 0]
cdValidIfng <- epCounts_filter2$cd[countOverlaps(epCounts_filter2, ifng) > 0]
cdValidCiff <- epCounts_filter2$cd[countOverlaps(epCounts_filter2, ciff) > 0]
## Visualize ABC scores among eQTL-validated subsets
plot(density(epCounts_filter2$cd), main = "CD score Distributions")
# lines(density(epCounts_filter2$cd[countOverlaps(epCounts_filter2, ifng) == 0]), main = "CD score Distributions")
lines(density(cdValidNaive), col="blue")
lines(density(cdValidIfng), col="darkgreen")
# lines(density(cdValidCiff), col="orange3")
legend(x = "topright",
       legend = c("All ABCD Scores",
                  "Naive eQTL-Validated ABCD Scores",
                  "Ifng eQTL-Validated ABCD Scores"
                  # "CRISPRi-FlowFISH Validated"
                  ),
       text.col = c("black", "blue", "darkgreen"), bty = 'n')
wilcox.test(cdValidCiff,epCounts_filter2$cd[countOverlaps(epCounts_filter2, ciff) == 0]) # 0.016
wilcox.test(cdValidNaive,epCounts_filter2$cd[countOverlaps(epCounts_filter2, naive) == 0]) # 0.016
wilcox.test(cdValidIfng,epCounts_filter2$cd[countOverlaps(epCounts_filter2, ifng) == 0]) # 0.0008348
plot(density(corValidIfng), main = "Correlation Distributions")
lines(density(data2$corr[countOverlaps(epCounts_filter2, ifng) == 0]), col="blue")

overpos <- findOverlaps(epCounts_filter2, ciff)
epCounts_filter2$ciff[overpos@from] = 1
fisher.test(table(epCounts_filter2$candidate2,epCounts_filter2$ciff))

overpos <- findOverlaps(epCounts_filter2, naive)
epCounts_filter2$naive[overpos@from] = 1
fisher.test(table(epCounts_filter2$candidate2,epCounts_filter2$naive)) # 0.112

overpos <- findOverlaps(epCounts_filter2, ifng)
epCounts_filter2$ifng[overpos@from] = 1
tab<-table(epCounts_filter2$candidate2,epCounts_filter2$ifng)
pos <- which(epCounts_filter2$candidate1==1 & epCounts_filter2$ifng==1)
order(data2$corr[pos],decreasing = TRUE)
fisher.test(tab) # 0.0007328
chisq.test(tab) # p = 0.0001837

## look at correlation * hicmax in eQTLs
cd2ValidNaive <- epCounts_filter2$cd2[countOverlaps(epCounts_filter2, naive) > 0]
cd2ValidIfng <- epCounts_filter2$cd2[countOverlaps(epCounts_filter2, ifng) > 0]
cd2ValidCiff <- epCounts_filter2$cd2[countOverlaps(epCounts_filter2, ciff) > 0]
## Visualize ABC scores among eQTL-validated subsets
plot(density(epCounts_filter2$cd2), main = "cd2 score Distributions")
# lines(density(epCounts_filter2$cd2[countOverlaps(epCounts_filter2, ifng) == 0]), main = "cd2 score Distributions")
lines(density(cd2ValidNaive), col="blue")
lines(density(cd2ValidIfng), col="darkgreen")
# lines(density(cd2ValidCiff), col="orange3")
legend(x = "topright",
       legend = c("All ABCD Scores",
                  "Naive eQTL-Validated ABCD Scores",
                  "Ifng eQTL-Validated ABCD Scores"
                  # "CRISPRi-FlowFISH Validated"
       ),
       text.col = c("black", "blue", "darkgreen"), bty = 'n')
wilcox.test(cd2ValidCiff,epCounts_filter2$cd2[countOverlaps(epCounts_filter2, ciff) == 0]) # 0.016
wilcox.test(cd2ValidNaive,epCounts_filter2$cd2[countOverlaps(epCounts_filter2, naive) == 0]) # 0.016
wilcox.test(cd2ValidIfng,epCounts_filter2$cd2[countOverlaps(epCounts_filter2, ifng) == 0]) # 0.0008348
plot(density(corValidIfng), main = "Correlation Distributions")
lines(density(data2$corr[countOverlaps(epCounts_filter2, ifng) == 0]), col="blue")

overpos <- findOverlaps(epCounts_filter2, ciff)
epCounts_filter2$ciff[overpos@from] = 1
fisher.test(table(epCounts_filter2$candidate2,epCounts_filter2$ciff))

overpos <- findOverlaps(epCounts_filter2, naive)
epCounts_filter2$naive[overpos@from] = 1
fisher.test(table(epCounts_filter2$candidate2,epCounts_filter2$naive)) # 0.112

overpos <- findOverlaps(epCounts_filter2, ifng)
epCounts_filter2$ifng[overpos@from] = 1
tab<-table(epCounts_filter2$candidate2.2,epCounts_filter2$ifng)
pos <- which(epCounts_filter2$candidate1==1 & epCounts_filter2$ifng==1)
order(data2$corr[pos],decreasing = TRUE)
fisher.test(tab) # 0.0007328
chisq.test(tab) # p = 0.0001837

## look at correlation * hic * peaks in eQTLs
abcdValidNaive <- epCounts_filter2$abcd[countOverlaps(epCounts_filter2, naive) > 0]
abcdValidIfng <- epCounts_filter2$abcd[countOverlaps(epCounts_filter2, ifng) > 0]
abcdValidCiff <- epCounts_filter2$abcd[countOverlaps(epCounts_filter2, ciff) > 0]
## Visualize ABC scores among eQTL-validated subsets
plot(density(epCounts_filter2$abcd), main = "Correlation Distributions")
# plot(density(epCounts_filter2$abcd[countOverlaps(epCounts_filter2, ifng) == 0]), main = "abcd score Distributions")
lines(density(abcdValidNaive), col="blue")
lines(density(abcdValidIfng), col="darkgreen")
# lines(density(abcdValidCiff), col="orange3")
legend(x = "topright",
       legend = c("All ABCD Scores",
                  "Naive eQTL-Validated ABCD Scores",
                  "Ifng eQTL-Validated ABCD Scores"
                  # "CRISPRi-FlowFISH Validated"
                  ),
       text.col = c("black", "blue", "darkgreen"), bty = 'n')
wilcox.test(abcdValidCiff,epCounts_filter2$abcd[countOverlaps(epCounts_filter2, ciff) == 0]) #
wilcox.test(abcdValidNaive,epCounts_filter2$abcd[countOverlaps(epCounts_filter2, naive) == 0]) #
wilcox.test(abcdValidIfng,epCounts_filter2$abcd[countOverlaps(epCounts_filter2, ifng) == 0]) # 0.001285

overpos <- findOverlaps(epCounts_filter2, ciff)
epCounts_filter2$ciff[overpos@from] = 1
fisher.test(table(epCounts_filter2$candidate3,epCounts_filter2$ciff))

overpos <- findOverlaps(epCounts_filter2, naive)
epCounts_filter2$naive[overpos@from] = 1
fisher.test(table(epCounts_filter2$candidate3,epCounts_filter2$naive)) # 0.112

overpos <- findOverlaps(epCounts_filter2, ifng)
epCounts_filter2$ifng[overpos@from] = 1
tab<-table(epCounts_filter2$candidate3,epCounts_filter2$ifng)
pos <- which(epCounts_filter2$candidate1==1 & epCounts_filter2$ifng==1)
order(data2$corr[pos],decreasing = TRUE)
fisher.test(tab) # 0.0007328
chisq.test(tab) # p = 0.00448

## look at correlation * hic * peaks in eQTLs
abcd2ValidNaive <- epCounts_filter2$abcd2[countOverlaps(epCounts_filter2, naive) > 0]
abcd2ValidIfng <- epCounts_filter2$abcd2[countOverlaps(epCounts_filter2, ifng) > 0]
abcd2ValidCiff <- epCounts_filter2$abcd2[countOverlaps(epCounts_filter2, ciff) > 0]
## Visualize ABC scores among eQTL-validated subsets
plot(density(epCounts_filter2$abcd2), main = "Correlation Distributions")
# plot(density(epCounts_filter2$abcd2[countOverlaps(epCounts_filter2, ifng) == 0]), main = "abcd2 score Distributions")
lines(density(abcd2ValidNaive), col="blue")
lines(density(abcd2ValidIfng), col="darkgreen")
lines(density(abcd2ValidCiff), col="orange3")
legend(x = "topright",
       legend = c("All abcd2 Scores",
                  "Naive eQTL-Validated abcd2 Scores",
                  "Ifng eQTL-Validated abcd2 Scores",
                  "CRISPRi-FlowFISH Validated"
       ),
       text.col = c("black", "blue", "darkgreen","orange2"), bty = 'n')
wilcox.test(abcd2ValidCiff,epCounts_filter2$abcd2[countOverlaps(epCounts_filter2, ciff) == 0]) #
wilcox.test(abcd2ValidNaive,epCounts_filter2$abcd2[countOverlaps(epCounts_filter2, naive) == 0]) #
wilcox.test(abcd2ValidIfng,epCounts_filter2$abcd2[countOverlaps(epCounts_filter2, ifng) == 0]) # 0.001285

overpos <- findOverlaps(epCounts_filter2, ciff)
epCounts_filter2$ciff[overpos@from] = 1
fisher.test(table(epCounts_filter2$candidate3.2,epCounts_filter2$ciff))

overpos <- findOverlaps(epCounts_filter2, naive)
epCounts_filter2$naive[overpos@from] = 1
fisher.test(table(epCounts_filter2$candidate3.2,epCounts_filter2$naive)) # 0.112

overpos <- findOverlaps(epCounts_filter2, ifng)
epCounts_filter2$ifng[overpos@from] = 1
tab<-table(epCounts_filter2$candidate3.2,epCounts_filter2$ifng)
pos <- which(epCounts_filter2$candidate1==1 & epCounts_filter2$ifng==1)
order(data2$corr[pos],decreasing = TRUE)
fisher.test(tab) # 0.0007328
chisq.test(tab) # p = 0.00448

## look at  hic * peaks in eQTLs
abcValidNaive <- epCounts_filter2$abc[countOverlaps(epCounts_filter2, naive) > 0]
abcValidIfng <- epCounts_filter2$abc[countOverlaps(epCounts_filter2, ifng) > 0]
abcValidCiff <- epCounts_filter2$abc[countOverlaps(epCounts_filter2, ciff) > 0]
## Visualize ABC scores among eQTL-validated subsets
plot(density(epCounts_filter2$abc), main = "abc Distributions")
# plot(density(epCounts_filter2$abc[countOverlaps(epCounts_filter2, ifng) == 0]), main = "abc score Distributions")
lines(density(abcValidNaive), col="blue")
lines(density(abcValidIfng), col="darkgreen")
lines(density(abcValidCiff), col="orange3")
legend(x = "topright",
       legend = c("All abc Scores",
                  "Naive eQTL-Validated abc Scores",
                  "Ifng eQTL-Validated abc Scores",
                  "CRISPRi-FlowFISH Validated"
       ),
       text.col = c("black", "blue", "darkgreen","orange2"), bty = 'n')
wilcox.test(abcValidCiff,epCounts_filter2$abc[countOverlaps(epCounts_filter2, ciff) == 0]) #
wilcox.test(abcValidNaive,epCounts_filter2$abc[countOverlaps(epCounts_filter2, naive) == 0]) #
wilcox.test(abcValidIfng,epCounts_filter2$abc[countOverlaps(epCounts_filter2, ifng) == 0]) # 0.001285

overpos <- findOverlaps(epCounts_filter2, ciff)
epCounts_filter2$ciff[overpos@from] = 1
fisher.test(table(epCounts_filter2$candidate3.2,epCounts_filter2$ciff))

overpos <- findOverlaps(epCounts_filter2, naive)
epCounts_filter2$naive[overpos@from] = 1
fisher.test(table(epCounts_filter2$candidate3.2,epCounts_filter2$naive)) # 0.112

overpos <- findOverlaps(epCounts_filter2, ifng)
epCounts_filter2$ifng[overpos@from] = 1
tab<-table(epCounts_filter2$candidate3.2,epCounts_filter2$ifng)
pos <- which(epCounts_filter2$candidate1==1 & epCounts_filter2$ifng==1)
order(data2$corr[pos],decreasing = TRUE)
fisher.test(tab) # 0.0007328
chisq.test(tab) # p = 0.00448

epindex <- unique(overpos@from)
pdf("plots/eQTL-scores.pdf",         # File name
    width = 7, height = 3, # Width and height in inches
    colormodel = "cmyk")
for (i in 1:length(epindex)) {
  index = epCounts_filter2[epindex[i]]
  gind <- epCounts_filter2[which(epCounts_filter2$anchor2.RNA_hgnc==index$anchor2.RNA_hgnc)]
  dat_plot <- data.frame(ranges = start(anchors(gind)$first),
                         abc = gind$abc,
                         cd = gind$cd,
                         abcd = gind$abcd) %>%
    tidyr::pivot_longer(abc:abcd, names_to = "score")
  p <- ggplot(dat_plot, aes(x=ranges, y=value,color = score)) +
    geom_point(size=0.8) +
    facet_wrap(~score) +
    xlab(index$anchor2.RNA_hgnc)
  print(p)
}
dev.off()


pdf("plots/Top-eQTL-pairs.pdf",         # File name
    width = 8, height = 5, # Width and height in inches
    colormodel = "cmyk")
t=c(0,0.5,1,1.5,2,4,6)
n=length(t)
jpeg(file="plots/Top-eQTL-pairs.jpg",width = 8, height = 8,units = "in",res=450)
colnames(Y.shift) <- paste(colnames(peak),colnames(gene),sep = "-")
par(mfrow = c(3, 2))
pair.names <- colnames(Y.shift2)[pos[order(data2$corr[pos],decreasing = TRUE)]]
# pair.names2 <- str_replace_all(pair.names,"\\.","-")
for(j in c(2:4,6:7,9)){
  i <- which(colnames(dat_filter) == pair.names[j])[1]
  bound <-quantile(c(Y.shift[,i],Y.shift[(n+1):(2*n),i]-res2[4,i]),probs = c(0,1))
  matplot(t,Y.shift[1:n,i],type="l", col="gray", lty=1,xlab="t(h)", ylab="Expression", ylim = bound,lwd=3)

  title(main=colnames(Y.shift)[i],
        sub = paste0("cor = ", round(ci[1,i],3), ";  ", " TLCC = ", round(shiftcorr[i],3), "; ", "a =", round(res2[2,i],3), ";  "
                      # "hic = ", round(epCounts_filter$LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic[i],3), "obs hic = ", round(epCounts_filter$obs_exp_CMB_omegaMap_inter_withNorms.hic[i],3)
        ), cex.main = 2, cex.sub = 1.5)

  points(t,Y.shift[1:n,i], pch=20, cex=2)

  lines(t,Y.shift[(n+1):(2*n),i],lty=1,col="blue",lwd = 3)
  points(t,Y.shift[(n+1):(2*n),i], pch=17, cex=2,col="blue")
  # lines(t,hic[,i]-mean(hic[,i]),lty=2,col="green")
  # points(t,hic[,i]-mean(hic[,i]), pch=20, cex=1,col="green")
  # lines(t,Y.shift[(n+1):(2*n),i]-res2[4,i],lty=2,col="blue")
  # points(t,Y.shift[(n+1):(2*n),i]-res2[4,i], pch=20, cex=1,col="blue")

  lines(t+res2[3,i],Y.shift[(n+1):(2*n),i]-res2[4,i],lty=2,col="red",lwd=3)
  points(t+res2[3,i],Y.shift[(n+1):(2*n),i]-res2[4,i], pch=21, cex=2,col="red")
  lines(t+shiftprior[i],Y.shift[(n+1):(2*n),i],lty=2,col="orange3",lwd=3)
  points(t+shiftprior[i],Y.shift[(n+1):(2*n),i], pch=24, cex=2,col="orange3")

  legend("bottomright", legend=c("Gene", "H3K27ac","MTSGP", "TLCC"),col=c("black","blue","red","orange3"), lty=c(1,1,2,2),pch=c(20,17,21,24), cex=0.8,text.font=2)

}
dev.off()


iscandidate <- epCounts_filter2$candidate1[overpos@from]
wilcox.test(data2$corr[overpos@from][!iscandidate],data2$corr[overpos@from][iscandidate], alternative = "less") # 0.003
plot(density(data2$corr[overpos@from]), main = "ABC Score Distributions")
lines(density(data2$corr[overpos@from][iscandidate]), col="blue")
plot(density(data2$corr[overpos@from][!iscandidate]), col="darkgreen")

simValidNaive <- data2$similarity[countOverlaps(epCounts_filter[filter2], naive) > 0] %>% log2()
simValidIfng <- data2$similarity[countOverlaps(epCounts_filter[filter2], ifng) > 0]%>% log2()
simValidCiff <- data2$similarity[countOverlaps(epCounts_filter[filter2], ciff) > 0]%>% log2()
## Visualize ABC scores among eQTL-validated subsets
plot(density(log2(data2$similarity)), main = "ABC Score Distributions")
lines(density(simValidNaive), col="blue")
lines(density(simValidIfng), col="darkgreen")
lines(density(simValidCiff), col="orange3")
legend(x = "topright",
       legend = c("All ABC Scores",
                  "Naive eQTL-Validated ABC Scores",
                  "Ifng eQTL-Validated ABC Scores",
                  "CRISPRi-FlowFISH Validated"),
       text.col = c("black", "blue", "darkgreen", "orange3"), bty = 'n')

wilcox.test(simValidNaive,data2$similarity[countOverlaps(epCounts_filter2, naive) == 0] %>% log2()) # 0.016

overpos <- findOverlaps(epCounts_filter[filter2], ifng)
wilcox.test(simValidIfng,data2$similarity[countOverlaps(epCounts_filter2, ifng) == 0] %>% log2()) # 0.08745
iscandidate <- epCounts_filter$candidate1[overpos@from]

wilcox.test(data2$similarity[overpos@from][!iscandidate],data2$similarity[overpos@from][iscandidate]) # 0.08745
plot(density(log2(data2$similarity[overpos@from])), main = "ABC Score Distributions")
lines(density(log2(data2$similarity[overpos@from][iscandidate])), col="blue")
lines(density(log2(data2$similarity[overpos@from][!iscandidate])), col="darkgreen")

overpos <- findOverlaps(epCounts_filter2, ciff)
iscandidate <- epCounts_filter2$candidate1[overpos@from]
wilcox.test(simValidCiff,data2$similarity[countOverlaps(epCounts_filter[filter2], ciff) == 0] %>% log2()) # 0.016

wilcox.test(data2$similarity[overpos@from][!iscandidate],data2$similarity[overpos@from][iscandidate]) # 0.08745
plot(density(log2(data2$similarity[overpos@from])), main = "ABC Score Distributions")
lines(density(log2(data2$similarity[overpos@from][iscandidate])), col="blue")
lines(density(log2(data2$similarity[overpos@from][!iscandidate])), col="darkgreen")
