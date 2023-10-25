library(tidyverse)

dat <- get_exprdat_from_code(GseName = "GSE39001")

dat1 <- dat$`GSE39001-GPL201`
dat2 <- dat$`GSE39001-GPL6244`


expr1 <- chip_ids_transform(dat1$exprdat, input.type = "GPL.Platforms.id", input = "GPL201", I_have_ids = F)
expr2 <- chip_ids_transform(dat2$exprdat, input.type = "GPL.Platforms.id", input = "GPL6244", I_have_ids = F)
it <- intersect(rownames(expr1), rownames(expr2))
expr1 <- expr1[it,]
expr2 <- expr2[it,]

pd1 <- dat1$pd %>% dplyr::select(geo = geo_accession, group = `clinic:ch1`, sample = `sample id:ch1`) %>% mutate(batch = 1)
pd2 <- dat2$pd %>% dplyr::select(geo = geo_accession, group = `clinic:ch1`, sample = `sample id:ch1`) %>% mutate(batch = 2)

library(sva)
expr <- cbind(expr1, expr2)
pd <- rbind(pd1, pd2)
combat_edata = ComBat(dat=expr, batch=pd$batch, mod=NULL, par.prior=TRUE)
save(combat_edata, file = "GSE39001.Rdata")

load("GSE39001.Rdata")

# 主成分分析
if(T){
  library("FactoMineR")
  library("factoextra") 
  dat=t(combat_edata) 
  dat=as.data.frame(dat)
  dat=cbind(dat,group_vs = pd$group) ##给表达矩阵加上分组信息
  dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE) #dat最后一列是group_list。pca分析需要一个纯数值矩阵，所以将dat最后一列去掉以后赋值给dat.pca
  fviz_pca_ind(dat.pca,geom.ind = "point",
               col.ind = dat$group_vs,
               palette = c("#00AFBB", "#E7B800", "#9955FF","#77DDFF"),
               addEllipses = TRUE,
               legend.title = "Groups"#,ellipse.type = "convex" 
  )
}

sur <- read.csv("survivol.csv", row.names = 1)
group <- group[group$sample %in% rownames(sur),]
sur <- sur[group$sample,]
sur$geo <- group$geo
cancer <- combat_edata[,sur$geo]
cancer <- exprdat_Log(cancer) %>% exprdat_normalize()

cancer <- t(cancer)
sur <- rename(sur, status = Status)
cancer <- cbind(sur[,c("time", "status")], cancer)
saveRDS(cancer, file = "clinical.cancer.rds")

cancer <- readRDS("clinical.cancer.rds")

c("ITGAE", "IKZF3", "LSP1", "NEDD9", "CLEC2D", "RBPJ", "TRBC2", 
  "OXNAD1") %in% colnames(cancer)



# ggirsk ------------------------------------------------------------------

library(rms)
library(ggrisk)
library(tidyverse)
library(survival)
library(survminer)

lis <- c("time", "status", c("ITGAE", "IKZF3", "LSP1", "NEDD9", "CLEC2D", "RBPJ", "TRBC2", "OXNAD1"))
lis <- c("time", "status", c("ITGAE", "IKZF3", "LSP1", "NEDD9", "CLEC2D", "RBPJ"))

load("GSE39001.Rdata")
load("HNSCC.Rdata")
dat <- cancer[,lis]
dat <- dat[,lis]

fit <- cph(Surv(time, status)~
             ITGAE+IKZF3+LSP1+NEDD9+CLEC2D+RBPJ, dat)

ggrisk:::model.data(fit)
ggrisk(fit,
       code.highrisk = 'High Risk',
       code.lowrisk = 'Low Risk', 
       title.A.ylab='Risk Score', 
       title.B.ylab='Survival Time(year)', 
       title.A.legend='Risk Group', 
       title.B.legend='Status',     
       title.C.legend='Expression', 
       relative_heights=c(0.1,0.05,0.005,0.15),    
       color.A=c(low='#71a2f3',high='#f56565'),
       color.B=c(code.0='#b2df8a',code.1='#fb9a99'), 
       color.C=c(low='black',median='white',high='red'), 
       vjust.A.ylab=1, 
       vjust.B.ylab=2  
)

ggrist_cluster <- data.frame(ggrist = fit$linear.predictors, sample = rownames(dat)) %>%
  arrange(ggrist) %>%
  remove_rownames() %>%
  column_to_rownames("sample") %>%
  mutate(ggrisk_cluster = case_when(ggrist >= 0 ~ "high",
                                    ggrist < 0 ~ "low"))

ggrist_cluster <- cbind(ggrist_cluster,
                        dat[rownames(ggrist_cluster), c("time", "status")])
colnames(ggrist_cluster)[2] <- "riskscore"

all(rownames(ggrist_cluster) %in% rownames(dat))
ggrist_cluster <- dat[rownames(ggrist_cluster),] %>% cbind(ggrist_cluster)
ggrist_cluster$ss <- ifelse(ggrist_cluster$status == 1, "dead", "live")
ggrist_cluster$rink <- 1:nrow(ggrist_cluster)
ggrist_cluster[c(1,2)] <- NULL
save(ggrist_cluster, file = "ggrisk_cluster.Rdata")

fit <- survfit(Surv(time, status) ~ riskscore, data = ggrist_cluster)
fit$call$data <- ggrist_cluster ## fu**
kk = ggsurvplot(
  fit,                     
  risk.table = TRUE,      
  pval = TRUE,             
  conf.int = F,
  xlim = c(0,2500),        
  break.time.by = 2000,   
  risk.table.y.text.col = T, 
  risk.table.y.text = FALSE,
  #add.all = TRUE,
  palette = "lancet" # "grey","npg","aaas","lancet","jco", "ucscgb","uchicago","simpsons" and "rickandmorty", "hue"
  # fun = "cumhaz"
); kk






library(timeROC)
library(survival)

time_roc_res <- timeROC(
  T = ggrist_cluster$time,
  delta = ggrist_cluster$status,
  marker = ggrist_cluster$ggrist,
  cause = 1,
  weighting="marginal",
  times = c(3 * 365, 5 * 365, 10 * 365),
  ROC = TRUE,
  iid = TRUE
)

time_ROC_df <- data.frame(
  TP_3year = time_roc_res$TP[, 1],
  FP_3year = time_roc_res$FP[, 1],
  TP_5year = time_roc_res$TP[, 2],
  FP_5year = time_roc_res$FP[, 2]
  #TP_10year = time_roc_res$TP[, 3],
  #FP_10year = time_roc_res$FP[, 3]
)
library(ggplot2)
ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#0072B5FF") +
  #geom_line(aes(x = FP_10year, y = TP_10year), size = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(face = "bold", size = 12, color = "black")
  )




