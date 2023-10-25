library(tidyverse)
expr <- data.table::fread("TCGA-CESC.htseq_counts.tsv.gz") %>% 
  as.data.frame() %>% 
  remove_rownames() %>% 
  column_to_rownames("Ensembl_ID")
head(expr)
expr <- RNAseq_ids_transform(data = expr, species = "humanGTF", 
                             logcount_to_count = T)
saveRDS(expr, file = "expr.count.rds")

d <- data.table::fread("TCGA-CESC.GDC_phenotype.tsv.gz", header = T)
sur <- data.table::fread("TCGA-CESC.survival.tsv", header = T)

pd <- data.table::fread("pd.sampleMap_CESC_clinicalMatrix", header = T)
sur <- data.table::fread("sur.txt", header = T)


get_exprdat_from_file <- function(input){
  file_lines <- readLines(input)
  data_lines <- file_lines[!grepl("^!", file_lines)] 
  temp_file <- tempfile() 
  writeLines(data_lines, temp_file)
  GPLs <- data.table::fread(temp_file) %>% 
    as.data.frame() %>% 
    remove_rownames() %>% 
    column_to_rownames("ID_REF")
  pd_lines <- file_lines[grepl("(Sample_title)|(Sample_geo_accession)|
                               (Sample_platform_id)|(Sample_source_name_ch1)", 
                               file_lines)]
  temp_file <- tempfile() 
  writeLines(pd_lines, temp_file)
  pd <- read.table(temp_file, sep = "\t", 
                   header = F, 
                   row.names = 1) %>% 
    t() %>% 
    as.data.frame() %>% 
    remove_rownames()
  colnames(pd) <- gsub("^!Sample_", "", colnames(pd))
  dat <- list(data = GPLs, pd = pd)
  dat
}

GSE6791 <- get_exprdat_from_file("GSE6791_series_matrix.txt.gz")
GSE6791_expr <- GSE6791$data
GSE6791_group <- GSE6791$pd
GSE6791_group <- GSE6791_group[str_detect(GSE6791_group$source_name_ch1, "Cervical"),]
GSE6791_expr <- GSE6791_expr[,GSE6791_group$geo_accession]
GSE6791_expr <- chip_ids_transform(exprdat = GSE6791_expr, 
                               input.type = "GPL.Platforms.id", 
                               input = "GPL570") %>% 
  exprdat_Log() %>% 
  exprdat_normalize()
neg <- c("GSM155655", "GSM155653", "GSM155664", "GSM155665", "GSM155666", 
         "GSM155667", "GSM155668", "GSM155669", "GSM155670", "GSM155671")
GSE6791_group$HPV <- ifelse(GSE6791_group$geo_accession %in% neg, "neg", "pos")
save(GSE6791_expr, GSE6791_group, file = "GSE6791.Rdata")



# COX分析 -------------------------------------------------------------------

library(tidyverse)
fpkm <- data.table::fread("TCGA-CESC.htseq_fpkm.tsv.gz", header = T) %>% 
  as.data.frame() %>% 
  remove_rownames() %>% 
  column_to_rownames("Ensembl_ID") %>% 
  RNAseq_ids_transform(species = "humanGTF", logcount_to_count = F)

sur <- read.table("TCGA-CESC.survival.tsv", sep = "\t", header = T, row.names = 1)
int <- intersect(colnames(fpkm), rownames(sur))

sur <- sur[int,]
sur <- dplyr::select(sur, status = OS, time = OS.time)
f <- fpkm[,int]

TCGA <- t(f) %>% as.data.frame() %>% cbind(sur, .)
save(TCGA, file = "TCGA.Rdata")

setwd("./20220321_cox_1.0/")
source("function.R")
modul <- read.csv("../../GSE171894_RAW/top40_hub_genes.csv", header = T)
gene_list <- split(modul$gene_name, modul$module)
Biu_UniCOX(data = "TCGA",
           title = "M1",
           multi.COX = "",
           gene.select = gene_list$`gamma-M1`,
           sur.plot.pval = 0.05, 
           plot.survivol = T
)

Biu_UniCOX(data = "TCGA",
           title = "M2",
           multi.COX = "",
           gene.select = gene_list$`gamma-M2`,
           sur.plot.pval = 0.05, 
           plot.survivol = T
)

Biu_UniCOX(data = "TCGA",
           title = "M4",
           multi.COX = "",
           gene.select = gene_list$`gamma-M4`,
           sur.plot.pval = 0.05, 
           plot.survivol = T
)

Biu_UniCOX(data = "TCGA",
           title = "M6",
           multi.COX = "",
           gene.select = gene_list$`gamma-M6`,
           sur.plot.pval = 0.05, 
           plot.survivol = T
)

library(ggthemes)
unicox <- read.csv("../unicox_M1_M2_M4_M6.csv", header = T)
unicox <- unicox %>% arrange(module, HR)
unicox$Characteristics <- fct_inorder(unicox$Characteristics)
ggplot(unicox, aes(HR, Characteristics, col = module)) + 
  geom_point(size=1.5, aes(col = module)) +
  geom_errorbarh(aes(xmax = Upper, xmin = Lower,col= module),size= 1,height = 0.1) + 
  scale_color_manual(values=c("#54278f", "#b30000", "#238b45", "#fec44f")) +
  ggnewscale::new_scale_color() + 
  geom_point(data = filter(unicox, module == "M1"), size=1.5, 
             mapping = aes(HR, Characteristics, col = Characteristics), 
             inherit.aes = F) +
  geom_errorbarh(data = filter(unicox, module == "M1"),
                 aes(xmax = Upper, xmin = Lower,col= Characteristics), size= 1, height = 0.1) +
  scale_color_manual(values= colorRampPalette(c("#bcbddc", "#54278f"), alpha = TRUE)(14))+
  ggnewscale::new_scale_color() +
  geom_point(data = filter(unicox, module == "M2"), size=1.5, 
             mapping = aes(HR, Characteristics, col = Characteristics), 
             inherit.aes = F) +
  geom_errorbarh(data = filter(unicox, module == "M2"),
                 aes(xmax = Upper, xmin = Lower,col= Characteristics), size= 1, height = 0.1) +
  scale_color_manual(values= colorRampPalette(c("#fdbb84", "#b30000"), alpha = TRUE)(8))+
  
  ggnewscale::new_scale_color() +
  geom_point(data = filter(unicox, module == "M4"), size=1.5, 
             mapping = aes(HR, Characteristics, col = Characteristics), 
             inherit.aes = F) +
  geom_errorbarh(data = filter(unicox, module == "M4"),
                 aes(xmax = Upper, xmin = Lower,col= Characteristics), size= 1, height = 0.1) +
  scale_color_manual(values= colorRampPalette(c("#ffffe5", "#238b45"), alpha = TRUE)(19))+
  
  ggnewscale::new_scale_color() +
  geom_point(data = filter(unicox, module == "M6"), size=1.5, 
             mapping = aes(HR, Characteristics, col = Characteristics), 
             inherit.aes = F) +
  geom_errorbarh(data = filter(unicox, module == "M6"),
                 aes(xmax = Upper, xmin = Lower,col= Characteristics), size= 1, height = 0.1) +
  scale_color_manual(values= colorRampPalette(c("#fff7bc", "#fec44f"), alpha = TRUE)(4))+
  
  
  xlab('HR') + ylab(" ")+
  theme_few()+
  theme(axis.text.x = element_text(size = 7, color = "black"))+
  theme(axis.text.y = element_text(size = 7, color = "black"))+
  geom_hline(aes(yintercept = 14.5), colour = "gray", 
             linetype = "dashed", size = 0.7)+
  theme(legend.background = element_rect(fill="#e4f6ff", size=0.25, linetype="solid"),
        legend.position = "NULL") +
  coord_flip() +
  Seurat::RotatedAxis()



# lasso -------------------------------------------------------------------

library(glmnet)
library(dplyr)
library(ggplot2)
library(ggsci)
text_gene = unicox$gene
load("TCGA.Rdata")
d <- TCGA
x2 = d %>% 
  dplyr::select(text_gene) %>%
  as.matrix()

y <- TCGA[,c("time", "status")] %>% 
  as.matrix()

model_lasso2 <- glmnet(x2, y, nlambda=100, alpha=1, family = "cox")
plot(model_lasso2, xvar = "lambda")
lambdas <- seq(from = 0, to = 0.5, length.out = 200)
cv_fit2 <- cv.glmnet(x=x2, y=y, family = "cox", set.seed(666),lambda=lambdas, alpha=1)
cv_fit2 <- cv.glmnet(x=x2, y=y, family = "cox", set.seed(666), alpha=1)

cvfit = cv_fit2
saveRDS(cvfit, file = "lasso_cvfit.rds")

c(cvfit$lambda.min,cvfit$lambda.1se)
fit1 <- glmnet(x=x2, y=y, nlambda=200, alpha = 1, family = "cox", lambda=cvfit$lambda.min)
fit2 <- glmnet(x=x2, y=y, nlambda=200, alpha = 1, family = "cox", lambda=cvfit$lambda.1se)
rownames(fit1$beta)[as.numeric(fit1$beta)!=0] %>% edit()


cv.Lasso <- cv_fit2
plot(cv.Lasso)

lambda_lse <- cv.Lasso$lambda.1se
lambda_lse
lambda_lse.coef <- coef(cv.Lasso$glmnet.fit, s = lambda_lse)
mat <- as.matrix(lambda_lse.coef)
rownames(mat)[mat != 0] 


xx <- data.frame(lambda=cv_fit2[["lambda"]],cvm=cv_fit2[["cvm"]],cvsd=cv_fit2[["cvsd"]],
                 cvup=cv_fit2[["cvup"]],cvlo=cv_fit2[["cvlo"]],nozezo=cv_fit2[["nzero"]])
xx$ll <- log(xx$lambda)
xx <- xx[xx$ll != -Inf,]
xx <- arrange(xx, ll)
xx$NZERO <- paste0('vars', 1:dim(xx)[1])
xx$NZERO <- factor(xx$NZERO, levels = paste0('vars', 1:dim(xx)[1]))
library(RColorBrewer)
ggplot(xx,aes(ll,cvm, color = NZERO))+
  geom_errorbar(aes(x=ll,ymin=cvlo,ymax=cvup),width=0.05,size=1)+
  geom_vline(xintercept = xx["s15","ll"],size=0.8,color='grey60',alpha=0.8,linetype=2)+
  geom_point(size=2)+
  xlab("Log Lambda")+ylab('Partial Likelihood Deviance')+
  theme_bw(base_rect_size = 1.5)+ 
  # scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(9, "PRGn")))(189)) +
  # scale_color_manual(values = colorRampPalette(rev(brewer.pal(9, "PRGn")))(199))+
  scale_color_manual(values = c(colorRampPalette(brewer.pal(9, "Pastel1"))(80)[10:50],
                                rev(colorRampPalette(brewer.pal(9, "YlOrBr"))(1500)[1:159])
  )
  )+
  scale_x_continuous(expand = c(0.02,0.02))+
  scale_y_continuous(expand = c(0.02,0.02))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15,color='black'),
        axis.text = element_text(size=12,color='black'),
        legend.title = element_blank(),
        #legend.text = element_text(size=12,color='black'),
        legend.position = "none")+
  annotate('text',x = -4.5,y=15.2,label='Optimal Lambda = 0.03',color='black')



# ggirsk ------------------------------------------------------------------

library(rms)
library(ggrisk)
library(tidyverse)
library(survival)
library(survminer)

lis <- c("time", "status", c("ITGAE", "IKZF3", "LSP1", "NEDD9", "CLEC2D", "RBPJ", "TRBC2", 
                             "OXNAD1"))

dat <- TCGA[,lis]

fit <- cph(Surv(time, status)~
             ITGAE+IKZF3+LSP1+NEDD9+CLEC2D+RBPJ+TRBC2+OXNAD1
           ,dat)

coxm_2 <- cph(Surv(time,status)~ITGAE+IKZF3+LSP1+NEDD9+CLEC2D+RBPJ+TRBC2+OXNAD1,
              data=dat,surv=T,x=T,y=T,time.inc = 2*365)
cal_2 <- calibrate(coxm_2, u=2*365, cmethod='KM', m=30, B=1000)
plot(cal_2,lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-year DFS',#便签
     ylab='Actual 2-year DFS(proportion)',#标签
     col=c(rgb(192,98,83,maxColorValue = 255)),#设置一个颜色
     xlim = c(0.5,1),ylim = c(0.2,1)) ##x轴和y轴范围

f <- coxph(Surv(time,status)~ITGAE+IKZF3+LSP1+NEDD9+CLEC2D+RBPJ+TRBC2+OXNAD1,data=dat)
sum.surv<-summary(f)
c_index<-sum.surv$concordance
c_index
# C-index  se(C) 
# 0.714    0.033 



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

library(ComplexHeatmap)
library(circlize)
pheatmap(ggrist_cluster[1:9] %>% t(), show_colnames = F, 
         scale = "row", 
         cluster_rows = F, 
         color = c("#54278f","#fefbf9","#b30000"), 
         treeheight_row = 0, cluster_cols = F)





fit <- survfit(Surv(time, status) ~ riskscore, data = ggrist_cluster)
fit$call$data <- ggrist_cluster ## fu**
kk = ggsurvplot(
  fit,                     
  risk.table = TRUE,      
  pval = TRUE,            
  conf.int = F,
  xlim = c(0,6500),        
  break.time.by = 2000,   
  risk.table.y.text.col = T, 
  risk.table.y.text = FALSE,
  #add.all = TRUE,
  palette = "lancet" # "grey","npg","aaas","lancet","jco", "ucscgb","uchicago","simpsons" and "rickandmorty", "hue"
  # fun = "cumhaz"
); kk


# ROC曲线 -------------------------------------------------------------------

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
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(face = "bold", size = 11, color = "black")
  )




# 决策曲线 --------------------------------------------------------------------

library(rmda)
library(ggDCA)
library(ggplot2)
library(rms)
library(caret)


cph1 <- cph(Surv(time, status)~ITGAE, dat)
cph2 <- cph(Surv(time, status)~ITGAE+IKZF3, dat)
cph3 <- cph(Surv(time, status)~ITGAE+IKZF3+LSP1, dat)
cph4 <- cph(Surv(time, status)~ITGAE+IKZF3+LSP1+NEDD9, dat)
cph5 <- cph(Surv(time, status)~ITGAE+IKZF3+LSP1+NEDD9+CLEC2D, dat)
cph6 <- cph(Surv(time, status)~ITGAE+IKZF3+LSP1+NEDD9+CLEC2D+RBPJ, dat)
cph7 <- cph(Surv(time, status)~ITGAE+IKZF3+LSP1+NEDD9+CLEC2D+RBPJ+TRBC2, dat)
cph8 <- cph(Surv(time, status)~ITGAE+IKZF3+LSP1+NEDD9+CLEC2D+RBPJ+TRBC2+OXNAD1, dat)

options(
  ggplot2.discrete.colour = ggsci::scale_colour_d3,
  ggplot2.discrete.fill = ggsci::scale_fill_d3
)

dca_cph <- dca(cph1, cph2, cph3, cph4, cph5, cph6, cph7, cph8, 
               model.names = c("ITGAE", 
                               "ITGAE+IKZF3", 
                               "ITGAE+IKZF3+LSP1",
                               "ITGAE+IKZF3+LSP1+NEDD9",
                               "ITGAE+IKZF3+LSP1+NEDD9+CLEC2D",
                               "ITGAE+IKZF3+LSP1+NEDD9+CLEC2D+RBPJ",
                               "ITGAE+IKZF3+LSP1+NEDD9+CLEC2D+RBPJ+TRBC2",
                               "ITGAE+IKZF3+LSP1+NEDD9+CLEC2D+RBPJ+TRBC2+OXNAD1"))

ggplot2::ggplot(dca_cph, lwd = 0.5)





# nomogram ----------------------------------------------------------------

library(rms)
example <- ggrist_cluster
ddist <- datadist(example); options(datadist='ddist')
ggrist_cluster %>% head()
ggrist_cluster$risk <- as.numeric(ifelse(ggrist_cluster$riskscore == "high", 1, 0))
gg <- dplyr::select(ggrist_cluster, status, risk)
fit<- lrm(vs~am,data = example)
fit<- lrm(status~ggrist,data = ggrist_cluster)
print(fit,latex = TRUE)
nom <- nomogram(fit,fun=function(x)1/(1+exp(-x)), funlabel = "Probability",fun.at=c(.01,.05,seq(.1,.9,by=.1),.95,.99))
plot(nom)





# 免疫浸润分析 ------------------------------------------------------------------


options(
  ggplot2.discrete.colour = ggsci::scale_colour_d3,
  ggplot2.discrete.fill = ggsci::scale_fill_d3
)


setwd("./免疫浸润分析/")
library(dplyr)
library(tidyverse)
source("CIBERSORT.R")
write.table(GSE6791_expr, file = "exp.txt", sep = "\t")
results = CIBERSORT("LM22.txt", "exp.txt", perm=1000, QN=TRUE)
write.csv(results,"Output.csv")

results <- read.csv("Output.csv",row.names = 1,header = T)
results <- results[GSE6791_group$geo_accession,]
group_list <- GSE6791_group$HPV
group_list$group_list <- group_list$group
fin_N_T <- results %>%
  cbind(Group = group_list) %>%
  as.data.frame() %>%
  mutate(sample = rownames(results)) %>%
  pivot_longer(cols = 1:22,
               names_to = "CellType",
               values_to = "Composition")

save(results, fin_N_T, results_T, file = "immune_1.Rdata")
# boxplot(nomal and cancer)
library(ggpubr)
library(ggthemes)
ggboxplot(
  fin_N_T,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "Group",
  xlab = "",
  ylab = "Cell composition",
  main = ""
) +
  stat_compare_means(
    label = "p.signif",
    method = "t.test",
    ref.group = ".all.",
    hide.ns = T
  ) +
  theme_base() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.3
  ))

library(GSVA)
library(circlize)
library(ComplexHeatmap)
exp <- GSE6791_expr[,GSE6791_group$geo_accession] %>% as.matrix()
# download the immune-cell-type-gene-sets from [doi: 10.1016/j.celrep.2016.12.019]:
geneset <- read.table("immune-cell-type-gene-sets.txt",header = T,sep = "\t")
genes <- unique(geneset$immune.Type) %>%
  lapply(.,function(name){
    as.vector(geneset[geneset$immune.Type == name,"MoleculeName"])
  }); names(genes) <- unique(geneset$immune.Type)
gsva_matrix <- gsva(exp,genes, 
                    method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
exp <- gsva_matrix %>% t() %>% scale() %>% t()
exp[exp > 2] = 2 ; exp[exp < -2] = -2
result <- hclust(dist(t(exp),method = "euclidean"),
                 method = 'complete') %>%
  cutree(k = 2) %>% str_replace_all("1", "high") %>% 
  str_replace_all("2", "low") %>%
  cbind(sub = colnames(exp), immune_cluster = .) %>%
  as.data.frame()
results <- as.factor(result$immune_cluster)
col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
top_annotation = HeatmapAnnotation(
  #results = as.factor(gene_cluster),
  cluster = anno_block(gp = gpar(fill = c("#ffb39c", "#c6c4c3")),
                       labels = c("N","P"),
                       labels_gp = gpar(col = "#443c3a", fontsize = 12)))

Heatmap(exp,name = " ",
        col = col_fun,
        top_annotation = top_annotation,
        column_split = GSE6791_group$HPV,
        show_heatmap_legend = T,
        border = T,
        show_column_names = F,
        show_row_names = T,
        column_title = NULL,show_row_dend = F)

cc <- as.data.frame(net$colors)
write.csv(cc, file = "module_color.csv")



library(IOBR)
estimate<-deconvo_tme(eset = GSE6791_expr, method = "estimate")
estimate<-deconvo_tme(eset = GSE6791_expr, method = "epic", arrays = TRUE)

estimate$group <- GSE6791_group$HPV
et <- estimate
estimate[,2:5] <- scale(as.matrix(estimate[,2:5]))
esti <- pivot_longer(estimate, -c(ID, group), names_to = "gg", values_to = "score")

library(DAAG)  
library(nlme)  
library(ggplot2)  
library(ggpubr)  
library(viridisLite) 
library(ggbeeswarm)
scidat <- esti


scidat$like <- scidat$score
scidat$like <- scidat$like + 
  rnorm(n = length(scidat$like),mean = mean(scidat$like), sd = sd(scidat$like))

head(scidat)
ggplot(data = scidat,aes(x = group, y = like, fill = group))+
  scale_fill_viridis_d(option = "D")+
  geom_violin(alpha=0.5, position = position_dodge(width = .75),size=1,color=NA) +
  geom_point(shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
  geom_boxplot(notch = T, outlier.size = -1, color="black",lwd=1, alpha = 1,show.legend = F, 
               width = 0.3, fill = "white")+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, 
                               dodge.width = .75, 
                               color="black",alpha=.5,
                               show.legend = F) +
  theme_minimal() +
  facet_grid(~gg) + 
  Seurat::RotatedAxis()

a <- esti$group == "pos" & esti$gg == "TumorPurity_estimate"
b <- esti$group == "neg" & esti$gg == "TumorPurity_estimate"
t.test(esti$score[a], esti$score[b])





pacman::p_load(msigdbr, GSVA, tidyverse, clusterProfiler)

geneset <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  subset(select = c("gs_name", "gene_symbol")) %>% 
  as.data.frame(); head(geneset)
geneset <- split(geneset$gene_symbol, geneset$gs_name)
geneset$gammaT <- c("XCL1", "FCER1G", "XCL2", "TYROBP", "KLRC1", "GNLY", "B3GNT7", "TRDC", "NCAM1", "LAT2")


gsva.res <- gsva(GSE6791_expr, geneset, method = "ssgsea") # genesets × annotations


pacman::p_load("ComplexHeatmap", "circlize")

pheatmap(gsva.res, show_colnames = T, scale = "row", angle_col = "45",
         color = c("#443c3a","#fefbf9","#e02e08"), treeheight_row = 0, cluster_cols = F)


library(pheatmap)
m <- gsva.res
m <- t(scale(t(m)))
all(colnames(gsva.res) == GSE6791_group$geo_accession)

GSE6791_group <- arrange(GSE6791_group, HPV)
GSE6791_expr <- GSE6791_expr[]
m <- m[,GSE6791_group$geo_accession]
write.csv(m, file = "gsva.res.csv")
ac = data.frame(g = GSE6791_group$HPV)
rownames(ac) = colnames(m)
rownames(m) <- str_split(rownames(m), "HALLMARK_", simplify = T)[,2]
pheatmap(m,
         show_colnames =F,
         color = colorRampPalette(c("#A555EC", "#EEEEEE", "#FF597B"))(100),
         show_rownames = T,
         cluster_cols = F,
         border_color = NA,
         cluster_rows = T,
         annotation_col = ac)


library(correlation)
library(see)
library(ggplot2)


edit(rownames(gsva.res))

c("HALLMARK_ADIPOGENESIS", "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_ANDROGEN_RESPONSE", 
  "HALLMARK_ANGIOGENESIS", "HALLMARK_APICAL_JUNCTION", "HALLMARK_APICAL_SURFACE", 
  "HALLMARK_APOPTOSIS", "HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_CHOLESTEROL_HOMEOSTASIS", 
  "HALLMARK_COAGULATION", "HALLMARK_COMPLEMENT", "HALLMARK_DNA_REPAIR", 
  "HALLMARK_E2F_TARGETS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", 
  "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE", 
  "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_G2M_CHECKPOINT", 
  "HALLMARK_GLYCOLYSIS", "HALLMARK_HEDGEHOG_SIGNALING", "HALLMARK_HEME_METABOLISM", 
  "HALLMARK_HYPOXIA", "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_IL6_JAK_STAT3_SIGNALING", 
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE", 
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_KRAS_SIGNALING_DN", 
  "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_MTORC1_SIGNALING", 
  "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_MYOGENESIS", 
  "HALLMARK_NOTCH_SIGNALING", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", 
  "HALLMARK_P53_PATHWAY", "HALLMARK_PANCREAS_BETA_CELLS", "HALLMARK_PEROXISOME", 
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_PROTEIN_SECRETION", 
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_SPERMATOGENESIS", 
  "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_UV_RESPONSE_DN", 
  "HALLMARK_UV_RESPONSE_UP", "HALLMARK_WNT_BETA_CATENIN_SIGNALING", 
  "HALLMARK_XENOBIOTIC_METABOLISM", "gammaT")
result <- cor_test(as.data.frame(t(gsva.res)), "gammaT", "HALLMARK_TGF_BETA_SIGNALING"); # plot(result)
plot(result,
     point = list(
       aes = list(color = "gammaT", size = "HALLMARK_TGF_BETA_SIGNALING"),
       alpha = 0.66
     ),
     smooth = list(color = "black", se = F)
) +
  see::theme_modern() +
  see::scale_color_material_c(palette = "rainbow", guide = "none") +
  scale_size_continuous(limits = c(0,15), guide = "none") 
FeaturePlot(NMF, features = "LCN2", slot = "data")




