options(scipen = 1000000000)
options(stringsAsFactors = F)
set.seed(1009)
library(readxl)
library(dplyr)
library(corrplot)
library(corrr)
library(RColorBrewer)
library(psych)
library(GGally)

library(CBCgrps)
library(tibble)
library(ggrepel)

library(ggplot2)
library(hrbrthemes)
library(ggpubr)

library(caret)
library(randomForest)
library(pROC)
library(e1071)
library(kknn)
library(adabag)
packageVersion("caret")

rm(list = ls())

setwd("C:/project/pro/onlypro")
#################phenotype preparation##############
#load transcriptomics phenotype 1930
pheno_pro_2132 <- read.csv("pheno_pro_2132.csv",check.names = FALSE)

#subjects filter and preparation for analysis

pheno_all <- pheno_pro_2132[, c("zz_nr",
                                # "lg_dnaAxiom_s4f4","u3n_meth850K_ff4",
                                #all IDs
                                "u3talteru", "u3csex", "u3tgewi", "u3tbmi",
                                "u3tgroe","u3tbmiwho","u3talkkon","u3ttumf","u3twhrat",
                                
                                "u3tdm110y15","u3tglukfast_n","u3tsysmm","u3tnuecht",
                                "u3tdiamm", "u3tcigreg_sf","u3tphys","u3l_hdln",
                                "u3l_ldln","u3l_trin","u3l_choln",
                                "u3h_crp", "u3tbfp_k")]
sum(is.na(pheno_all))
clinic <- na.omit(pheno_all)

clinic$u3csex <- factor(clinic$u3csex,
                        levels = c(1,2),
                        labels = c("male", "female"))
clinic$u3tbmiwho <- factor(clinic$u3tbmiwho,
                           levels = c(1,2,3,4,5,6),
                           labels = c("underweight", "normal weight", "pre-obese","obese grade 1", "obese grade 2","obese grade 3"))
clinic$u3tdm110y15 <- factor(clinic$u3tdm110y15,
                             levels = c(0,1,2,3,4,5,6,7,8,9),
                             labels = c("normal","IFG","IGT","IFG/IGT","new diagnosed", "known Typ 2 Diabetes", "Typ 1 Diabetes",  "known drug",  "unknown", "no validation"))
table(clinic$u3tdm110y15)
#"1 = Regular smokers , 2 = Occasional smokers , 3 = Ex-smokers , 4 = Never smokers"
clinic$u3tcigreg_sf <- factor(clinic$u3tcigreg_sf,
                              levels = c(1,2,3,4),
                              labels = c("smoker","smoker", "ex-smoker", "never-smoker"))
clinic$u3tphys <- factor(clinic$u3tphys, levels = c(1,2),
                         labels = c("active", "inactive"))

#firstly filter underweight (BMI < 15 kg/m2) or all missing covariate values
summary(clinic$u3tbmiwho)

clinic <- filter(clinic, clinic$u3tbmiwho != "underweight")
summary(clinic)

#filer didn’t have at least 8 h of fasting 
table(clinic$u3tnuecht)
clinic <- filter(clinic, clinic$u3tnuecht != 2)

#filer out diabetes
summary(clinic$u3tdm110y15)
clinic$u3tdm110y15 <- factor(clinic$u3tdm110y15,
                             levels = c("normal","IFG","IGT","IFG/IGT","new diagnosed", "known Typ 2 Diabetes", "Typ 1 Diabetes",  "known drug",  "unknown", "no validation"),
                             labels = c("Non_T2D", "Non_T2D", "Non_T2D", "Non_T2D", "T2D", "T2D", "Non_T2D", "Non_T2D", "Non_T2D", "Non_T2D"))
clinic$u3tbmiwho <- factor(clinic$u3tbmiwho,
                           levels = c("normal weight", "pre-obese", "obese grade 1", "obese grade 2", "obese grade 3"),
                           labels = c("Non_Obese", "Non_Obese", "Obese", "Obese", "Obese"))
summary(clinic[ , 2:23])
write.csv(clinic, file = "phenotype_final.csv", row.names = FALSE)

###########################stastical analysis####################
clinic <- read.csv("phenotype_final.csv", header = 1, row.names = 1)
continous <- clinic[,c("u3tbmiwho","u3talteru", "u3tgewi", "u3tgroe","u3talkkon",
                       "u3ttumf","u3twhrat", "u3tglukfast_n","u3tsysmm",
                       "u3tdiamm", "u3l_hdln", "u3l_ldln","u3l_trin",
                       "u3l_choln", "u3h_crp", "u3tbfp_k")]
category <- clinic[,c("u3tbmiwho", "u3csex", "u3tcigreg_sf", "u3tphys",
                      "u3tdm110y15")]
table1 <- twogrps(continous, gvar = "u3tbmiwho", norm.rd = 1, sk.rd =1 , cat.rd =1 ,
                  p.rd =3,ShowStatistic = T,  pnormtest = 0)
print(table1$Table, quote = T)
table2 <- twogrps(category, gvar = "u3tbmiwho", norm.rd = 1, sk.rd =1 , cat.rd =1 ,
                  p.rd =3,ShowStatistic = F,  pnormtest = 0)
print(table2$Table, quote = T)

write.csv(table1$Table, file = "continous_var_stastic.csv")
write.csv(table2$Table, file = "category_var_stastic.csv")

########################prepare data############################
rm(list = ls())
clinic <- read.csv("phenotype_final.csv", header = 1)
proteomic_2132 <- read.csv("proteomic_2132_uniprot.csv", header = 1)
#adjust/match protein names and numbers
offical_name <- read_excel("offical_name.xlsx")
pro_name <- offical_name[, c("koraID", "Genenames", "Entry", "Protein names")]
pro_name <- pro_name[complete.cases(pro_name), ]
#write.csv(pro_name, file = "pro_name_809.csv")
proteomic_3 <- proteomic_2132[, 1:3]
selected_columns <- names(proteomic_2132[,4:845]) %in% pro_name$koraID
proteomic_809 <- proteomic_2132[, selected_columns]
colnames(proteomic_809)[4:812] <- pro_name$Genenames

com_data <- merge(clinic, proteomic_809, by = "zz_nr")
write.csv(com_data, file = "pheno_pro_2045.csv", row.names = FALSE)
write.csv(proteomic_809, file = "proteomic_809.csv", row.names = FALSE)


#############"SAA2-SAA4"####################
##############PCA##############
rm(list = ls())
data <- read.csv(file = "pheno_11_pro_809.csv", row.names = 1)
data$BMI_grade <- factor(data$BMI_grade,
                         levels = c(0,1),
                         labels = c("Non_Obese", "Obese"))
poulet <- data[, -c(1, 3:11)]
colnames(poulet)[1] <- "BMI"
poulet$BMI <- as.factor(poulet$BMI)

pro_matrix <- as.matrix(poulet[,-1])
pca2 = pca(pro_matrix, method = "svd", scale = "uv", nPcs = 3)
summary(pca2)
Groups <- poulet$BMI
# Create the plot
png("plot/pca.png", width=20, height=20, units = "cm", res=600) 
ggplot(as.data.frame(pca2@scores), aes(x = PC1, y = PC2, color = Groups)) +
  geom_point(size = 0.5)+
  scale_color_manual(values = c("#1597A5","#FFC24B","#FEB3AE") ) +
  scale_fill_manual(values = c("#1597A5","#FFC24B"))+
  theme_bw() +
  labs(color = "Groups", fill = NULL) +
  ggtitle("PCA plot,  non−obesity vs. obesity") +
  guides(color=guide_legend(title= NULL))+
  xlab(paste("PC1",round( pca2@R2[1], digits = 3) * 100, "% variance explained")) +
  ylab(paste("PC2",round( pca2@R2[2], digits = 3) * 100, "% variance explained")) +
  stat_ellipse(data=as.data.frame(pca2@scores),
               geom = "polygon",level = 0.95,type = "norm",
               linetype = 2,size=0.5,
               aes(fill=Groups),
               alpha=0.2,
               show.legend = T)
dev.off()
#PLS-DA
otu_raw <- read.csv(file = "pheno_11_pro_809.csv", row.names = 1)
otu_raw$BMI_grade <- factor(otu_raw$BMI_grade,
                            levels = c(0,1),
                            labels = c("Non_Obese", "Obese"))
group <- data.frame(otu_raw[, "BMI_grade" ])
colnames(group) <- 'bmi'
otu <- otu_raw[, 12:820]

df1_plsda <- opls(otu, group$bmi, orthoI = 0)

data <- as.data.frame(df1_plsda@scoreMN)
data$group = group$bmi
data$samples = rownames(data)
#提取解释度
x_lab <- df1_plsda@modelDF[1, "R2X"] * 100
y_lab <- df1_plsda@modelDF[2, "R2X"] * 100
png("plot/plsda1.png", width=20, height=20, units = "cm", res=600) 
ggplot(data,aes(x=p1,y=p2,color=group))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(size= 0.5)+#绘制点图并设定大小
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed",color="red")+
  geom_hline(yintercept = 0,lty="dashed",color="red")+#图中虚线
  guides(color=guide_legend(title= NULL))+#去除图例标题
  labs(x=paste0("P1 (",x_lab,"%)"),
       y=paste0("to1"), fill = NULL)+ #将x、y轴标题改为贡献度
  stat_ellipse(data=data,type = "norm",
               geom = "polygon",level = 0.95,
               linetype = 2,size=0.5,
               aes(fill= group),
               alpha=0.2,
               show.legend = T)+
  scale_color_manual(values = col) +#点的颜色设置
  scale_fill_manual(values = c("#1597A5","#FFC24B"))+
  theme(axis.title.x=element_text(size=12),#修改X轴标题文本
        axis.title.y=element_text(size=12,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=10),#修改x轴刻度标签文本
        axis.text.x=element_text(size=10),#修改y轴刻度标签文本
        panel.grid=element_blank())#隐藏网格线
dev.off()
data_VIP <- df1_plsda@vipVn
data_VIP_select <- data_VIP[data_VIP > 1] #阈值通常设为1
as.data.frame(data_VIP_select)
#将VIP值与原始数据合并
otu1 <- as.data.frame(t(otu_raw))
data_VIP_select <- cbind(otu1[names(data_VIP_select), ], data_VIP_select)
data_VIP_select$data_VIP_select <- "VIP"
#排序
data_VIP_select <- data_VIP_select[order(data_VIP_select$data_VIP_select, decreasing = TRUE), ]
# plot(df1_plsda, typeVc = "x-loading") #展示前10个
head(data_VIP_select)

################logistical regression model######################
rm(list = ls())
com_data <- read.csv("pheno_pro_2045.csv",check.names = FALSE)
com_data[,c("SAA2-SAA4")]

com_data$u3tbmiwho <- ifelse(com_data$u3tbmiwho == "Obese", 1, 0)
clinic_logistic <- data.frame(row.names = com_data$zz_nr,
                              BMI = com_data$u3tbmi,
                              BMI_grade= com_data$u3tbmiwho,
                              age = com_data$u3talteru,
                              sex = com_data$u3csex,
                              activities = com_data$u3tphys,
                              SBP = com_data$u3tsysmm,
                              TG = log(com_data$u3l_trin),
                              HDL = com_data$u3l_hdln,
                              smoking = com_data$u3tcigreg_sf,
                              FPG = com_data$u3tglukfast_n,
                              T2D = com_data$u3tdm110y15)
sum(is.na(clinic_logistic))
co_data <- cbind(clinic_logistic, com_data[, 26:834])
#write.csv(co_data, file = "pheno_11_pro_809.csv")
glm(BMI ~ + age + sex + activities +
      SBP + TG + smoking + HDL + FPG + T2D, data = clinic_logistic) 
#basic analysis
set.seed(123)
outcome1 <- c()
for (i in 12:820){
  model_logistic <- glm(BMI_grade ~ co_data[,i] + age + sex, data = co_data,
                        family = "binomial") 
  outcome1 <- rbind(outcome1, c(colnames(co_data)[i], coef(summary(model_logistic))[2, c(1,2,4)]))
}

outcome1 <- data.frame(outcome1)
names(outcome1)[4]  <- "pvalue"
names(outcome1)[1] <- "Proteins"
for (i in 2:4) {outcome1[,i] <- as.numeric(outcome1[,i])}
outcome1 <- outcome1[order(outcome1$pvalue),]
outcome1$P.adjust <- p.adjust(outcome1$pvalue,method="bonferroni",n=length(outcome1$pvalue))

outcome1[,2:3] <- round(outcome1[,2:3], digits = 3)
outcome1$pvalue <- signif(outcome1$pvalue, digits = 3)
outcome1$P.adjust <- signif(outcome1$P.adjust, digits = 3)
outcome1$log10_adj_pval = -log10(outcome1$P.adjust)
outcome1$log10_pval = -log10(outcome1$pvalue)

#OR 95%CI
outcome1$OR <- round(exp(outcome1$Estimate), digits = 3)
for (i in 1:nrow(outcome1)) {
  model_logistic <- glm(BMI_grade ~ co_data[, i + 11] + age + sex, data = co_data, family = "binomial")
  conf_int <- confint(model_logistic)
  lower_bound <- exp(conf_int[2, 1])
  upper_bound <- exp(conf_int[2, 2])
  
  #calculate 95%CI
  CI <- paste("(", round(lower_bound, 3), "; ", round(upper_bound, 3), ")", sep = "")
  outcome1[i, "CI"] <- CI
}
outcome1$OR_CI <- paste(outcome1$OR, paste(outcome1$CI), sep = " ")

outcome1 <- outcome1[, -c(8:9)]
#beta 95%CI
#outcome1$conf_interval <- NA
#for (i in 1:nrow(outcome1)) {
#  model <- glm(BMI_grade ~ co_data[, i + 11] + age + sex, data = co_data, family = "binomial")
#  conf_int <- confint(model)
#  outcome1$conf_interval[i] <- paste0("(", round(conf_int[2,1], 3), "; ", round(conf_int[2,2], 3), ")")}
#outcome1$Estimate_conf <- paste(outcome1$Estimate, outcome1$conf_interval, sep = " ")

#table(outcome1[,5] < 0.05)
outcome1_sig <- outcome1[outcome1$P.adjust < 0.05, ]

#save it
write.csv(outcome1, file = "logistic_basic.csv", quote = FALSE, row.names = FALSE)
write.csv(outcome1_sig, file = "logistic_basic_sig.csv", quote = FALSE, row.names = FALSE)

#sensitive analysis
#+age,sex(BMI)", "+BP", "+FPG,T2D", "+lipids", "+lifestyle"
set.seed(123)
outcome <- c()
for (i in 12:820){
  model_logistic <- glm(BMI_grade ~ co_data[,i] + age + sex + SBP + FPG + T2D
                        +TG +  HDL+ activities +smoking, data = co_data,
                        family = "binomial") 
  outcome <- rbind(outcome, c(colnames(co_data)[i], coef(summary(model_logistic))[2, c(1,2,4)]))
}

outcome
outcome <- data.frame(outcome)
names(outcome)[4]  <- "pvalue"
names(outcome)[1] <- "Proteins"
for (i in 2:4) {outcome[,i] <- as.numeric(outcome[,i])}
outcome <- outcome[order(outcome$pvalue),]
outcome$P.adjust <- p.adjust(outcome$pvalue,method="bonferroni",n=length(outcome$pvalue))

# export the outcome
outcome[,2:3] <- round(outcome[,2:3], digits = 3)
outcome$pvalue <- signif(outcome$pvalue, digits = 3)
outcome$P.adjust <- signif(outcome$P.adjust, digits = 3)
outcome$log10_adj_pval = -log10(outcome$P.adjust)
outcome$log10_pval = -log10(outcome$pvalue)

#OR 95%CI
outcome$OR <- round(exp(outcome$Estimate), digits = 3)
for (i in 1:nrow(outcome)) {
  model_logistic <- glm(BMI_grade ~ co_data[, i + 11]+ age + sex + SBP + FPG + T2D
                        +TG +  HDL+ activities +smoking, data = co_data, family = "binomial")
  conf_int <- confint(model_logistic)
  lower_bound <- exp(conf_int[2, 1])
  upper_bound <- exp(conf_int[2, 2])
  
  #calculate 95%CI
  CI <- paste("(", round(lower_bound, 3), "; ", round(upper_bound, 3), ")", sep = "")
  outcome[i, "CI"] <- CI
}
outcome$OR_CI <- paste(outcome$OR, paste(outcome$CI), sep = " ")

outcome <- outcome[, -c(8:9)]
#table(outcome[,5] < 0.05)
outcome_sig <- outcome[outcome$P.adjust < 0.05, ]

#save it
write.csv(outcome, file = "logistic_full.csv", quote = FALSE, row.names = FALSE)
write.csv(outcome_sig, file = "logistic_full_sig.csv", quote = FALSE, row.names = FALSE)

indices1 <- match(outcome1_sig$Proteins, pro_name$Genenames)
name1 <- pro_name[indices1, ]

indices <- match(outcome_sig$Proteins, pro_name$Genenames)
name <- pro_name[indices, ]

write.csv(name1, file = "sig_name1.csv", row.names = F)
write.csv(name, file = "sig_name.csv", row.names = F)
#plot
outcome <- read.csv("logistic_basic.csv", header=TRUE, row.names=1,check.names = FALSE)
outcome$OR <- exp(outcome$Estimate)
table(outcome[,4] < 0.05)
#filtered_outcome <- outcome1[outcome1$P.adjust < 0.05, ]
plot_log = data.frame(row.names = rownames(outcome),
                      Estimate = outcome$OR, 
                      log10_pval = outcome$log10_pval)
#filter upregulated and down regulated proteins and use ggrepel function to put it in volcano plot
sig_logistic <- outcome[outcome[,4] < 0.05, ]
sig_logistic$labels <- rownames(sig_logistic)
sig_logistic$Estimate <- sig_logistic$OR
head(sig_logistic,20)
relative_size <- (outcome$log10_pval)/3

#x:c(-1.1, 1.1),y: c(-1,19)
logistic_basic <- ggplot(plot_log, aes(x = Estimate, y=log10_pval,
                                       color = ifelse(outcome$P.adjust < 0.05, "Association",
                                                      "Unassociation")))+
  geom_point(shape = 16, alpha = 0.7, size = 3.5) + 
  scale_x_continuous(limits=c(0, 2.5), expand=c(0,0)) +
  scale_y_continuous(limits=c(-1,29), expand=c(0,0)) + 
  labs(x= "Odds ratio", y=bquote('P ('*-Log[10]*')'), #Beta coefficient
       color = "Association Status") +
  ggtitle("Obesity association: Basic Model") +
  geom_hline(yintercept = -log10(0.05/809), linetype="dashed")+
  theme_classic() +
  scale_color_manual(values=c("Association"="red", "Unassociation" = "grey")) + 
  geom_text_repel(data = sig_logistic,
                  aes(label = labels), color = "black",
                  size = 3, family = "Arial", force = 1,
                  segment.color = "black", segment.size = 0.8, segment.alpha = 1,
                  hjust = 0.5, show.legend=F)+ 
  theme(legend.position = "none")#delete legend
logistic_basic
warnings()

outcome1 <- read.csv("logistic_full.csv", header=TRUE, row.names=1,check.names = FALSE)
outcome1$OR <- exp(outcome1$Estimate)
table(outcome1[,4] < 0.05)
#filtered_outcome <- outcome1[outcome1$P.adjust < 0.05, ]
plot_log1 = data.frame(row.names = rownames(outcome1),
                       Estimate = outcome1$OR, 
                       log10_pval = outcome1$log10_pval)
#filter upregulated and down regulated proteins and use ggrepel function to put it in volcano plot
sig_logistic1 <- outcome1[outcome1[,4] < 0.05, ]
sig_logistic1$labels <- rownames(sig_logistic1)
sig_logistic1$Estimate <- sig_logistic1$OR
head(sig_logistic1,20)
relative_size1 <- (outcome1$log10_pval)/3

#x:c(-1.1, 1.1),y: c(-1,19)
logistic_full <- ggplot(plot_log1, aes(x = Estimate, y=log10_pval,
                                       color = ifelse(outcome1$P.adjust < 0.05, "Association",
                                                      "Unassociation")))+
  geom_point(shape = 16, alpha = 0.7, size = 3.5) + 
  scale_x_continuous(limits=c(0.4, 1.8), expand=c(0,0)) +
  scale_y_continuous(limits=c(-1,15), expand=c(0,0)) + 
  labs(x= "Odds ratio", y=bquote('P ('*-Log[10]*')'), #Beta coefficient
       color = "Association Status") +
  ggtitle("Obesity association: Full Model") +
  geom_hline(yintercept = -log10(0.05/809), linetype="dashed")+
  theme_classic() +
  scale_color_manual(values=c("Association"="red", "Unassociation" = "grey")) + 
  geom_text_repel(data = sig_logistic1,
                  aes(label = labels), color = "black",
                  size = 3, family = "Arial", force = 1,
                  segment.color = "black", segment.size = 0.8, segment.alpha = 1,
                  hjust = 0.5, show.legend=F)+ 
  theme(legend.position = "none")#delete legend
logistic_full
warnings()
logis <- ggarrange(logistic_basic, logistic_full,
                   labels = c("A", "B"), font.label = list(size = 20),label.x = -0.01,
                   label.y = 1,
                   ncol = 2, nrow = 1, common.legend = FALSE,widths = c(4, 4, 4),
                   heights = c(4, 4, 4), align = "hv") 

logis

ggsave("plot/logistic_model_OR.png", logis, dpi=600, width= 30, height= 15, unit="cm")
#ggsave("plot/logistic_model.png", logis, dpi=600, width= 30, height= 15, unit="cm")

######################t test plot parts########################################
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
library(ggsignif)
library(ggpubr)
data <- read.csv(file = "pheno_11_pro_809.csv", row.names = 1,check.names = FALSE)
data$BMI_grade <- ifelse(data$BMI_grade == 1, "Obese", "Non-Obese")
spro <- data[,c("BMI_grade", "CRP", "APOD", "LGALS3BP", "SAA2-SAA4")]
data1 <- spro
colnames(data1)[1] <- "BMI"
colnames(data1)[5] <- "SAA2"
data1$BMI <- factor(data1$BMI)
summary(data1)
#load function
theme_niwot <- function(){
  theme_bw() +
    theme(text = element_text(family = "Helvetica Light"),
          axis.text = element_text(size = 16), 
          axis.title = element_text(size = 18),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 18, vjust = 1, hjust = 0),
          legend.text = element_text(size = 12),          
          legend.title = element_blank(),                              
          legend.position = c(0.95, 0.15), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))
}
#plot
colnames(data1)
#"CRP"      "APOD"     "LGALS3BP"  SAA2-SAA4
p1 <- ggplot(data1, aes(x = BMI, y = CRP, fill = BMI)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.5) +
  geom_point(aes(y = SAA2, color = BMI), 
             position = position_jitter(width = 0.15), size = 1, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5) +
  labs(y = "CRP", x = NULL) +
  guides(fill = FALSE, color = FALSE) +
  scale_y_continuous(limits = c(-7, 5)) +
  scale_fill_manual(values = c("#FA7070", "#C4E4FF")) +
  scale_colour_manual(values = c("#FA7070", "#C4E4FF")) +
  theme_niwot()+
  ggtitle("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        axis.title.x = element_text(size = 18),  # 调整 X 轴标签字体大小
        axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, vjust = 0),
        plot.caption = element_text(hjust = 0.5, vjust = 0)) +
  stat_compare_means(method = "t.test") #, aes(label = ..p.signif..)

p2 <- ggplot(data1, aes(x = BMI, y = APOD, fill = BMI)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.5) +
  geom_point(aes(y = SAA2, color = BMI), 
             position = position_jitter(width = 0.15), size = 1, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5) +
  labs(y = "APOD", x = NULL) +
  guides(fill = FALSE, color = FALSE) +
  scale_y_continuous(limits = c(-7, 5)) +
  scale_fill_manual(values = c("#FA7070", "#C4E4FF")) +
  scale_colour_manual(values = c("#FA7070", "#C4E4FF")) +
  theme_niwot()+
  ggtitle("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        axis.title.x = element_text(size = 18),  # 调整 X 轴标签字体大小
        axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, vjust = 0),
        plot.caption = element_text(hjust = 0.5, vjust = 0)) +
  stat_compare_means(method = "t.test") #, aes(label = ..p.signif..)
p3 <- ggplot(data1, aes(x = BMI, y = LGALS3BP, fill = BMI)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.5) +
  geom_point(aes(y = SAA2, color = BMI), 
             position = position_jitter(width = 0.15), size = 1, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5) +
  labs(y = "LGALS3BP", x = NULL) +
  guides(fill = FALSE, color = FALSE) +
  scale_y_continuous(limits = c(-7, 5)) +
  scale_fill_manual(values = c("#FA7070", "#C4E4FF")) +
  scale_colour_manual(values = c("#FA7070", "#C4E4FF")) +
  theme_niwot()+
  ggtitle("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        axis.title.x = element_text(size = 18),  # 调整 X 轴标签字体大小
        axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, vjust = 0),
        plot.caption = element_text(hjust = 0.5, vjust = 0)) +
  stat_compare_means(method = "t.test") #, aes(label = ..p.signif..)

p4 <- ggplot(data1, aes(x = BMI, y = SAA2, fill = BMI)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.5) +
  geom_point(aes(y = SAA2, color = BMI), 
             position = position_jitter(width = 0.15), size = 1, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5) +
  labs(y = "SAA2-SAA4", x = NULL) +
  guides(fill = FALSE, color = FALSE) +
  scale_y_continuous(limits = c(-7, 5)) +
  scale_fill_manual(values = c("#FA7070", "#C4E4FF")) +
  scale_colour_manual(values = c("#FA7070", "#C4E4FF")) +
  theme_niwot()+
  ggtitle("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        axis.title.x = element_text(size = 18),  # 调整 X 轴标签字体大小
        axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, vjust = 0),
        plot.caption = element_text(hjust = 0.5, vjust = 0)) +
  stat_compare_means(method = "t.test") #, aes(label = ..p.signif..)



s <- ggarrange(p1, p2, p3, p4,labels = c( "C"),
               font.label = list(size = 20),
               ncol = 2, nrow = 2, common.legend = FALSE,widths = c(4, 4, 4),
               heights = c(4, 4, 4), align = "hv")
s
final <- ggarrange(logis, s,#labels = c("A", "B", "C", "D"),
                   font.label = list(size = 20),
                   ncol = 1, nrow = 2, common.legend = FALSE,widths = c(4, 4, 4),
                   heights = c(4, 8, 4), align = "hv")
ggsave("plot/logis_final.png", final, dpi=600, width= 30, height= 45, units="cm", bg = "white")

ggsave("plot/proteins.png", s, dpi=600, width= 30, height= 30, units="cm", bg = "white")
###################model type ################################
a <- data.frame(
  x = factor(c("+age,sex(BMI)", "+BP", "+FPG,T2D", "+lipids", "+lifestyle"), 
             levels = c("+age,sex(BMI)", "+BP", "+FPG,T2D", "+lipids", "+lifestyle")),
  y = as.numeric(c(25, 24, 21, 12, 11))
)

model_type <- ggplot(a, aes(x = x, y = y, group = 1)) +
  geom_line( color="#FF8F8F",size= 2, alpha = 1) +
  geom_point(shape=21, color="black", fill="#FF8F8F", size=5) +
  theme_ipsum()+
  labs(x = "Model Type", y = "Number of Significant Proteins")+
  theme_bw() +
  geom_text(aes(label = y), vjust = -0.8, hjust = -0, size = 4)+ 
  theme(axis.title.x = element_text(size = 14),  # 调整 X 轴标签字体大小
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),  # 调整 X 轴标签字体大小
        axis.text.y = element_text(size = 12)) 
model_type
ggsave("plot/model_type.png", model_type, dpi=600, width= 15, height= 15, unit="cm")

#####################priority lasso#####################
library(prioritylasso)
library(pROC)
library(dplyr)

rm(list = ls())
logistic_basic_sig <- read.csv("logistic_basic_sig.csv", header=TRUE, row.names=1)
pheno_pro_1557 <- read.csv("pheno_pro_2045.csv")
pro <- pheno_pro_1557[,26:834]
clinic<- data.frame(row.names = pheno_pro_1557$zz_nr,
                    BMI = pheno_pro_1557$u3tbmiwho,
                    age = pheno_pro_1557$u3talteru,
                    sex = pheno_pro_1557$u3csex,
                    activities = pheno_pro_1557$u3tphys,
                    SBP = pheno_pro_1557$u3tsysmm,
                    TG = log(pheno_pro_1557$u3l_trin),
                    HDL = pheno_pro_1557$u3l_hdln,
                    smoking = pheno_pro_1557$u3tcigreg_sf,
                    FPG = pheno_pro_1557$u3tglukfast_n,
                    T2D = pheno_pro_1557$u3tdm110y15)

clinic<- data.frame(row.names = pheno_pro_1557$zz_nr,
                    BMI = pheno_pro_1557$u3tbmiwho,
                    age = pheno_pro_1557$u3talteru,
                    sex = pheno_pro_1557$u3csex)
data<- cbind(clinic, pro[,colnames(pro) %in% rownames(logistic_basic_sig)])
# change "male"and "female"to numeric
data$sex <- ifelse(data$sex == "male", 1, 0)
data$BMI <- ifelse(data$BMI == "Obese", 1, 0)

data <- as.matrix(data)
pl_out <- as.numeric(data[,1])  # outcome
pl_pred <- data[,2:27]  # age+sex+ proteins
blocks <- list(bp1=1:2, bp2=3:26) #in 2:25, 1:2 is age+sex, 3:24 is proteins

#left:lambda.min, right:lambda.1se
pl1<- prioritylasso(X = pl_pred, Y = pl_out, family = "gaussian",
                    type.measure = "mse",  block1.penalization = FALSE, 
                    blocks = blocks, max.coef = c(Inf,Inf),standardize = FALSE,
                    nfolds = 10,lambda.type = "lambda.1se",
                    cvoffset = FALSE)
pl1$block1unpen

pl1$nzero
#the proteins results
pl1$coefficients
coeff1 <- pl1$coefficients
coeff1 <- coeff1[coeff1 != 0]
print(round(coeff1, 4))
names(coeff1)

selected_proteins_list <- list()
# repeat 1000 times
set.seed(1009)
for (i in 1:1000) {
  pl <- prioritylasso(X = pl_pred, Y = pl_out, family = "gaussian",
                      type.measure = "mse",  block1.penalization = FALSE, 
                      blocks = blocks, max.coef = c(Inf,Inf),standardize = FALSE,
                      nfolds = 10, lambda.type = "lambda.1se")
  
  #extract proteins in each model with non-zero cofficident
  selected_proteins <- names(pl$coefficients[pl$coefficients != 0])
  
  # add the extracted proteins name in the list
  selected_proteins_list[[i]] <- selected_proteins
}

# extracted the proteins that frequency >20%
frequent_pro <- table(unlist(selected_proteins_list))
a <- data.frame(frequent_pro)
write.csv(a, file = "lasso-frequency.csv")
names(frequent_pro[frequent_pro >= 200])

data<- cbind(clinic, pro[,colnames(pro) %in% rownames(logistic_basic_sig)])
select_pro_agesex<- data[ ,c("BMI", "age", "sex", "A2M", "ADIPOQ", "AFM","APCS",
                             "APOD", "APOF", "APOM","AZGP1","C4B", "CFH",
                             "CRP", "LGALS3BP","GPX3", "SERPINF1", "SERPINA6", "S100A9")]
write.csv(select_pro_agesex, file = "select_pro_agesex.csv")

####################MC#####################
######load pro data#############
rm(list = ls())

#prepare and devide data (7:3)
setwd("C:/project/pro/onlypro")
data <- read.csv(file = "select_pro_agesex.csv", row.names = 1) #don't , check.names = FALSE here
data <- data[,-c(2,3)] #delete sex and age colums
# change "male"and "female"to numeric
#data$sex <- ifelse(data$sex == "male", 1, 0)
data$BMI <- as.factor(data$BMI)
set.seed(1009)
ind<-sample(2,nrow(data),replace = TRUE,prob = c(0.7,0.3))
traindata<-data[ind == 1,]
testdata<-data[ind == 2,]

########################load phenotype data#################
rm(list = ls())
data <- read.csv(file = "pheno_11_pro_809.csv", row.names = 1)
data <- data[,-c(1,12:820)] #delete sex and age colums
# change "male"and "female"to numeric
#data$sex <- ifelse(data$sex == "male", 1, 0)
data$BMI_grade <- as.factor(data$BMI)
colnames(data)[1] <- "BMI"
set.seed(1009)
ind<-sample(2,nrow(data),replace = TRUE,prob = c(0.7,0.3))
traindata<-data[ind == 1,]
testdata<-data[ind == 2,]

################10 fold validation####################
# include ntree, mtry(waiting to be choosen), run two times
library(caret)
grid <- expand.grid(mtry = c(2, 3, 4,5,6,7,8,9, 10))
# select by train()
names(getModelInfo()) #https://topepo.github.io/caret/train-models-by-tag.html#Neural_Network
set.seed(1009)
model <- train(BMI ~ . , #model design
               data = traindata,method = "rf",# method RF
               trControl = trainControl(method = "repeatedcv",  # method
                                        number = 10, repeats = 10) , 
               tuneGrid = grid) 
model[["results"]][["Accuracy"]]
print(model$bestTune)

set.seed(1009)
param_grid <- expand.grid(sigma = c(1/15, 1, 10), C = c(0.1, 1, 10))#15 = protein numbers
model <- train(BMI ~ ., #model design
               data = traindata, method = "svmRadial",# method svm
               tuneGrid = param_grid,
               trControl = trainControl(method = "repeatedcv",  # method
                                        number = 10)) 
model[["results"]][["Accuracy"]]
print(model$bestTune)

set.seed(1009)
model <- train(BMI ~ ., #model design
               data = traindata, method = "kknn",# method knn
               trControl = trainControl(method = "repeatedcv",  # method
                                        number = 10)) 
model[["results"]][["Accuracy"]]
print(model$bestTune)


set.seed(1009)
#traindata$smoking <- ifelse(traindata$smoking == "never-smoker", 0, 1)
model <- train(BMI ~ ., #model design
               data = traindata[,1:10],method = "AdaBag",# method RF
               trControl = trainControl(method = "repeatedcv",  # method
                                        number = 10)) 
model[["results"]][["Accuracy"]]
print(model$bestTune)

#####################proteins MODEL################
#random forest
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9789589/
#use 10 fold cross validation to choose the best parameters
library(randomForest)
library(pROC)

set.seed(1009)
#default is 500 trees
randomForest_rf<-randomForest(BMI ~.,data = traindata, #BMI ~ . - age - sex
                              ntree = 500,
                              mtry = 5,
                              proximity = TRUE,
                              importance = TRUE)
table(predict(randomForest_rf),traindata$BMI)
randomForest_rf
#Plot error rates based on different numbers in the random forest
#plot(randomForest_rf)
#importance score
importance(randomForest_rf,type = 1)
#importance plot
varImpPlot(randomForest_rf)

#model validation
set.seed(1009)
irispred <-predict(randomForest_rf,newdata = testdata)
freq1 <- table(irispred,testdata$BMI)
#plot(margin(randomForest_rf,testdata$BMI))
#mixed matrix
caret::confusionMatrix(testdata$BMI, irispred) #AUC=0.7841

irispred1 <- predict(randomForest_rf,newdata = testdata, type = "pro")
rocc1 <- roc(testdata$BMI, as.matrix(irispred1[ , 1]),ci=T,auc=T, direction=">")
rocc1

png("plot/RF-pro.png", width=20, height=20, units = "cm", res=600) 
plot(rocc1, print.auc = TRUE,
     print.auc.adj = (0.05), 
     print.auc.x = 0.63,
     print.auc.y = 0.3,
     print.auc.col = "#E72929",
     print.auc.cex =2, 
     auc.polygon = TRUE,
     max.auc.polygon = TRUE,
     auc.polygon.col = "#F0F3FF",
     grid = c(0.1, 0.2),
     cex.lab = 2, ci.col = "red",
     grid.col = c("green", "red"))
dev.off()

#SVM
#chose parameters

library(e1071)
set.seed(1009)
svm <- svm(BMI ~ . , 
           data = traindata,
           kernel = 'radial',C = 1, 
           sigma = 0.06, probability = TRUE) #type ="C-classification"
# prediction
pre_svm <- predict(svm, newdata = testdata, probability = TRUE)
table(pre_svm,testdata$BMI)
#mixed matrix
caret::confusionMatrix(testdata$BMI, pre_svm) 
#roc
#predict can output different typrs, such as "decision.values",
#"probabilities"
pre_svm1 <- attr(pre_svm, "probabilities")
rocc2 <- roc(testdata$BMI,as.matrix(pre_svm1[ , 1]),ci=T,auc=T)
rocc2
auc(rocc2)

png("plot/SVM-pro.png", width=20, height=20, units = "cm", res=600) 
plot(rocc2, print.auc = TRUE,
     print.auc.adj = (0.05), 
     print.auc.x = 0.63,
     print.auc.y = 0.3,
     print.auc.col = "#E72929",
     print.auc.cex =2, 
     auc.polygon = TRUE,
     max.auc.polygon = TRUE,
     auc.polygon.col = "#F0F3FF",
     grid = c(0.1, 0.2),
     cex.lab = 2, ci.col = "red",
     grid.col = c("green", "red"))
dev.off()
#KNN
library(caret)
library(kknn)
set.seed(1009)
pre <- kknn(BMI ~ . ,traindata, testdata,
            distance = 2, k =9, 
            kernel= "triangular")
summary(pre)
fit <- predict(pre,newdata = testdata)
t<-table(testdata$BMI, fit)
sum(diag(t))/sum(t)
caret::confusionMatrix(testdata$BMI, fit)
fit1 <- predict(pre,newdata = testdata, type = "prob")
head(fit1)
rocc3 <- roc(testdata$BMI, as.matrix(fit1[ , 1]),ci=T,auc=T, direction=">")
rocc3
auc(rocc3)

png("plot/KNN-pro.png", width=20, height=20, units = "cm", res=600) 
plot(rocc3, print.auc = TRUE,
     print.auc.adj = (0.05), 
     print.auc.x = 0.63,
     print.auc.y = 0.3,
     print.auc.col = "#E72929",
     print.auc.cex =2, 
     auc.polygon = TRUE,
     max.auc.polygon = TRUE,
     auc.polygon.col = "#F0F3FF",
     grid = c(0.1, 0.2),
     cex.lab = 2, ci.col = "red",
     grid.col = c("green", "red"))
dev.off()

#adaboost
set.seed(1009)
library(adabag)
ada <-boosting(BMI ~.,data= traindata,mfinal=150, 
               control=rpart.control(maxdepth=3))
pred<-predict.boosting(ada,newdata=testdata, type = "class")
pred$confusion
ad <- as.factor(pred$class)
#mixed matrix
caret::confusionMatrix(testdata$BMI, ad)
rocc4 <- roc(testdata$BMI, data.frame(pred$prob)[,-1],ci=T,auc=T, direction="<")
rocc4
png("plot/ada-pro.png", width=20, height=20, units = "cm", res=600) 
plot(rocc4, print.auc = TRUE,
     print.auc.adj = (0.05), 
     print.auc.x = 0.63,
     print.auc.y = 0.3,
     print.auc.col = "#E72929",
     print.auc.cex =2, 
     auc.polygon = TRUE,
     max.auc.polygon = TRUE,
     auc.polygon.col = "#F0F3FF",
     grid = c(0.1, 0.2),
     cex.lab = 2, ci.col = "red",
     grid.col = c("green", "red")) #,print.thres = TRUE
dev.off()

#output result
# Create a vector of prediction results
predictions <- list("irispred", "pre_svm", "fit", "ad")

# Create an empty dataframe to store results
result <- data.frame(
  Value = numeric(),
  Model = character()
)
#AccuracyNull is  NoInformationRate
# Loop through each prediction
for (pred in predictions) {
  # Get confusion matrix
  conf_matrix <- confusionMatrix(testdata$BMI, get(pred))
  
  # Extract metrics
  accuracy <- conf_matrix$overall['Accuracy']
  kappa <- conf_matrix$overall['Kappa']
  sensitivity <- conf_matrix$byClass['Sensitivity']
  specificity <- conf_matrix$byClass['Specificity']
  ppv <- conf_matrix$byClass['Pos Pred Value']
  npv <- conf_matrix$byClass['Neg Pred Value']
  f1 <- conf_matrix$byClass['F1']
  NoInformationRate <- conf_matrix$overall['AccuracyNull']
  
  # Add results to dataframe
  result <- rbind(result, data.frame(
    Value = c( sensitivity, specificity, NoInformationRate, accuracy, kappa,ppv, npv, f1),
    Model = rep(pred, 8)
  ))
}

result$Value <- signif(result$Value, digits = 3)
# Transpose the result dataframe
result <- t(result)

write.csv(result, file = "result/mc_pro_result.csv")

#extract AUC
res <- data.frame(
  Model = c("rocc1", "rocc2", "rocc3", "rocc4"),#"rocc5", "rocc6", "rocc7", "rocc8"
  #"rocc1", "rocc2", "rocc3", "rocc4"
  AUC = numeric(4),
  CI_lower = numeric(4),
  CI_upper = numeric(4)
)

rocc_list <- list(rocc1, rocc2, rocc3, rocc4)#rocc5, rocc6, rocc7, rocc8
#rocc1, rocc2, rocc3, rocc4
for (i in 1:length(rocc_list)) {
  auc_value <- auc(rocc_list[[i]])
  ci <- ci(rocc_list[[i]])
  res$AUC[i] <- auc_value
  res$CI_lower[i] <- ci[1]
  res$CI_upper[i] <- ci[3]
}

res$AUC <- sprintf("%.3f", res$AUC)
res$CI_lower <- sprintf("%.3f", res$CI_lower)
res$CI_upper <- sprintf("%.3f", res$CI_upper)
res_with_format <- paste(res$AUC, " (", res$CI_lower, "-", res$CI_upper, ")", sep="")

res_table <- data.frame(Model = res$Model, AUC_CI = res_with_format)
write.csv(res_table, file = "result/mc_pro_auc.csv")
##################phenotype machine learning###################
#random forest
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9789589/
#use 10 fold cross validation to choose the best parameters
library(caret)
library(randomForest)
library(pROC)
set.seed(1009)
#default is 500 trees
randomForest_rf<-randomForest(BMI ~.,data = traindata, #BMI ~ . - age - sex
                              ntree = 500,
                              mtry = 2,
                              proximity = TRUE,
                              importance = TRUE)
table(predict(randomForest_rf),traindata$BMI)
randomForest_rf
#Plot error rates based on different numbers in the random forest
#plot(randomForest_rf)
#importance score
importance(randomForest_rf,type = 1)
#importance plot
varImpPlot(randomForest_rf)

#model validation
irispred<-predict(randomForest_rf,newdata = testdata)
freq1 <- table(irispred,testdata$BMI)
#plot(margin(randomForest_rf,testdata$BMI))
#mixed matrix
caret::confusionMatrix(testdata$BMI, irispred) #AUC=0.7841
set.seed(1009)
irispred1 <- predict(randomForest_rf,newdata = testdata, type = "pro")
rocc5 <- roc(testdata$BMI, as.matrix(irispred1[ , 1]),ci=T,auc=T, direction=">")
rocc5

png("plot/RF-pheno.png", width=20, height=20, units = "cm", res=600) 
plot(rocc5, print.auc = TRUE,
     print.auc.adj = (0.05), 
     print.auc.x = 0.63,
     print.auc.y = 0.3,
     print.auc.col = "#E72929",
     print.auc.cex =2, 
     auc.polygon = TRUE,
     max.auc.polygon = TRUE,
     auc.polygon.col = "#F0F3FF",
     grid = c(0.1, 0.2),
     cex.lab = 2, ci.col = "red",
     grid.col = c("green", "red"))
dev.off()

#SVM
#chose parameters
library(e1071)
set.seed(1009)
svm <- svm(BMI ~ . , 
           data = traindata,
           kernel = 'radial',C = 1, 
           sigma = 0.111, probability = TRUE) #type ="C-classification"
# prediction
pre_svm <- predict(svm, newdata = testdata, probability = TRUE)
table(pre_svm,testdata$BMI)
#mixed matrix
caret::confusionMatrix(testdata$BMI, pre_svm) 
#roc
#predict can output different typrs, such as "decision.values",
#"probabilities"
pre_svm1 <- attr(pre_svm, "probabilities")
rocc6 <- roc(testdata$BMI,as.matrix(pre_svm1[ , 1]),ci=T,auc=T)
rocc6
auc(rocc6)

png("plot/SVM-pheno.png", width=20, height=20, units = "cm", res=600) 
plot(rocc6, print.auc = TRUE,
     print.auc.adj = (0.05), 
     print.auc.x = 0.63,
     print.auc.y = 0.3,
     print.auc.col = "#E72929",
     print.auc.cex =2, 
     auc.polygon = TRUE,
     max.auc.polygon = TRUE,
     auc.polygon.col = "#F0F3FF",
     grid = c(0.1, 0.2),
     cex.lab = 2, ci.col = "red",
     grid.col = c("green", "red"))
dev.off()
#KNN
library(caret)
library(kknn)
set.seed(1009)
pre <- kknn(BMI ~ . ,traindata, testdata,
            distance = 2, k =9, 
            kernel= "triangular")
summary(pre)
fit <- predict(pre,newdata = testdata)
t<-table(testdata$BMI, fit)
sum(diag(t))/sum(t)
caret::confusionMatrix(testdata$BMI, fit)
fit1 <- predict(pre,newdata = testdata, type = "prob")
head(fit1)
rocc7 <- roc(testdata$BMI, as.matrix(fit1[ , 1]),ci=T,auc=T, direction=">")
rocc7
auc(rocc7)

png("plot/KNN-pheno.png", width=20, height=20, units = "cm", res=600) 
plot(rocc7, print.auc = TRUE,
     print.auc.adj = (0.05), 
     print.auc.x = 0.63,
     print.auc.y = 0.3,
     print.auc.col = "#E72929",
     print.auc.cex =2, 
     auc.polygon = TRUE,
     max.auc.polygon = TRUE,
     auc.polygon.col = "#F0F3FF",
     grid = c(0.1, 0.2),
     cex.lab = 2, ci.col = "red",
     grid.col = c("green", "red"))
dev.off()

#adaboost
set.seed(1009)
library(adabag)
ada <-boosting(BMI ~.,data= traindata,mfinal=100, 
               control=rpart.control(maxdepth=3))
pred<-predict.boosting(ada,newdata=testdata, type = "class")
pred$confusion

ad <- as.factor(pred$class)

#mixed matrix
caret::confusionMatrix(testdata$BMI, ad)
rocc8 <- roc(testdata$BMI, data.frame(pred$prob)[,-1],ci=T,auc=T, direction="<")
rocc8
png("plot/ada-pheno.png", width=20, height=20, units = "cm", res=600) 
plot(rocc8, print.auc = TRUE,
     print.auc.adj = (0.05), 
     print.auc.x = 0.63,
     print.auc.y = 0.3,
     print.auc.col = "#E72929",
     print.auc.cex =2, 
     auc.polygon = TRUE,
     max.auc.polygon = TRUE,
     auc.polygon.col = "#F0F3FF",
     grid = c(0.1, 0.2),
     cex.lab = 2, ci.col = "red",
     grid.col = c("green", "red")) #,print.thres = TRUE
dev.off()
#output result
# Create a vector of prediction results
predictions <- list("irispred", "pre_svm", "fit", "ad")

# Create an empty dataframe to store results
result <- data.frame(
  Value = numeric(),
  Model = character()
)
#AccuracyNull is  NoInformationRate
# Loop through each prediction
for (pred in predictions) {
  # Get confusion matrix
  conf_matrix <- confusionMatrix(testdata$BMI, get(pred))
  
  # Extract metrics
  accuracy <- conf_matrix$overall['Accuracy']
  kappa <- conf_matrix$overall['Kappa']
  sensitivity <- conf_matrix$byClass['Sensitivity']
  specificity <- conf_matrix$byClass['Specificity']
  ppv <- conf_matrix$byClass['Pos Pred Value']
  npv <- conf_matrix$byClass['Neg Pred Value']
  f1 <- conf_matrix$byClass['F1']
  NoInformationRate <- conf_matrix$overall['AccuracyNull']
  
  # Add results to dataframe
  result <- rbind(result, data.frame(
    Value = c( sensitivity, specificity, NoInformationRate, accuracy, kappa,ppv, npv, f1),
    Model = rep(pred, 8)
  ))
}

result$Value <- signif(result$Value, digits = 3)
# Transpose the result dataframe
result <- t(result)

write.csv(result, file = "result/mc_pheno_result.csv")

#extract AUC
res <- data.frame(
  Model = c("rocc5", "rocc6", "rocc7", "rocc8"),#"rocc5", "rocc6", "rocc7", "rocc8"
  #"rocc1", "rocc2", "rocc3", "rocc4"
  AUC = numeric(4),
  CI_lower = numeric(4),
  CI_upper = numeric(4)
)

rocc_list <- list(rocc5, rocc6, rocc7, rocc8)#rocc5, rocc6, rocc7, rocc8
#rocc1, rocc2, rocc3, rocc4
for (i in 1:length(rocc_list)) {
  auc_value <- auc(rocc_list[[i]])
  ci <- ci(rocc_list[[i]])
  res$AUC[i] <- auc_value
  res$CI_lower[i] <- ci[1]
  res$CI_upper[i] <- ci[3]
}

res$AUC <- sprintf("%.3f", res$AUC)
res$CI_lower <- sprintf("%.3f", res$CI_lower)
res$CI_upper <- sprintf("%.3f", res$CI_upper)
res_with_format <- paste(res$AUC, " (", res$CI_lower, "-", res$CI_upper, ")", sep="")

res_table <- data.frame(Model = res$Model, AUC_CI = res_with_format)
write.csv(res_table, file = "result/mc_pheno_auc.csv")

######################plot mc##############
dev.new()
#png("plot/mc_auc.png", width=60, height=30, units = "cm", res=600)
par(mfrow = c(2, 4), cex.lab = 1.5)  

plot(rocc8, print.auc = TRUE,
     print.auc.adj = (0.05), 
     print.auc.x = 0.63,
     print.auc.y = 0.3,
     print.auc.col = "#E72929",
     print.auc.cex =2, 
     auc.polygon = TRUE,
     max.auc.polygon = TRUE,
     auc.polygon.col = "#F0F3FF",
     grid = c(0.1, 0.2),
     cex.lab = 2, ci.col = "red",
     grid.col = c("green", "red"))


dev.off()
###################correationship between photype and proteins########
rm(list = ls())
data <- read.csv(file = "pheno_11_pro_809.csv", row.names = 1, check.names = FALSE)
com_data <- read.csv("pheno_pro_2045.csv", check.names = FALSE)
data1 <- data.frame(row.names = com_data$zz_nr,
                    BMI = com_data$u3tbmi,
                    WHR = com_data$u3twhrat,
                    age = com_data$u3talteru,
                    SBP = com_data$u3tsysmm,
                    DBP = com_data$u3tdiamm,
                    alcohol = com_data$u3talkkon,
                    FPG = com_data$u3tglukfast_n,
                    TG =  com_data$u3l_trin,
                    HDL = com_data$u3l_hdln,
                    LDL = com_data$u3l_ldln,
                    TCHO = com_data$u3l_choln) # CRP = com_data$u3h_crp
pro <- read.csv(file = "select_pro_agesex.csv", row.names = 1, check.names = FALSE)
data2 <- pro[, -c(1:3)]

cor.result <- corr.test(data2, data1, use = "pairwise.complete.obs", method = "spearman")
r <- cor.result$r
p <- cor.result$p

#write.csv(r, file = "result/pheno_pro_correlation_R.csv")
#write.csv(p, file = "result/pheno_pro_correlation_pvalue.csv")

p[ p >= (0.05)] <- NA
rowSums(!is.na(p))
colSums(!is.na(p))
x <- data.frame(r[p < 0.05])

r <- cor.result$r
p <- cor.result$p
a <- data.frame(
  x = factor(c(rownames(p)), 
             levels = c(rownames(p))),
  y = as.numeric(c(11 ,9,  8, 11,  7, 10, 9,  4,  6, 10, 9, 9,  8, 9, 10, 9))
)
b <- data.frame(
  x = factor(c(colnames(p)), 
             levels = c(colnames(p))),
  y = as.numeric(c(16,15,15,11,11,11,15,14,15, 9, 7))
)
dev.new()
model_type <- ggplot(a, aes(x = x, y = y, group = 1)) +
  geom_line( color="#FF8F8F",size= 2, alpha = 1) +
  geom_point(shape=21, color="black", fill="#FF8F8F", size=5) +
  theme_ipsum()+
  labs(x = "Proteins", y = "Number of Significant Risk factors")+
  theme_bw() +
  #geom_text(aes(label = y), vjust = -1.2, hjust = 0.45, size = 5)+ 
  theme(axis.title.x = element_text(size = 20),  # 调整 X 轴标签字体大小
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 18),  # 调整 X 轴标签字体大小
        axis.text.y = element_text(size = 18)) +
  scale_y_continuous(breaks = seq(0, 11, 1), limits = c(3, 11))
model_type
model_type1 <- ggplot(b, aes(x = x, y = y, group = 1)) +
  geom_line( color="#FF8F8F",size= 2, alpha = 1) +
  geom_point(shape=21, color="black", fill="#FF8F8F", size=5) +
  theme_ipsum()+
  labs(x = "Risk factors", y = "Number of Significant Proteins")+
  theme_bw() +
  #geom_text(aes(label = y), vjust = -1.2, hjust = 0.45, size = 5)+ 
  theme(axis.title.x = element_text(size = 20),  # 调整 X 轴标签字体大小
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 18),  # 调整 X 轴标签字体大小
        axis.text.y = element_text(size = 18)) +
  scale_y_continuous(breaks = seq(0, 16, 1), limits = c(3, 16))
model_type1

library(ggpubr)
num <- ggarrange(model_type1, model_type,
                 #labels = c("A", "B"), font.label = list(size = 20),
                 ncol = 1, nrow = 2, common.legend = FALSE,widths = c(4, 4, 4),
                 heights = c(4, 4, 4), align = "hv") 

num
ggsave("plot/pro_cor_num.png", num, dpi=600, width= 60, height= 45, unit="cm")
dev.off()
#col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))(10) #"#EE9988", "#FFFFFF",
col <- colorRampPalette(c("blue", "white", "red"))(10)
#col <- colorRampPalette(c("darkred","white","midnightblue"))(10)
dev.new()
png("plot/pheno_pro_cor_nocrp.png", width=27, height=21, units = "cm", res=600) 
p <- cor.result$p
corrplot(r, method = "square", p.mat = p, sig.level = 0.05,
         insig = "blank", tl.col = "black", col =col, tl.cex = 1.5)
dev.off()

########################enrichment#####################
library(org.Hs.eg.db)
#https://mp.weixin.qq.com/s/5tPFV_cX-QjFNMWEiK0SBg
#library("R.utils")
#R.utils::setOption("clusterProfiler.download.method", "auto")
#getOption("clusterProfiler.download.method")
#or
#devtools::install_github("YuLab-SMU/DOSE")
#devtools::install_github("YuLab-SMU/HDO.db")   
#devtools::install_github('YuLab-SMU/clusterProfiler')  
library(clusterProfiler)
library(HDO.db)
library(DOSE)
#install.packages("R.utils")
#library("R.utils")  #this part is more useful
#R.utils::setOption("clusterProfiler.download.method",'auto') 
#getOption("clusterProfiler.download.method")
library(openxlsx)#读取.xlsx文件
library(ggplot2)#柱状图和点状图
library(stringr)#transfer gene ID
library(enrichplot)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
library(stats)
##GO analysis of specified species, species abbreviation index table details
#http: //bioconductor.org/packages/release/BiocViews.html#___OrgDb
keytypes(org.Hs.eg.db)
#KEGG analysis of specified species, species abbreviation index table details
#http//www.genome.jp/kegg/catalog/org_list.html

library(msigdbr)
library(reactome.db)
library(ReactomePA)
library(graphite)


library(tidyverse)
library(clusterProfiler)
library(msigdbr)  #install.packages("msigdbr")
library(GSVA) 
library(GSEABase)
library(pheatmap)
library(limma)
library(BiocParallel)
library(tidyverse)  # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
library(ggthemes)
library(ggprism)
library(ggplot2)
library(stringr)

rm(list = ls())
setwd("C:/project/pro/onlypro")
data <- read.csv(file = "select_pro_agesex.csv", row.names = 1,check.names = FALSE) #don't , check.names = FALSE here
data <- colnames(data)[4:19] #delete sex and age colums
data
gene_sig <- bitr(data,fromType = 'SYMBOL',toType = 'ENTREZID',
                 OrgDb = org.Hs.eg.db)
#write.csv(gene_sig, file = "C:/project/3/gene_list.csv")
#GO 
GO <- enrichGO(gene_sig$ENTREZID,#GO analysis
               OrgDb = 'org.Hs.eg.db',
               keyType = "ENTREZID",#read gene ID
               ont = "ALL",#(ont is ALL means includes Biological Process,Cellular Component, Mollecular Function）
               pvalueCutoff = 0.05,#assign p value threshold
               qvalueCutoff = 0.05,#specify q value threshold
               pAdjustMethod = "bonferroni", 
               readable = T)
GO1 <- data.frame(GO)
GO0 <- clusterProfiler::simplify(GO)
GO00 <- data.frame(GO0)
#write.csv(GO00, file = "result/GO.csv")
#library(devtools)
#install_github("ievaKer/aPEAR")
library(aPEAR)
set.seed(1009)
enrichmentNetwork(GO@result)
#write.csv(GO00, file = "result/GO.csv")
barplot(GO0, split="ONTOLOGY",showCategory = 5)+facet_grid(ONTOLOGY~., scale="free")
dotplot(GO0, split="ONTOLOGY",showCategory = 5)+facet_grid(ONTOLOGY~., scale="free")#dotplot


#############KEEG#############
search_kegg_organism('hsa', by='kegg_code')
# p adjust method is "BH" and filter q value
KEGG <- enrichKEGG(gene_sig$ENTREZID,#KEGG analysis
                   organism = 'hsa',
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "bonferroni", 
                   qvalueCutoff = 0.05, 
                   use_internal_data = FALSE)

KEGG1 <- as.data.frame(KEGG)
#write.csv(KEGG1, file = "result/KEGG.csv")

barplot(KEGG,showCategory = 40,title = 'KEGG Pathway')
dotplot(KEGG)

enrichplot::cnetplot(GO,circular=FALSE,colorEdge = TRUE)#基因-通路关联网络图
enrichplot::cnetplot(KEGG,circular=FALSE,colorEdge = TRUE)#circluar为指定是否环化，基因过多时建议设置为FALSE
enrichplot::heatplot(GO,showCategory = 16)#基因-通路关联热图
enrichplot::heatplot(KEGG,showCategory = 10)
GO2 <- pairwise_termsim(GO0)
enrichplot::emapplot(GO2,showCategory = 16, color = "p.adjust", layout = "kk",
                     group_category = T, cex_category= , min_edge = 0.3,
                     nClustrer = NULL)


#tree pathways
edox2 <- pairwise_termsim(GO0)
png("plot/GO_tree.png", width=40, height=20, units = "cm", res=600) 
#treeplot(edox2)
tree <- treeplot(edox2, split = "ONTOLOGY", cluster.params = list(label_words_n = 3, n = 3,method ="ward.D"), 
                 showCategory = , hilight = T) 
#"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
dev.off()

#GO bubble
logis <- read.csv("logistic_basic_sig.csv",check.names = FALSE, row.names = 1)
log <- logis[rownames(logis) %in% data,]
genedata<-data.frame(ID=gene_sig$SYMBOL,logFC=log2(log$Estimate))

GO_Term<-data.frame(Category = GO00$ONTOLOGY,
                    ID = GO00$ID,Term = GO00$Description, 
                    Genes = gsub("/", ", ", GO00$geneID), 
                    adj_pval = GO00$p.adjust)
circ <- circle_dat(GO_Term, genedata)#改成goplot输入格式
png("plot/GObubble_id.png", width=20, height=20, units = "cm", res=600) 
bubble <- GOBubble(circ, labels = 5, ID = T, bg.col = T, table.legend = F, table.col = T, display = "single")
dev.off()

png("plot/GO.png", width=40, height=40, units = "cm", res=600) 
ggarrange(tree, bubble, labels = c("A", "B"),
          font.label = list(size = 20),
          ncol = 1, nrow = 2, common.legend = FALSE,widths = c(4, 4, 4),
          heights = c(4, 6, 4), align = "hv")
dev.off()

heatplot(GO0,  #前面的logFC向量
         showCategory = 60, #显示前8
         symbol = "dot", # 图形，可选"dot"
         pvalue = NULL,
         label_format = 30)

library(DOSE)
x <- enrichDO(gene          = gene_sig$ENTREZID,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = gene_sig$ENTREZID,
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)
data(geneList, package="DOSE")
geneList <- gene_sig$ENTREZID
library(enrichplot)
de <-  gene_sig$ENTREZID
edo <- enrichDO(de, ont           = "DO",
                pvalueCutoff  = 0.05,
                pAdjustMethod = "BH",
                
                minGSSize     = 10,
                maxGSSize     = 500,
                qvalueCutoff  = 0.05,
                readable      = FALSE)
dotplot(edo, showCategory=10) 

###################MR########################

