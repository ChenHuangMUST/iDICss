# Grouping samples using consistent clustering
library(ConsensusClusterPlus)
d <- t(scale(t(iDC_matrix),center = T,scale = T))
range(d)
set.seed(1)
results<-ConsensusClusterPlus(as.matrix(d),maxK=8,reps=1000,
                              pItem=0.80,pFeature=1,
                              innerLinkage="ward.D",
                              finalLinkage="ward.D",
                              title='consensus',
                              clusterAlg="pam",
                              distance="pearson",
                              plot="pdf")

classes<-results[[2]]$consensusClass
sample_c1<-names(classes)[which(classes==1)]#Samples classified into class1
sample_c2<-names(classes)[which(classes==2)]#Samples classified into class2
sur1$con_class <- ifelse(rownames(sur1)%in%sample_c1,"class1",
                         ifelse(rownames(sur1)%in%sample_c2,"class2",NA))

# Filtering for differential genes using limma
{library(tidyr)
  library(dplyr)
  library(limma)}
group_list <- as.factor(classes)
design <- model.matrix(~0+factor(group_list))
colnames(design) <- paste("class",levels(group_list),sep = "")
rownames(design) <- names(classes)
# This contrast matrix statement
contrast.matrix <- makeContrasts("class1-class2",
                                 levels = design)
# step1 
fit <- lmFit(exp2,design)
# step2 
fit2 <- contrasts.fit(fit, contrast.matrix) 
# step3  
fit2 <- eBayes(fit2) 
# step4 
tempOutput <- topTable(fit2,adjust="BH",sort.by="B",n=Inf)
DEG <- tempOutput[which(abs(tempOutput$logFC)>1&
                          tempOutput$adj.P.Val<0.05),]
exp_deg <- exp[rownames(DEG),]


# Filtering for important genes using Boruta
library(Boruta)
data_Boruta <- cbind.data.frame(classes,as.data.frame(t(exp_deg)))
set.seed(1)
res_Boruta <- Boruta(classes ~ ., data = data_Boruta, 
                     doTrace = 2, ntree = 1500)
# Median Z is used to adjust the classification of unrecognized factors (important or unimportant) 
zmed <- TentativeRoughFix(res_Boruta)
#getConfirmedFormula(res_Boruta)
sig_gene1 <- attStats(zmed)
sig_gene2 <- rownames(sig_gene1)[which(sig_gene1$decision=="Confirmed")]


# Filtering of signature genes using lasso-cox
DEG2 <- as.data.frame(t(exp_deg))
sur2 <- cbind(sur1[,c("status","time")],
              DEG2[match(rownames(sur1),rownames(DEG2)),sig_gene2])
library(glmnet)
set.seed(2)
lasso_x <- as.matrix(sur2[,-c(1,2)])
lasso_y <- as.matrix(sur2[,c("status","time")])
fit <- glmnet(lasso_x,lasso_y,family = "cox")
pdf("lasso coefficients.pdf",onefile = F,width = 7,height = 7)
plot(fit, label=TRUE, xvar="lambda")
dev.off()
lasso.cv <- cv.glmnet(lasso_x,lasso_y,family="cox",alpha=1,
                      type.measure="deviance")
pdf("lasso Partial Likelihood Deviance.pdf",onefile = F,width = 7,height = 7)
plot(lasso.cv)
dev.off()
lasso.cv$lambda.min
lasso.cv$lambda.1se
l.coef1<-coef(lasso.cv$glmnet.fit,s=lasso.cv$lambda.min,exact = F)
l.coef2<-coef(lasso.cv$glmnet.fit,s=lasso.cv$lambda.1se,exact = F)
coef_lasso <- as.data.frame(as.matrix(coef(lasso.cv,s=lasso.cv$lambda.min)))
gene_lasso <- rownames(coef_lasso)[which(coef_lasso$`1`!=0)]

coef_lasso1 <- coef_lasso[which(coef_lasso$`1`!=0),1]
names(coef_lasso1) <- gene_lasso

# Constructing a scoring system
exp_gene_lasso <- exp_deg[gene_lasso,]
risk_score <- apply(exp_gene_lasso,2,function(x){x*coef_lasso1})
my_coef <- coef_lasso1
risk_score1 <- data.frame(risk_score = colSums(risk_score))
risk_score2 <- cbind.data.frame(risk_score1,sur1[match(rownames(risk_score1),rownames(sur1)),])
risk_score2$risk_class <- ifelse(risk_score2$risk_score >= median(risk_score2$risk_score),"high","low")

