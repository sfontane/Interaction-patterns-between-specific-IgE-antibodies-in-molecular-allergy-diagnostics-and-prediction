#####################################################################################
#Machine learning to identify pairwise interactions between specific IgE antibodies and their association with asthma: a
#cross-sectional analysis within a population-based birth cohort
#Sara Fontanella, Clement Frainay, Clare S Murray, Angela Simpson, Adnan Custovic

#### CODE to run full analysis 
###  As seeds were not set in the primary analysis, results might vary due to the random initialisation of the classification models.

##########
#Load packages
##########
library(ggplot2)
library(igraph)
library(Pigengene)
library(NbClust)

##########
#Load data
##########
rm(list=ls())

MyData<-read.csv(file="Filtered_data11yrs.csv")
sIgE<-MyData[,1:44]
asthma<-MyData[,45]
###########################################
#Statistical grouping of allergen components and their connectivity structure: component clusters
###########################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#build adjacency on distance correlation matrix
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AdjMatIgE<-dcor.matrix(sIgE)
P<-matrix(NA,ncol(sIgE),ncol(sIgE))

for (i in 1:ncol(sIgE)){
  for (j in 1:ncol(sIgE)){
    if (j>i){
      P[i,j]<-dcor.ttest(sIgE[,i],sIgE[,j])$p.value
      P[j,i]<-P[i,j]
    }
  }
}
alpha = .05   #level of significance
AdjMatIgE[P>alpha]<-0
diag(AdjMatIgE)=0
#~~~~~~~~~~~~~~~
#Cluster
#~~~~~~~~~~~~~~~
h=0.4   #maximum dissimilarity between two element in the same cluster
fit1<-hclust(as.dist(1-AdjMatIgE), method="average")
groups <- cutree(fit1, h=h)

sort(table(groups))   #number of allergen components per cluster
groups   #cluster assignment of allergen components
###########################################
#Patterns of sensitisation among study participants: sensitisation clusters
###########################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#create bipartite network 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
bin.sIgE<-matrix(0,nrow(sIgE),ncol(sIgE))
t = .3   # sIgE response threshold to be considered as positive
bin.sIgE[sIgE>t]<-1

minC = 2   # minimal number of clusters
maxC = 6   # maximal number of clusters
NB<-NbClust(bin.sIgE,diss=NULL,distance="binary",min.nc=minC,max.nc=maxC,method="ward.D2",index="ch")
Clusters<-as.factor(NB$Best.partition)


####################################################################################
# Differential sIgE co-expression patterns in asthma & CLASSIFICATION (random seed)
####################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~JDINAC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(caret)
library(e1071)
library(pROC)

log_sIgE<-as.data.frame(cbind(log(sIgE+1),asthma))

#  unbalanced dataset correction
weights <-asthma
weights[weights==0] <-.3
weights[weights==1] <-.7

# extract edges list
AdjMatIgE[lower.tri(AdjMatIgE,diag=T)] <-0
EDGE <- which(AdjMatIgE!=0, arr.ind=T) 

source("jdinac_weighted.R")

difnet <- jdinac_weighted(EDGE=EDGE,classLabel=asthma,DataFit=log_sIgE,weight=weights,nsplit=50,nfolds=10)

Prob.difnet<-difnet$yPre
Hard.cl.difnet<-NULL
Hard.cl.difnet[Prob.difnet>.5]<-1
Hard.cl.difnet[Prob.difnet<=.5]<-0

#confusion matrix
Conf.matrix.JDINAC<-caret::confusionMatrix(as.factor(Hard.cl.difnet),as.factor(asthma),"1")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~LOGISTIC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(glmnet)

log_sIgE$asthma[log_sIgE$asthma==0]<-"no"
log_sIgE$asthma[log_sIgE$asthma==1]<-"yes"

# define training control
train_control<- trainControl(method= "repeatedcv", number=50,repeats=10,classProbs = TRUE,summaryFunction=twoClassSummary,savePred =T)
# train the model 
tuneGrid=expand.grid(alpha=1,lambda=0.02) #parameters chosen with cv.
model<- train(asthma~., data=log_sIgE, trControl=train_control, method="glmnet",family="binomial",
              metric = "ROC",tuneGrid =tuneGrid,weights = weights)

#posterior probability computation
fitpred <- model$pred
Prob.log.reg=matrix(NA,213,1)
for (i in 1:213){
  Temp=model$pred[model$pred$rowIndex==i,]
  Prob.log.reg[i,]=mean(Temp[,5])
}

Hard.cl.logreg<-NULL
Hard.cl.logreg[Prob.log.reg>.5]<-1
Hard.cl.logreg[Prob.log.reg<=.5]<-0

#confusion matrix
Conf.matrix.logreg<-confusionMatrix(as.factor(Hard.cl.logreg),as.factor(asthma),"1")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Performance measures~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(hmeasure)

Perf.JDinac<-HMeasure(asthma, Prob.difnet, threshold=0.5, level=0.95)$metrics
Perf.logreg<-HMeasure(asthma, Prob.log.reg, threshold=0.5, level=0.95)$metrics

