## RF Methyl classification VUmc
library(caret)
library(randomForest)
library(doMC)
library(e1071)
library(minfi)

###PART 1: make new probe lists
##1 co-del vs non co-del probes in VUmc
# load VUmc data
load("/projects/verbun/VUMC/Rdata/gRatioSet.RData")
# select only patients with IDH mutation 
gRatioSet.MT <- gRatioSet[,gRatioSet$IDH == "MT"]
#gRatioSet.MT <- gRatioSet[,gRatioSet$Patient == 4]

VUmc.MT <- getBeta(gRatioSet.MT)
colnames(VUmc.MT) <- gRatioSet.MT$Sample_Name
# make new probe list and dataset with B values from 1308 probes VUmc
IDH1308Probes <- read.csv("/projects/verbun/BACKUP/JAX/Data/Other_data/Probes/1308IDHmutantProbes.csv", sep="/", header=FALSE, stringsAsFactors=FALSE)
keep1308 <- (rownames(VUmc.MT) %in% IDH1308Probes[[1]])
VUmc.MT.probes.df <- VUmc.MT[keep1308,]
VUmc.MT.probes.df <- t(VUmc.MT.probes.df)
VUmc.MT.probes.names <- colnames(VUmc.MT.probes.df)

##2 GIMP low vs GIMP high probes in VUmc
# make make new probe list and dataset with B values from 163 probes VUmc
GIMP163Probes <- read.csv("/projects/verbun/BACKUP/JAX/Data/Other_data/Probes/163GIMPProbes.csv", sep="/", header=FALSE, stringsAsFactors=FALSE)
keep163 <- (rownames(VUmc.MT) %in% GIMP163Probes[[1]])
VUmc.GIMP.probes.df <- VUmc.MT[keep163,]
VUmc.GIMP.probes.df <- t(VUmc.GIMP.probes.df)
VUmc.GIMP.probes.names <- colnames(VUmc.GIMP.probes.df)

##3 WT probes in VUmc
# select only patients with IDH wild type 
gRatioSet.WT <- gRatioSet[,gRatioSet$Patient == "WT"]
VUmc.WT <- getBeta(gRatioSet.WT)
colnames(VUmc.WT) <- gRatioSet.WT$Sample_Name
# make make new probe list and dataset with B values from 914 probes VUmc
WT914Probes <- read.csv("/projects/verbun/BACKUP/JAX/Data/Other_data/Probes/914IDHwtProbes.csv", sep="/", header=FALSE, stringsAsFactors=FALSE)
keep914 <- (rownames(VUmc.WT) %in% WT914Probes[[1]])
VUmc.WT.probes.df <- VUmc.WT[keep914,]
VUmc.WT.probes.df <- t(VUmc.WT.probes.df)
VUmc.WT.probes.names <- colnames(VUmc.WT.probes.df)
  
### PART 2: make new trainingsets 
##1 MT set
load("/projects/verbun/BACKUP/JAX/Data/Rdata/RF/trainingdata.1308.RData")
keep.MT <- (colnames(trainingdata.1308) %in% VUmc.MT.probes.names)
training.MT <- trainingdata.1308[,keep.MT]
training.MT$cartoon <- trainingdata.1308$cartoon
##2 GIMP set
load("/projects/verbun/BACKUP/JAX/Data/Rdata/RF/trainingdata.163.RData")
keep.GIMP <- (colnames(trainingdata.163) %in% VUmc.GIMP.probes.names)
training.GIMP <- trainingdata.163[,keep.GIMP]
training.GIMP$cartoon <- trainingdata.163$cartoon
##3 WT set
load("/projects/verbun/BACKUP/JAX/Data/Rdata/RF/trainingdata.914.RData")
keep.WT <- (colnames(trainingdata.914) %in% VUmc.WT.probes.names)
training.WT <- trainingdata.914[,keep.WT]
training.WT$cartoon <- trainingdata.914$cartoon

### PART 3: train random forest with new TCGA trainingsets and predict VUmc data 
## A) MT
trainingdata <- training.MT
# register cores for doMC (8)
registerDoMC(cores = 1) #Number of CPUs
# set up k-fold cross validation
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)

# Set your seed so your work is repeatable
set.seed(42)
# Create a subset of your data to train your model on.  This makes sure you have
# equal representation of the 'cartoon' groups in your training set
#inTraining <- createDataPartition(trainingdata$cartoon, p=0.8, list=FALSE, times=1) #using 80% of data to train and the remaining 20% for testing
# Training Set
#myTrain <- trainingdata[inTraining, ]
# Testing Set
#myTest <- trainingdata[-inTraining, ]
# Confirm seed is set
#set.seed(210)
# set values for mtry
# mtry is the "Number of variables randomly sampled as candidates at each split"
# traditionally for classification you use the sqrt of the number of variables
# but here we try a range of mtry values to find the best parameters for our model
mtryVals <- floor(c(seq(100, 2000, by=100),
                    sqrt(ncol(trainingdata))))
mtryGrid <- data.frame(.mtry=mtryVals)
# Confirm seed again
set.seed(420)
# Set number of cores (8)
registerDoMC(cores = 1) #Number of CPUs
# Run Training
TCGA.MT.RF <- train(cartoon ~ ., # variable to be trained on
                    data = trainingdata, # Data we are using
                    method = "rf", # Method we are using
                    trControl = fitControl, # How we validate
                    # We created this object above
                   ntree = 5000, # number of trees
                    # is dependent on training data size
                    importance = TRUE, # calculate varible importance
                    # can be omitted to speed up calc
                    tuneGrid = mtryGrid, # set mtrys
#                    subset = inTraining # define training set
                    allowParallel = TRUE # run parallel cores 
)

save(TCGA.MT.RF, file ="/projects/verbun/VUMC/RF/TCGA.MT.RF.RData")
# make sure you have your object that you wish to reclassify
# this object must have only the features that you have trained on, and must have
# all the features that are trained on.

# predict the probability with which each sample of the test set is classified
#TCGA.MT.RF.pred <- predict(TCGA.MT.RF, myTest, type="prob")
# record the best classification in the test set
#myTest$RF.classification <- predict(TCGA.MT.RF, myTest)
# show the confusion matrix
#confusionMatrix(data = myTest$RF.classification, reference = myTest$cartoon)

# predict the classification of the new data set (new gliomas in this case)
# this data contains the same common set of probes
TCGA.MT.RF.pred <- predict(TCGA.MT.RF, VUmc.MT.probes.df, type="prob")
# record the best classification
VUmc.MT.probes.df$RF.classification <- predict(TCGA.MT.RF, VUmc.MT.probes.df)

# show how the new data falls into your classification groups
result <- table(VUmc.MT.probes.df$RF.classification)

# save everything
save(VUmc.MT.probes.df, file="/projects/verbun/VUMC/RF/VUmc.MT.output.RData")
save(result, file="/projects/verbun/VUMC/RF/VUmc.MT.result.RData")
save(TCGA.MT.RF.pred, file="/projects/verbun/VUMC/RF/TCGA.MT.RF.pred.RData")

#######################

##B) GIMP
trainingdata <- training.GIMP
# register cores for doMC (8)
registerDoMC(cores = 1) #Number of CPUs
# set up k-fold cross validation
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)

# Set your seed so your work is repeatable
set.seed(42)
# Create a subset of your data to train your model on.  This makes sure you have
# equal representation of the 'cartoon' groups in your training set
#inTraining <- createDataPartition(trainingdata$cartoon, p=0.8, list=FALSE, times=1) #using 80% of data to train and the remaining 20% for testing
# Training Set
#myTrain <- trainingdata[inTraining, ]
# Testing Set
#myTest <- trainingdata[-inTraining, ]
# Confirm seed is set
set.seed(210)
# set values for mtry
# mtry is the "Number of variables randomly sampled as candidates at each split"
# traditionally for classification you use the sqrt of the number of variables
# but here we try a range of mtry values to find the best parameters for our model
mtryVals <- floor(c(seq(100, 2000, by=100),
                    sqrt(ncol(trainingdata))))
mtryGrid <- data.frame(.mtry=mtryVals)
# Confirm seed again
set.seed(420)
# Set number of cores (8)
registerDoMC(cores = 1) #Number of CPUs
# Run Training
TCGA.GIMP.RF <- train(cartoon ~ ., # variable to be trained on
                    data = trainingdata, # Data we are using
                    method = "rf", # Method we are using
                    trControl = fitControl, # How we validate
                    # We created this object above
                    ntree = 5000, # number of trees
                    # is dependent on training data size
                    importance = TRUE, # calculate varible importance
                    # can be omitted to speed up calc
                    tuneGrid = mtryGrid # set mtrys
#                    subset = inTraining # define training set
)

save(TCGA.GIMP.RF, file ="/projects/verbun/VUMC/RF/TCGA.GIMP.RF.RData")
# make sure you have your object that you wish to reclassify
# this object must have only the features that you have trained on, and must have
# all the features that are trained on.

# predict the probability with which each sample of the test set is classified
#TCGA.GIMP.RF.pred <- predict(TCGA.GIMP.RF, myTest, type="prob")
# record the best classification in the test set
#myTest$RF.classification <- predict(TCGA.GIMP.RF, myTest)
# show the confusion matrix
#confusionMatrix(data = myTest$RF.classification, reference = myTest$cartoon)

# predict the classification of the new data set (new gliomas in this case)
# this data contains the same common set of probes
TCGA.GIMP.RF.pred <- predict(TCGA.GIMP.RF, VUmc.GIMP.probes.df, type="prob")
# record the best classification
VUmc.GIMP.probes.df$RF.classification <- predict(TCGA.GIMP.RF, VUmc.GIMP.probes.df)

# show how the new data falls into your classification groups
result <- table(VUmc.GIMP.probes.df$RF.classification)

# save everything
save(VUmc.GIMP.probes.df, file="/projects/verbun/VUMC/RF/VUmc.GIMP.output.RData")
save(result, file="/projects/verbun/VUMC/RF/VUmc.GIMP.result.RData")
save(TCGA.GIMP.RF.pred, file="/projects/verbun/VUMC/RF/TCGA.GIMP.RF.pred.RData")

#################

##C) WT
trainingdata <- training.WT
# register cores for doMC (8)
registerDoMC(cores = 1) #Number of CPUs
# set up k-fold cross validation
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)

# Set your seed so your work is repeatable
set.seed(42)
# Create a subset of your data to train your model on.  This makes sure you have
# equal representation of the 'cartoon' groups in your training set
#inTraining <- createDataPartition(trainingdata$cartoon, p=0.8, list=FALSE, times=1) #using 80% of data to train and the remaining 20% for testing
# Training Set
#myTrain <- trainingdata[inTraining, ]
# Testing Set
#myTest <- trainingdata[-inTraining, ]
# Confirm seed is set
set.seed(210)
# set values for mtry
# mtry is the "Number of variables randomly sampled as candidates at each split"
# traditionally for classification you use the sqrt of the number of variables
# but here we try a range of mtry values to find the best parameters for our model
mtryVals <- floor(c(seq(100, 2000, by=100),
                    sqrt(ncol(trainingdata))))
mtryGrid <- data.frame(.mtry=mtryVals)
# Confirm seed again
set.seed(420)
# Set number of cores (8)
registerDoMC(cores = 1) #Number of CPUs
# Run Training
TCGA.WT.RF <- train(cartoon ~ ., # variable to be trained on
                      data = trainingdata, # Data we are using
                      method = "rf", # Method we are using
                      trControl = fitControl, # How we validate
                      # We created this object above
                      ntree = 5000, # number of trees
                      # is dependent on training data size
                      importance = TRUE, # calculate varible importance
                      # can be omitted to speed up calc
                      tuneGrid = mtryGrid # set mtrys
                      #subset = inTraining # define training set
)

save(TCGA.WT.RF, file ="/projects/verbun/VUMC/RF/TCGA.WT.RF.RData")
# make sure you have your object that you wish to reclassify
# this object must have only the features that you have trained on, and must have
# all the features that are trained on.

# predict the probability with which each sample of the test set is classified
#TCGA.WT.RF.pred <- predict(TCGA.WT.RF, myTest, type="prob")
# record the best classification in the test set
#myTest$RF.classification <- predict(TCGA.WT.RF, myTest)
# show the confusion matrix
#confusionMatrix(data = myTest$RF.classification, reference = myTest$cartoon)

# predict the classification of the new data set (new gliomas in this case)
# this data contains the same common set of probes
TCGA.WT.RF.pred <- predict(TCGA.WT.RF, VUmc.WT.probes.df, type="prob")
# record the best classification
VUmc.WT.probes.df$RF.classification <- predict(TCGA.WT.RF, VUmc.WT.probes.df)

# save everything
save(VUmc.WT.probes.df, file="/projects/verbun/VUMC/RF/VUmc.WT.output.RData")
save(result, file="/projects/verbun/VUMC/RF/VUmc.WT.result.RData")
save(TCGA.WT.RF.pred, file="/projects/verbun/VUMC/RF/TCGA.WT.RF.pred.RData")



##### patient wise GIMP vs codel
# this data contains the same common set of probes
#load("/projects/verbun/VUMC/RF/TCGA.MT.RF.RData")
#TCGA.P04.RF.pred <- predict(TCGA.MT.RF, VUmc.MT.probes.df, type="prob")
# record the best classification
#VUmc.MT.probes.df$RF.classification <- predict(TCGA.MT.RF, VUmc.MT.probes.df)

##### patient wise GIMP high vs low
# this data contains the same common set of probes
#load("/projects/verbun/VUMC/RF/TCGA.GIMP.RF.RData")
#TCGA.P04.GIMP.RF.pred <- predict(TCGA.GIMP.RF, VUmc.GIMP.probes.df, type="prob")
# record the best classification
#VUmc.GIMP.probes.df$RF.classification <- predict(TCGA.GIMP.RF, VUmc.GIMP.probes.df)


### patient wise WT
# predict the classification of the new data set (new gliomas in this case)
# this data contains the same common set of probes
#load("/projects/verbun/VUMC/RF/TCGA.WT.RF.RData")
#TCGA.WT.RF.pred <- predict(TCGA.WT.RF, VUmc.WT.probes.df, type="prob")
# record the best classification
#VUmc.WT.probes.df$RF.classification <- predict(TCGA.WT.RF, VUmc.WT.probes.df)


# show how the new data falls into your classification groups
#result <- table(VUmc.WT.probes.df$RF.classification)

# predict the classification of the new data set (new gliomas in this case)
# this data contains the same common set of probes
#TCGA.WT.RF.pred <- predict(TCGA.WT.RF, VUmc.WT.probes.df, type="prob")
# record the best classification
#VUmc.WT.probes.df$RF.classification <- predict(TCGA.WT.RF, VUmc.WT.probes.df)


# save everything
#save(VUmc.WT.probes.df, file="/projects/verbun/VUMC/RF/VUmc.WT.output.RData")
#save(result, file="/projects/verbun/VUMC/RF/VUmc.WT.result.RData")

#########

#save(list=ls(), file="/projects/verbun/VUMC/RF/RF.all.Rdata")


# show how the new data falls into your classification groups
#result$Class.163 <- VUmc.163$RF.classification
#result$Class <- ifelse((result$Class.1308 %in% c("IDHmut-K1","IDHmut-K2") & result$Class.163 == "G-CIMP-low"), "G-CIMP-low",
#                       ifelse((result$Class.1308 %in% c("IDHmut-K1","IDHmut-K2") & result$Class.163 == "G-CIMP-high"), "G-CIMP-high",
#                              "Codel"))
#result[nrow(result),]$Class <- NA
#setwd("/projects/verbun/VUMC/MethClass/")
#save(result, file ="VUmc.result.RData")
#gRatioSet$LGm.RF <- result$Class
#setwd("/projects/verbun/VUMC/MethClass/")
#save(gRatioSet, file="gRatioSet.RData")

