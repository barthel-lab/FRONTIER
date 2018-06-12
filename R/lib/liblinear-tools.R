
library(LiblineaR)
library(caret)
library(ROCR)

run_pred <- function(data_matrix, data_labels, k = 10) {
  
  stopifnot(k > 1)
  
  folds = createFolds(data_labels, k=k)
  
  kfoldcrossvalidation = lapply(names(folds), function(fold_name) {
    message(fold_name)
    fold = folds[[fold_name]]
    
    train_data  = data_matrix[-fold, , drop=FALSE]
    train_targ  = factor(data_labels[-fold], levels = unique(data_labels))
    
    test_data   = data_matrix[fold, , drop=FALSE]
    test_targ   = factor(data_labels[fold], levels = unique(data_labels))
    
    # Re-train best model with best cost value.
    m = LiblineaR(data=train_data, target=train_targ, type=0, cost=1000, bias=T, verbose=F)
    
    # features = as.data.frame(t(m$W)) %>% 
    #  mutate(region = colnames(m$W), fold = fold_name) %>% 
    #  select(region, fold, everything())
    
    # %>%
    #  gather(key = "target", value = "weight", -fold, -region)
    
    # Make prediction
    p = predict(m, test_data, proba=T, decisionValues=T)
    
    decisions = as.data.frame(p$decisionValues[,unique(data_labels)]) %>%
      mutate(true_label = as.character(test_targ), pred_label = as.character(p$predictions), fold = fold_name, sample = rownames(test_data)) %>%
      select(sample, fold, true_label, pred_label, everything()) 
    
    #%>%
    #  gather(everything(), -pred_label, -true_label, -fold, -sample, -truth, key = "target", value = "decision")
    
    return(decisions)
  })
  
  decisions = data.table::rbindlist(kfoldcrossvalidation) %>% as.data.frame()
  decisions$truth = decisions$true_label == decisions$pred_label
  decisions$decision_value = as.numeric(apply(decisions, 1, function(x) x[x['pred_label']]))
  
  return(decisions)
}


multi_roc <- function(decision_values, categ) {
  if(length(categ) > 1) {
    pred = prediction(as.numeric(sapply(1:length(categ), function(i) decision_values[i, categ[i]])), decision_values$true_label == categ)
    n = length(categ)
  }
  else {
    pred = prediction(as.numeric(decision_values[,categ]), decision_values$true_label == categ)
    n = sum(decision_values$true_label == categ)
  }
  perf = performance(pred, measure = "tpr", x.measure = "fpr")
  auc = performance(pred, measure = "auc")
  return(data.frame(fpr = perf@x.values[[1]], tpr = perf@y.values[[1]], categ = ifelse(length(categ) > 1, "Overall", categ), n = n, auc = unlist(auc@y.values)))
}


# # Find the best model with the best cost parameter via 10-fold cross-validations
# tryTypes = c(0)
# tryCosts = c(1000, 1, 0.001)
# bestCost = NA
# bestAcc = 0
# bestType = NA
# 
# for(ty in tryTypes) {
#   for(co in tryCosts) {
#     
#     acc=LiblineaR(data = tcgameth, target = target, type = ty, cost = co, bias = 1, cross = 5, verbose = FALSE)
#     message("Results for C=",co," : ",acc," accuracy.\n",sep="")
#     
#     if(acc>bestAcc) {
#       bestCost=co
#       bestAcc=acc
#       bestType=ty
#     }
#   }
# }
# 
# cat("Best model type is:",bestType,"\n")
# cat("Best cost is:",bestCost,"\n")
# cat("Best accuracy is:",bestAcc,"\n")