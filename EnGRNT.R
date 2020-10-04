# EnGRNT to GRN Inferece 

################ load necessary libraries ############################################################
library(igraph)
library(readr)
library(e1071) 
library(ebmc)
library(pROC)

################ read input data #####################################################################

# read gene expression profiles in a variety of knockdown, knockout, and multifactorial 
expression_path = read_tsv("./Data/Ecoli-1_knockouts.tsv",col_names = TRUE)
expression_data = as.matrix(expression_path) 

# read known regulatory interactions between TF and target gene
regulation_path = read_tsv("./Data/Ecoli-1_goldstandard.tsv",col_names = FALSE)
regulation_data = as.matrix(regulation_path)

#----------------------------------------------------------------------------------
#Prepare Regulation Data

total_regulation = matrix(0,nrow(regulation_data),3)

for(i in 1:ncol(expression_data)){
  for(j in 1:nrow(regulation_data)){
    
    if(colnames(expression_data)[i]==regulation_data[j,1]){
      total_regulation[j,1] = i
    }
    if(colnames(expression_data)[i]==regulation_data[j,2]){
      total_regulation[j,2] = i 
    }
    
  }
}

total_regulation[,3] = regulation_data[,3]
#------------------------------------------------------------------------------------
#Required Variables

exp_size = ncol(expression_data)
no_of_exp = nrow(expression_data)

gene_number = matrix(1:exp_size)
expression_trans = t(expression_data)
expression_trans = cbind(gene_number,expression_trans)
expression_trans = as.data.frame(expression_trans)
colnames(expression_trans)[1] = "g.Number"
no_of_exp =no_of_exp+1

test_size= round(0.3* no_of_exp)

expressed_regulation = subset(total_regulation,total_regulation[,3]>0)
TF_list = unique(expressed_regulation[,1])

list_lables = matrix(0,nrow(expression_trans),length(TF_list))

test_predictions = matrix(0,test_size,length(TF_list))

AUROC = matrix(0,1,length(TF_list))
#AUPR = matrix(0,1,length(TF_list))

index=list()

for(k in 1:length(TF_list)){
  row2 = which(expressed_regulation[,1]==TF_list[k],arr.ind = TRUE)
  #print(row2)
  index[[k]]=as.list(row2)
  #lables
  labs = matrix(0,nrow = exp_size,ncol = 1,byrow = TRUE)
  
  for(m in 1:length(index[[k]])){
    row_labs = as.numeric(expressed_regulation[index[[k]][[m]],2])
    labs[row_labs,1]=1
    
  }
  #local_list = expression_trans
  local_list = cbind(expression_trans,labs)
  #print(labs)
  local_list <- local_list[sample(1:nrow(local_list)), ]
  
  list_lables[,k] = local_list[,-1:-no_of_exp]
  
  training_set = local_list[(test_size+1):exp_size,1:no_of_exp]
  test_set = local_list[1:test_size,1:no_of_exp]
  
  nrow_train = nrow(training_set)
  
  test_set = as.data.frame(test_set)
  test_lable = list_lables[1:test_size,k]
  
  #prepare train set and test set
  svm_lable = list_lables[(test_size+1):exp_size,k]
  
  w = cbind(training_set,svm_lable,deparse.level = 1)
  w = as.data.frame(w)
  ncol_w = ncol(w)
  colnames(w)[ncol_w]= "lable1"
  w$lable1 <- factor(w$lable1, levels = c("0", "1"), labels = c("0", "1"))
  data1 = w[1:nrow_train,]
  
  
  #generate model for svm learn
  svm.model <- ub(lable1 ~ ., data = data1, size = 20,ir = 1, alg = "svm",
                  svm.ker = "radial")
  
  pred =  predict(svm.model,newdata = test_set, type = "class",
                  decision.values = TRUE)
  
  #prediction
  pred= factor(pred,levels = c("0", "1"), labels = c("0", "1"))
  mat <- matrix( as.numeric(as.character(unlist(pred))), nrow=length(pred) )
  test_predictions[,k]= mat[1:test_size,]
  
  #Confusion matrix result of prediction
  #print(table(pred,test_lable))
  
  AUROC[,k] = auc(pred, test_lable)
  
  #AUROC[is.nan(AUROC)] <- 0
}


mean_AUROC = sum(AUROC)/length(TF_list)

mean_AUROC = round(mean_AUROC,digits = 4)

print(mean_AUROC)














