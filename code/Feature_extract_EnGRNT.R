#Feature Extraction From GRN

library(igraph)
library(readr)

################ read input data ###################################################################

# read gene expression profiles in a variety of knockdown, knockout and multifactorial 
expression_path = read_tsv("./Data/Ecoli-1_knockout.tsv",col_names = TRUE)
expression_data = as.matrix(expression_path) 

# read known regulatory interactions between TF and target gene
regulation_path = read_tsv("./Data/Ecoli-1_goldstandard.tsv",col_names = FALSE)
regulation_data = as.matrix(regulation_path)

#---------------------------------------------------------------------------------------------------
#Prepare Regulation Data

total_regulation = matrix(0,nrow(regulation_data),3)

for(i in 1:ncol(expression_data)){
  for(j in 1:nrow(regulation_data)){
    
    if(colnames(expression_data)[i] == regulation_data[j,1]){
      total_regulation[j,1] = i
    }
    if(colnames(expression_data)[i] == regulation_data[j,2]){
      total_regulation[j,2] = i 
    }
    
  }
}

total_regulation[,3] = regulation_data[,3]

exp_size = ncol(expression_data)
no_of_exp = nrow(expression_data)

#test_size= round(0.3* no_of_exp)
#expression_trans = t(expression_data)

expressed_regulation = subset(total_regulation,total_regulation[,3]>0)
TF_list = unique(expressed_regulation[,1])

expressed_reg = subset(regulation_data,regulation_data[,3]>0)

grn_net = graph_from_data_frame(expressed_reg,directed = TRUE)

nodes_name <- data.frame(name = V(grn_net)$name)

#Degree---------------------------------------------------------------------------------------------

out_degree = centralization.degree(grn_net, mode="out")

out_degVal = out_degree$res
out_degVal = as.data.frame(out_degVal)

#normalize
out_degVal_nor = out_degVal$out_degVal/(exp_size-1)
out_degVal_nor = round(out_degVal_nor, 4)

out_degVal_nor = cbind(nodes_name,out_degVal_nor)
write_tsv(out_degVal_nor,"./outputs/out_degree_50.tsv",append = FALSE)

#Betweenness-----------------------------------------------------------------
net_betweenness = centralization.betweenness(grn_net,directed = TRUE, nobigint = TRUE, normalized = TRUE)

net_betVal = net_betweenness$res
net_betVal = as.data.frame(net_betVal)

#normalize
net_betVal = net_betVal$net_betVal/((exp_size-1)*(exp_size-2)/2)
net_betVal = round(net_betVal, 4)

net_betVal = cbind(nodes_name,net_betVal)
write_tsv(net_betVal,"./outputs/net_betVal_50.tsv",append = FALSE)

#Clustering Coefficient through Wu method----------------------------------------------------------------

#Find triangles in GRN
tiangl_grn = count_triangles(grn_net, vids = V(grn_net))
tiangl_grn = as.data.frame(tiangl_grn)

clust_coeVal = (2*tiangl_grn)/(out_degVal*(out_degVal-1))

clust_coeVal = as.numeric(unlist(clust_coeVal))
clust_coeVal[which(!is.finite(clust_coeVal))] <- 0
clust_coeVal = round(clust_coeVal, 4)

clust_coeVal = cbind(nodes_name,clust_coeVal)
write_tsv(clust_coeVal,"./outputs/net_clusVal_50.tsv",append = FALSE)

#---------------------------------------------------------------------------------------------------





