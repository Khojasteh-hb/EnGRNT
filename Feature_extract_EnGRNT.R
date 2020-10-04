#Feature Extraction phase in EnGRNT 

################ load necessary libraries ############################################################
library(igraph)
library(readr)
library(ggplot2)
library(intergraph)

################ read input data #####################################################################

# read gene expression profiles in a variety of knockdown, knockout, and multifactorial 
expression_path = read_tsv("./Data/Ecoli-1_knockouts.tsv",col_names = TRUE)
expression_data = as.matrix(expression_path) 

# read known regulatory interactions between TF and target gene 
regulation_path = read_tsv("./Data/Ecoli-1_goldstandard.tsv",col_names = FALSE)
regulation_data = as.matrix(regulation_path)

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

exp_size = ncol(expression_data)
no_of_exp = nrow(expression_data)

#test_size= round(0.3* no_of_exp)
#expression_trans = t(expression_data)

expressed_regulation = subset(total_regulation,total_regulation[,3]>0)
TF_list = unique(expressed_regulation[,1])

expressed_reg = subset(regulation_data,regulation_data[,3]>0)

g = as.data.frame(expressed_reg)
net = graph.data.frame(g,directed = TRUE)

nodes_name <- data.frame(name = V(net)$name)

#adjc_net = as.matrix(net,matrix.type = "edgelist")

#Degree-----------------------------------------------------------------
in_degree = centralization.degree(net, mode="in")

in_degVal = in_degree$res
in_degVal = as.data.frame(in_degVal)
in_degVal = cbind(nodes_name,in_degVal)
#write_tsv(in_degVal,"D:/Feature-Extraction-GRN/in_degree_values/in_degree_100.tsv",append = FALSE)

out_degree = centralization.degree(net, mode="out")

out_degVal = out_degree$res
out_degVal = as.data.frame(out_degVal)
out_degVal = cbind(nodes_name,out_degVal)
#write_tsv(out_degVal,"D:/Feature-Extraction-GRN/out_degree_values/out_degree_100.tsv",append = FALSE)
#------------------------------------------------------------------

#Clustering Coefficient.----------------------------------------------------------------
clust_coeVal=matrix(0,1,ncol(expression_data))
clust_coe = transitivity(net, type= "local", weights=NULL, isolates= "zero")

clust_coeVal = clust_coe
clust_coeVal = as.data.frame(clust_coeVal)
clust_coeVal = cbind(nodes_name,clust_coeVal)
#write_tsv(clust_coeVal,"D:/Feature-Extraction-GRN/clustering-coef_values/clust_coeVal_100.tsv",append = FALSE)
#------------------------------------------------------------------

#Betweenness-----------------------------------------------------------------
net_betweenness = centralization.betweenness(net,directed = TRUE, nobigint = TRUE,
                                             normalized = TRUE)
net_betVal = net_betweenness$res
net_betVal = as.data.frame(net_betVal)
net_betVal = cbind(nodes_name,net_betVal)
#write_tsv(net_betVal,"D:/Feature-Extraction-GRN/between_values/net_betVal_100.tsv",append = FALSE)
#------------------------------------------------------------------

#Degree Distribution-----------------------------------------------------------------
p=degree.distribution(net)
#plot(p,type = 'o')
#------------------------------------------------------------------

#Plot Gene Regulatory Network-----------------------------------------------------------------
#tkplot(net)
#------------------------------------------------------------------


