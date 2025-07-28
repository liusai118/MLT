library(data.table)
library(preprocessCore)
library(cmapR)
library(Seurat)
#############################
matrix <- my_ds@mat
id2list <- fread("data/GSE92742_Broad_LINCS_inst_info.txt.gz")
a <- id2list[which( id2list$inst_id %in% colnames(matrix)),]
table(a$cell_id)
gene <- fread("data/META/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz")
data <- unique(allcell@meta.data$celltype)
mat <- as.data.frame(matrix)
mat$pr_gene_id <- rownames(mat)
gene <- as.data.frame(gene)
mat <- merge(gene,mat,by="pr_gene_id")
rownames(mat) <- mat[,2]
mat <- mat[,-c(1:5)]
##########
common <- intersect(rownames(mat), rownames(allcell@assays$RNA@counts))
a2 <- as.data.frame(table(allcell@meta.data$celltype))
a1 <- unique(allcell@meta.data$celltype)
##############
Idents(allcell) <- "celltype"
allcell@active.ident <- as.factor(allcell@active.ident)
for(i in 5:length(file_list)){
    my_ds = parse_gctx(paste0("data/output_gctx_files/",file_list[i]))
    mat <- as.data.frame(my_ds@mat)
    mat$pr_gene_id <- rownames(mat)
    mat <- merge(gene,mat,by="pr_gene_id")
    rownames(mat) <- mat[,2]
    mat <- mat[,-c(1:5)]
    for(j in 1:length(a1)){
      class(allcell)
      cell1 <- subset(allcell,celltype == paste0(a1[j]))
      scmatrix <- cell1@assays$RNA$data
      matrix_a1 <- cbind(mat[common,], scmatrix[common,])
      matrix_a1 <- as.matrix(matrix_a1)
      print("normalization")
      matrix_a1 <- normalize.quantiles(matrix_a1)
      L1000 <- matrix_a1[, 1:ncol(mat)]
      sc_cell <- matrix_a1[, -c(1:ncol(mat))]
      print("caculating")
      L1000_c <- calculate_centroid(L1000,10)######参数
      sc_cell_c <- calculate_centroid(sc_cell,10)
      cor_peason[j,i] <- cor(L1000_c[[2]], sc_cell_c[[2]])
      eu_dis[j,i] <- euclidean_distance(L1000_c[[2]], sc_cell_c[[2]])
      cos_dis[j,i]   <- cosine_similarity(L1000_c[[2]], sc_cell_c[[2]])
      #cov_matrix <- cov(t(matrix_a2))
      #X_ms <- mahalanobis_distance(L1000_c[[2]],sc_cell_c[[2]],cov_matrix,tol=1e-40)
      #cor_peason[j,i] <- matrix(0,length(file_list),length(a1))
      # eu_dis[j,i] <- matrix(0,length(file_list),length(a1))
      #ms_dis[j,i] <- matrix(0,length(file_list),length(a1))
      #cos_dis[j,i] <- matrix(0,length(file_list),length(a1))
      print(paste0(a1[j],file_list[i],"has been done",i,j))
    }
   rm(cell1)
   rm(scmatrix)
   rm(matrix_a1)
   rm(L1000 )
   rm(sc_cell)
   gc()
   write.csv(cos_dis,paste0("csv/",i,"+",j,"cos_dis.csv"))
   write.csv(eu_dis,paste0("csv/",i,"+",j,"eu_dis.csv"))
   write.csv(cor_peason,paste0("csv/",i,"+",j,"cor_peason.csv"))
}

