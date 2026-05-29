#https://github.com/edavis71/scEVA
#function to calculate eva statistics
eva_stats <- function(df, dist, pheno) {
  phenotypes = as.factor(pheno)
  samplesC1 <- colnames(dist[,phenotypes==0])
  samplesC2 <- colnames(dist[,phenotypes==1])
  
  dist1 <- dist[samplesC1,samplesC1]
  n1 <- length(samplesC1)
  E1 <- sum(dist1)/(n1*(n1-1))
  Dis1P2 <- sum(dist1^2)
  E1P2 <- Dis1P2/(n1*(n1-1))
  E1DCross <- (sum(apply(dist1,1,sum)^2)-Dis1P2)/(n1*(n1-1)*(n1-2))
  VarD1 <- E1P2 - E1^2
  CovE1Cross <- E1DCross - E1^2
  VarEta1 <- 4*((n1*(n1-1))/2*VarD1 + n1*(n1-1)*(n1-2)*CovE1Cross)/((n1*(n1-2))^2)
  
  dist2 <- dist[samplesC2,samplesC2]
  n2 <- length(samplesC2)
  E2 <- sum(dist2)/(n2*(n2-1))
  Dis2P2 <- sum(dist2^2)
  E2P2 <- Dis2P2/(n2*(n2-1))
  E2DCross <- (sum(apply(dist2,1,sum)^2)-Dis2P2)/(n2*(n2-1)*(n2-2))
  VarD2 <- E2P2 - E2^2
  CovE2Cross <- E2DCross - E2^2
  VarEta2 <- 4*((n2*(n2-1))/2*VarD2 + n2*(n2-1)*(n2-2)*CovE2Cross)/((n2*(n2-2))^2)
  
  vartotal <- VarEta2 + VarEta1;
  zscore <- (E1-E2)/sqrt(abs(vartotal+0.0000001));
  
  pvalue=2*(1-pnorm(abs(zscore)))
  
  nAll = n1 + n2
  EAll <- sum(dist)/(nAll*(nAll-1))
  DisAllP2 <- sum(dist^2)
  EAllP2 <- DisAllP2/(nAll*(nAll-1))
  EAllDCross <- (sum(apply(dist,1,sum)^2)-DisAllP2)/(nAll*(nAll-1)*(nAll-2))
  CovEAllCross <- EAllDCross - EAll^2
  lambda = n1/(n1+n2)
  
  myvartotal = 4*(1/lambda+1/(1-lambda))*CovEAllCross/nAll
  
  return(list(E1=E1,E2=E2,
              zscore=zscore,     
              VarEta1=VarEta1,VarEta2=VarEta2, #EVA 1
              sdtotal=sqrt(abs(vartotal)),#EVA 1 sd 
              pvalue=pvalue))
 
}

############################################################################################
#Read in Regulons and Load genes in each regulon
suppressPackageStartupMessages({
   library(SCopeLoomR)
   library(SCENIC)
})    
wd="D:/16p_11/SCIENIC_result"
setwd(wd)

input_file = "Deletions" #Controls Deletions
loom <- open_loom(paste(input_file,"auc_mtx.loom",sep="."), mode="r") 

#Regulons:
regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons') # as incid matrix
regulons_motif <- regulonsToGeneLists(regulons_incidMat) # convert to list
lengths(regulons_motif)

regulons_subset <- regulons_motif[which(lengths(regulons_motif) > 10)]
regulons_subset <- regulons_subset[which(lengths(regulons_subset) < 200)]
length(regulons_subset)

load("D:/16p_11/regulon_filtering_25th_Jan.RData") 
#Ctrl_regulons_union, Ctrl_regulons_intersect, Del_regulons_union, Del_regulons_intersect
regulons_subset <- regulons_subset[Del_regulons_intersect]

############################################################################################
library(GSBenchMark)
library(ggplot2)
library(reshape2)
library(GSReg)
library(philentropy)
library(dplyr) 
library(Seurat) 
library(monocle3) 
#source("multiplot.R")

file_dir <- "D:/16p_11"
setwd(file_dir)

iPSC_object <- readRDS("./filtered_object_19th_Jan.rds")
DefaultAssay(iPSC_object) <- "RNA"
iPSC_object@meta.data$Sample <- unlist(lapply(rownames(iPSC_object@meta.data),
                                              function(x){strsplit(x, split = "_")[[1]][2]}))
iPSC_object@meta.data$Replicate <- paste("iPSC",iPSC_object@meta.data$Sample,sep='_')
iPSC_object@meta.data$Dataset <- "iPSC"
iPSC_object@meta.data$cellType <- iPSC_object@meta.data$broad_Type
iPSC_object@meta.data$clustered_prev <- iPSC_object@meta.data$seurat_clusters
iPSC_object <- SetIdent(object = iPSC_object, value = "seurat_clusters")

iPSC_object
iPSC_object@meta.data$cellType_tmp <- paste(iPSC_object@meta.data$broad_Type, iPSC_object@meta.data$Replicate, iPSC_object@meta.data$Conditions, sep="_")
table(iPSC_object@meta.data$cellType_tmp)



##############################################################################
#Running begin here
##############################################################################
#Real dataset
cds <- SeuratWrappers::as.cell_data_set(iPSC_object)
levels(as.factor(pData(cds)$Replicate)) #"C1" "C2" "C3" "D1" "D2" "D3"
pData(cds)$cellType_tmp <- paste(pData(cds)$broad_Type, pData(cds)$Conditions, sep="_")
levels(as.factor(pData(cds)$cellType_tmp)) 

#data(diracpathways)
hallmark_pathways = regulons_subset #diracpathways
#EVA anlaysis for all celltypes and samples/patients/conditions
#https://github.com/edavis71/scEVA/blob/3a4b7a625021dc876763880b81db7b335487cdee/breastcancerBC9_11.R
do_pathways2 <- function(pathway, sample1, sample2, celltype) {
  # subset cds by Replicate ID
  test<- cds[,(pData(cds)$Replicate==sample1 | pData(cds)$Replicate==sample2)]
  # subset cancer cells
  testc <- test[,(pData(test)$cellType_tmp==celltype)]
  pData(testc)$compare <- ifelse(pData(testc)$Replicate == sample1, 0,1)
  #celltype1 =0; celltype2 =1 in pData(testc)$compare
  # get expression matrix for pt cancer cells
  exprsdata <- as.matrix(exprs(testc))
  #rownames(exprsdata) <- fData(testc)$gene_short_name
  #classic eva
  VarAnKendallV = GSReg::GSReg.GeneSets.EVA(geneexpres=exprsdata, pathways=pathway, phenotypes=as.factor(pData(testc)$compare)) 
  
  pvalustat = sapply(VarAnKendallV,function(x) x$pvalue);
  kendall <- lapply(VarAnKendallV, function (x) x[1:20])
  dd  <-  as.data.frame(matrix(unlist(kendall), nrow=length(unlist(kendall[1]))))
  rownames(dd) <- make.unique(names(kendall[[1]]))
  colnames(dd) <- names(kendall)
  dd <- as.data.frame(t(dd))
  dd$pathway_name <- rownames(dd)
  
  dd$metric <- "kendall_tau"
  dd$patient1 <- sample1
  dd$patient2 <- sample2
  dd$celltype <- celltype
  return(dd)
}

#Using Replicate instead of sample; Using cellType_tmp instead of celltype
a <- do_pathways2(pathway = hallmark_pathways, sample1 = "iPSC_C1", sample2 = "iPSC_C2", celltype = "Progenitor_Controls")
b <- do_pathways2(pathway = hallmark_pathways, sample1 = "iPSC_C2", sample2 = "iPSC_C3", celltype = "Progenitor_Controls")
c <- do_pathways2(pathway = hallmark_pathways, sample1 = "iPSC_C3", sample2 = "iPSC_C1", celltype = "Progenitor_Controls")

d <- do_pathways2(pathway = hallmark_pathways, sample1 = "iPSC_C1", sample2 = "iPSC_C2", celltype = "Precursor_Controls")
e <- do_pathways2(pathway = hallmark_pathways, sample1 = "iPSC_C2", sample2 = "iPSC_C3", celltype = "Precursor_Controls")
f <- do_pathways2(pathway = hallmark_pathways, sample1 = "iPSC_C3", sample2 = "iPSC_C1", celltype = "Precursor_Controls")

g <- do_pathways2(pathway = hallmark_pathways, sample1 = "iPSC_D1", sample2 = "iPSC_D2", celltype = "Progenitor_Deletions")
h <- do_pathways2(pathway = hallmark_pathways, sample1 = "iPSC_D2", sample2 = "iPSC_D3", celltype = "Progenitor_Deletions")
i <- do_pathways2(pathway = hallmark_pathways, sample1 = "iPSC_D3", sample2 = "iPSC_D1", celltype = "Progenitor_Deletions")

j <- do_pathways2(pathway = hallmark_pathways, sample1 = "iPSC_D1", sample2 = "iPSC_D2", celltype = "Precursor_Deletions")
k <- do_pathways2(pathway = hallmark_pathways, sample1 = "iPSC_D2", sample2 = "iPSC_D3", celltype = "Precursor_Deletions")
l <- do_pathways2(pathway = hallmark_pathways, sample1 = "iPSC_D3", sample2 = "iPSC_D1", celltype = "Precursor_Deletions")

all <- rbind(a,b,c,d,e,f,g,h,i,j,k,l)
all$pvaluadj <- p.adjust(all$pvalue,method='BH')

common_cols <- c("estat", "pathway", "patient", "UniqueID")
all$UniqueID <- paste(all$patient1, all$celltype, sep="_")
table(all$UniqueID)
E1 <- all[c("E1","pathway_name","patient1","UniqueID")] #all[c(1,21,23,27)]
colnames(E1) <- common_cols

mat <- acast(E1, pathway ~ UniqueID, value.var='estat') #col  celltype fun.aggregate = sum
mats <- t(apply(mat,1,scale)) #the rows will be the name of pathway/geneset 
colnames(mats) <- colnames(mat)

#############
#Draw heatmap
library(ComplexHeatmap)
#############
#add heatmap annotation
df = data.frame(patient = sapply(colnames(mats), function(x){strsplit(x, split = "_")[[1]][2]}), #12 columns
                celltype = sapply(colnames(mats), function(x){strsplit(x, split = "_")[[1]][3]}))
#Green as highlight
greenfocus = c("#41AB5D", "#252525", "#525252", "#737373", "#969696", "#BDBDBD", "#D9D9D9", "#F0F0F0")
top = HeatmapAnnotation(df = df, col = list(patient = c("C1" = "#619CFF", "C2" = "#F8766D", "C3" = "#00BA38",
                                                        "D1" = "#8B7355", "D2" = "#CABA83", "D3" = "#EEE5AD"), 
                                            celltype = c("Precursor" = "#252525", 
                                                         "Progenitor" = "#525252")))

Heatmap(mats, name = "EVA Statistic", c("steelblue3", "khaki1", "red1"), top_annotation = top, row_names_side = "left", 
        row_dend_side = "right", width = unit(80, "mm"), row_names_gp = gpar(fontsize = 8))

pdf(file = paste("heatmap_Del_intersect_", Sys.Date(), "EVA_Statistic.pdf", sep = ""), width = 8, height = 12, useDingbats = F)
Heatmap(mats, name = "EVA Statistic", c("steelblue3", "khaki1", "red1"), top_annotation = top, row_names_side = "left", 
        row_dend_side = "right", width = unit(80, "mm"), row_names_gp = gpar(fontsize = 8))
dev.off(which = dev.cur())

save(all,mats, file='EVA_Ctrl_Union_Statistic.Rdata')

#####################################################################
setwd("~/Downloads")
load("EVA_Ctrl_Union_Statistic.Rdata")
mats <- as.data.frame(mats)
dataframe <- data.frame(matrix(nrow = nrow(mats), ncol = 18))
colnames(dataframe) <- c("D_N_Prog","p_Prog","adj_p_Prog","Prog_D2","Prog_D1","Prog_D3","Prog_C2","Prog_C3","Prog_C1",
                         "D_N_Prec","p_Prec","adj_p_Prec","Prec_D2","Prec_D1","Prec_D3","Prec_C1","Prec_C2","Prec_C3")
rownames(dataframe) <- rownames(mats)

D_N_Prog <- read.csv("Controls_Progenitor_regulon_summary_25th_Jan.csv", header=T, row.names=1)
D_N_Prog <- D_N_Prog[rownames(dataframe),]
D_N_Prec <- read.csv("Controls_Precursor_regulon_summary_25th_Jan.csv", header=T, row.names=1)
D_N_Prec <- D_N_Prec[rownames(dataframe),]
	
dataframe$D_N_Prog <- D_N_Prog$Progenitor_D_N
dataframe$Prog_D2 <- mats$iPSC_D2_Progenitor_Deletions
dataframe$Prog_D1 <- mats$iPSC_D1_Progenitor_Deletions
dataframe$Prog_D3 <- mats$iPSC_D3_Progenitor_Deletions
dataframe$Prog_C2 <- mats$iPSC_C2_Progenitor_Controls
dataframe$Prog_C3 <- mats$iPSC_C3_Progenitor_Controls
dataframe$Prog_C1 <- mats$iPSC_C1_Progenitor_Controls

for(i in 1:nrow(dataframe)){
data_a = unname(unlist(dataframe[i, 4:6]))
data_b = unname(unlist(dataframe[i, 7:9]))
dataframe[i,2] <- t.test(data_a, data_b)$p.value
}
dataframe$adj_p_Prog <- p.adjust(dataframe$p_Prog, method = "BH")


dataframe$D_N_Prec <- D_N_Prec$Precursor_D_N
dataframe$Prec_D2 <- mats$iPSC_D2_Precursor_Deletions
dataframe$Prec_D1 <- mats$iPSC_D1_Precursor_Deletions
dataframe$Prec_D3 <- mats$iPSC_D3_Precursor_Deletions
dataframe$Prec_C1 <- mats$iPSC_C1_Precursor_Controls
dataframe$Prec_C2 <- mats$iPSC_C2_Precursor_Controls
dataframe$Prec_C3 <- mats$iPSC_C3_Precursor_Controls
		
for(i in 1:nrow(dataframe)){
data_a = unname(unlist(dataframe[i, 13:15]))
data_b = unname(unlist(dataframe[i, 16:18]))
dataframe[i,11] <- t.test(data_a, data_b)$p.value
}
dataframe$adj_p_Prec <- p.adjust(dataframe$p_Prec, method = "BH")

write.csv(dataframe, "./EVA_Ctrl_union_Statistic_08th_Feb.csv")

##########################################################################################################################
setwd("~/Downloads")
load("EVA_Del_union_Statistic.Rdata")
mats <- as.data.frame(mats)
dataframe <- data.frame(matrix(nrow = nrow(mats), ncol = 18))
colnames(dataframe) <- c("D_N_Prog","p_Prog","adj_p_Prog","Prog_D2","Prog_D1","Prog_D3","Prog_C2","Prog_C3","Prog_C1",
                         "D_N_Prec","p_Prec","adj_p_Prec","Prec_D2","Prec_D1","Prec_D3","Prec_C1","Prec_C2","Prec_C3")
rownames(dataframe) <- rownames(mats)

D_N_Prog <- read.csv("Deletions_Progenitor_regulon_summary_25th_Jan.csv", header=T, row.names=1)
D_N_Prog <- D_N_Prog[rownames(dataframe),]
D_N_Prec <- read.csv("Deletions_Precursor_regulon_summary_25th_Jan.csv", header=T, row.names=1)
D_N_Prec <- D_N_Prec[rownames(dataframe),]
	
dataframe$D_N_Prog <- D_N_Prog$Progenitor_D_N
dataframe$Prog_D2 <- mats$iPSC_D2_Progenitor_Deletions
dataframe$Prog_D1 <- mats$iPSC_D1_Progenitor_Deletions
dataframe$Prog_D3 <- mats$iPSC_D3_Progenitor_Deletions
dataframe$Prog_C2 <- mats$iPSC_C2_Progenitor_Controls
dataframe$Prog_C3 <- mats$iPSC_C3_Progenitor_Controls
dataframe$Prog_C1 <- mats$iPSC_C1_Progenitor_Controls

for(i in 1:nrow(dataframe)){
data_a = unname(unlist(dataframe[i, 4:6]))
data_b = unname(unlist(dataframe[i, 7:9]))
dataframe[i,2] <- t.test(data_a, data_b)$p.value
}
dataframe$adj_p_Prog <- p.adjust(dataframe$p_Prog, method = "BH")


dataframe$D_N_Prec <- D_N_Prec$Precursor_D_N
dataframe$Prec_D2 <- mats$iPSC_D2_Precursor_Deletions
dataframe$Prec_D1 <- mats$iPSC_D1_Precursor_Deletions
dataframe$Prec_D3 <- mats$iPSC_D3_Precursor_Deletions
dataframe$Prec_C1 <- mats$iPSC_C1_Precursor_Controls
dataframe$Prec_C2 <- mats$iPSC_C2_Precursor_Controls
dataframe$Prec_C3 <- mats$iPSC_C3_Precursor_Controls
		
for(i in 1:nrow(dataframe)){
data_a = unname(unlist(dataframe[i, 13:15]))
data_b = unname(unlist(dataframe[i, 16:18]))
dataframe[i,11] <- t.test(data_a, data_b)$p.value
}
dataframe$adj_p_Prec <- p.adjust(dataframe$p_Prec, method = "BH")

write.csv(dataframe, "EVA_Del_union_Statistic_08th_Feb.csv")
