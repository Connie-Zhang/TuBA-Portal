load("up_cancer_list.RData")
load("down_cancer_list.RData")
load("up_data_list.RData")
load("down_data_list.RData")
load("survival_data_list.RData")
load("down_BP_list.RData")
load("up_BP_list.RData")
load("up_pval_list.RData")
load("down_pval_list.RData")
load("down_BP_full.RData")
load("up_BP_full.RData")
load("up_foo_list.RData")
load("down_foo_list.RData")

cancer_names <- names(up_data_list)
full_cancer_names <- c("Bladder urothelial carcinoma (BLCA)" = "BLCA", "Breast invasive carcinoma (BRCA)" = "BRCA", "Cervical squamous cell carcinoma and endocervical adenocarcinoma (CESC)" = "CESC", "Colon and Rectum adenocarcinoma (COADREAD)" = "COADREAD",
                       "Esophageal carcinoma (ESCA)" = "ESCA", "Glioblastoma multiforme (GBM)" = "GBM", "Glioblastoma multiforme and Brain lower grade glioma (GBMLGG)" = "GBMLGG", "Head and Neck squamous cell carcinoma (HNSC)" = "HNSC", 
                       "Kidney renal clear cell carcinoma (KIRC)" = "KIRC", "Kidney renal papillary cell carcinoma (KIRP)" = "KIRP", "Brain lower grade glioma (LGG)" = "LGG", "Liver hepatocellular carcinoma (LIHC)" = "LIHC", "Lung adenocarcinoma (LUAD)" = "LUAD",
                       "Lung squamous cell carcinoma (LUSC)" = "LUSC", "Ovarian serous cystadenocarcinoma (OV)" = "OV", "Pancreatic adenocarcinoma (PAAD)" = "PAAD", "Pheochromocytoma and Paraganglioma (PCPG)" = "PCPG", "Prostate adenocarcinoma (PRAD)" = "PRAD",
                       "Sarcoma (SARC)" = "SARC", "Skin cutaneous melanoma (SKCM)" = "SKCM", "Stomach adenocarcinoma (STAD)" = "STAD", "Thyroid carcinoma (THCA)" = "THCA")
chrom <-  c("chr19","Missing","chr11","chr16","chr17","chr1","chr2","chr12","chr7","chr10",  
            "chr14","chr9","chrX","chr8","chr22","chr20","chr5","chr15","chr3","chr6","chr4","chr13","chr18","chr21","all")
myColors <- c("#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20","#808080")
names(myColors) <- c("high_amplification","amplification","no_change","1_copy_del","2_copy_del","Missing")

names(survival_data_list) <- names(up_data_list)
variable_names <- c("OS"="Overall Survival","DSS" = "Disease-specific Survival","DFI"="Disease-free Interval","PFI"="Progression-free Interval")

full_data_list <- list("Up"=up_data_list,"Down"=down_data_list)
full_pval_list <- list("Up" = up_pval_list,"Down" = down_pval_list)
full_cancer_list <- list("Up" = up_cancer_list,"Down" = down_cancer_list)
full_BP_list <- list("Up"=up_BP_full,"Down"=down_BP_full)
full_foo_list <- list("Up"=up_foo_list,"Down"=down_foo_list)

myColors[["1_copy_del"]] <- "#F21A00"
myColors[["2_copy_del"]] <- "#972D15"

param_high <- list("BLCA"=c(0.05,0.3),"BRCA"=c(0.05,0.3),"CESC"=c("0.075,0.3"),
                   "COADREAD"=c(0.075,0.4),"ESCA"=c(0.1,0.3),"GBM"=c(0.1,0.5),
                   "GBMLGG"=c(0.05,0.4),"HNSC"=c(0.05,0.4),"KIRC"=c(0.05,0.4),
                   "KIRP"=c(0.075,0.4),"LGG"=c(0.05,0.5),"LIHC"=c(0.075,0.3),
                   "LUAD"=c(0.05,0.3),"LUSC"=c(0.05,0.3),"OV"=c(0.075,0.3),
                   "PAAD"=c(0.1,0.5),"PCPG"=c(0.1,0.4),"PRAD"=c(0.05,0.4),
                   "SARC"=c(0.075,0.3),"SKCM"=c(0.075,0.4),"STAD"=c(0.075,0.4),
                   "THCA"=c(0.05,0.4))
param_low <- list("BLCA"=c(0.05,0.3),"BRCA"=c(0.05,0.3),"CESC"=c("0.075,0.3"),
                  "COADREAD"=c(0.075,0.4),"ESCA"=c(0.1,0.3),"GBM"=c(0.1,0.5),
                  "GBMLGG"=c(0.05,0.4),"HNSC"=c(0.05,0.4),"KIRC"=c(0.05,0.4),
                  "KIRP"=c(0.075,0.4),"LGG"=c(0.05,0.5),"LIHC"=c(0.075,0.3),
                  "LUAD"=c(0.05,0.3),"LUSC"=c(0.05,0.3),"OV"=c(0.075,0.3),
                  "PAAD"=c(0.1,0.5),"PCPG"=c(0.1,0.4),"PRAD"=c(0.05,0.4),
                  "SARC"=c(0.075,0.3),"SKCM"=c(0.075,0.4),"STAD"=c(0.075,0.4),
                  "THCA"=c(0.05,0.4))
TuBA_params <- list("Up"=param_high,"Down"=param_low)
locus_dat <- read.delim("gene_with_protein_product.txt")
locus_dat <- locus_dat %>% select(symbol, location)
