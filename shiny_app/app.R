library(shiny)
library(shinyWidgets)
library(ggplot2)
library(dplyr)
library(patchwork)
library(rsconnect)
library(base)
library(lubridate)
library(shinyjs)
library(shinydashboard)
library(plotly)
library(htmlwidgets)
library(rlang)
library(forcats)
library(tidyr)
library(rlang)
library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(survminer)
library(rlist)
library(DT)
library(tidyverse)
library(stringr)

cancer_names <- names(up_data_list)
chrom <-  c("chr19","Missing","chr11","chr16","chr17","chr1","chr2","chr12","chr7","chr10",  
"chr14","chr9","chrX","chr8","chr22","chr20","chr5","chr15","chr3","chr6","chr4","chr13","chr18","chr21","all")
myColors <- c("#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20","#808080")
names(myColors) <- c("high_amplification","amplification","no_change","1_copy_del","2_copy_del","Missing")

load("up_cancer_list.RData")
load("down_cancer_list.RData")
load("up_data_list.RData")
load("down_data_list.RData")
load("survival_data_list.RData")
load("up_foo_list.RData")
load("down_foo_list.RData")
load("down_BP_list.RData")
load("up_BP_list.RData")

names(survival_data_list) <- names(up_data_list)

BLCA_data_H <- up_data_list["BLCA"][[1]]
BLCA_data_filtered_H <- BLCA_data_H %>% filter(Samples.In.Bicluster > 20)
BLCA_data_L <- down_data_list["BLCA"][[1]]
BLCA_data_filtered_L <- BLCA_data_L %>% filter(Samples.In.Bicluster > 20)

BLCA_genes_H <- BLCA_data_filtered_H$Gene.ID
BLCA_genes_L <- BLCA_data_filtered_L$Gene.ID

# code for BLCA Go Terms
# BLCA_bp_H <-
# BLCA_bp_L <- 

BLCA_gene_list_H <- list()
for (i in unique(BLCA_data_filtered_H$Bicluster.No)){
    new <- BLCA_data_filtered_H %>% filter(Bicluster.No==i)
    genes <- new$Gene.ID
    BLCA_gene_list_H <- list.append(BLCA_gene_list_H,genes)
}
names(BLCA_gene_list_H) <- unique(BLCA_data_filtered_H$Bicluster.No)

BLCA_gene_list_L <- list()
for (i in unique(BLCA_data_filtered_L$Bicluster.No)){
    new <- BLCA_data_filtered_L %>% filter(Bicluster.No==i)
    genes <- new$Gene.ID
    BLCA_gene_list_L <- list.append(BLCA_gene_list_L,genes)
}
names(BLCA_gene_list_L) <- unique(BLCA_data_filtered_L$Bicluster.No)

BLCA_BP_list_H <- list()
for (i in unique(BLCA_data_filtered_H$Bicluster.No)){
    new <- up_BP_list["BLCA"][[1]] %>% filter(Bicluster.No==i)
    BPs <- c(new$GoTerm1, new$GoTerm2, new$GoTerm3, new$GoTerm4, new$GoTerm5)
    BLCA_BP_list_H <- list.append(BLCA_BP_list_H,BPs)
}
names(BLCA_BP_list_H) <- unique(BLCA_data_filtered_H$Bicluster.No)

BLCA_BP_list_L <- list()
for (i in unique(BLCA_data_filtered_L$Bicluster.No)){
    new <- down_BP_list["BLCA"][[1]] %>% filter(Bicluster.No==i)
    BPs <- c(new$GoTerm1, new$GoTerm2, new$GoTerm3, new$GoTerm4, new$GoTerm5)
    BLCA_BP_list_L <- list.append(BLCA_BP_list_L,BPs)
}
names(BLCA_BP_list_L) <- unique(BLCA_data_filtered_L$Bicluster.No)

bladder_hBP <- unique(c(up_BP_list["BLCA"][[1]]$GoTerm1,up_BP_list["BLCA"][[1]]$GoTerm2,up_BP_list["BLCA"][[1]]$GoTerm3,up_BP_list["BLCA"][[1]]$GoTerm4,up_BP_list["BLCA"][[1]]$GoTerm5))
bladder_hBP <- c(bladder_hBP, "all")
bladder_lBP <- unique(c(down_BP_list["BLCA"][[1]]$GoTerm1,down_BP_list["BLCA"][[1]]$GoTerm2,down_BP_list["BLCA"][[1]]$GoTerm3,down_BP_list["BLCA"][[1]]$GoTerm4,down_BP_list["BLCA"][[1]]$GoTerm5))
bladder_lBP <- c(bladder_lBP, "all")

# Define UI for application that draws a histogram
ui <- dashboardPage (
    
    dashboardHeader(title="Bicluster Visualizations"),
    
    dashboardSidebar(
        sidebarMenu(
            menuItem("Introduction", tabName = "intro"),
            menuItem("Biological Pathways", tabName = "biopath"),
            menuItem("Copy Number", tabName="copynum"),
            menuItem("Contact", tabName = "contact"))),
    
    dashboardBody(
        tabItems(
            tabItem("intro",h4("to be added")),
            tabItem("biopath",
                    useShinyjs(),
                    fluidRow(
                        column(3,selectInput(
                            inputId = "type",
                            label = "Choose cancer of interest:",
                            selected = "BLCA",
                            choices = cancer_names
                        )),
                        column(2, selectInput(
                            inputId = "reg",
                            label = "Choose up or down regulated gene expression:",
                            selected="Up",
                            choices = c("Up","Down"))),
                        column(3,selectInput(
                            inputId = "variable",
                            label = "Select variable of interest",
                            choices = c("Overall survival", "Disease-specific survival", "Disease-free interval", "Progression-free interval"))),
                        # be able to input multiple genes
                        column(2, selectizeInput(
                            inputId = "gene",
                            label = "Gene of interest",
                            choices = BLCA_genes_H,
                            multiple = TRUE)),
                        column(3, selectizeInput(
                            inputId = "path",
                            label = "Biological Pathway of interest",
                            choices = bladder_hBP,
                            multiple=TRUE)),
                        # change bic input choices to dynamic
                        column(2,selectInput(
                            inputId = "bic",
                            label = "Select the bicluster of interest:",
                            choices = unique(BLCA_data_filtered_H$Bicluster.No)))),
                    sidebarPanel(
                        textOutput(outputId="bicinfo"),width = 3),
                    mainPanel(
                        tabsetPanel(
                            tabPanel("Visualization",plotOutput(outputId = "survivalvis"),width = 9),
                            tabPanel("Gene Information", dataTableOutput(outputId = "survivaltable"),style = "height:500px; overflow-y: scroll;overflow-x: scroll;",width = 9)))),
            tabItem("copynum",
                    useShinyjs(),
                    fluidRow(
                        column(2, selectInput(
                            inputId = "type_copy",
                            label = "Choose cancer of interest:",
                            selected = "BLCA",
                            choices = cancer_names)),
                        column(3, selectInput(
                            inputId = "reg_copy",
                            label = "Choose up or down regulated gene expression:",
                            selected="Up",
                            choices = c("Up","Down"))),
                        column(3,selectInput(
                            inputId = "num",
                            label = "Select the bicluster of interest:",
                            choices = unique(BLCA_data_filtered_H$Bicluster.No))),
                        column(3,selectInput(
                            inputId = "withcopy",
                            label = "Fill with copy number?",
                            choices = c("Yes","No"))),
                        column(3,selectizeInput(
                            inputId = "chrom",
                            label = "Select the chromosome of interest:",
                            choices = chrom, multiple=TRUE
                        ))),
                    sidebarPanel(
                        textOutput(outputId= "info"), width = 3),
                    mainPanel(
                        tabsetPanel(
                            tabPanel("Visualization",plotlyOutput(outputId = "mapvis"),width = 9),
                            tabPanel("Gene Information",dataTableOutput(outputId = "table"),style = "height:500px; overflow-y: scroll;overflow-x: scroll;",width = 9)))),
            tabItem("contact",h4("to be added")))))

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    # for biological pathway page
    observeEvent(input$reg, 
                 if (input$reg=="Up"){
                     data_list <- up_data_list
                     observeEvent(input$type, 
                                  if (!is.null(input$type)){
                                    dat_filtered <- data_list[input$type][[1]] %>% filter(Samples.In.Bicluster > 20)
                                    updateSelectInput(session,"bic", choices = dat_filtered$Bicluster.No)
                                    updateSelectInput(session,"gene", choices = c(unique(dat_filtered$Gene.ID),"all"))
                                    gene_list <- list()
                                    for (i in unique(dat_filtered$Bicluster.No)){
                                        new <- dat_filtered %>% filter(Bicluster.No==i)
                                        genes <- new$Gene.ID
                                        gene_list <- list.append(gene_list,genes)
                                    }
                                    names(gene_list) <- unique(dat_filtered$Bicluster.No)
                                    observeEvent(input$gene,
                                                 if (input$gene=="all") {
                                                     updateSelectInput(session,"bic",choices = unique(dat_filtered$Bicluster.No))
                                                 }
                                                 else {
                                                     updated_bic<-c()
                                                     for (i in 1:length(gene_list)){
                                                         if (input$gene %in% gene_list[[i]]) {
                                                             updated_bic <- append(updated_bic, names(gene_list[i]))
                                                         }
                                                     }
                                                     updateSelectInput(session,"bic", choices = updated_bic)})
                                    
                                    bp_list <- list()
                                    for (i in unique(dat_filtered$Bicluster.No)){
                                        new <- up_BP_list[input$type][[1]] %>% filter(Bicluster.No==i)
                                        BPs <- c(new$GoTerm1, new$GoTerm2, new$GoTerm3, new$GoTerm4, new$GoTerm5)
                                        bp_list <- list.append(bp_list,BPs)
                                    }
                                    names(bp_list) <- unique(dat_filtered$Bicluster.No)
                                    
                                    bp <- unique(c(up_BP_list[input$type][[1]]$GoTerm1,up_BP_list[input$type][[1]]$GoTerm2,up_BP_list[input$type][[1]]$GoTerm3,up_BP_list[input$type][[1]]$GoTerm4,up_BP_list[input$type][[1]]$GoTerm5))
                                    updateSelectInput(session,"path",choices = c(bp,"all"))
                                    
                                    observeEvent(input$path,
                                                 if (input$path=="all"){
                                                     updateSelectInput(session,"bic",choices = unique(dat_filtered$Bicluster.No))
                                                 }
                                                 else {
                                                     updated_bic<-c()
                                                     for (i in 1:length(bp_list)){
                                                         if (input$path %in% bp_list[[i]]) {
                                                             updated_bic <- append(updated_bic, names(bp_list[i]))
                                                         }
                                                     }
                                                     updateSelectInput(session,"bic",choices = updated_bic)
                                                 })
                                    
                                    observeEvent(input$variable,
                                                 if (input$variable=="Progression-free interval"){
                                                         output$survivalvis <- renderPlot({
                                                             samples <- colnames(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]][2:length(colnames(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]]))])
                                                             samples <- str_replace_all(samples,"\\.","-")
                                                             survival <- survival_data_list[input$type][[1]]
                                                             biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                                                             km_fit <- survfit(Surv(PFI.time,PFI)~bicluster, data=biclist)
                                                             autoplot(km_fit) + 
                                                                 labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
                                                                      title = " Progression-free Interval of \n Cancer Patients \n") + 
                                                                 theme(plot.title = element_text(hjust = 0.5), 
                                                                       axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                                                                       axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                                                                       legend.title = element_text(face="bold", size = 10))})
                                                         output$bicinfo <- renderText({
                                                             samples <- colnames(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]][2:length(colnames(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]]))])
                                                             samples <- str_replace_all(samples,"\\.","-")
                                                             survival <- survival_data_list[input$type][[1]]
                                                             biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                                                             km_fit <- survfit(Surv(PFI.time,PFI)~bicluster, data=biclist)
                                                             print(paste0("This bicluster contains: \n ",toString(nrow(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[1]]))," genes, \n ",toString(unique(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[1]]$Samples.In.Bicluster))," samples. \n","The KM-analysis outputs a p-value of ",round(surv_pvalue(km_fit,data=biclist)$pval,3)))
                                                         })}
                                                 else if (input$variable=="Overall survival"){
                                                     output$survivalvis <- renderPlot({
                                                         samples <- colnames(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]][2:length(colnames(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]]))])
                                                         samples <- str_replace_all(samples,"\\.","-")
                                                         survival <- survival_data_list[input$type][[1]]
                                                         biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                                                         km_fit <- survfit(Surv(OS.time,OS)~bicluster, data=biclist)
                                                         autoplot(km_fit) + 
                                                             labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
                                                                  title = " Overall Survival of \n Cancer Patients \n") + 
                                                             theme(plot.title = element_text(hjust = 0.5), 
                                                                   axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                                                                   axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                                                                   legend.title = element_text(face="bold", size = 10))})
                                                     output$bicinfo <- renderText({
                                                         samples <- colnames(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]][2:length(colnames(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]]))])
                                                         samples <- str_replace_all(samples,"\\.","-")
                                                         survival <- survival_data_list[input$type][[1]]
                                                         biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                                                         km_fit <- survfit(Surv(OS.time,OS)~bicluster, data=biclist)
                                                         print(paste0("This bicluster contains: \n ",toString(nrow(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[1]]))," genes, \n ",toString(unique(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[1]]$Samples.In.Bicluster))," samples. \n","The KM-analysis outputs a p-value of ",round(surv_pvalue(km_fit,data=biclist)$pval,3)))
                                                     })
                                                 }
                                                 else if (input$variable=="Disease-free interval"){
                                                     output$survivalvis <- renderPlot({
                                                         samples <- colnames(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]][2:length(colnames(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]]))])
                                                         samples <- str_replace_all(samples,"\\.","-")
                                                         survival <- survival_data_list[input$type][[1]]
                                                         biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                                                         km_fit <- survfit(Surv(DFI.time,DFI)~bicluster, data=biclist)
                                                         autoplot(km_fit) + 
                                                             labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
                                                                  title = " Disease-free Interval of \n Cancer Patients \n") + 
                                                             theme(plot.title = element_text(hjust = 0.5), 
                                                                   axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                                                                   axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                                                                   legend.title = element_text(face="bold", size = 10))})
                                                     output$bicinfo <- renderText({
                                                         samples <- colnames(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]][2:length(colnames(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]]))])
                                                         samples <- str_replace_all(samples,"\\.","-")
                                                         survival <- survival_data_list[input$type][[1]]
                                                         biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                                                         km_fit <- survfit(Surv(DFI.time,DFI)~bicluster, data=biclist)
                                                         print(paste0("This bicluster contains: \n ",toString(nrow(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[1]]))," genes, \n ",toString(unique(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[1]]$Samples.In.Bicluster))," samples. \n","The KM-analysis outputs a p-value of ",round(surv_pvalue(km_fit,data=biclist)$pval,3)))
                                                     })
                                                 }
                                                 else if (input$variable=="Disease-specific survival"){
                                                     output$survivalvis <- renderPlot({
                                                         samples <- colnames(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]][2:length(colnames(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]]))])
                                                         samples <- str_replace_all(samples,"\\.","-")
                                                         survival <- survival_data_list[input$type][[1]]
                                                         biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                                                         km_fit <- survfit(Surv(DSS.time,DSS)~bicluster, data=biclist)
                                                         autoplot(km_fit) + 
                                                             labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
                                                                  title = " Disease-specific Survival of \n Cancer Patients \n") + 
                                                             theme(plot.title = element_text(hjust = 0.5), 
                                                                   axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                                                                   axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                                                                   legend.title = element_text(face="bold", size = 10))})
                                                     output$bicinfo <- renderText({
                                                         samples <- colnames(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]][2:length(colnames(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]]))])
                                                         samples <- str_replace_all(samples,"\\.","-")
                                                         survival <- survival_data_list[input$type][[1]]
                                                         biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                                                         km_fit <- survfit(Surv(DSS.time,DSS)~bicluster, data=biclist)
                                                         print(paste0("This bicluster contains: \n ",toString(nrow(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[1]]))," genes, \n ",toString(unique(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[1]]$Samples.In.Bicluster))," samples. \n","The KM-analysis outputs a p-value of ",round(surv_pvalue(km_fit,data=biclist)$pval,3)))
                                                     })})
                                    output$survivaltable <- renderDataTable({
                                        datatable(up_cancer_list[input$type][[1]][[as.integer(input$bic)]][[1]]%>% select(Gene.ID,chrom),options = list(paging=FALSE))
                                    })
                                    # updateSelectInput(session,"path", choices = )
                                    })}
                 else {
                     data_list <- down_data_list
                     observeEvent(input$type, 
                                  if (!is.null(input$type)){
                                      dat_filtered <- data_list[input$type][[1]] %>% filter(Samples.In.Bicluster > 20)
                                      updateSelectInput(session,"bic", choices = dat_filtered$Bicluster.No)
                                      updateSelectInput(session,"gene", choices = c(unique(dat_filtered$Gene.ID),"all"))
                                      gene_list <- list()
                                      for (i in unique(dat_filtered$Bicluster.No)){
                                          new <- dat_filtered %>% filter(Bicluster.No==i)
                                          genes <- new$Gene.ID
                                          gene_list <- list.append(gene_list,genes)
                                      }
                                      names(gene_list) <- unique(dat_filtered$Bicluster.No)
                                      observeEvent(input$gene,
                                                   if (input$gene=="all") {
                                                       updateSelectInput(session,"bic",choices = unique(dat_filtered$Bicluster.No))
                                                   }
                                                   else {
                                                       updated_bic<-c()
                                                       for (i in 1:length(gene_list)){
                                                           if (input$gene %in% gene_list[[i]]) {
                                                               updated_bic <- append(updated_bic, names(gene_list[i]))
                                                           }
                                                       }
                                                       updateSelectInput(session,"bic", choices = updated_bic)})
                                      
                                      bp_list <- list()
                                      for (i in unique(dat_filtered$Bicluster.No)){
                                          new <- down_BP_list[input$type][[1]] %>% filter(Bicluster.No==i)
                                          BPs <- c(new$GoTerm1, new$GoTerm2, new$GoTerm3, new$GoTerm4, new$GoTerm5)
                                          bp_list <- list.append(bp_list,BPs)
                                      }
                                      names(bp_list) <- unique(dat_filtered$Bicluster.No)
                                      
                                      bp <- unique(c(down_BP_list[input$type][[1]]$GoTerm1,down_BP_list[input$type][[1]]$GoTerm2,down_BP_list[input$type][[1]]$GoTerm3,down_BP_list[input$type][[1]]$GoTerm4,down_BP_list[input$type][[1]]$GoTerm5))
                                      updateSelectInput(session,"path",choices = c(bp,"all"))
                                      
                                      observeEvent(input$path,
                                                   if (input$path=="all"){
                                                       updateSelectInput(session,"bic",choices = unique(dat_filtered$Bicluster.No))
                                                   }
                                                   else {
                                                       updated_bic<-c()
                                                       for (i in 1:length(bp_list)){
                                                           if (input$path %in% bp_list[[i]]) {
                                                               updated_bic <- append(updated_bic, names(bp_list[i]))
                                                           }
                                                       }
                                                       updateSelectInput(session,"bic",choices = updated_bic)
                                                   })
                                      
                                      observeEvent(input$variable,
                                                   if (input$variable=="Progression-free interval"){
                                                       output$survivalvis <- renderPlot({
                                                           samples <- colnames(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]][2:length(colnames(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]]))])
                                                           samples <- str_replace_all(samples,"\\.","-")
                                                           survival <- survival_data_list[input$type][[1]]
                                                           biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                                                           km_fit <- survfit(Surv(PFI.time,PFI)~bicluster, data=biclist)
                                                           autoplot(km_fit) + 
                                                               labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
                                                                    title = " Progression-free Interval of \n Cancer Patients \n") + 
                                                               theme(plot.title = element_text(hjust = 0.5), 
                                                                     axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                                                                     axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                                                                     legend.title = element_text(face="bold", size = 10))})
                                                       output$bicinfo <- renderText({
                                                           samples <- colnames(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]][2:length(colnames(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]]))])
                                                           samples <- str_replace_all(samples,"\\.","-")
                                                           survival <- survival_data_list[input$type][[1]]
                                                           biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                                                           km_fit <- survfit(Surv(PFI.time,PFI)~bicluster, data=biclist)
                                                           print(paste0("This bicluster contains: \n ",toString(nrow(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[1]]))," genes, \n ",toString(unique(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[1]]$Samples.In.Bicluster))," samples. \n","The KM-analysis outputs a p-value of ",round(surv_pvalue(km_fit,data=biclist)$pval,3)))
                                                       })}
                                                   else if (input$variable=="Overall survival"){
                                                       output$survivalvis <- renderPlot({
                                                           samples <- colnames(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]][2:length(colnames(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]]))])
                                                           samples <- str_replace_all(samples,"\\.","-")
                                                           survival <- survival_data_list[input$type][[1]]
                                                           biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                                                           km_fit <- survfit(Surv(OS.time,OS)~bicluster, data=biclist)
                                                           autoplot(km_fit) + 
                                                               labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
                                                                    title = " Overall Survival of \n Cancer Patients \n") + 
                                                               theme(plot.title = element_text(hjust = 0.5), 
                                                                     axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                                                                     axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                                                                     legend.title = element_text(face="bold", size = 10))})
                                                       output$bicinfo <- renderText({
                                                           samples <- colnames(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]][2:length(colnames(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]]))])
                                                           samples <- str_replace_all(samples,"\\.","-")
                                                           survival <- survival_data_list[input$type][[1]]
                                                           biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                                                           km_fit <- survfit(Surv(OS.time,OS)~bicluster, data=biclist)
                                                           print(paste0("This bicluster contains: \n ",toString(nrow(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[1]]))," genes, \n ",toString(unique(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[1]]$Samples.In.Bicluster))," samples. \n","The KM-analysis outputs a p-value of ",round(surv_pvalue(km_fit,data=biclist)$pval,3)))
                                                       })
                                                   }
                                                   else if (input$variable=="Disease-free interval"){
                                                       output$survivalvis <- renderPlot({
                                                           samples <- colnames(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]][2:length(colnames(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]]))])
                                                           samples <- str_replace_all(samples,"\\.","-")
                                                           survival <- survival_data_list[input$type][[1]]
                                                           biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                                                           km_fit <- survfit(Surv(DFI.time,DFI)~bicluster, data=biclist)
                                                           autoplot(km_fit) + 
                                                               labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
                                                                    title = " Disease-free Interval of \n Cancer Patients \n") + 
                                                               theme(plot.title = element_text(hjust = 0.5), 
                                                                     axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                                                                     axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                                                                     legend.title = element_text(face="bold", size = 10))})
                                                       output$bicinfo <- renderText({
                                                           samples <- colnames(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]][2:length(colnames(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]]))])
                                                           samples <- str_replace_all(samples,"\\.","-")
                                                           survival <- survival_data_list[input$type][[1]]
                                                           biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                                                           km_fit <- survfit(Surv(DFI.time,DFI)~bicluster, data=biclist)
                                                           print(paste0("This bicluster contains: \n ",toString(nrow(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[1]]))," genes, \n ",toString(unique(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[1]]$Samples.In.Bicluster))," samples. \n","The KM-analysis outputs a p-value of ",round(surv_pvalue(km_fit,data=biclist)$pval,3)))
                                                       })
                                                   }
                                                   else if (input$variable=="Disease-specific survival"){
                                                       output$survivalvis <- renderPlot({
                                                           samples <- colnames(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]][2:length(colnames(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]]))])
                                                           samples <- str_replace_all(samples,"\\.","-")
                                                           survival <- survival_data_list[input$type][[1]]
                                                           biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                                                           km_fit <- survfit(Surv(DSS.time,DSS)~bicluster, data=biclist)
                                                           autoplot(km_fit) + 
                                                               labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
                                                                    title = " Disease-specific Survival of \n Cancer Patients \n") + 
                                                               theme(plot.title = element_text(hjust = 0.5), 
                                                                     axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                                                                     axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                                                                     legend.title = element_text(face="bold", size = 10))})
                                                       output$bicinfo <- renderText({
                                                           samples <- colnames(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]][2:length(colnames(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[3]]))])
                                                           samples <- str_replace_all(samples,"\\.","-")
                                                           survival <- survival_data_list[input$type][[1]]
                                                           biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                                                           km_fit <- survfit(Surv(DSS.time,DSS)~bicluster, data=biclist)
                                                           print(paste0("This bicluster contains: \n ",toString(nrow(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[1]]))," genes, \n ",toString(unique(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[1]]$Samples.In.Bicluster))," samples. \n","The KM-analysis outputs a p-value of ",round(surv_pvalue(km_fit,data=biclist)$pval,3)))
                                                       })})
                                      # updateSelectInput(session,"path", choices = )
                                      output$survivaltable <- renderDataTable({
                                          datatable(down_cancer_list[input$type][[1]][[as.integer(input$bic)]][[1]]%>% select(Gene.ID,chrom),options = list(paging=FALSE))
                                      })
                                      })})
    
    # for copy number page
    observeEvent(input$reg_copy,
                 if (input$reg_copy == "Up"){
                     data_list <- up_data_list
                     observeEvent(input$type_copy,
                                  if(!is.null(input$type_copy)){
                                      dat_filtered <- data_list[input$type_copy][[1]] %>% filter(Samples.In.Bicluster > 20)
                                      updateSelectInput(session,"num", choices = dat_filtered$Bicluster.No)
                                      
                                      hbic_chrom_list <- list()
                                      for (i in unique(dat_filtered$Bicluster.No)){
                                          new <- up_cancer_list[input$type_copy][[1]][[i]][[2]]
                                          bic_chrom <- new$prop
                                          names(bic_chrom) <- new$chrom
                                          hbic_chrom_list <- list.append(hbic_chrom_list,names(which(bic_chrom>0.8)))
                                      }
                                      names(hbic_chrom_list) <- unique(dat_filtered$Bicluster.No)
                                      
                                      observeEvent(input$chrom,
                                                   if (input$chrom=="all"){
                                                       updateSelectInput(session,"num",choices = dat_filtered$Bicluster.No)
                                                   }
                                                   else {
                                                       updateSelectInput(session,"num",choices = names(which(hbic_chrom_list==input$chrom)))})
                                      
                                      observeEvent(input$withcopy,
                                                   if (input$withcopy=="No"){
                                                        output$mapvis <- renderPlotly({
                                                        dat <- up_cancer_list[input$type_copy][[1]][[as.integer(input$num)]][[2]]
                                                        ggplotly(dat %>% ggplot(aes(x=chrom,y=prop,fill="#46ACC8"))+geom_bar(stat='identity')+
                                                       labs(title = "") +
                                                       theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "light grey"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "light grey"),legend.position = "none",plot.title = element_text(size=10)) + 
                                                       coord_cartesian(ylim=c(0,NA)) +
                                                       xlab("Chromosome") + 
                                                       ylab("Percentage") + 
                                                       scale_y_continuous(labels=scales::percent_format()),tooltip="text")})}
                                                   else{
                                                       output$mapvis <- renderPlotly({
                                                           foo <- up_foo_list[input$type_copy][[1]][[as.integer(input$num)]]
                                                           ggplotly(foo %>% ggplot(aes(x=chrom,fill=copynumber))+geom_bar()+theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "light grey"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "light grey"),plot.title = element_text(size=10)) + 
                                                                        coord_cartesian(ylim=c(0,NA)) + 
                                                                        xlab("Chromosome") + scale_y_continuous(labels = function(x) x/10)+
                                                                        ylab("Percentage")+ scale_fill_manual(values = myColors),tooltip="text")})})
                                      output$info <- renderText({
                                          biclist <- up_cancer_list[input$type_copy][[1]][[as.integer(input$num)]][[1]]
                                          print(paste0("This bicluster contains: \n ",toString(nrow(biclist))," genes, \n ",toString(biclist$Samples.In.Bicluster[1])," samples."))
                                      })
                                      output$table <- renderDataTable({
                                          biclist <- up_cancer_list[input$type_copy][[1]][[as.integer(input$num)]][[1]]
                                          datatable(biclist %>% select(Gene.ID,chrom),options = list(paging=FALSE))
                                          
                                      })
                                      })}
                else {
                    data_list <- down_data_list
                    observeEvent(input$type_copy,
                                 if(!is.null(input$type_copy)){
                                     dat_filtered <- data_list[input$type_copy][[1]] %>% filter(Samples.In.Bicluster > 20)
                                     updateSelectInput(session,"num", choices = dat_filtered$Bicluster.No)
                                     lbic_chrom_list <- list()
                                     for (i in unique(dat_filtered$Bicluster.No)){
                                         new <- down_cancer_list[input$type_copy][[1]][[i]][[2]]
                                         bic_chrom <- new$prop
                                         names(bic_chrom) <- new$chrom
                                         lbic_chrom_list <- list.append(lbic_chrom_list,names(which(bic_chrom>0.8)))
                                     }
                                     names(lbic_chrom_list) <- unique(dat_filtered$Bicluster.No)
                                     
                                     observeEvent(input$chrom,
                                                  if (input$chrom=="all"){
                                                      updateSelectInput(session,"num",choices = dat_filtered$Bicluster.No)
                                                  }
                                                  else {
                                                      updateSelectInput(session,"num",choices = names(which(lbic_chrom_list==input$chrom)))})
                                     observeEvent(input$withcopy,
                                                  if (input$withcopy=="No"){
                                                      output$mapvis <- renderPlotly({
                                                          dat <- down_cancer_list[input$type_copy][[1]][[as.integer(input$num)]][[2]]
                                                          ggplotly(dat %>% ggplot(aes(x=chrom,y=prop,fill="#46ACC8"))+geom_bar(stat='identity')+
                                                                       labs(title = "") +
                                                                       theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "light grey"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "light grey"),legend.position = "none",plot.title = element_text(size=10)) + 
                                                                       coord_cartesian(ylim=c(0,NA)) +
                                                                       xlab("Chromosome") + 
                                                                       ylab("Percentage") + 
                                                                       scale_y_continuous(labels=scales::percent_format()),tooltip="text")})}
                                                  else{
                                                      output$mapvis <- renderPlotly({
                                                          foo <- down_foo_list[input$type_copy][[1]][[as.integer(input$num)]]
                                                          ggplotly(foo %>% ggplot(aes(x=chrom,fill=copynumber))+geom_bar()+theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "light grey"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "light grey"),plot.title = element_text(size=10)) + 
                                                                       coord_cartesian(ylim=c(0,NA)) + 
                                                                       xlab("Chromosome") + scale_y_continuous(labels = function(x) x/10)+
                                                                       ylab("Percentage")+ scale_fill_manual(values = myColors),tooltip="text")})})
                                     output$info <- renderText({
                                         biclist <- down_cancer_list[input$type_copy][[1]][[as.integer(input$num)]][[1]]
                                         print(paste0("This bicluster contains: \n ",toString(nrow(biclist))," genes, \n ",toString(biclist$Samples.In.Bicluster[1])," samples."))
                                     })
                                     output$table <- renderDataTable({
                                         biclist <- down_cancer_list[input$type_copy][[1]][[as.integer(input$num)]][[1]]
                                         datatable(biclist %>% select(Gene.ID,chrom),options = list(paging=FALSE))
                                         
                                     })})})}

# Run the application 
shinyApp(ui = ui, server = server)
