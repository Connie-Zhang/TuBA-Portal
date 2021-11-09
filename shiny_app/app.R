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

load("biclusterlist2.RData")
load("survival.RData")
load("foolist2.RData")
load("foolist3.RData")
load("BLCA_PrimaryTumors_Cleaned_H0_05_JcdInd0_3_MinGenes2_MinSamples2_GenesInBiclusters_EdgeBasedSampleEnrichment.RData")
load("BLCA_PrimaryTumors_Cleaned_L0_05_JcdInd0_3_MinGenes2_MinSamples2_GenesInBiclusters_EdgeBasedSampleEnrichment.RData")
load("locus.RData")

round_preserve_sum <- function(x, digits = 0) {
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}

chrom <- biclusterlist2[[1]][[2]]
chrom <- chrom$chrom # where did "All" go
chrom <- c(as.character(chrom),"all")

myColors <- c("#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20","#808080")
names(myColors) <- c("high_amplification","amplification","no_change","1_copy_del","2_copy_del","Missing")

hbic <- BLCA_PrimaryTumors_Cleaned_H0_05_JcdInd0_3_MinGenes2_MinSamples2_GenesInBiclusters_EdgeBasedSampleEnrichment
hbic_filtered <- hbic %>% filter(Samples.In.Bicluster>20)

lbic <- BLCA_PrimaryTumors_Cleaned_L0_05_JcdInd0_3_MinGenes2_MinSamples2_GenesInBiclusters_EdgeBasedSampleEnrichment
lbic_filtered <- lbic %>% filter(Samples.In.Bicluster>20)

bladder_hgenes <- unique(BLCA_PrimaryTumors_Cleaned_H0_05_JcdInd0_3_MinGenes2_MinSamples2_GenesInBiclusters_EdgeBasedSampleEnrichment$Gene.ID)
bladder_hgenes <- c(bladder_hgenes, "all")
bladder_lgenes <- unique(BLCA_PrimaryTumors_Cleaned_L0_05_JcdInd0_3_MinGenes2_MinSamples2_GenesInBiclusters_EdgeBasedSampleEnrichment$Gene.ID)
bladder_lgenes <- c(bladder_lgenes, "all")
bladder_hBP <- unique(c(BP_info_H$GoTerm1,BP_info_H$GoTerm2,BP_info_H$GoTerm3,BP_info_H$GoTerm4,BP_info_H$GoTerm5))
bladder_hBP <- c(bladder_hBP, "all")
bladder_lBP <- unique(c(BP_info_L$GoTerm1,BP_info_L$GoTerm2,BP_info_L$GoTerm3,BP_info_L$GoTerm4,BP_info_L$GoTerm5))
bladder_lBP <- c(bladder_lBP, "all")

chrom_list <- unique(map$chrom)

hbic_gene_list <- list()

for (i in unique(hbic_filtered$Bicluster.No)){
  new <- hbic_filtered %>% filter(Bicluster.No==i)
  genes <- new$Gene.ID
  hbic_gene_list <- list.append(hbic_gene_list,genes)
}
names(hbic_gene_list) <- unique(hbic_filtered$Bicluster.No)

lbic_gene_list <- list()

for (i in unique(lbic_filtered$Bicluster.No)){
  new <- lbic_filtered %>% filter(Bicluster.No==i)
  genes <- new$Gene.ID
  lbic_gene_list <- list.append(lbic_gene_list,genes)
}
names(lbic_gene_list) <- unique(lbic_filtered$Bicluster.No)

hbic_BP_list <- list()
for (i in unique(hbic_filtered$Bicluster.No)){
  new <- BP_info_H %>% filter(Bicluster.No==i)
  BPs <- c(new$GoTerm1, new$GoTerm2, new$GoTerm3, new$GoTerm4, new$GoTerm5)
  hbic_BP_list <- list.append(hbic_BP_list,BPs)
}
names(hbic_BP_list) <- unique(hbic_filtered$Bicluster.No)

lbic_BP_list<- list()
for (i in unique(lbic_filtered$Bicluster.No)){
  new <- BP_info_L %>% filter(Bicluster.No==i)
  BPs <- c(new$GoTerm1, new$GoTerm2, new$GoTerm3, new$GoTerm4, new$GoTerm5)
  lbic_BP_list <- list.append(lbic_BP_list,BPs)
}
names(lbic_BP_list) <- unique(lbic_filtered$Bicluster.No)

hbic_chrom_list <- list()
for (i in unique(hbic_filtered$Bicluster.No)){
  new <- biclusterlist2[[i]][[2]]
  bic_chrom <- new$prop
  names(bic_chrom) <- new$chrom
  hbic_chrom_list <- list.append(hbic_chrom_list,names(which(bic_chrom>0.8)))
}
names(hbic_chrom_list) <- unique(hbic_filtered$Bicluster.No)

lbic_chrom_list <- list()
for (i in unique(lbic_filtered$Bicluster.No)){
  new <- biclusterlist3[[i]][[2]]
  bic_chrom <- new$prop
  names(bic_chrom) <- new$chrom
  lbic_chrom_list <- list.append(lbic_chrom_list,names(which(bic_chrom>0.8)))
}
names(lbic_chrom_list) <- unique(lbic_filtered$Bicluster.No)


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
                column(3, selectInput(
                  inputId = "reg",
                  label = "Choose up or down regulated gene expression:",
                  selected="Up",
                  choices = c("Up","Down"))),
                column(3,selectInput(
                  inputId = "variable",
                  label = "Select variable of interest",
                  choices = c("Overall survival", "Disease-specific survival", "Disease-free interval", "Progression-free interval"))),
                # be able to input multiple genes
                column(3, selectizeInput(
                  inputId = "gene",
                  label = "Gene of interest",
                  choices = bladder_hgenes,
                  multiple = TRUE)),
                column(3, selectizeInput(
                  inputId = "path",
                  label = "Biological Pathway of interest",
                  choices = bladder_hBP,
                  multiple=TRUE)),
                # change bic input choices to dynamic
                column(3,selectInput(
                  inputId = "bic",
                  label = "Select the bicluster of interest:",
                  choices = unique(hbic_filtered$Bicluster.No)))),
              sidebarPanel(
                textOutput(outputId="bicinfo"),width = 3),
              mainPanel(
                tabsetPanel(
                  tabPanel("Visualization",plotOutput(outputId = "survivalvis"),width = 9),
                  tabPanel("Gene Information", dataTableOutput(outputId = "survivaltable"),style = "height:500px; overflow-y: scroll;overflow-x: scroll;",width = 9)))),
      tabItem("copynum",
              useShinyjs(),
              fluidRow(
                column(3, selectInput(
                  inputId = "reg_copy",
                  label = "Choose up or down regulated gene expression:",
                  selected="Up",
                  choices = c("Up","Down"))),
                column(3,selectInput(
                  inputId = "num",
                  label = "Select the bicluster of interest:",
                  choices = unique(hbic_filtered$Bicluster.No))),
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

  
  server <- function(input, output, session) {
    gene_list <- list()
    observeEvent(input$reg,{
      if(input$reg=="Up"){
        biclist <- biclusterlist2
        gene_list <- hbic_gene_list
        path_list <- hbic_BP_list
        updateSelectizeInput(session,"gene",choices = bladder_hgenes)
        updateSelectizeInput(session,"path",choices = bladder_hBP)
        # if the input gene is not in any of the biclusters, the bicluster list won't update
        observeEvent(input$gene,
          if (input$gene=="all") {
            updateSelectInput(session,"bic",choices = unique(hbic_filtered$Bicluster.No))
          }
          else {
          updated_bic<-c()
          for (i in 1:length(gene_list)){
            if (input$gene %in% gene_list[[i]]) {
              updated_bic <- append(updated_bic, names(gene_list[i]))
            }
          }
          updateSelectInput(session,"bic",choices = updated_bic)
        })
        observeEvent(input$path,
          if (input$path=="all"){
            updateSelectInput(session,"bic",choices = unique(hbic_filtered$Bicluster.No))
          }
          else {
          updated_bic<-c()
          for (i in 1:length(path_list)){
            if (input$path %in% path_list[[i]]) {
              updated_bic <- append(updated_bic, names(path_list[i]))
            }
          }
          updateSelectInput(session,"bic",choices = updated_bic)
        })
        output$bicinfo <- renderText({
          samples <- colnames(biclist[[as.integer(input$bic)]][[5]])[2:length(colnames(biclist[[as.integer(input$bic)]][[5]]))]
          surv_dat <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
          km_fit <- survfit(Surv(PFI.time,PFI)~bicluster, data=surv_dat)
          print(paste0("This bicluster contains: \n ",toString(nrow(biclist[[as.integer(input$bic)]][[5]]))," genes, \n ",toString(unique(biclist[[as.integer(input$bic)]][[1]]$Samples.In.Bicluster))," samples. \n","The KM-analysis outputs a p-value of ",round(surv_pvalue(km_fit,data=surv_dat)$pval,3)))
        })
        observeEvent(input$variable,
                     if (input$variable=="Progression-free interval"){
                       output$survivalvis <- 
                         output$survivalvis <- renderPlot({
                           samples <- colnames(biclist[[as.integer(input$bic)]][[5]])[2:length(colnames(biclist[[as.integer(input$bic)]][[5]]))]
                           biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                           km_fit <- survfit(Surv(PFI.time,PFI)~bicluster, data=biclist)
                           autoplot(km_fit) + 
                             labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
                                  title = " Progression-free Interval Of \n Bladder Cancer Patients \n") + 
                             theme(plot.title = element_text(hjust = 0.5), 
                                   axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                                   axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                                   legend.title = element_text(face="bold", size = 10))})
                     }
                     else if (input$variable=="Overall survival"){
                       output$survivalvis <- renderPlot({
                         samples <- colnames(biclist[[as.integer(input$bic)]][[5]])[2:length(colnames(biclist[[as.integer(input$bic)]][[5]]))]
                         biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                         km_fit <- survfit(Surv(OS.time,OS)~bicluster, data=biclist)
                         autoplot(km_fit) + 
                           labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
                                title = "Overall Survival Times Of \n Bladder Cancer Patients \n") + 
                           theme(plot.title = element_text(hjust = 0.5), 
                                 axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                                 axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                                 legend.title = element_text(face="bold", size = 10))})
                     }
                     else if (input$variable=="Disease-specific survival"){
                       output$survivalvis <- renderPlot({
                         samples <- colnames(biclist[[as.integer(input$bic)]][[5]])[2:length(colnames(biclist[[as.integer(input$bic)]][[5]]))]
                         biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                         km_fit <- survfit(Surv(DSS.time,DSS)~bicluster, data=biclist)
                         autoplot(km_fit) + 
                           labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
                                title = "Disease-specific Survival Times Of \n Bladder Cancer Patients \n") + 
                           theme(plot.title = element_text(hjust = 0.5), 
                                 axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                                 axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                                 legend.title = element_text(face="bold", size = 10))})
                     }
                     else if (input$variable=="Disease-free interval"){
                       output$survivalvis <- renderPlot({
                         samples <- colnames(biclist[[as.integer(input$bic)]][[5]])[2:length(colnames(biclist[[as.integer(input$bic)]][[5]]))]
                         survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0)) %>% ggplot(aes(x=DFI.time,fill=as.factor(bicluster)))+geom_boxplot(alpha=0.5)+
                           scale_fill_manual(values = c("#78B7C5","#F2AD00"),name = paste0("Bicluster ",as_string(input$bic)), labels = c("No","Yes"),guide = guide_legend(reverse=TRUE)) +
                           theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "light grey"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "light grey"),plot.title = element_text(size=10))})
                       output$survivalvis <- renderPlot({
                         samples <- colnames(biclist[[as.integer(input$bic)]][[5]])[2:length(colnames(biclist[[as.integer(input$bic)]][[5]]))]
                         biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                         km_fit <- survfit(Surv(DFI.time,DFI)~bicluster, data=biclist)
                         autoplot(km_fit) + 
                           labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
                                title = "Disease-free Interval Times Of \n Bladder Cancer Patients \n") + 
                           theme(plot.title = element_text(hjust = 0.5), 
                                 axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                                 axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                                 legend.title = element_text(face="bold", size = 10))})
                     }
        )
        output$survivaltable <- renderDataTable({
          datatable(biclist[[as.integer(input$bic)]][[1]]%>% select(Gene.ID,chrom),options = list(paging=FALSE))
        })
        }
      else if(input$reg=="Down") {
        biclist <- biclusterlist3
        gene_list <- lbic_gene_list
        path_list <- lbic_BP_list
        updateSelectizeInput(session,"gene",choices = bladder_lgenes)
        updateSelectizeInput(session,"path",choices = bladder_lBP)
        updateSelectInput(session,"bic", choices = unique(lbic_filtered$Bicluster.No))
        observeEvent(input$gene,
          if (input$gene=="all") {
              updateSelectInput(session,"bic",choices = unique(lbic_filtered$Bicluster.No))
              }
          else {
          updated_bic<-c()
          for (i in 1:length(gene_list)){
            if (input$gene %in% gene_list[[i]]) {
              updated_bic <- append(updated_bic, names(gene_list[i]))
            }
          }
          updateSelectInput(session,"bic",choices = updated_bic)})
        observeEvent(input$path,
          if (input$path=="all") {
              updateSelectInput(session,"bic",choices = unique(lbic_filtered$Bicluster.No))}
          else {
          updated_bic<-c()
          for (i in 1:length(path_list)){
            if (input$path %in% path_list[[i]]) {
              updated_bic <- append(updated_bic, names(path_list[i]))
            }
          }
          updateSelectInput(session,"bic",choices = updated_bic)
        })
        }
      output$bicinfo <- renderText({
        # biclist <- biclist[[as.integer(input$bic)]][[1]]
        samples <- colnames(biclist[[as.integer(input$bic)]][[5]])[2:length(colnames(biclist[[as.integer(input$bic)]][[5]]))]
        surv_dat <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
        km_fit <- survfit(Surv(PFI.time,PFI)~bicluster, data=surv_dat)
        print(paste0("This bicluster contains: \n ",toString(nrow(biclist[[as.integer(input$bic)]][[5]]))," genes, \n ",toString(unique(biclist[[as.integer(input$bic)]][[1]]$Samples.In.Bicluster))," samples. \n","The KM-analysis outputs a p-value of ",round(surv_pvalue(km_fit,data=surv_dat)$pval,3)))
      })
      observeEvent(input$variable,
                   if (input$variable=="Progression-free interval"){
                     output$survivalvis <- 
                       output$survivalvis <- renderPlot({
                         samples <- colnames(biclist[[as.integer(input$bic)]][[5]])[2:length(colnames(biclist[[as.integer(input$bic)]][[5]]))]
                         biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                         km_fit <- survfit(Surv(PFI.time,PFI)~bicluster, data=biclist)
                         autoplot(km_fit) + 
                           labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
                                title = " Progression-free Interval Of \n Bladder Cancer Patients \n") + 
                           theme(plot.title = element_text(hjust = 0.5), 
                                 axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                                 axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                                 legend.title = element_text(face="bold", size = 10))})
                   }
                   else if (input$variable=="Overall survival"){
                     output$survivalvis <- renderPlot({
                       samples <- colnames(biclist[[as.integer(input$bic)]][[5]])[2:length(colnames(biclist[[as.integer(input$bic)]][[5]]))]
                       biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                       km_fit <- survfit(Surv(OS.time,OS)~bicluster, data=biclist)
                       autoplot(km_fit) + 
                         labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
                              title = "Overall Survival Times Of \n Bladder Cancer Patients \n") + 
                         theme(plot.title = element_text(hjust = 0.5), 
                               axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                               axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                               legend.title = element_text(face="bold", size = 10))})
                   }
                   else if (input$variable=="Disease-specific survival"){
                     output$survivalvis <- renderPlot({
                       samples <- colnames(biclist[[as.integer(input$bic)]][[5]])[2:length(colnames(biclist[[as.integer(input$bic)]][[5]]))]
                       biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                       km_fit <- survfit(Surv(DSS.time,DSS)~bicluster, data=biclist)
                       autoplot(km_fit) + 
                         labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
                              title = "Disease-specific Survival Times Of \n Bladder Cancer Patients \n") + 
                         theme(plot.title = element_text(hjust = 0.5), 
                               axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                               axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                               legend.title = element_text(face="bold", size = 10))})
                   }
                   else if (input$variable=="Disease-free interval"){
                     output$survivalvis <- renderPlot({
                       samples <- colnames(biclist[[as.integer(input$bic)]][[5]])[2:length(colnames(biclist[[as.integer(input$bic)]][[5]]))]
                       survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0)) %>% ggplot(aes(x=DFI.time,fill=as.factor(bicluster)))+geom_boxplot(alpha=0.5)+
                         scale_fill_manual(values = c("#78B7C5","#F2AD00"),name = paste0("Bicluster ",as_string(input$bic)), labels = c("No","Yes"),guide = guide_legend(reverse=TRUE)) +
                         theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "light grey"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "light grey"),plot.title = element_text(size=10))})
                     output$survivalvis <- renderPlot({
                       samples <- colnames(biclist[[as.integer(input$bic)]][[5]])[2:length(colnames(biclist[[as.integer(input$bic)]][[5]]))]
                       biclist <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
                       km_fit <- survfit(Surv(DFI.time,DFI)~bicluster, data=biclist)
                       autoplot(km_fit) + 
                         labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
                              title = "Disease-free Interval Times Of \n Bladder Cancer Patients \n") + 
                         theme(plot.title = element_text(hjust = 0.5), 
                               axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                               axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                               legend.title = element_text(face="bold", size = 10))})
                   }
      )
      output$survivaltable <- renderDataTable({
        datatable(biclist[[as.integer(input$bic)]][[1]] %>% select(Gene.ID,chrom),options = list(paging=FALSE))
      })})
    
    # find biclusters witth more than 80% from the same chromosome. 
    # Give each bicluster like this a chromosome label
    
    observeEvent(input$reg_copy,
                 if (input$reg_copy=="Up"){
                   biclist <- biclusterlist2
                   observeEvent(input$withcopy,
                                if (input$withcopy=="No"){
                                  output$mapvis <- renderPlotly({
                                    biclist <- biclist[[as.integer(input$num)]][[2]]
                                    ggplotly(biclist %>% ggplot(aes(x=chrom,y=prop,fill=col))+geom_bar(stat='identity')+
                                               labs(title = "") +
                                               theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "light grey"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "light grey"),legend.position = "none",plot.title = element_text(size=10)) + 
                                               coord_cartesian(ylim=c(0,NA)) + scale_fill_manual(values = levels(as.factor(biclist$col)))+
                                               xlab("Chromosome") + 
                                               ylab("Percentage") + 
                                               scale_y_continuous(labels=scales::percent_format()),tooltip="text")})
                                }
                                else {
                                  output$mapvis <- renderPlotly({
                                    foo <- foolist2[[as.integer(input$num)]]
                                    ggplotly(foo %>% ggplot(aes(x=chrom,fill=copynumber,text=paste(copynumber,": \n",round(copyprop*100,1),"%")))+geom_bar()+theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "light grey"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "light grey"),plot.title = element_text(size=10)) + 
                                               coord_cartesian(ylim=c(0,NA)) + 
                                               xlab("Chromosome") + scale_y_continuous(labels = function(x) x/10)+
                                               ylab("Percentage")+ scale_fill_manual(values = myColors),tooltip="text")})})
                   observeEvent(input$chrom,
                                if (input$chrom=="all"){
                                  updateSelectInput(session,"num",choices = hbic_filtered$Bicluster.No)
                                }
                                else {
                                updateSelectInput(session,"num",choices = names(which(hbic_chrom_list==input$chrom)))})
                   output$info <- renderText({
                     biclist <- biclist[[as.integer(input$num)]][[1]]
                     print(paste0("This bicluster contains: \n ",toString(nrow(biclist))," genes, \n ",toString(biclist$Samples.In.Bicluster[1])," samples."))
                   })
                   output$table <- renderDataTable({
                     biclist <- biclist[[as.integer(input$num)]][[1]]
                     datatable(biclist %>% select(Gene.ID,chrom),options = list(paging=FALSE))
                     
                   })}
                 else if (input$reg_copy=="Down"){
                   biclist <- biclusterlist3
                   updateSelectInput(session,"num", choices = unique(lbic_filtered$Bicluster.No))
                   observeEvent(input$withcopy,
                                if (input$withcopy=="No"){
                                  output$mapvis <- renderPlotly({
                                    biclist <- biclist[[as.integer(input$num)]][[2]]
                                    ggplotly(biclist %>% ggplot(aes(x=chrom,y=prop,fill=col))+geom_bar(stat='identity')+
                                               labs(title = "") +
                                               theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "light grey"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "light grey"),legend.position = "none",plot.title = element_text(size=10)) + 
                                               coord_cartesian(ylim=c(0,NA)) + scale_fill_manual(values = levels(as.factor(biclist$col)))+
                                               xlab("Chromosome") + 
                                               ylab("Percentage") + 
                                               scale_y_continuous(labels=scales::percent_format()),tooltip="text")})
                                }
                                else {
                                  output$mapvis <- renderPlotly({
                                    foo <- foolist3[[as.integer(input$num)]]
                                    ggplotly(foo %>% ggplot(aes(x=chrom,fill=copynumber,text=paste(copynumber,": \n",round(copyprop*100,1),"%")))+geom_bar()+theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "light grey"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "light grey"),plot.title = element_text(size=10)) + 
                                               coord_cartesian(ylim=c(0,NA)) + 
                                               xlab("Chromosome") + scale_y_continuous(labels = function(x) x/10)+
                                               ylab("Percentage")+ scale_fill_manual(values = myColors),tooltip="text")})})
                   observeEvent(input$chrom,
                                if (input$chrom=="all"){
                                  updateSelectInput(session,"num",choices = lbic_filtered$Bicluster.No)
                                }
                                else {
                                  updateSelectInput(session,"num",choices = names(which(lbic_chrom_list==input$chrom)))})
                   output$info <- renderText({
                     biclist <- biclist[[as.integer(input$num)]][[1]]
                     print(paste0("This bicluster contains: \n ",toString(nrow(biclist))," genes, \n ",toString(biclist$Samples.In.Bicluster[1])," samples."))
                   })
                   output$table <- renderDataTable({
                     biclist <- biclist[[as.integer(input$num)]][[1]]
                     datatable(biclist %>% select(Gene.ID,chrom),options = list(paging=FALSE))
                   })})

  }
  
  
  shinyApp(ui = ui, server = server)
  
  
  