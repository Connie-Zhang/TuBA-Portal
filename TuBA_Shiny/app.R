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

source("load_data.R")

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
              shinyjs::useShinyjs(),
              fluidRow(
                shinyjs::useShinyjs(),
                column(3,selectInput(
                  inputId = "type",
                  label = "Choose cancer of interest:",
                  selected = "BLCA",
                  choices = full_cancer_names
                )),
                column(3, selectInput(
                  inputId = "reg",
                  label = "Up or down regulated gene expression:",
                  selected="Up",
                  choices = c("Up","Down"))),
                column(3,selectInput(
                  inputId = "variable",
                  label = "Select variable of interest",
                  choices = c("Overall survival"="OS", "Disease-specific survival"="DSS", "Disease-free interval"="DFI", "Progression-free interval"="PFI"),
                  selected = "Overall survival")),
                column(3,selectInput(
                  inputId = "pathsig",
                  label = "Filter by significance of biological pathway:",
                  choices = c("No"=20,"0.05 significance"=0.05,"0.10 significance"=0.10),
                  selected = "No"
                )),
                column(3, selectizeInput(
                  inputId = "gene",
                  label = "Gene of interest",
                  choices = c(unique(up_data_list[["BLCA"]] %>% pull(Gene.ID))),
                  selected = NULL,
                  multiple = TRUE)),
                column(3, selectizeInput(
                  inputId = "path",
                  label = "Biological Pathway of interest",
                  choices = c(unique(up_BP_full[["BLCA"]] %>% pull(GO_term))),
                  selected = NULL,
                  multiple = TRUE)),
                # change bic input choices to dynamic
                column(3,selectInput(
                  inputId = "bic",
                  label = "Select the bicluster of interest:",
                  choices = unique(up_BP_full[["BLCA"]] %>% pull(bic)),
                  selected = 1)),
                column(3,selectInput(
                  inputId = "sig",
                  label = "Filter by significance of survival analysis:",
                  choices = c("No"=20,"0.05 significance"=0.05,"0.10 significance"=0.10),
                  selected = "No"))),
              sidebarPanel(
                textOutput(outputId="bicinfo"),width = 3),
              mainPanel(
                tabsetPanel(
                  tabPanel("Visualization",plotOutput(outputId = "survivalvis"),width = 9),
                  tabPanel("Gene Information", dataTableOutput(outputId = "survivaltable"),textOutput(outputId = "bic_genes"),style = "height:500px; overflow-y: scroll;overflow-x: scroll;",width = 9)))),
      tabItem("copynum",
              useShinyjs(),
              fluidRow(
                column(2, selectInput(
                  inputId = "type_copy",
                  label = "Choose cancer of interest:",
                  selected = "BLCA",
                  choices = full_cancer_names)),
                column(3, selectInput(
                  inputId = "reg_copy",
                  label = "Up or down regulated gene expression:",
                  selected="Up",
                  choices = c("Up","Down"))),
                column(3,selectInput(
                  inputId = "num",
                  label = "Select the bicluster of interest:",
                  choices = unique(up_BP_full[["BLCA"]] %>% pull(bic)))),
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
                  tabPanel("Gene Information",dataTableOutput(outputId = "table"),textOutput(outputId="bic_genes2"), style = "height:500px; overflow-y: scroll;overflow-x: scroll;",width = 9)))),
      tabItem("contact",h4("to be added")))))

# Define server logic required to draw a histogram
server <- function(input, output,session) {
  toListen <- reactive({
    list(input$reg,input$type,input$variable,input$sig,input$gene,input$path,input$pathsig)
  })
  observeEvent(toListen(),{
    if(!is.null(input$type)){
      dat_filtered <- full_data_list[[input$reg]][input$type_copy][[1]] %>% filter(Samples.In.Bicluster > 20)
      updateSelectInput(session,"bic", choices = dat_filtered$Bicluster.No)}
    
      pval_filtered <- full_pval_list[[input$reg]][[input$type]] %>% filter((!!sym(input$variable))<input$sig) %>% pull(bic)
      if(is.null(input$path)){
        path_filtered <- NULL
      }
      else{
        path_filtered <- full_BP_list[[input$reg]][[input$type]] %>% filter(GO_term==input$path) %>% filter(p_val < input$pathsig) %>% pull(bic)}
      
      if(is.null(input$gene)){
        gene_filtered <- NULL
      }
      else {
        gene_filtered <- full_data_list[[input$reg]][[input$type]] %>% filter(Gene.ID %in% input$gene) %>% count(Bicluster.No) %>% filter(n==length(input$gene)) %>% pull(Bicluster.No)}
      if (length(gene_filtered)==0){gene_filtered <- NULL}
      if(length(path_filtered)==0){path_filtered <- NULL}
      if(length(pval_filtered)==0){pval_filtered <- NULL}
      mylist <- list(pval_filtered,path_filtered,gene_filtered)
      updated_bic <- Reduce(intersect,mylist[vapply(mylist, Negate(is.null),NA)])
      if(length(updated_bic)!=0){
        updateSelectInput(session,"bic",choices=updated_bic)}
      else if (length(updated_bic)==0){
        updateSelectInput(session,"bic",choices = NULL)}
      if(is.null(input$gene)&is.null(input$path)&is.null(input$sig)&is.null(input$pathsig)){
        updateSelectInput(session,"bic",choices=full_data_list[[input$reg]][[input$type]] %>% pull(Bicluster.No))
        
      }
  })

  
  output$survivalvis <- renderPlot({
    samples <- colnames(full_cancer_list[[input$reg]][[input$type]][[as.integer(input$bic)]][[3]])[-1]
    samples <- str_replace_all(samples,"\\.","-")
    survival <- survival_data_list[[input$type]]
    biclist <- survival %>% mutate(bicluster = as.factor(ifelse(sample %in% samples,1,0)))

    form <- as.formula(paste0("Surv(",input$variable,".time,",input$variable,")~bicluster"))
    km_fit <- survfit(form, data=biclist)
    # fortify(km_fit)
    ggplot2::autoplot(km_fit) +
      labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n",
           title = paste0(variable_names[[input$variable]]," of \n Cancer Patients \n"), colour = "Samples") +
      theme(plot.title = element_text(hjust = 0.5),
            axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
            axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
            legend.title = element_text(face="bold", size = 10),
            panel.background = element_rect(fill="white"),
            panel.grid = element_line(colour="grey", size = 0.25)) +
      guides(fill=FALSE) +
      scale_color_manual(labels = c("not in bicluster", "in bicluster"), values = c(2,1))
    })
  output$bicinfo <- renderText({
    dat <- full_cancer_list[[input$reg]][[input$type]][[as.integer(input$bic)]][[3]]
    samples <- colnames(dat)[-1]
    samples <- str_replace_all(samples,"\\.","-")
    survival <- survival_data_list[[input$type]]
    biclist <- survival %>% mutate(bicluster = as.factor(ifelse(sample %in% samples,1,0)))

    form <- as.formula(paste0("Surv(",input$variable,".time,",input$variable,")~bicluster"))
    km_fit <- surv_fit(form, data=biclist)
    print(paste0("This bicluster contains: \n ",toString(nrow(full_cancer_list[[input$reg]][[input$type]][[as.integer(input$bic)]][[1]]))," genes, \n ",toString(unique(full_cancer_list[[input$reg]][input$type][[1]][[as.integer(input$bic)]][[1]]$Samples.In.Bicluster))," samples. \n","The KM-analysis outputs a p-value of ",round(surv_pvalue(km_fit,data=biclist)$pval,3)))
  })
  
  observeEvent(input$type_copy,
               if(!is.null(input$type_copy)){
                 dat_filtered <- full_data_list[[input$reg_copy]][input$type_copy][[1]] %>% filter(Samples.In.Bicluster > 20)
                 updateSelectInput(session,"num", choices = dat_filtered$Bicluster.No)
                 
                 hbic_chrom_list <- list()
                 for (i in unique(dat_filtered$Bicluster.No)){
                   new <- full_cancer_list[[input$reg_copy]][input$type_copy][[1]][[i]][[2]]
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
                                  dat <- full_cancer_list[[input$reg_copy]][input$type_copy][[1]][[as.integer(input$num)]][[2]]
                                  ggplotly(dat %>% ggplot(aes(x=chrom,y=prop,fill="#46ACC8"))+geom_bar(stat='identity')+
                                             labs(title = "") +
                                             theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "light grey"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "light grey"),legend.position = "none",plot.title = element_text(size=10)) + 
                                             coord_cartesian(ylim=c(0,NA)) +
                                             xlab("Chromosome") + 
                                             ylab("Percentage") + 
                                             scale_y_continuous(labels=scales::percent_format()),tooltip="text")})}
                              else{
                                output$mapvis <- renderPlotly({
                                  foo <- full_foo_list[[input$reg_copy]][input$type_copy][[1]][[as.integer(input$num)]]
                                  ggplotly(foo %>% ggplot(aes(x=chrom,fill=copynumber))+geom_bar()+theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "light grey"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "light grey"),plot.title = element_text(size=10)) + 
                                             coord_cartesian(ylim=c(0,NA)) + 
                                             xlab("Chromosome") + scale_y_continuous(labels = function(x) x/10)+
                                             ylab("Percentage")+ scale_fill_manual(values = myColors),tooltip="text")})})
                 output$info <- renderText({
                   # biclist <- full_cancer_list[[input$reg_copy]][input$type_copy][[1]][[as.integer(input$num)]][[1]]
                   # print(paste0("This bicluster contains: \n ",toString(nrow(biclist))," genes, \n ",toString(biclist$Samples.In.Bicluster[1])," samples."))
                   dat <- full_cancer_list[[input$reg_copy]][[input$type_copy]][[as.integer(input$num)]][[3]]
                   samples <- colnames(dat)[-1]
                   samples <- str_replace_all(samples,"\\.","-")
                   survival <- survival_data_list[[input$type_copy]]
                   biclist <- survival %>% mutate(bicluster = as.factor(ifelse(sample %in% samples,1,0)))
                   
                   form <- as.formula(paste0("Surv(",input$variable,".time,",input$variable,")~bicluster"))
                   km_fit <- surv_fit(form, data=biclist)
                   print(paste0("This bicluster contains: \n ",toString(nrow(full_cancer_list[[input$reg_copy]][[input$type_copy]][[as.integer(input$num)]][[1]]))," genes, \n ",toString(unique(full_cancer_list[[input$reg_copy]][input$type_copy][[1]][[as.integer(input$num)]][[1]]$Samples.In.Bicluster))," samples. \n","The KM-analysis outputs a p-value of ",round(surv_pvalue(km_fit,data=biclist)$pval,3)))
                 })
                 output$table <- renderDataTable({
                   biclist <- full_cancer_list[[input$reg_copy]][input$type_copy][[1]][[as.integer(input$num)]][[1]]
                   datatable(biclist %>% select(Gene.ID,chrom),options = list(paging=FALSE))
                 })
                 output$bic_genes2 <- renderText({
                   paste0('\"',paste(full_cancer_list[[input$reg_copy]][input$type][[1]][[as.integer(input$bic)]][[1]]%>% pull(Gene.ID),collapse='","'),'\"')
                 })

})}

# Run the application 
shinyApp(ui = ui, server = server)
