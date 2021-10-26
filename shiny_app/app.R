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

load("biclusterlist2.RData")
load("survival.RData")
load("foolist2.RData")

round_preserve_sum <- function(x, digits = 0) {
    up <- 10 ^ digits
    x <- x * up
    y <- floor(x)
    indices <- tail(order(x-y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    y / up
}

myColors <- c("#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20","#808080")
names(myColors) <- c("high_amplification","amplification","no_change","1_copy_del","2_copy_del","Missing")

ui <- dashboardPage (
    
    dashboardHeader(title="Bicluster Visualizations"),
    
    dashboardSidebar(
        sidebarMenu(
            menuItem("Introduction", tabName = "intro"),
            menuItem("Survival Visualizations", tabName = "survival"),
            menuItem("Probe Map Visualizations", tabName="map"),
            menuItem("Contact", tabName = "contact")
        )
    ),
    
    dashboardBody(
        tabItems(
            tabItem("intro",h4("to be added")),
            tabItem("survival",
                    useShinyjs(),
                    fluidRow(
                            column(4,selectInput(
                            inputId = "bic",
                            label = "Select the bicluster of interest:",
                            choices = c(1:663))),
                            column(4,selectInput(
                                inputId = "variable",
                                label = "Select variable of interest",
                                choices = c("Overall survival", "Disease-specific survival", "Disease-free interval", "Progression-free interval")))),
                    sidebarPanel(
                        textOutput(outputId="bicinfo"),width = 3),
                    mainPanel(
                        tabsetPanel(
                            tabPanel("Visualization",plotOutput(outputId = "survivalvis"),width = 9),
                            tabPanel("Gene Information", tableOutput(outputId = "survivaltable"),width = 9)))),
            tabItem("map",
                    useShinyjs(),
                    fluidRow(
                        column(4,selectInput(
                            inputId = "num",
                            label = "Select the bicluster of interest:",
                            choices = c(1:663))),
                        column(4,selectInput(
                            inputId = "withcopy",
                            label = "Fill with copy number?",
                            choices = c("Yes","No")))),
                    sidebarPanel(
                        textOutput(outputId= "info"), width = 3),
                    mainPanel(
                        tabsetPanel(
                            tabPanel("Visualization",plotlyOutput(outputId = "mapvis"),width = 9),
                            tabPanel("Gene Information",tableOutput(outputId = "table"),width = 9)))),
            tabItem("contact",h4("to be added"))
    )))

server <- function(input, output, session) {
    
    output$bicinfo <- renderText({
        dat <- biclusterlist2[[as.integer(input$bic)]][[1]]
        samples <- colnames(biclusterlist2[[as.integer(input$bic)]][[5]])[2:length(colnames(biclusterlist2[[as.integer(input$bic)]][[5]]))]
        surv_dat <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
        km_fit <- survfit(Surv(PFI.time,PFI)~bicluster, data=surv_dat)
        print(paste0("This bicluster contains: \n ",toString(nrow(dat))," genes, \n ",toString(dat$Samples.In.Bicluster[1])," samples. \n","The KM-analysis outputs a p-value of ",round(surv_pvalue(km_fit,data=surv_dat)$pval,3)))
        # print(paste0("The KM-analysis outputs a p-value of ",round(surv_pvalue(km_fit,data=surv_dat)$pval,3)))
    })
    observeEvent(input$variable,
        if (input$variable=="Progression-free interval"){
    output$survivalvis <- 
        # renderPlot({
        # samples <- colnames(biclusterlist2[[as.integer(input$bic)]][[5]])[2:length(colnames(biclusterlist2[[as.integer(input$bic)]][[5]]))]
        # survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0)) %>% ggplot(aes(x=PFI.time,fill=as.factor(bicluster)))+geom_boxplot(alpha=0.5)+
        #     scale_fill_manual(values = c("#78B7C5","#F2AD00"), name = paste0("Bicluster ",as_string(input$bic)), labels = c("No","Yes"),guide = guide_legend(reverse=TRUE)) +
        #     theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "light grey"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "light grey"),plot.title = element_text(size=10))})
    output$survivalvis <- renderPlot({
        samples <- colnames(biclusterlist2[[as.integer(input$bic)]][[5]])[2:length(colnames(biclusterlist2[[as.integer(input$bic)]][[5]]))]
        dat <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
        km_fit <- survfit(Surv(PFI.time,PFI)~bicluster, data=dat)
        autoplot(km_fit) + 
            labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
                 title = " Progression-free Interval Of \n Bladder Cancer Patients \n") + 
            theme(plot.title = element_text(hjust = 0.5), 
                  axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                  axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                  legend.title = element_text(face="bold", size = 10))})
    }
        else if (input$variable=="Overall survival"){
        # output$survivalvis <- renderPlot({
        #     samples <- colnames(biclusterlist2[[as.integer(input$bic)]][[5]])[2:length(colnames(biclusterlist2[[as.integer(input$bic)]][[5]]))]
        #     survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0)) %>% ggplot(aes(x=OS.time,fill=as.factor(bicluster)))+geom_boxplot(alpha=0.5)+
        #         scale_fill_manual(values = c("#78B7C5","#F2AD00"),name = paste0("Bicluster ",as_string(input$bic)), labels = c("No","Yes"),guide = guide_legend(reverse=TRUE)) +
        #         theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "light grey"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "light grey"),plot.title = element_text(size=10))})
        output$survivalvis <- renderPlot({
            samples <- colnames(biclusterlist2[[as.integer(input$bic)]][[5]])[2:length(colnames(biclusterlist2[[as.integer(input$bic)]][[5]]))]
            dat <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
            km_fit <- survfit(Surv(OS.time,OS)~bicluster, data=dat)
            autoplot(km_fit) + 
                labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
                     title = "Overall Survival Times Of \n Bladder Cancer Patients \n") + 
                theme(plot.title = element_text(hjust = 0.5), 
                      axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                      axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                      legend.title = element_text(face="bold", size = 10))})
    }
        else if (input$variable=="Disease-specific survival"){
        # output$survivalvis <- renderPlot({
        #     samples <- colnames(biclusterlist2[[as.integer(input$bic)]][[5]])[2:length(colnames(biclusterlist2[[as.integer(input$bic)]][[5]]))]
        #     survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0)) %>% ggplot(aes(x=DSS.time,fill=as.factor(bicluster)))+geom_boxplot(alpha=0.5)+
        #         scale_fill_manual(values = c("#78B7C5","#F2AD00"),name = paste0("Bicluster ",as_string(input$bic)), labels = c("No","Yes"),guide = guide_legend(reverse=TRUE)) +
        #         theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "light grey"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "light grey"),plot.title = element_text(size=10))})
        output$survivalvis <- renderPlot({
            samples <- colnames(biclusterlist2[[as.integer(input$bic)]][[5]])[2:length(colnames(biclusterlist2[[as.integer(input$bic)]][[5]]))]
            dat <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
            km_fit <- survfit(Surv(DSS.time,DSS)~bicluster, data=dat)
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
            samples <- colnames(biclusterlist2[[as.integer(input$bic)]][[5]])[2:length(colnames(biclusterlist2[[as.integer(input$bic)]][[5]]))]
            survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0)) %>% ggplot(aes(x=DFI.time,fill=as.factor(bicluster)))+geom_boxplot(alpha=0.5)+
                scale_fill_manual(values = c("#78B7C5","#F2AD00"),name = paste0("Bicluster ",as_string(input$bic)), labels = c("No","Yes"),guide = guide_legend(reverse=TRUE)) +
                theme(panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "light grey"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "light grey"),plot.title = element_text(size=10))})
        output$survivalvis <- renderPlot({
            samples <- colnames(biclusterlist2[[as.integer(input$bic)]][[5]])[2:length(colnames(biclusterlist2[[as.integer(input$bic)]][[5]]))]
            dat <- survival %>% mutate(bicluster = ifelse(sample %in% samples,1,0))
            km_fit <- survfit(Surv(DFI.time,DFI)~bicluster, data=dat)
            autoplot(km_fit) + 
                labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
                     title = "Disease-free Interval Times Of \n Bladder Cancer Patients \n") + 
                theme(plot.title = element_text(hjust = 0.5), 
                      axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                      axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                      legend.title = element_text(face="bold", size = 10))})
    }
    )
    output$survivaltable <- renderTable({
        biclusterlist2[[as.integer(input$bic)]][[1]] %>% select(Gene.ID,chrom)
    })
    
    observeEvent(input$withcopy,
                 if (input$withcopy=="No"){
                     output$mapvis <- renderPlotly({
                         dat <- biclusterlist2[[as.integer(input$num)]][[2]]
                         ggplotly(dat %>% ggplot(aes(x=chrom,y=prop,fill=col))+geom_bar(stat='identity')+
                             labs(title = "") +
                             theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.background = element_rect(fill = "white",colour = "white",size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "light grey"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "light grey"),legend.position = "none",plot.title = element_text(size=10)) + 
                             coord_cartesian(ylim=c(0,NA)) + scale_fill_manual(values = levels(as.factor(dat$col)))+
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
    output$info <- renderText({
        dat <- biclusterlist2[[as.integer(input$num)]][[1]]
        print(paste0("This bicluster contains: \n ",toString(nrow(dat))," genes, \n ",toString(dat$Samples.In.Bicluster[1])," samples."))
    })
    output$table <- renderTable({
        dat <- biclusterlist2[[as.integer(input$num)]][[1]]
        dat %>% select(Gene.ID,chrom)
    })
}


shinyApp(ui = ui, server = server)
