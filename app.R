#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

##The code in this script is adapted from the skeleton application provided in class

library(shiny)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(pheatmap)
library(shinyjs)
library(EnhancedVolcano)

#load the microarray expression data
load("expression.Rdata")
#load the differential expression data
load("differential_expression.Rdata")



#turn off the RStudio graphics
graphics.off()

# Define UI for application that draws the heatmap and volcano plot
ui <- navbarPage("Assignment Shiny App",
                 
                 # Create tab for heatmap
                 tabPanel("Heatmap",
                          
                          # Heatmap application title
                          titlePanel("Microarray sample expression heatmap"),
                          
                          # Sidebar with a slider input, checkboxes and radio buttons for various properties of the plot 
                          sidebarLayout(
                              sidebarPanel("Customisation options",
                                  sliderInput("font_row",
                                              "Font size row:",
                                              min = 6,
                                              max = 14,
                                              value = 10),
                                  sliderInput("font_col",
                                              "Font size col:",
                                              min = 6,
                                              max = 14,
                                              value = 10),
                                  sliderInput("cell_height",
                                              "Cell Height:",
                                              min = 4,
                                              max = 11,
                                              value = 90),
                                  sliderInput("cell_width",
                                              "Cell Width:",
                                              min = 10,
                                              max = 100,
                                              value = 10,
                                              step = 5),
                                  checkboxInput("srownames", "Show Row Names", FALSE),
                                  checkboxInput("logtansform", "Log transform values", FALSE),
                                  radioButtons("norm", "Scale by", choices=c("none","row","column"))
                                  
                              ),
                              
                              # Show heatmap
                              mainPanel(
                                  plotOutput("distPlot")
                              )
                          )
                 ),
                 
                 # Create new tab for volcano plot
                 tabPanel("Volcano plot",
    
    # Application title
    titlePanel("Volcano plot of differential expression in DNMT-Tg mice"),
    
    # Sidebar with a slider input and checkboxes for properties of the plot 
    sidebarLayout(
        sidebarPanel("Customisation options",
            sliderInput("font_axis",
                        "Font size axis labels:",
                        min = 6,
                        max = 20,
                        value = 10,
                        step = 1),
            sliderInput("pval_cutoff",
                        "-log10 P-value cutoff:",
                        min = 0,
                        max = 10,
                        value = 2,
                        step = 0.1),
            sliderInput("FC_cutoff",
                        "Log2 fold change cutoff:",
                        min = 0,
                        max = 7,
                        value = 2,
                        step = 0.5),
            sliderInput("gene_label",
                        "Gene label size:",
                        min = 0,
                        max = 5,
                        value = 4,
                        step = 1),
            sliderInput("point_size",
                        "Gene point size:",
                        min = 0,
                        max = 5,
                        value = 1,
                        step = 1),
            
            
            checkboxInput("sgenes", "Show Genes", FALSE),
            checkboxInput("line_connectors", "Show line connectors", FALSE)
            
        ),
        
        # Show volcano plot
        mainPanel(
            plotOutput("volPlot", height = "600px", width = "800px")
        )
    )
)
)

# Define server logic required to draw plots
server <- function(input, output,session) {
    
    output$distPlot <- renderPlot({
        
        #if the log transform checkbox is selected, log transform the expression value
        if(input$logtansform){
            expression <- log2(expression + 1)
        }
        if(is.null(input$select)){
            mysel<-NULL
        }else if(input$select[1]=="none"){
            mysel<-NULL
        }else if(length(input$select)==1){
            #if the data frame has one column it converts to a factor
            #force the type to be a data frame and restore row and column names
            mysel <-as.data.frame(experiment[,input$select[1]])
            rownames(mysel) <-rownames(experiment)
            colnames(mysel) <-input$select[1]
        }else{
            mysel<-experiment[,input$select]
        }
        #plot heatmap with chosen properties
        pheatmap(expression,
                 fontsize_row = input$font_row,
                 fontsize_col=input$font_col,
                 cellheight=input$cell_height,
                 cellwidth=input$cell_width,
                 show_rownames=input$srownames,
                 scale=input$norm,
                 annotation_col=mysel)
    }, execOnResize = F,height = 700)
    
    observeEvent(input$refresh, {
        session$invalidate
    })
    
    output$volPlot <- renderPlot({
        #if the show genes checkbox is selected, then the gene symbols are assigned as labels
        if(input$sgenes){
             gnames <- myresults$Symbol
        }else{gnames <- NA
        }
        #plot volcano plot
EnhancedVolcano(myresults, lab = gnames, x = 'logFC', y= 'adj.P.Val', pCutoff=10^-(input$pval_cutoff), FCcutoff=input$FC_cutoff, axisLabSize = input$font_axis, labSize = input$gene_label, pointSize = input$point_size, drawConnectors= input$line_connectors, border = "full", xlim = c(-9,9), ylim = c(0,8))
    }, execOnResize = F )
    
    observeEvent(input$refresh, {
        session$invalidate
    })
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
