library(BiocManager)
options(repos = BiocManager::repositories())
library(shiny)
library(shinydashboard)
library(knitr)
library(markdown)
library(rmarkdown)
library(goSorensen)
library(ggplot2)
source("info/simfuncsPar.R")
library("org.Ag.eg.db")
library("org.At.tair.db")
library("org.Bt.eg.db")
library("org.Ce.eg.db")
library("org.Cf.eg.db")
library("org.Dm.eg.db")
library("org.Dr.eg.db")
library("org.EcK12.eg.db")
library("org.EcSakai.eg.db")
library("org.Gg.eg.db")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("org.Mmu.eg.db")
library("org.Mxanthus.db")
library("org.Pt.eg.db")
library("org.Rn.eg.db")
library("org.Sc.sgd.db")
library("org.Ss.eg.db")
library("org.Xl.eg.db")

dashboardPage(skin = "blue", title = "goSorensenApp",
  dashboardHeader(title = div('goSorensen: App to detect equivalence between feature lists'), 
                            titleWidth = "1000"),
  
  dashboardSidebar(sidebarMenu(
    menuItem(div(hr(), img(src = "goSorensen.png", height = "100px", style="text-align: center;"))),
    menuItem('How to use this app',
             tabName = 'how',
             icon = icon("info")),
    menuItem("Use example data", 
             icon = icon("table"),
             menuSubItem("From a contingecy table",
                         tabName = "excont"),
             menuSubItem("From feature lists",
                         tabName = "exfeat")),
    menuItem('Enter your data',
             icon = icon("dna"),
             menuSubItem("Enter a contingecy table",
                            tabName = "contin"),
             menuSubItem("Enter feature lists",
                            tabName = "feature"))
  )),
  
  dashboardBody(
    tabItems(
      #tabItem(tabName = "how", 
      #        includeHTML(render("info/userinfo.html"))
      #),
      tabItem(tabName = "how",
              fluidRow(box(htmlOutput("inc"), width =12))),
      tabItem(tabName = "excont",
        fluidPage(tabsetPanel(
          tabPanel("Contingency table", 
                   sidebarLayout(
                     sidebarPanel(
                       wellPanel(
                         p("n11 = 56"),
                         p("n10 = 30"),
                         p("n01 = 1"),
                         p("n00 = 471")
                       )
                     ),
                     mainPanel(
                       verbatimTextOutput("ctable1")
                     )
                   )),
          tabPanel("Equivalence test",
                   sidebarLayout(
                     sidebarPanel( 
                                  titlePanel(h3("Enter your parameters")),
                       wellPanel(
                         numericInput(inputId = "d0", 
                                      label = "Irrelevance Limit d0", 
                                      value = 0.2857, min = 0, step = 0.0001),
                         numericInput(inputId = "conf", 
                                      label = "Confidence Level", 
                                      value = 0.95, min = 0.5, max = 1, step = 0.01),
                         radioButtons(inputId = "samp", 
                                      label = "Sample distribution", 
                                      choices = c("Normal" = F, "Bootstrap" = T), 
                                      selected = F),
                         numericInput(inputId = "nSim", 
                                      label = "Number of simulations", 
                                      value = 10000, min = 1000, max = 1000000, 
                                      step = 1)
                       )          
                     ),
                     mainPanel(
                         verbatimTextOutput(outputId = "equiv"),
                         plotOutput(outputId = "graph", width = "400px", height = "250px")
                     )
                   ))
        ))
      ),
      tabItem(tabName = "exfeat",
              fluidPage(tabsetPanel(
                tabPanel("Select Feature Lists",
                         h5("This list of genes belongs to human organisms and Genome wide annotation for Human (org.Hs.eg.db) is used"),
                         splitLayout(
                           wellPanel(
                             selectInput(inputId = "gl1", 
                                         label = "Select Gene List 1", 
                                         choices = c("atlas", "sanger", "cangenes", "cis", "miscellaneous", "Vogelstein", "waldman"),
                                         selected = "atlas"),
                             verbatimTextOutput("list1")
                           ),
                           wellPanel(
                             selectInput(inputId = "gl2", 
                                         label = "Select Gene List 2", 
                                         choices = c("atlas", "sanger", "cangenes", "cis", "miscellaneous", "Vogelstein", "waldman"),
                                         selected = "sanger"),
                             verbatimTextOutput("list2")
                           )
                         )),
                tabPanel("Contingency table",
                         sidebarLayout(
                           sidebarPanel(
                             radioButtons("onto", "Select Ontology", 
                                          c("BP", "MF", "CC")),
                             radioButtons("golevel", "Select GO level", 
                                          3:10),
                             actionButton("go", "Compute")
                           ),
                           mainPanel(
                             verbatimTextOutput("ctonco")
                           ))),
                tabPanel("Equivalence test",
                         sidebarLayout(
                           sidebarPanel(
                             titlePanel(h3("Enter your parameters")),
                             wellPanel(
                               numericInput(inputId = "d0_1", 
                                            label = "Irrelevance Limit d0", 
                                            value = 0.2857, min = 0, step = 0.0001),
                               numericInput(inputId = "conf_1", 
                                            label = "Confidence Level", 
                                            value = 0.95, min = 0.5, max = 1, step = 0.01),
                               radioButtons(inputId = "samp_1", 
                                            label = "Sample distribution", 
                                            choices = c("Normal" = F, "Bootstrap" = T), 
                                            selected = F),
                               numericInput(inputId = "nSim_1", 
                                            label = "Number of simulations", 
                                            value = 10000, min = 1000, max = 1000000, 
                                            step = 1),
                               actionButton("go_1", "Compute")
                             )
                           ),
                           mainPanel(
                             verbatimTextOutput(outputId = "equiv_1"),
                             plotOutput(outputId = "graph_1", width = "400px", height = "250px")
                           )
                         ))
              ))
              
      ),
      tabItem(tabName = "contin",
              tabsetPanel(
                tabPanel("Enter contingency table",
                         sidebarLayout(
                           sidebarPanel(
                             wellPanel(
                             numericInput("n11", "n11", NA, 0, step = 1, 
                                          width = "100px"),
                             numericInput("n10", "n10", NA, 0, step = 1, 
                                          width = "100px"),
                             numericInput("n01", "n01", NA, 0, step = 1, 
                                          width = "100px"),
                             numericInput("n00", "n00", NA, 0, step = 1, 
                                          width = "100px")),
                             wellPanel(
                               textInput("ngl1", "Name - Gene List 1", 
                                         placeholder = "Enter name of gene list 1"),
                               textInput("ngl2", "Name - Gene List 2", 
                                         placeholder = "Enter name of gene list 2")
                             )
                           ),
                           mainPanel(
                             verbatimTextOutput("ctcomp")
                           )
                         )),
                tabPanel("Results for equivalence test",
                         sidebarLayout(
                           sidebarPanel(
                             titlePanel(h3("Enter your parameters")),
                             wellPanel(
                               numericInput(inputId = "d0_2", 
                                            label = "Irrelevance Limit d0", 
                                            value = 0.2857, min = 0, step = 0.0001),
                               numericInput(inputId = "conf_2", 
                                            label = "Confidence Level", 
                                            value = 0.95, min = 0.5, max = 1, step = 0.01),
                               radioButtons(inputId = "samp_2", 
                                            label = "Sample distribution", 
                                            choices = c("Normal" = F, "Bootstrap" = T), 
                                            selected = F),
                               numericInput(inputId = "nSim_2", 
                                            label = "Number of simulations", 
                                            value = 10000, min = 1000, max = 1000000, 
                                            step = 1)
                             )
                           ),
                           mainPanel(
                             verbatimTextOutput(outputId = "equiv_2"),
                             plotOutput(outputId = "graph_2", width = "400px", height = "250px")
                           )
                         ))
              )
      ),
      tabItem(tabName = "feature",
              selectInput("genome", "Using Genome Wide Annotation for:", 
                          c("Anopheles" = "org.Ag.eg.db",
                            "Arabidopsis" = "org.At.tair.db",
                            "Bovine" = "org.Bt.eg.db",
                            "Worm" = "org.Ce.eg.db",
                            "Canine" = "org.Cf.eg.db",
                            "Fly" = "org.Dm.eg.db",
                            "Zebrafish" = "org.Dr.eg.db",
                            "E coli strain K12" = "org.EcK12.eg.db",
                            "E coli strain Sakai" = "org.EcSakai.eg.db",
                            "Chicken" = "org.Gg.eg.db",
                            "Human" = "org.Hs.eg.db",
                            "Mouse" = "org.Mm.eg.db",
                            "Rhesus" = "org.Mmu.eg.db",
                            "Myxococcus xanthus DK 1622" = "org.Mxanthus.db",
                            "Chimp" = "org.Pt.eg.db",
                            "Rat" = "org.Rn.eg.db",
                            "Yeast" = "org.Sc.sgd.db",
                            "Pig" = "org.Ss.eg.db",
                            "Xenopus" = "org.Xl.eg.db"), 
                          selected = "org.Hs.eg.db"),
              tabsetPanel(
                tabPanel("Enter your feature lists",
                         splitLayout(
                           wellPanel(
                             textInput("imname1", "Enter name Feature List 1"),
                             fileInput("file1", "Upload Feature List 1 (txt File)"),
                             verbatimTextOutput("fl1")
                           ),
                           wellPanel(
                             textInput("imname2", "Enter name Feature List 2"),
                             fileInput("file2", "Upload Feature List 2 (txt File)"),
                             verbatimTextOutput("fl2")
                           )
                         )),
                tabPanel("Contingency table",
                         sidebarLayout(
                           sidebarPanel(
                             radioButtons("imponto", "Select Ontology", 
                                          c("BP", "MF", "CC")),
                             radioButtons("impgolevel", "Select GO level", 
                                          3:10),
                             actionButton("goimp", "Compute")
                           ),
                           mainPanel(
                             verbatimTextOutput("impctonco")
                           ))),
                tabPanel("Equivalence test", 
                         sidebarLayout(
                           sidebarPanel(
                             titlePanel(h3("Enter your parameters")),
                             wellPanel(
                               numericInput(inputId = "d0_3", 
                                            label = "Irrelevance Limit d0", 
                                            value = 0.2857, min = 0, step = 0.0001),
                               numericInput(inputId = "conf_3", 
                                            label = "Confidence Level", 
                                            value = 0.95, min = 0.5, max = 1, step = 0.01),
                               radioButtons(inputId = "samp_3", 
                                            label = "Sample distribution", 
                                            choices = c("Normal" = F, "Bootstrap" = T), 
                                            selected = F),
                               numericInput(inputId = "nSim_3", 
                                            label = "Number of simulations", 
                                            value = 10000, min = 1000, max = 1000000, 
                                            step = 1),
                               actionButton("go3", "Compute")
                             )
                           ),
                           mainPanel(
                             verbatimTextOutput(outputId = "equiv_3"),
                             plotOutput(outputId = "graph_3", width = "400px", height = "250px")
                           )
                         ))
              )
              
      )
    )
  ) 
)

