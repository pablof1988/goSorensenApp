server <- function(input, output){
  output$inc <- renderUI({
    includeHTML("info/userinfo.html")
  })

  Gene.List1_Vs_Gene.List2_BP_3 <- reactive({
    tab <- as.table(matrix(c(56, 1, 30, 471), 2, 2, 
                           dimnames = list(c(T, F), c(T, F))))
    
    names(dimnames(tab)) <- c("Enriched in Gene List 1", "Enriched in Gene List 2")
    tab
  })
  
  oncoLists <- reactive({
    data("allOncoGeneLists")
    allOncoGeneLists
  })
  
  output$ctable1 <- renderPrint({
    print("CONTINGENCY TABLE FOR GENE LIST 1 Vs GENE LIST 2")
    Gene.List1_Vs_Gene.List2_BP_3()
  })
  
  output$equiv <- renderPrint({
    equivTestSorensen( Gene.List1_Vs_Gene.List2_BP_3(), d0 = input$d0,
                      conf.level = input$conf, boot = input$samp, 
                      nboot = input$nSim)
  })
  
  output$graph <- renderPlot({
    nsims <- input$nSim
    n <- sum( Gene.List1_Vs_Gene.List2_BP_3())
    p01 <-  Gene.List1_Vs_Gene.List2_BP_3()[2, 1] / n
    p10 <-  Gene.List1_Vs_Gene.List2_BP_3()[1, 2] / n
    d0 <- input$d0
    p11 <- dSorensen2p11(d0, p01, p10) 
    
    d <- dSorensen(c(p11, p01, p10, 1 - (p11 + p01 + p10)))
    se <- seSorensen(c(p11, p01, p10, 1 - (p11 + p01 + p10)) * n)
    
    set.seed(1234)
    tabs <- rmultinom(nsims, n, c(p11, p01, p10, 1 - (p11 + p01 + p10)))
    
    stat <- function(tab, d) {
      return((dSorensen(tab) - d) / seSorensen(tab))
    }
    
    statVals <- apply(tabs, 2, function(tab, d) {
      return(stat(as.numeric(tab), d))
    }, d = d)
    
    
    iRan <- sample.int(length(statVals), 1)
    oneRandomTab <- tabs[,iRan]
    boot.tabs <- rmultinom(nsims, n, oneRandomTab / n)
    
    boot.statVals <- apply(boot.tabs, 2, function(tab, d) {
      return(stat(as.numeric(tab), d))
    }, d = dSorensen(oneRandomTab))
    
    Values <- c(statVals, boot.statVals)
    Density <- c(rep("True distribution", length(statVals)), 
                 rep("Bootstrap", length(boot.statVals)))
    
    data <- data.frame(Values, Density)
    
    ggplot(data, aes(x = Values, color = Density)) +
      geom_density(linewidth = 0.75) +
      stat_function(fun = function(x) dnorm(x, mean = 0, sd = 1),
                    linetype = "longdash", linewidth = 0.75,
                    aes(colour = 'Normal')) +
      scale_colour_manual("", values = c('deepskyblue1', 'deeppink', 'black')) +
      guides(color = guide_legend(
        override.aes = list(
          #fill  = c("white", "white", 'white'), 
          linetype = c("solid", 'longdash', "solid")))) +
      theme(legend.position = c(0.15, 0.875), 
            legend.background = element_blank()) +
      labs(x = "", y = "")
  })
  
  output$list1 <- renderPrint({
    oncoLists()[[input$gl1]]
  })
  
  output$list2 <- renderPrint({
    oncoLists()[[input$gl2]]
  })
  
  allOncoCT <- reactive({
    buildEnrichTable(oncoLists()[[input$gl1]], oncoLists()[[input$gl2]], 
                     onto = input$onto, GOLevel = as.numeric(input$golevel), 
                     orgPackg = "org.Hs.eg.db", 
                     listNames = c(input$gl1, input$gl2))
  })
  
  output$ctonco <- renderPrint({
    input$go
    isolate(allOncoCT())
  })
  
  output$equiv_1 <- renderPrint({
    input$go_1
    isolate(res1 <- equivTestSorensen(allOncoCT(), d0 = input$d0_1,
                       conf.level = input$conf_1, boot = input$samp_1, 
                       nboot = input$nSim_1))
    res1$data.name <- paste(input$gl1, "Vs", input$gl2, input$onto, input$golevel)
    res1
  })
  
  output$graph_1 <- renderPlot({
    input$go_1
    nsims <- input$nSim_1
    n <- sum(allOncoCT())
    p01 <-  allOncoCT()[2, 1] / n
    p10 <-  allOncoCT()[1, 2] / n
    d0 <- input$d0_1
    p11 <- dSorensen2p11(d0, p01, p10) 
    
    d <- dSorensen(c(p11, p01, p10, 1 - (p11 + p01 + p10)))
    se <- seSorensen(c(p11, p01, p10, 1 - (p11 + p01 + p10)) * n)
    
    set.seed(1234)
    tabs <- rmultinom(nsims, n, c(p11, p01, p10, 1 - (p11 + p01 + p10)))
    
    stat <- function(tab, d) {
      return((dSorensen(tab) - d) / seSorensen(tab))
    }
    
    statVals <- apply(tabs, 2, function(tab, d) {
      return(stat(as.numeric(tab), d))
    }, d = d)
    
    
    iRan <- sample.int(length(statVals), 1)
    oneRandomTab <- tabs[,iRan]
    boot.tabs <- rmultinom(nsims, n, oneRandomTab / n)
    
    boot.statVals <- apply(boot.tabs, 2, function(tab, d) {
      return(stat(as.numeric(tab), d))
    }, d = dSorensen(oneRandomTab))
    
    Values <- c(statVals, boot.statVals)
    Density <- c(rep("True distribution", length(statVals)), 
                 rep("Bootstrap", length(boot.statVals)))
    
    data <- data.frame(Values, Density)
    
    rg1 <- ggplot(data, aes(x = Values, color = Density)) +
      geom_density(linewidth = 0.75) +
      stat_function(fun = function(x) dnorm(x, mean = 0, sd = 1),
                    linetype = "longdash", linewidth = 0.75,
                    aes(colour = 'Normal')) +
      scale_colour_manual("", values = c('deepskyblue1', 'deeppink', 'black')) +
      guides(color = guide_legend(
        override.aes = list(
          #fill  = c("white", "white", 'white'), 
          linetype = c("solid", 'longdash', "solid")))) +
      theme(legend.position = c(0.15, 0.875), 
            legend.background = element_blank()) +
      labs(x = "", y = "")
    isolate(rg1)
  })
  
  com_contable <- reactive({
    tab <- as.table(matrix(c(input$n11, input$n01, input$n10, input$n00), 2, 2, 
                           dimnames = list(c(T, F), c(T, F))))
    
    names(dimnames(tab)) <- paste("Enrich in", c(input$ngl1, input$ngl2))
    tab
  })
  
  output$ctcomp <- renderPrint({
    print(paste("CONTINGENCY TABLE FOR", input$ngl1, "VS", input$ngl2))
    com_contable()
  })
  
  output$equiv_2 <- renderPrint({
    resequiv_2 <- equivTestSorensen(com_contable(), d0 = input$d0_2,
                       conf.level = input$conf_2, boot = input$samp_2, 
                       nboot = input$nSim_2)
    resequiv_2$data.name <- paste(input$ngl1, "Vs", input$ngl2)
    resequiv_2
    
  })

  
  output$graph_2 <- renderPlot({
    nsims <- input$nSim_2
    n <- sum(com_contable())
    p01 <-  com_contable()[2, 1] / n
    p10 <-  com_contable()[1, 2] / n
    d0 <- input$d0_2
    p11 <- dSorensen2p11(d0, p01, p10) 
    
    d <- dSorensen(c(p11, p01, p10, 1 - (p11 + p01 + p10)))
    se <- seSorensen(c(p11, p01, p10, 1 - (p11 + p01 + p10)) * n)
    
    set.seed(1234)
    tabs <- rmultinom(nsims, n, c(p11, p01, p10, 1 - (p11 + p01 + p10)))
    
    stat <- function(tab, d) {
      return((dSorensen(tab) - d) / seSorensen(tab))
    }
    
    statVals <- apply(tabs, 2, function(tab, d) {
      return(stat(as.numeric(tab), d))
    }, d = d)
    
    
    iRan <- sample.int(length(statVals), 1)
    oneRandomTab <- tabs[,iRan]
    boot.tabs <- rmultinom(nsims, n, oneRandomTab / n)
    
    boot.statVals <- apply(boot.tabs, 2, function(tab, d) {
      return(stat(as.numeric(tab), d))
    }, d = dSorensen(oneRandomTab))
    
    Values <- c(statVals, boot.statVals)
    Density <- c(rep("True distribution", length(statVals)), 
                 rep("Bootstrap", length(boot.statVals)))
    
    data <- data.frame(Values, Density)
    
    ggplot(data, aes(x = Values, color = Density)) +
      geom_density(linewidth = 0.75) +
      stat_function(fun = function(x) dnorm(x, mean = 0, sd = 1),
                    linetype = "longdash", linewidth = 0.75,
                    aes(colour = 'Normal')) +
      scale_colour_manual("", values = c('deepskyblue1', 'deeppink', 'black')) +
      guides(color = guide_legend(
        override.aes = list(
          #fill  = c("white", "white", 'white'), 
          linetype = c("solid", 'longdash', "solid")))) +
      theme(legend.position = c(0.15, 0.875), 
            legend.background = element_blank()) +
      labs(x = "", y = "")
  })
  
  imp_list1 <- reactive({
    dir <- input$file1
    dat <- read.table(dir$datapath)
    as.character(as.vector(unlist(dat)))
  })
  
  imp_list2 <- reactive({
    dir <- input$file2
    dat <- read.table(dir$datapath)
    as.character(as.vector(unlist(dat)))
  })
  
  output$fl1 <- renderPrint({
    imp_list1()
  })
  
  output$fl2 <- renderPrint({
    imp_list2()
  })
  
  impallOncoCT <- reactive({
    buildEnrichTable(imp_list1(), imp_list2(), 
                     onto = input$imponto, GOLevel = as.numeric(input$impgolevel), 
                     orgPackg = input$genome, 
                     listNames = c(input$imname1, input$imname2))
  })
  
  output$impctonco <- renderPrint({
    input$goimp
    isolate(print(paste("CONTINGENCY TABLE FOR", input$imname1, "VS", input$imname2, 
                "ONTO:", input$imponto, "GO_LEVEL:", input$impgolevel)))
    isolate(impallOncoCT())
  })
  
  output$equiv_3 <- renderPrint({
    input$go3
    isolate(resequiv_3 <- equivTestSorensen(impallOncoCT(), d0 = input$d0_3,
                                    conf.level = input$conf_3, boot = input$samp_3, 
                                    nboot = input$nSim_3))
    resequiv_3$data.name <- paste0(input$imname1, "_Vs_", input$imname2, 
                                   "_GO:", input$imponto, "_LEVEL", input$impgolevel)
    resequiv_3
  })
  
  
  output$graph_3 <- renderPlot({
    input$go3
    nsims <- input$nSim_3
    n <- sum(impallOncoCT())
    p01 <-  impallOncoCT()[2, 1] / n
    p10 <-  impallOncoCT()[1, 2] / n
    d0 <- input$d0_3
    p11 <- dSorensen2p11(d0, p01, p10) 
    
    d <- dSorensen(c(p11, p01, p10, 1 - (p11 + p01 + p10)))
    se <- seSorensen(c(p11, p01, p10, 1 - (p11 + p01 + p10)) * n)
    
    set.seed(1234)
    tabs <- rmultinom(nsims, n, c(p11, p01, p10, 1 - (p11 + p01 + p10)))
    
    stat <- function(tab, d) {
      return((dSorensen(tab) - d) / seSorensen(tab))
    }
    
    statVals <- apply(tabs, 2, function(tab, d) {
      return(stat(as.numeric(tab), d))
    }, d = d)
    
    
    iRan <- sample.int(length(statVals), 1)
    oneRandomTab <- tabs[,iRan]
    boot.tabs <- rmultinom(nsims, n, oneRandomTab / n)
    
    boot.statVals <- apply(boot.tabs, 2, function(tab, d) {
      return(stat(as.numeric(tab), d))
    }, d = dSorensen(oneRandomTab))
    
    Values <- c(statVals, boot.statVals)
    Density <- c(rep("True distribution", length(statVals)), 
                 rep("Bootstrap", length(boot.statVals)))
    
    data <- data.frame(Values, Density)
    
    isolate(ggplot(data, aes(x = Values, color = Density)) +
      geom_density(linewidth = 0.75) +
      stat_function(fun = function(x) dnorm(x, mean = 0, sd = 1),
                    linetype = "longdash", linewidth = 0.75,
                    aes(colour = 'Normal')) +
      scale_colour_manual("", values = c('deepskyblue1', 'deeppink', 'black')) +
      guides(color = guide_legend(
        override.aes = list(
          #fill  = c("white", "white", 'white'), 
          linetype = c("solid", 'longdash', "solid")))) +
      theme(legend.position = c(0.15, 0.875), 
            legend.background = element_blank()) +
      labs(x = "", y = ""))
  })
}