library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(colourpicker)
library(shinythemes)
library(thematic)
library(cli)
library('RColorBrewer')
library(shinycssloaders) 
library(scales)
library(stringr)
library(shinyWidgets)
library(tidyr)
library(stats)
library(forcats)
library('ggbeeswarm')
library(DESeq2)
library(DT)
library(tibble)
library(readr)
library(gplots)
library(shinyjs)

css<-'#gc+div div a {color: black;}'
ui<-fluidPage(
  #bslib package allows me to use a theme from Bootswatch v5 and allows for further customization
  #I based my app off of the 'lux' style available on bootswatch.com
  theme=bslib::bs_theme(bootswatch = "lux",bg='black',fg='white',primary='white',secondary='white'),
  #I define some inline css stylings for various aspects of my application, mainly the slider bars and navbar
  tags$head(tags$style(css)),
  chooseSliderSkin('Square'),
  tags$style(".irs--square .irs-line {background-color: #333}
             .irs--square .irs-bar {background-color: white}
             .irs--square .irs-handle {border: 3px solid white; background-color:black; transform: rotate(90deg);}
             .progress-bar {background-color:#4bbf73}
             .thead, .th {color:white}
             .multi-wrapper .non-selected-wrapper, .multi-wrapper .selected-wrapper{width:100%;height:75px}
             .navbar {font-size: 1.5rem;}
             .navbar-header {padding-top:1rem; padding-bottom:0rem}
             .navbar-nav {display: flex !important;justify-content: center !important;width: 100%;}"),
  #I use a navbar instead of nested tabs, since it looks a bit cleaner
  navbarPage(div(h3('Final Project'),h5('Jackson Faulx',style='color:#999; text-align:center; font-size: 16px')),
    tabPanel('Samples',
      sidebarLayout(
        sidebarPanel(
          #file input only accepts csv and tsv files, wont let you submit other types of files
          fileInput('file1',p('LOAD COUNTS MATRIX',style='color:white; text-align:center'),accept=c('.csv','.tsv')),
          #these conditional panels will only output the UI inside if a certain tab is selected
          conditionalPanel(condition="input.mt=='SUMMARY'",
            h2('Summary Table',style='color:white; text-align:center'),
            p('Provides a breakdown of each column in the metadata specifying the column name, data type, and distinct values')
                           ),
          conditionalPanel(condition="input.mt=='SAMPLE INFORMATION'",
            h2('Sample Information Table',style='color:white; text-align:center'),
            p('A sortable table containing all of the metadata provided for each sample')
                           )),
        mainPanel(
          tabsetPanel(id='mt',
            tabPanel('SUMMARY',
              tableOutput('meta'),
                    ),
            tabPanel('SAMPLE INFORMATION',
              dataTableOutput('si') 
                        )
                      )
                    )
                  )
                ),
     tabPanel('Counts',
      sidebarLayout(
        sidebarPanel(
          #File inputs will be on every analysis page, but the data also carries throughout the app
          fileInput('file2',p('LOAD COUNTS MATRIX',style='color:white; text-align:center'),accept=c('.csv','.tsv')),
          #more conditional description panels
          conditionalPanel(condition="input.ct=='FILTER INFO'",
            h2('Filter Statistics',style='color:white; text-align:center'),
            p('A breakdown of the total number of genes and samples in the counts matrix, 
                             as well as the number and percentage of genes that pass the two thresholds set by the sliders below')
                          ),
          conditionalPanel(condition="input.ct=='FILTER PLOTS'",
            h2('Diagnostic Scatter Plots',style='color:white; text-align:center'),
            h5('Median Count vs Variance',style='color:white; text-align:center'),
            p('A plot showing median count vs log10 transformed variance for each gene, with genes passing the thresholds
                             set below highlighted in green'),
            h5('Median Count vs Number of Zeros',style='color:white; text-align:center'),
            p('A plot showing median count vs number of samples containing zero reads for each gene, 
                             with genes passing the thresholds set below highlighted in green')
                          ),
          conditionalPanel(condition="input.ct=='HEATMAP'",
            h2('Clustered Heatmap',style='color:white; text-align:center'),
            p('A clustered heatmap depicting log10-transformed and normalized read counts of genes that pass the threshold set below')
                          ),
          conditionalPanel(condition="input.ct=='PCA'",
            h2('Principal Component Analysis',style='color:white; text-align:center'),
            h5('Scatter Plot',style='color:white; text-align:center'),
            p('A scatter plot showing variance contribution of the two chosen principal components for genes passing 
                             the thresholds set below'),
            h5('Beeswarm Plot',style='color:white; text-align:center'),
            p('A beeswarm plot of the top N principal components and their contribution to overall variance for genes passing
                             the thresholds set below')
                          ),
          #slider bars to set filter thresholds
          sliderInput('vslide',min=0,max=100,
                      label=p('MINIMUM PERCENTILE OF VARIANCE',style='color:white; text-align:center'),
                      value=50,step=1),
          sliderInput('nzslide',min=0,max=9,label=p('MINIMUM NON-ZERO SAMPLES',style='color:white; text-align:center'),
                      value=5,step=1)),
        mainPanel(
          tabsetPanel(id='ct',
            tabPanel('FILTER INFO',
               #I use column()s throughout the UI to center and format different elements of my UI
               column(align='center', width=12,
                 #I use br() in order to add a line break between each text line
                 h2('Total Number of Samples'),
                 br(),
                 #the UI output for the filter stats is dependent on calculations that are performed on the server side
                 #and will be updated as the slider bars are changed
                 h3(uiOutput('filtsum1',class='text-info')),
                 br(),
                 h2('Total Number of Genes'),
                 br(),
                 h3(uiOutput('filtsum2',class='text-info')),
                 br(),
                 h2('Genes Passing Filter'),
                 br(),
                 h3(uiOutput('filtsum3',class='text-success')),
                 br(),
                 h2('Genes Failing Fitler'),
                 br(),
                 h3(uiOutput('filtsum4',class='text-danger'))
                        )
                       ),
            tabPanel('FILTER PLOTS',
              fluidRow(
                #split layout is used to put two graphs side by side in the UI
                splitLayout(cellWidths = c('50%','50%'),
                  #i use WithSpinner() throughout the app in order to add animated spinners for when plots or other UI aspects are loading          
                  withSpinner(plotOutput('medvarplot'),type=8,color='white'),
                  withSpinner(plotOutput('mednzplot'),type=8,color='white')
                             )
                            )
                           ),
            tabPanel('HEATMAP',
               withSpinner(plotOutput('heatmap',width='90%',height='750px'),type=8,color='white')),
            tabPanel('PCA',
              fluidRow(
                splitLayout(cellWidths = c('50%','50%'),
                  withSpinner(plotOutput('pca'),type=8,color='white'),
                  withSpinner(plotOutput('bees'),type=8,color='white')
                            )
                          ),
              fluidRow(
                #I attempt to place each of graph customizer inputs below their respective plots
                #(emphasis on attempt)
                column(width=3, offset=1,
                  uiOutput('pcx')),
                column(width=3,
                  uiOutput('pcy')
                        ),
                column(width=4, 
                  uiOutput('pcs1')
      
                        )
                      )
                    )
                  )
                )
              )
            ),
     tabPanel('DE',
      sidebarLayout(
        sidebarPanel(
          fileInput('file3',p('LOAD COUNTS MATRIX',style='color:white; text-align:center'),accept=c('.csv','.tsv')),
          conditionalPanel(condition="input.det=='RESULTS'",
            h2('Results Table',style='color:white; text-align:center'),
            p('A sortable and searchable table containing the results of a differential expression analysis performed on the given
               counts matrix, comparing gene expression between adult (Ad) and neonatal (P) mouse myocytes')
                          ),
          conditionalPanel(condition="input.det=='VOLCANO PLOT'",
            h2('Volcano Plot',style='color:white; text-align:center'),
            p('A custom scatter plot of two result statistics, with color customization for genes that pass a specified adjusted
              p-value threshold')
                          ),
                    ),
        mainPanel(
          tabsetPanel(id='det',
            tabPanel('RESULTS',
              withSpinner(DT::dataTableOutput('detab'),type=8,color='white')
                      ),
            tabPanel('VOLCANO PLOT',
              sidebarLayout(
                sidebarPanel(
                  #I use radio buttons with a pre selected option
                  radioButtons('radx',p('CHOOSE STATISTIC FOR THE X-AXIS',style='color:white; text-align:center'),
                              c('baseMean','log2FoldChange','lfcSE','stat','pvalue','padj'),
                                selected='log2FoldChange'),
                  radioButtons('rady',p('CHOOSE STATISTIC FOR THE Y-AXIS',style='color:white; text-align:center'),
                              c('baseMean','log2FoldChange','lfcSE','stat','pvalue','padj'),
                                selected='padj'),
                  #color picker allows customization of the colors on the graph
                  colourpicker::colourInput('color1',p('BASE POINT COLOR',style='color:white; text-align:center'),value='#860038'),
                  colourpicker::colourInput('color2',p('HIGHLIGHT POINT COLOR',style='color:white; text-align:center'),value='#FDBB30'),
                  sliderInput('slider',min=-300,max=0,
                              label=p('SELECT THE -LOG10 P-ADJUSTED THRESHOLD',style='color:white; text-align:center'),
                              value=-150,step=1),
                  #all inputs should wait to process until the action button is activated
                  actionButton('button','Plot',width='100%',class='btn-primary')
                             ),
                mainPanel(
                  withSpinner(plotOutput('volcano'),type=8,color='white')
                            )
                          )     
                        )
                      )
                    )
                  )
                ),
     tabPanel('Gene Expression Visualizer',
      sidebarLayout(
       sidebarPanel(
        fileInput('file4',p('LOAD COUNTS MATRIX',style='color:white; text-align:center'),accept=c('.csv','.tsv')),
        h2('Individual Gene Analysis',style='color:white; text-align:center'),
        p('A custom plot of normalized counts for a chosen gene, split by a chosen metadata variable'),
        #all of the UI here will only show up after a file has been uploaded
        withSpinner(uiOutput('info'),type=8,color='white'),
        uiOutput('genes'),
        uiOutput('plts'),
        actionButton('button2','Plot',width='100%',class='btn-primary')
                    ),
       mainPanel(
        withSpinner(plotOutput('graph'),type=8,color='white')
                    )
                  )
                ) 
              )
            )

server<-function(input,output,session){
  #the options() function increases the file size limit for upload
  options(shiny.maxRequestSize=30*1024^2)
  #thematic_shiny() styles plots in accordance with the selected 'lux' style from bootswatch
  thematic::thematic_shiny()
  
  #DATA PROCESSING - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  #Filter with CPM normalization via DESeq2
  filter1 <- function(verse_counts,variance,nz) {
    verse_counts1<-data.frame(verse_counts)%>%column_to_rownames(var='gene')
    vari<-apply(verse_counts1,1,var)
    v<-verse_counts1 %>% filter(vari>quantile(vari,na.rm=TRUE,probs=(variance/100)),rowSums(verse_counts1>0)>nz)
    v1<-deseq_normalize(v)
    return(v1)
  }
  
  #filter with no normalization (for running DESeq2)
  filternn <- function(verse_counts,variance,nz) {
    verse_counts1<-data.frame(verse_counts)%>%column_to_rownames(var='gene')
    vari<-apply(verse_counts1,1,var)
    v<-verse_counts1 %>% filter(vari>quantile(vari,na.rm=TRUE,probs=(variance/100)),rowSums(verse_counts1>0)>nz)
    return(v)
  }
  
  #normalization function for running DESeq2
  deseq_normalize <- function(count_data) {
    cm<-as.matrix(count_data)
    des<-DESeqDataSetFromMatrix(countData = cm,
                                colData = tibble(sample_name=colnames(count_data)),design=~1)
    des<-estimateSizeFactors(des)
    ndseq<-data.frame(counts(des,normalized=TRUE))
    row.names(ndseq)<-row.names(count_data)
    return(ndseq)
  }
  
  #DESeq2 function
  run_deseq <- function(count_dataframe) {
    meta<-metainfo(count_dataframe)
    meta$timepoint<-factor(meta$timepoint)
    meta$timepoint<-relevel(meta$timepoint, ref = "P0")
    rcounts<-data.frame(count_dataframe) %>% dplyr::select(starts_with(c('vP0', 'vAd')))
    metac<-data.frame(as_tibble(meta) %>% dplyr::select(sample,timepoint)%>%
                        filter(timepoint %in% c('P0', 'Ad')))
    dds <- DESeqDataSetFromMatrix(rcounts,colData=metac,design=~timepoint) 
    dds <- DESeq(dds)
    res <- results(dds)
    return(res)
  }
  
  #REACTIVES - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  #Load data
  load_data <- reactive({
    #This req() fucntion allows the file to be input on any tab and will automatically load in data across the whole app
    #conditional statements for each possible file upload input
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    if (isTruthy(input$file1)){
      f<-read_delim(input$file1$datapath)
      return(f)
    }
    if (isTruthy(input$file2)){
      f<-read_delim(input$file2$datapath)
      return(f)
    }
    if (isTruthy(input$file3)){
      f<-read_delim(input$file3$datapath)
      return(f)
    }
    if (isTruthy(input$file4)){
      f<-read_delim(input$file4$datapath)
      return(f)
    }
  })
  
  #Filtered data
  filtered<-reactive({
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    if (isTruthy(input$file1)){
      f<-read_delim(input$file1$datapath)
      filtered<-filter1(f,input$vslide,input$nzslide)
      return(filtered)
    }
    if (isTruthy(input$file2)){
      f<-read_delim(input$file2$datapath)
      filtered<-filter1(f,input$vslide,input$nzslide)
      return(filtered)
    }
    if (isTruthy(input$file3)){
      f<-read_delim(input$file3$datapath)
      filtered<-filter1(f,input$vslide,input$nzslide)
      return(filtered)
    }
    if (isTruthy(input$file4)){
      f<-read_delim(input$file4$datapath)
      filtered<-filter1(f,input$vslide,input$nzslide)
      return(filtered)
    }
  })
  
  #Filtered data with no normalization (for graphing)
  filteredNN<-reactive({
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    if (isTruthy(input$file1)){
      f<-read_delim(input$file1$datapath)
      filtered<-filternn(f,input$vslide,input$nzslide)
      return(filtered)
    }
    if (isTruthy(input$file2)){
      f<-read_delim(input$file2$datapath)
      filtered<-filternn(f,input$vslide,input$nzslide)
      return(filtered)
    }
    if (isTruthy(input$file3)){
      f<-read_delim(input$file3$datapath)
      filtered<-filternn(f,input$vslide,input$nzslide)
      return(filtered)
    }
    if (isTruthy(input$file4)){
      f<-read_delim(input$file4$datapath)
      filtered<-filternn(f,input$vslide,input$nzslide)
      return(filtered)
    }
  })
  
  #DESeq2 analysis results
  DE<-reactive({
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    if (isTruthy(input$file1)){
      f<-read_delim(input$file1$datapath)
      filtered<-filternn(f,input$vslide,input$nzslide)
      ds<-run_deseq(filtered)
      return(data.frame(ds))
    }
    if (isTruthy(input$file2)){
      f<-read_delim(input$file2$datapath)
      filtered<-filternn(f,input$vslide,input$nzslide)
      ds<-run_deseq(filtered)
      return(data.frame(ds))
    }
    if (isTruthy(input$file3)){
      f<-read_delim(input$file3$datapath)
      filtered<-filternn(f,input$vslide,input$nzslide)
      ds<-run_deseq(filtered)
      return(data.frame(ds))
    }
    if (isTruthy(input$file4)){
      f<-read_delim(input$file4$datapath)
      filtered<-filternn(f,input$vslide,input$nzslide)
      ds<-run_deseq(filtered)
      return(data.frame(ds))
    }
  })
  
  #METADATA TAB - - - - - - - - - - - - -
  
  #function to create metadata table from sample names in counts matrix
  metainfo <- function(counts) {
    sample_names<-colnames(counts)
    data.frame(sample=sample_names)%>%
  mutate(
        timepoint=substr(sample,2,3),
        replicate=substr(sample,5,5)
      )
  }
  
  #outputs the summary of the metadata table created in the metainfo() function
  output$meta<- renderTable({
  data<-load_data()[,-1]
   metadata<-metainfo(data)
   me<-tibble('Column'=colnames(metadata), 'Type'=c(class(metadata$sample),class(metadata$timepoint),class(metadata$replicate)),
          'Distinct Values'=c(paste(unique(metadata$sample),collapse=', '),
                              paste(unique(metadata$timepoint),collapse=', '),
                              paste(unique(metadata$replicate),collapse=', ')))
   return(me)
   #I add this to add some minor styling to the UI table
  },width = "100%",bordered=TRUE,hover=TRUE)
  
  #Output interactive metadata table using DT
  output$si<-DT::renderDataTable({
    data<-load_data()[,-1]
  metadata<-metainfo(data)
  DT::datatable(metadata, class='table-light',style='bootstrap5', rownames = FALSE,
                options = list(paging = FALSE, searching = TRUE,
                               fixedColumns = TRUE, autoWidth = TRUE,
                               ordering = TRUE, dom = 'Bfrtip'))
  })
  
  
  #COUNTS TAB- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  #Output filter pass/fail statistics
  #Number of samples
  output$filtsum1<-renderUI({
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    return(length(load_data()))
    })
  #number of genes
  output$filtsum2<-renderUI({
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    return(nrow(load_data()))
    })
  #number/percent passing
  output$filtsum3<-renderUI({
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    filtg<-nrow(filtered())
   return(
    list(filtg,paste(' [',round((filtg/nrow(load_data()))*100,2),'%','] '))
    )
    })
  #number/percent failing filter
  output$filtsum4<-renderUI({
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    filtg<-nrow(filtered())
    datag<-nrow(load_data())
    return(
      list((datag-filtg),paste(' [',round(((datag-filtg)/nrow(load_data()))*100,2),'%','] '))
    )
    })
  
  #PLOTS - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  #output plot of median count vs variance, colored by passing/failing the filter
  output$medvarplot<-renderPlot({
    data<-load_data()[,-1]%>%mutate(med=apply(.,1,median),vari=apply(.,1,var),nonzero=rowSums(.>0))
    plot1<-ggplot(data,mapping=aes(x=vari,y=med))+
      geom_point(mapping=aes(color=((vari>quantile(vari,na.rm=TRUE,probs=(input$vslide/100)))&(nonzero>=input$nzslide))),na.rm=TRUE)+
      labs(title='Median Count vs Variance',x='Variance',y='Median Count',
           color=paste0('Over ',input$vslide,'% Variance \nand ',input$nzslide,' Non-Zero\n Samples'))+
      theme(plot.title=element_text(hjust=0.5))+
      scale_color_manual(values=c("FALSE"='#d9534f',"TRUE"="#4bbf73"))+
      scale_y_continuous(trans='log10')+
      scale_x_continuous(trans='log10')
    return(plot1)
    })
  
  #output non zero counts vs variance plot, colored by passing/failing the filter
  output$mednzplot<-renderPlot({
    data<-load_data()[,-1]%>%mutate(med=apply(.,1,median),vari=apply(.,1,var),nonzero=rowSums(.>0))
    plot1<-ggplot(data,mapping=aes(x=(length(data[1,])-(nonzero+2)),y=med))+
      geom_point(mapping=aes(color=((vari>quantile(vari,na.rm=TRUE,probs=(input$vslide/100)))&(nonzero>=input$nzslide))),na.rm=TRUE)+
      labs(title='Median Count vs Number of Zeros',x='Number of Zeros',y='Median Count',
           color=paste0('Over ',input$vslide,'% Variance \nand ',input$nzslide,' Non-Zero\n Samples'))+
      theme(plot.title=element_text(hjust=0.5))+
      scale_color_manual(values=c("FALSE"='#d9534f',"TRUE"="#4bbf73"))+
      scale_y_continuous(trans='log10')+
      xlim(0,10)
    return(plot1)
  })
  
  #Output heatmap using heatmap.2 function with color bar
  output$heatmap <- renderPlot({
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    int_matrix<-as.matrix(dplyr::select(as.data.frame(filtered()),c(colnames(filtered()))))
    int_matrix<-log10(int_matrix+1)
    pal<-colorRampPalette(brewer.pal(11, 'RdBu'))(25)
    #I use heatmap.2 since it clusters better and has a continuous color legend
    heatmap.2(int_matrix,col = pal,density.info='none',scale='none',trace='none',keysize=1,key.xlab='Log10 Read Counts')
  })
  
  #Output PCA plot
  #runs PCA, then principal components are used for graphing
  output$pca<- renderPlot({
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    s<-as.matrix(filtered())
    expr_mat<-as.data.frame(t(s))
    pca<-prcomp(expr_mat,center=TRUE,scale=FALSE)
    pca_var<-tibble(PC=factor(str_c("PC",1:length(colnames(pca$x))),str_c("PC",1:length(colnames(pca$x)))),
                    Variance=pca$sdev**2,vex=Variance/sum(Variance)*100)
    xv<-as.integer(str_sub(input$PCX,start=3))
    yv<-as.integer(str_sub(input$PCY,start=3))
    t<-as_tibble(pca$x)%>%mutate(timepoint=str_sub(rownames(pca$x),2,3))
    ggplot(t,aes(x=!!sym(input$PCX),y=!!sym(input$PCY),color=timepoint))+
      #labs show the percent variance for each PC being graphed
      geom_point()+labs(title='PCA',x=paste0(input$PCX,' ',as.character(round(pca_var$vex[xv])),'% variance'),
                                   y=paste0(input$PCY,' ',as.character(round(pca_var$vex[yv])),'% variance'))+
      theme(plot.title=element_text(hjust=0.5))
  })
  
  #Output bee swarm plot of PCA results
  output$bees<- renderPlot({
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    s<-as.matrix(filtered())
    expr_mat<-as.data.frame(t(s))
    pca<-prcomp(expr_mat,center=TRUE,scale=FALSE)
    bee<-as_tibble(pca$x) %>% select(1:input$pcs)%>%
      pivot_longer(everything(),names_to="PC",values_to="projection") %>%
      mutate(PC=forcats::fct_relevel(PC,str_c("PC",1:input$pcs))) %>%
      ggplot(aes(x=PC,y=projection,color=PC)) +
      geom_beeswarm() + labs(title="PCA Projection Plot") +
      theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0.5),plot.title=element_text(hjust=0.5))
    return(bee)
  })
  
  #Output UI for showing top N principal components for beeswarm plot
  #will only appear if the data is loaded in
  output$pcs1<-renderUI({
    if (!is.null(load_data())){
      sliderInput('pcs',min=1,max=8,
                  label='Show Top N Principle Components',
                  value=5,step=1)
    }
  })
  
  #Output UI for selecting which principal components are graphed in PCA
  #The plot customization UIs will render only once data has been uploaded
  output$pcx<- renderUI({
    if (!is.null(load_data())){
      selectizeInput('PCX',label='PC X',choices=str_c("PC",1:8),selected='PC1',width='50%')
    }
  })
  output$pcy<- renderUI({
    if (!is.null(load_data())){
      selectizeInput('PCY',label='PC Y',choices=str_c("PC",1:8),selected='PC2',width='50%')
    }
  })
  
#DIFFERENTIAL EXPRESSION TAB - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  #function for creating a volcano plot
  volcano_plot <-function(dataf, x_name, y_name, slider, color1, color2) {
    #If statement will only transform the y-axis to log10 if 'padj' is the y-axis variable
    if (y_name=='padj'){
      dataf$padj<-log10(dataf$padj)
      TF<-slider
    }
    else {
      TF<-10^slider
    }
    vp<-ggplot(dataf)+geom_point(aes(x=!!sym(x_name),y=!!sym(y_name),color=padj<TF,fill=NA))+
      scale_color_manual(values=c('FALSE'=color1,'TRUE'=color2))+
      labs(title=paste(y_name,'vs',x_name),x=x_name,y=y_name,color=paste0('padj<10^',TF))+
      theme(plot.title=element_text(hjust=0.5))
    #will output reverse y scale if the y value is padj
    ifelse(y_name=='padj',return(vp+scale_y_reverse()),return(vp))
  }
  
  #Output interactive DESeq2 results datatable 
  output$detab<-DT::renderDataTable({
    des<-data.frame(DE())%>%rownames_to_column(var='gene')
    DT::datatable(des, class='table-light',style='bootstrap5', width='auto', rownames=FALSE,
                  options = list(paging = TRUE, searching = TRUE, scrollX=TRUE,
                                 fixedColumns = TRUE, autoWidth = TRUE,
                                 ordering = TRUE, dom = 'Bfrtip'))
 
     })
  #Output for volcano plot 
  output$volcano <- renderPlot({
    input$button
    #Isolate function makes it so that the code will not run until the button is clicked
    isolate({
      volcano_plot(DE(),input$radx,input$rady,input$slider,input$color1,input$color2)
    })
  })
  



#INDIVIDUAL GENE EXPRESSION TAB - - - - - - - - - - - - - - - - - - - - - - - - -

  #Output UI for choosing which variable to graph
  output$info<-renderUI({
  req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
  meta1<-metainfo(filtered())
  return(selectizeInput('cat',p('CHOOSE A VARIABLE',style='color:white; text-align:center'),choices=colnames(meta1)))
  })

  #Output UI for choosing what type of graph to plot
  output$plts<-renderUI({
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    radioButtons('plots', p('CHOOSE PLOT TYPE',style='color:white; text-align:center'),
                 c('Bar Plot','Box Plot','Violin Plot','Beeswarm Plot')
    )
  })
  
  #Output UI for selecting which gene to graph
  #multiInput allows the user to select a gene from a list of all genes, with search functionality
  output$genes<-renderUI({
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    multiInput('gc',p('CHOOSE A GENE',style='color:white; text-align:center'),choices=rownames(filteredNN()),options=list(limit=1))
    
  })
  
  #output the chosen graph with the chosen variables
  output$graph<-renderPlot({
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    req(input$gc)
    input$button2
    isolate({
    td<-data.frame(t(filtered()))%>%rownames_to_column(var='sample')
    met<-metainfo(filtered())
    #combine the metadata and counts data, then pivots longer to make it graphable
    comb<-left_join(met,td,by=c('sample'='sample'))
    #conditionals based on the type of graph chosen
    if (input$plots=='Bar Plot'){
      sel<-comb%>%select(input$cat,starts_with('E'))%>%pivot_longer(starts_with('E'),names_to = 'gene',values_to='value')%>%
        filter(gene==input$gc)
      plot<-ggplot(sel,aes(!!sym(input$cat),value))+geom_bar(stat='identity',aes(fill=!!sym(input$cat)))+
        labs(title=paste('Normalized Counts by',input$cat),y='Number of Reads')+theme(plot.title=element_text(hjust=0.5))
      return(plot)
    }
    if (input$plots=='Box Plot'){
      sel<-comb%>%select(input$cat,starts_with('E'))%>%pivot_longer(starts_with('E'),names_to = 'gene',values_to='value')%>%
        filter(gene==input$gc)
      plot<-ggplot(sel,aes(!!sym(input$cat),value))+geom_boxplot(aes(fill=!!sym(input$cat)))+
        labs(title=paste('Normalized Counts by',input$cat),y='Number of Reads')+theme(plot.title=element_text(hjust=0.5))
      return(plot)
    }
    if (input$plots=='Violin Plot'){
      sel<-comb%>%select(input$cat,starts_with('E'))%>%pivot_longer(starts_with('E'),names_to = 'gene',values_to='value')%>%
        filter(gene==input$gc)
      plot<-ggplot(sel,aes(!!sym(input$cat),value))+geom_violin(aes(fill=!!sym(input$cat)))+
        labs(title=paste('Normalized Counts by',input$cat),y='Number of Reads')+theme(plot.title=element_text(hjust=0.5))
      return(plot)
    }
    if (input$plots=='Beeswarm Plot'){
      sel<-comb%>%select(input$cat,starts_with('E'))%>%pivot_longer(starts_with('E'),names_to = 'gene',values_to='value')%>%
        filter(gene==input$gc)
      plot<-ggplot(sel,aes(!!sym(input$cat),value))+geom_beeswarm(stat='identity',aes(color=!!sym(input$cat)))+
        labs(title=paste('Normalized Counts by',input$cat),y='Number of Reads')+theme(plot.title=element_text(hjust=0.5))
      return(plot)
    }
    })
  })
}

shinyApp(ui=ui,server=server)