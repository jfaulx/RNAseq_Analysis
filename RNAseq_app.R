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
library(rsconnect)
library(biomaRt)

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
  #I use a navbar to organize the apps into 4 different sections
  navbarPage(div(h3('RNA-seq Analysis'),h5('Jackson Faulx',style='color:#999; text-align:center; font-size: 16px')),
    #Samples tab - - - - - - - - - - - - - - 
    tabPanel('Samples',
      sidebarLayout(
        sidebarPanel(
          #file input only accepts csv and tsv files, wont let you submit other types of files
          fileInput('file1',p('LOAD COUNTS MATRIX',style='color:white; text-align:center'),accept=c('.csv','.tsv')),
          fileInput('meta1',p('LOAD METADATA',style='color:white; text-align:center'),accept=c('.csv','.tsv')),
          #these conditional panels will only output the text inside to the UI if a certain tab is selected
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
     #Counts tab - - - - - - - - - - - - - - - -
     tabPanel('Counts',
      sidebarLayout(
        sidebarPanel(
          #File inputs will be on every analysis page, but the data also carries throughout the app
          fileInput('file2',p('LOAD COUNTS MATRIX',style='color:white; text-align:center'),accept=c('.csv','.tsv')),
          fileInput('meta2',p('LOAD METADATA',style='color:white; text-align:center'),accept=c('.csv','.tsv')),
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
          #slider bars for the user to set filter thresholds
          sliderInput('vslide',min=0,max=100,
                      label=p('MINIMUM PERCENTILE OF VARIANCE',style='color:white; text-align:center'),
                      value=50,step=1),
          #This points to the number of non-zero samples slider, which depends on the number of samples in the counts file
          uiOutput('nonz')),
        mainPanel(
          #I create separate tabs within the 'Counts' tab
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
                column(width=3, offset=1,
                  uiOutput('pcx')),
                column(width=3,
                  uiOutput('pcy')
                        ),
                column(width=4, 
                  uiOutput('pcs1')
                        )
                      ),
              fluidRow(
                column(width=6, offset=1,
                  uiOutput('pcgroup')
                        )
                      )
                    )
                  )
                )
              )
            ),
     #Differential Expression tab - - - - - - - - - - - - - - - - - - - - - - -
     tabPanel('DE',
      sidebarLayout(
        sidebarPanel(
          fileInput('file3',p('LOAD COUNTS MATRIX',style='color:white; text-align:center'),accept=c('.csv','.tsv')),
          fileInput('meta3',p('LOAD METADATA',style='color:white; text-align:center'),accept=c('.csv','.tsv')),
          conditionalPanel(condition="input.det=='RESULTS'",
            h2('Results Table',style='color:white; text-align:center'),
            p('A sortable and searchable table containing the results of a differential expression analysis performed on the given
               counts matrix, comparing gene expression between the two groups of samples specified below')
                          ),
          conditionalPanel(condition="input.det=='VOLCANO PLOT'",
            h2('Volcano Plot',style='color:white; text-align:center'),
            p('A custom scatter plot of two result statistics, with color customization for genes that pass a specified adjusted
              p-value threshold')
                          ),
          #Allows the user to select design factors for running DESeq2
          textInput('des',p('ENTER A DESIGN',style='color:white; text-align:center'),value='~timepoint',placeholder='ex. ~timepoint + batch',),
          withSpinner(uiOutput('con'),type=8,color='white'),
          withSpinner(uiOutput('vs'),type=8,color='white'),
          actionButton('butt','RUN DESEQ2',width='100%',class='btn-primary')
                    ),
        mainPanel(
          tabsetPanel(id='det',
            tabPanel('RESULTS',
              withSpinner(DT::dataTableOutput('detab'),type=8,color='white')
                      ),
            tabPanel('VOLCANO PLOT',
              sidebarLayout(
                sidebarPanel(
                  #color picker allows customization of the colors on the graph
                  colourpicker::colourInput('color1',p('FAILING COLOR',style='color:white; text-align:center'),value='#860038'),
                  colourpicker::colourInput('color2',p('PASSING COLOR',style='color:white; text-align:center'),value='#FDBB30'),
                  colourpicker::colourInput('color3',p('HIGHLIGHT COLOR',style='color:white; text-align:center'),value='white'),
                  sliderInput('slider',min=0,max=15,
                              label=p('SELECT THE LOG2 FOLD CHANGE THRESHOLD',style='color:white; text-align:center'),
                              value=2,step=0.1)
                             ),
                mainPanel(
                  withSpinner(plotOutput('volcano'),type=8,color='white'),
                  br(),
                  #This part prints out the top 10 differentially expressed genes, allowing the user to select one and get additional info
                  fluidRow(
                    column(width=6,
                      withSpinner(uiOutput('topg'),type=8,color='white')
                                ),
                    column(width=6,
                      uiOutput('hg'),     
                      withSpinner(uiOutput('gid'),type=8,color='white'),
                      h4('DESCRIPTION:',style='color:white; text-align:center'),
                      withSpinner(uiOutput('description'),type=8,color='white')
                                )
                              )
                            )
                          )     
                        )
                      )
                    )
                  )
                ),
     #Individual Gene Expression tab - - - - - - - - - - - - - - - - - - - - - - 
     tabPanel('Gene Expression Visualizer',
      sidebarLayout(
       sidebarPanel(
        fileInput('file4',p('LOAD COUNTS MATRIX',style='color:white; text-align:center'),accept=c('.csv','.tsv')),
        fileInput('meta4',p('LOAD METADATA',style='color:white; text-align:center'),accept=c('.csv','.tsv')),
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
  
  #DATA PROCESSING FUNCTIONS - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  #Filter with CPM normalization via DESeq2
  filter1 <- function(verse_counts,variance,nz) {
    verse_counts1<-data.frame(verse_counts)%>%column_to_rownames(var='gene')
    vari<-apply(verse_counts1,1,var)
    #I filter the data using user-provided values for percentile of variance and number of non-zero samples
    v<-verse_counts1 %>% dplyr::filter(vari>quantile(vari,na.rm=TRUE,probs=(variance/100)))
    v<-v%>%dplyr::filter(rowSums(.>0)>nz)
    #I normalize using a function based off of DESeq2 normalization methods (defined below)
    v1<-deseq_normalize(v)
    return(v1)
  }
  
  #filter with no normalization (for running DESeq2)
  #Same as the above function, but doesn't include normalization (running DESeq2 will normalize the data)
  filternn <- function(verse_counts,variance,nz) {
    verse_counts1<-data.frame(verse_counts)%>%column_to_rownames(var='gene')
    vari<-apply(verse_counts1,1,var)
    v<-verse_counts1 %>% dplyr::filter(vari>quantile(vari,na.rm=TRUE,probs=(variance/100)))
    v<-v%>%dplyr::filter(rowSums(.>0)>nz)
    return(v)
  }
  
  #normalization function 
  deseq_normalize <- function(count_data) {
    cm<-as.matrix(count_data)
    #I turn the counts data into a matrix and turn it into a summarized experiment object
    des<-DESeqDataSetFromMatrix(countData = cm,
                                colData = tibble(sample_name=colnames(count_data)),design=~1)
    #I then normalize the counts and return them
    des<-estimateSizeFactors(des)
    ndseq<-data.frame(counts(des,normalized=TRUE))
    row.names(ndseq)<-row.names(count_data)
    return(ndseq)
  }
  
  #DESeq2 function
  run_deseq <- function(count_dataframe) {
    meta<-meta_data()
    #I create a summarized experiment object with the counts and associated metadata file, using the user-defined design
    dds <- DESeqDataSetFromMatrix(count_dataframe,colData=meta,design=as.formula(input$des)) 
    #I then run DESeq2 and return the results generated from the user-defined contrast
    dds <- DESeq(dds)
    res <- results(dds, contrast=c(input$contrast,input$num1,input$num2))
    return(res)
  }
  
  #REACTIVES - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #Load in Metadata
  #This req() fucntion allows the file to be input on any tab and will automatically load in data across the whole app
  #conditional statements for each possible file upload input
  meta_data <- reactive({
    req(isTruthy(input$meta1) || isTruthy(input$meta2) || isTruthy(input$meta3) || isTruthy(input$meta4))
    if (isTruthy(input$meta1)){
      f<-read_delim(input$meta1$datapath)
      return(f)
    }
    if (isTruthy(input$meta2)){
      f<-read_delim(input$meta2$datapath)
      return(f)
    }
    if (isTruthy(input$meta3)){
      f<-read_delim(input$meta3$datapath)
      return(f)
    }
    if (isTruthy(input$meta4)){
      f<-read_delim(input$meta4$datapath)
      return(f)
    }
  })
  
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
  
  #Filtered data is loaded in
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
  
  #outputs the summary of the metadata reactive
  output$meta<- renderTable({
   metadata<-meta_data()
   #I create a tibble containing the column name and type of data
   me<-tibble('Column'=colnames(metadata), 'Type'=sapply(metadata,class))
   #I create a fucntion that prints dicrete values for non-numeric columns, and mean and standard deviations for numeric columns
   distfunc<-function(x){
     case_when(x[['Type']]!='numeric' ~ paste(unique(metadata[[x[['Column']]]]),collapse=', '),
               x[['Type']]=='numeric' ~ paste0(mean(metadata[[x[['Column']]]]), ' (',sd(metadata[[x[['Column']]]]),')'))
   }
   #I apply the function to the tibble of columns that I created and add it to the table
   dsv<-apply(me,1,distfunc)
   me<-mutate(me, 'Distinct Values/Mean (SD)'= dsv)
   return(me)
   #I add this to add some minor styling to the UI table
  },width = "100%",bordered=TRUE,hover=TRUE)
  
  #Output interactive metadata table using DT
  output$si<-DT::renderDataTable({
  metadata<-meta_data()
  #I use DT to create an interactive data table of the user-uploaded metadata
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
    #I find the number of columns in the count data
    return(length(load_data()))
    })
  #number of genes
  output$filtsum2<-renderUI({
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    #I find the number of rows in the count data (number of genes)
    return(nrow(load_data()))
    })
  #number/percent passing
  output$filtsum3<-renderUI({
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    #I find the number of genes in the filtered data set and divide it by the number of genes in the unfiltered to get the percentage
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
    #I find the difference in the number of genes in the filtered data set vs the number of genes in the unfiltered set, dividing by the total to get percentage
    return(
      list((datag-filtg),paste(' [',round(((datag-filtg)/nrow(load_data()))*100,2),'%','] '))
    )
    })
  
  #This will be output in the UI to allow the user to select the filter threshold for minimum number of non-zero samples
  output$nonz<-renderUI({
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    sliderInput('nzslide',min=0,max=length(load_data()),label=p('MINIMUM NON-ZERO SAMPLES',style='color:white; text-align:center'),
                value=5,step=1)
  })
  
  #PLOTS - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  #output plot of median count vs variance, colored by passing/failing the filter
  output$medvarplot<-renderPlot({
    #I create columns that show the median count and variance for each row (gene)
    data<-load_data()[,-1]%>%mutate(med=apply(.,1,median),vari=apply(.,1,var),nonzero=rowSums(.>0))
    #I then plot the data, coloring based on if the genes would pass the user-chosen filter threshold
    plot1<-ggplot(data,mapping=aes(x=vari,y=med))+
      geom_point(mapping=aes(color=((vari>quantile(vari,na.rm=TRUE,probs=(input$vslide/100)))&(nonzero>=input$nzslide))),na.rm=TRUE)+
      labs(title='Median Count vs Variance',x='Variance',y='Median Count',
           color=paste0('Over ',input$vslide,'% Variance \nand ',input$nzslide,' Non-Zero\n Samples'))+
      theme(plot.title=element_text(hjust=0.5))+
      scale_color_manual(values=c("FALSE"='#d9534f',"TRUE"="#4bbf73"))+
      #I make both axes graphed on a log10 scale to improve readability
      scale_y_continuous(trans='log10')+
      scale_x_continuous(trans='log10')
    return(plot1)
    })
  
  #output non zero counts vs variance plot, colored by passing/failing the filter
  output$mednzplot<-renderPlot({
    #I create columns that show the median count,variance, and number of nonzero samples for each row (gene)
    data<-load_data()[,-1]%>%mutate(med=apply(.,1,median),vari=apply(.,1,var),nonzero=rowSums(.>0))
    #I then plot the data, coloring based on if the genes would pass the user-chosen filter threshold
    plot1<-ggplot(data,mapping=aes(x=(length(data[1,])-(nonzero+2)),y=med))+
      geom_point(mapping=aes(color=((vari>quantile(vari,na.rm=TRUE,probs=(input$vslide/100)))&(nonzero>=input$nzslide))),na.rm=TRUE)+
      labs(title='Median Count vs Number of Zeros',x='Number of Zeros',y='Median Count',
           color=paste0('Over ',input$vslide,'% Variance \nand ',input$nzslide,' Non-Zero\n Samples'))+
      theme(plot.title=element_text(hjust=0.5))+
      scale_color_manual(values=c("FALSE"='#d9534f',"TRUE"="#4bbf73"))+
      #I log10 transform the y axis to improve readability
      scale_y_continuous(trans='log10')+
      scale_x_continuous(breaks=seq(0,10,by=1))
    return(plot1)
  })
  
  #Output heatmap using heatmap.2 function with color bar
  output$heatmap <- renderPlot({
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    #I turn the filtered counts data into a matrix adn log transform the counts
    int_matrix<-as.matrix(dplyr::select(as.data.frame(filtered()),c(colnames(filtered()))))
    int_matrix<-log10(int_matrix+1)
    #I define the color palette that I want to use
    pal<-colorRampPalette(brewer.pal(11, 'RdBu'))(25)
    #I use heatmap.2 since it clusters better and has a continuous color legend
    heatmap.2(int_matrix,col = pal,density.info='none',scale='none',trace='none',keysize=1,key.xlab='Log10 Read Counts')
  })
  
  #Output PCA plot - - - - - - - - - -
  #runs PCA, then principal components are used for graphing
  output$pca<- renderPlot({
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    s<-as.matrix(filtered())
    #I first transpose the counts matrix and run the pca function
    expr_mat<-as.data.frame(t(s))
    pca<-prcomp(expr_mat,center=TRUE,scale=FALSE)
    #I create a tibble with all of the principal components, their variance, and their proportion of the total variance
    pca_var<-tibble(PC=factor(str_c("PC",1:length(colnames(pca$x))),str_c("PC",1:length(colnames(pca$x)))),
                    Variance=pca$sdev**2,vex=Variance/sum(Variance)*100)
    #I then create variables for the PCs that the user chooses to graph
    xv<-as.integer(str_sub(input$PCX,start=3))
    yv<-as.integer(str_sub(input$PCY,start=3))
    #I then add the metadata category that will group the points by color
    t<-as_tibble(pca$x)%>%mutate(cat=as.character(meta_data()[[input$GB]]))
    #I then plot the PCA variances according to which PCs the user wants graphed
    ggplot(t,aes(x=!!sym(input$PCX),y=!!sym(input$PCY),color=cat))+
      #labels show the percent variance for each PC being graphed
      geom_point()+labs(title='PCA',x=paste0(input$PCX,' ',as.character(round(pca_var$vex[xv])),'% variance'),
                                   y=paste0(input$PCY,' ',as.character(round(pca_var$vex[yv])),'% variance'),
                        color=input$GB)+
      theme(plot.title=element_text(hjust=0.5))
  })
  
  #Output bee swarm plot of PCA results
  output$bees<- renderPlot({
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    s<-as.matrix(filtered())
    #I first transpose the counts matrix and run the pca function
    expr_mat<-as.data.frame(t(s))
    pca<-prcomp(expr_mat,center=TRUE,scale=FALSE)
    #I then graph the number of PCs that the user specifies on a beeswarm plot
    bee<-as_tibble(pca$x) %>% dplyr::select(1:input$pcs)%>%
      #I pivot longer to put the tibble in a format for graphing
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
  output$pcgroup<- renderUI({
    if (!is.null(load_data()) & !is.null(meta_data())){
      selectInput('GB',label='GROUP BY',choices=colnames(meta_data()),width='50%')
    }
  })
  
#DIFFERENTIAL EXPRESSION TAB - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  output$con<-renderUI({
    #once data is loaded in, allows user to select which metadata group they want to run DESeq2 on
    req(isTruthy(input$meta1) || isTruthy(input$meta2) || isTruthy(input$meta3) || isTruthy(input$meta4))
    selectInput('contrast',p('SELECT CONTRAST',style='color:white; text-align:center'),choices=colnames(meta_data()[,-1]))
  })
  
  output$vs<-renderUI({
    #allows the user to select which two variables are to be compared for differential expression in DESeq2
    req(isTruthy(input$meta1) || isTruthy(input$meta2) || isTruthy(input$meta3) || isTruthy(input$meta4))
    fluidRow(
     column(4, selectInput('num1','',choices=unique(meta_data()[input$contrast]))),
     column(4, p('VS',style='color:white; text-align:center')),
     column(4, selectInput('num2','',choices=unique(meta_data()[input$contrast])))
     )
  })
  
  #function for creating a volcano plot for differential expression data
  volcano_plot <-function(dataf, slider, color1, color2, color3) {
    filt<-dataf%>%rownames_to_column(var='gene')
    filt$gene<-gsub("\\..*","",filt$gene)
    #this variable will be used to highlight one of the top 10 genes selected by the user
    filt<-filt%>%dplyr::filter(gene==input$tg)
    #I then graph the log2fold change vs the p-adjusted value for each gene, coloring based on a user-set log2 fold change threshold
    vp<-ggplot(dataf)+geom_point(aes(x=dataf$log2FoldChange,y=log10(dataf$padj),color=abs(dataf$log2FoldChange)>slider,fill=NA))+
      scale_color_manual(values=c('FALSE'=color1,'TRUE'=color2))+
      labs(title='Volcano Plot',x='Log2 Fold Change',y='P-Adjusted',color=paste('Log2FC >',slider))+
      theme(plot.title=element_text(hjust=0.5))+
      #I add vertical lines to denote the log2 fold change threshold
      geom_vline(xintercept = -slider, linetype="dashed", color='white')+
      geom_vline(xintercept = slider, linetype="dashed", color='white')+
      scale_y_reverse()+
      geom_point(data=filt,mapping=aes(x=filt$log2FoldChange,y=log10(filt$padj)),
                 pch=21, fill=NA, size=4, color=color3, stroke=1)
    return(vp)
 }
  
  #Output interactive DESeq2 results datatable 
  output$detab<-DT::renderDataTable({
    req(isTruthy(input$meta1) || isTruthy(input$meta2) || isTruthy(input$meta3) || isTruthy(input$meta4))
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    req(input$des)
    #input wont run unless the button is pressed
    input$butt
    isolate({
    #I use DTable to turn the DESeq2 results tibble into a searchable and sortable data table
    des<-data.frame(DE())%>%rownames_to_column(var='gene')
    DT::datatable(des, class='table-light',style='bootstrap5', width='auto', rownames=FALSE,
                  options = list(paging = TRUE, searching = TRUE, scrollX=TRUE,
                                 fixedColumns = TRUE, autoWidth = TRUE,
                                 ordering = TRUE, dom = 'Bfrtip'))
      })
     })
  
  #Output for volcano plot 
  output$volcano <- renderPlot({
      volcano_plot(DE(),input$slider,input$color1,input$color2, input$color3)
  })
  
  #Function to order genes in DESeq2 results by padj value
  best<-function(dat){
    ret<-data.frame(dat)%>%rownames_to_column(var='gene')
    #apply log2 fold change threshold filters
    ret<-ret%>%filter(abs(log2FoldChange)>input$slider)
    ret$gene<-gsub("\\..*","",ret$gene)
    return(ret[order(ret$padj),])
  }
  
  #Outputs top 10 DE genes by padj value
   output$topg<-renderUI({
     req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
     #I return the top 10 more differentially expressed genes
     bd<-head(best(DE()),10)
     return(radioButtons('tg',h3('TOP 10 DIFFERENTIALLY EXPRESSED GENES',style='color:white; text-align:center'),choices=bd$gene))
   })
   
  #function to run biomart and retrieve gene info
   gi<-function(data){
     #I pull otu annotation information for the top 10 genes found above
     bg<-best(data)
     ens<-bg$gene
     #If it is a mouse dataset, it will look for annotations in the mouse genome
     if (substr(ens[1],1,7)=='ENSMUSG'){
       mart<-useMart("ensembl", dataset = "mmusculus_gene_ensembl")
       ids<-getBM(attributes = c("ensembl_gene_id", "mgi_symbol",'mgi_description'), values = ens, mart = mart)
     }
     #If it is a human dataset, it will look for annotations in the human genome
     if (substr(ens[1],1,4)=='ENSG'){
       mart<-useMart("ensembl", dataset = "hsapiens_gene_ensembl")
       ids<-getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", 'description'), values = ens, mart = mart)
     }
     #I then combine the list of genes with the list of annotations and output the selected gene
     bm<-left_join(bg,ids,by=join_by(gene==ensembl_gene_id))
     bms<-bm%>%dplyr::filter(gene==input$tg)
     return(bms)
   }
   
   #UI outputs for gene name and description
   output$gid<-renderUI({
     g<-gi(DE())
     #outputs gene name
     h4(g[,8],style='color:grey; text-align:center')
   })
   
  output$description<-renderUI({
    g<-gi(DE())
    #outputs gene description
    h4(g[,9],style='color:grey; text-align:center')
  })

  output$hg<-renderUI({
    req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
    h3('HIGHLIGHTED GENE:',style='color:white; text-align:center')
  })
  
#INDIVIDUAL GENE EXPRESSION TAB - - - - - - - - - - - - - - - - - - - - - - - - -

  #Output UI for choosing which metadata variable to graph
  output$info<-renderUI({
  req(isTruthy(input$file1) || isTruthy(input$file2) || isTruthy(input$file3) || isTruthy(input$file4))
  meta1<-meta_data()
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
    #requires data to be loaded and a gene to be selected
    input$button2
    isolate({
    td<-data.frame(t(filtered()))%>%rownames_to_column(var='sample')
    met<-meta_data()
    #combine the metadata and counts data, then pivots longer to make it graphable
    comb<-left_join(met,td,by=c('sample'='sample'))
    #conditionals based on the type of graph chosen
    if (input$plots=='Bar Plot'){
      sel<-data.frame(comb)%>%dplyr::select(input$cat,starts_with('E'))%>%pivot_longer(starts_with('E'),names_to = 'gene',values_to='value')%>%
        filter(gene==input$gc)
      plot<-ggplot(sel,aes(!!sym(input$cat),value))+geom_bar(stat='identity',aes(fill=!!sym(input$cat)))+
        labs(title=paste('Normalized Counts by',input$cat),y='Number of Reads')+theme(plot.title=element_text(hjust=0.5))
      return(plot)
    }
    if (input$plots=='Box Plot'){
      sel<-data.frame(comb)%>%dplyr::select(input$cat,starts_with('E'))%>%pivot_longer(starts_with('E'),names_to = 'gene',values_to='value')%>%
        filter(gene==input$gc)
      plot<-ggplot(sel,aes(!!sym(input$cat),value))+geom_boxplot(aes(fill=!!sym(input$cat)))+
        labs(title=paste('Normalized Counts by',input$cat),y='Number of Reads')+theme(plot.title=element_text(hjust=0.5))
      return(plot)
    }
    if (input$plots=='Violin Plot'){
      sel<-data.frame(comb)%>%dplyr::select(input$cat,starts_with('E'))%>%pivot_longer(starts_with('E'),names_to = 'gene',values_to='value')%>%
        filter(gene==input$gc)
      plot<-ggplot(sel,aes(!!sym(input$cat),value))+geom_violin(aes(fill=!!sym(input$cat)))+
        labs(title=paste('Normalized Counts by',input$cat),y='Number of Reads')+theme(plot.title=element_text(hjust=0.5))
      return(plot)
    }
    if (input$plots=='Beeswarm Plot'){
      sel<-data.frame(comb)%>%dplyr::select(input$cat,starts_with('E'))%>%pivot_longer(starts_with('E'),names_to = 'gene',values_to='value')%>%
        filter(gene==input$gc)
      plot<-ggplot(sel,aes(!!sym(input$cat),value))+geom_beeswarm(stat='identity',aes(color=!!sym(input$cat)))+
        labs(title=paste('Normalized Counts by',input$cat),y='Number of Reads')+theme(plot.title=element_text(hjust=0.5))
      return(plot)
    }
    })
  })
}

shinyApp(ui=ui,server=server)