#Modified for Pathway App, in pathwayApp_files/ subfolder...

output$upload.files.ui <- renderUI({
  tagList(tags$div(
    tags$p(str_c("Prepare your own comparison data in Excel, save as a csv file and upload here. ", 
                 "The comparison data should have the following columns,")), 
    tags$ul(
      tags$li("UniqueID: gene symbol or other common identifiers like NCBI geneID, Ensembl ID, Unigene, UniProtKB ID, etc."),
      tags$li("test: name of the comparison. The file can include multiple comparisons, each with its own name listed in the test column."),
      tags$li("logFC: log2 fold change."),
      tags$li("P.Value: P value."),
      tags$li("Adj.P.Value: P value after mulitple test correction.")
      #tags$li("logExp: (optional), log value of average expression. Can be used to filter out low expressors.")
    ),
    tags$p(em(str_c("The column order doesn't matter. The file can be compresed as .gz or .zip file."))), 
    tags$a(href="human_brain_comparison_data.csv", "Download an example csv file (AD and BP brain proteome data, ~ 6K genes, 3 comparisons)"),
    tags$hr()
  ),
  
  radioButtons("upload_species", label= "Species",choices=c("human","mouse", "rat"), inline = TRUE, selected = "human"),
  checkboxInput("autogene", "Add gene symbols automatically (Uncheck if your file has gene symbols in the UniqueID column)",  TRUE, width="90%"),
  fileInput("F_comp", "Comparison Data File", width=500),
  actionButton("uploadData", "Submit Data")
  )
})


observeEvent(input$uploadData, {  
  cleanup_empty<-function(df) {
    df.empty=(is.na(df) | df=="")
    selCol=!(colSums( df.empty)==nrow(df))
    selRow=!(rowSums( df.empty)==ncol(df))
    return(df[selRow, selCol])
  }
  
  Project_name=format(Sys.time(), "%Y.%m.%d_%H.%M")
  ProjectID=str_c("PA.Upload_",  Project_name)
  if (nchar(ProjectID)>45)  {ProjectID=substr(ProjectID, 1, 45)}
  ProjectID=str_c(ProjectID,"_", stri_rand_strings(1,6) )
  #cat(ProjectID, "created.\n")
  #get expression data
  withProgress(message = 'Processing...', value = 0, {
    comp_file=input$F_comp$datapath
    if (str_detect(comp_file, "gz$") ) {comp_file=gzfile(comp_file, "rt")}
    if (str_detect(comp_file, "zip$") ) {
      fnames = as.character(unzip(comp_file, list = TRUE)$Name)
      comp_file=unz(comp_file, fnames[1])
    }
    results_long=read.csv(comp_file, header=T, check.names=F)
    #browser() #debug
    required_headers=c("UniqueID", "test", "logFC", "P.Value", "Adj.P.Value")
    if (!all(required_headers %in% colnames(results_long) )){
      showModal(modalDialog(
        title = "File Format Error",
        "Please include all required columns in the file. Please update the file and upload again.",
        easyClose = FALSE,
        footer = modalButton("OK")
      ))
      shiny::validate("Please include all required columns in the file.")
    }
   
    
    results_long=cleanup_empty(results_long)
    species=input$upload_species
    if (input$autogene){
      lookup_file=str_c("pathwayApp_files/", species, "_lookup.rds")
      lookup=readRDS(lookup_file)
      IDs=unique(results_long$UniqueID)
      ensembl_Per<-sum(str_detect(IDs, "ENSG|ENSMUSG|ENSRNOG"))/length(IDs)
      if (ensembl_Per>0.5) { #if more than 50% Ensembl gene IDs, remove version number
        results_long<-results_long%>%mutate(UniqueID=str_replace(UniqueID, "\\.\\d+$", "") )
      }
      IDs=unique(results_long$UniqueID)
      ProteinGeneName<-data.frame(id=0:(length(IDs)-1), UniqueID=IDs, Protein.ID=NA)%>%left_join(lookup%>%transmute(UniqueID=ID, Gene.Name=GeneName))
    } else {
      IDs=unique(results_long$UniqueID)
      ProteinGeneName<-data.frame(id=0:(length(IDs)-1), UniqueID=IDs, Protein.ID=NA, Gene.Name=IDs)
    }
    tests=unique(results_long$test)
    strOut=str_c("unlisted/", ProjectID, ".RData")
    save(ProteinGeneName, results_long, tests, file=strOut)
    unlisted_project<-data.frame(Name=ProjectID,  ShortName=ProjectID, 
                                 ProjectID, Species=species)
    write.csv(unlisted_project, str_c("unlisted/", ProjectID, ".csv") )
    unlisted_project=read.csv(str_c("unlisted/", ProjectID, ".csv"))
    ProjectInfo$ProjectID=ProjectID
    ProjectInfo$Name=unlisted_project$Name
    ProjectInfo$Species=unlisted_project$Species
    ProjectInfo$ShortName=unlisted_project$ShortName
    ProjectInfo$file1= paste("unlisted/",  ProjectID, ".RData", sep = "")  #data file
    ProjectInfo$file2= NULL #Correlation results
    showModal(modalDialog(
      title = "File Uploaded Successfully",
      "Now you can review the comparison data here, and perform pathway analysis by going to the Gene Set Enrichment Tab.",
      easyClose = FALSE,
      footer = modalButton("OK")
    ))
    updateTabsetPanel(session, "Tables", selected = "Comparisons")
    Reactome_iframe_visible(FALSE)
    ora_normal(ora_normal()+1)
    
  })
})

observeEvent(input$uploadList, {
   ora_list(ora_list()+1)
   updateTabsetPanel(session, "menu", selected = "Gene Set Enrichment")
   #updateTabsetPanel(session, "geneset_tabset", selected = "Over-Representation Analysis (ORA)")
})



observe({
  DataIn = DataReactive()
  tests = DataIn$tests
  updateSelectizeInput(session,'deg_test',choices=tests, selected=tests[1])
})


deg_counts_data <-reactive ({
  DataIn = DataReactive()
  results_long = DataIn$results_long
  FCcut = log2(as.numeric(input$deg_fccut))
  pvalcut = as.numeric(input$deg_pvalcut)

  if (input$deg_psel == "Padj") {
    tmpdat = results_long %>% dplyr::filter( Adj.P.Value < pvalcut & abs(logFC) > FCcut) 
  } else {
    tmpdat = results_long %>% dplyr::filter( P.Value < pvalcut & abs(logFC) > FCcut) 
  }
  deg_stat<-tmpdat%>%group_by(test)%>%dplyr::summarize(DEG=n(), Up=sum(logFC>0), Down=sum(logFC<0))%>%ungroup
  names(deg_stat)[1]="Comparison"
  more_comp=setdiff(unique(results_long$test), deg_stat$Comparison)
  if (length(more_comp)>0) {deg_stat<-rbind(deg_stat, data.frame(Comparison=more_comp, DEG=0, Up=0, Down=0))}
  deg_stat<-deg_stat%>%arrange(Comparison)%>%filter(!is.na(Comparison))

  return(deg_stat)
})

output$deg_counts <- DT::renderDT(server=FALSE,{
  DT::datatable(deg_counts_data(),extensions = 'Buttons',  selection = 'none', class = 'cell-border strip hover',
    options = list(dom = 'ft', pageLength = 20), 
    rownames= FALSE) %>% formatStyle(1, cursor = 'pointer',color='blue')
})
observeEvent(input$deg_counts_cell_clicked, {
  info = input$deg_counts_cell_clicked
  #browser() #debug
  if (is.null(info$value) || info$col != 0) return()
  updateSelectizeInput(session, 'deg_test', selected = info$value)
})


comp_data <-reactive ({
  DataIn = DataReactive()
  headers=c("UniqueID", "test", "logFC", "P.Value", "Adj.P.Value", "Gene.Name", "logExp")
  selected_results <-DataIn$results_long%>%dplyr::select(any_of(headers))%>%dplyr::filter(test==input$deg_test)
  if (input$subset=="Only DEGs"){
    FCcut = log2(as.numeric(input$deg_fccut))
    pvalcut = as.numeric(input$deg_pvalcut)
    if (input$deg_psel == "Padj") {
      selected_results<-selected_results%>%dplyr::filter( Adj.P.Value < pvalcut & abs(logFC) > FCcut) 
    } else {
      selected_results<-selected_results%>% dplyr::filter( P.Value < pvalcut & abs(logFC) > FCcut) 
    }
  }
  return(selected_results)
})

output$comp_table <- DT::renderDT(server=FALSE,{
  data=comp_data()
  DT::datatable(comp_data(),extensions = 'Buttons',  selection = 'none', class = 'cell-border strip hover',
    options = list(
      dom = 'lBfrtip', pageLength = 20,
      buttons = list(
        list(extend = "csv", text = "Download Page", filename = "Page_results",
             exportOptions = list(modifier = list(page = "current"))),
        list(extend = "csv", text = "Download All", filename = "Comparison_Results",
             exportOptions = list(modifier = list(page = "all")))
      )
    ), 
    rownames= FALSE)
})




