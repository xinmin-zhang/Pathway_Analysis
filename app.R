###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: ui.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 9/6/2019
##@version 3.0
###########################################################################################################
source("pathwayApp_files/global.R",local = TRUE)$value

ui <- fluidPage(
          theme = shinytheme("cerulean"),
          windowTitle = "Pathway Analysis Tool",
          tagList(tags$head(tags$style(type = 'text/css','.navbar-brand{display:none;}')),
          tags$head(tags$style(".shiny-notification {position: fixed; top: 50% ;left: 30%; width: 50%")),
          useShinyjs(),
          titlePanel(
            fluidRow(
              column(8, h2(strong("Pathway Analysis", align = 'left')))
             # column(8,  h2(strong(textOutput('project')), align = 'left'))
            ),
              windowTitle = "Pathway Analysis" ),


navbarPage(title ="", id="menu",

##########################################################################################################
## Select Dataset
##########################################################################################################
tabPanel("Load Data",
  fluidRow(
    column(12,
    tabsetPanel(id="Tables",
       #tabPanel(title="Introduction",htmlOutput('intro')),
       #tabPanel(title="Project Table", DT::dataTableOutput('projecttable')),
       tabPanel(title="Upload Data",
         tags$p(tags$b("Web Tool for Functional Enrichment Analysis and Pathway Visualization")),
         tags$div(
           tags$img(src = "pathways.JPG", style = "height: 200px; display: block; margin-left: 10px; margin-right: auto;")
         ),
        tags$p(tags$i("Metabase and Reactome pathway view are for human gene set only.")),
        tags$hr(style="border-color: RoyalBlue;"),
         radioButtons("select_dataset",label="Data Source", choices=c("Upload Comparison Data (csv)", "Gene List", "Demo Data"),
            inline = T, selected= "Upload Comparison Data (csv)"),
         tags$hr(),
         conditionalPanel("input.select_dataset=='Demo Data'",
              selectInput("sel_project", label="Available Dataset",
                choices=c("", projects), selected=NULL)),
         conditionalPanel("input.select_dataset=='Upload Comparison Data (csv)'",
              uiOutput('upload.files.ui'), #see process_uploaded_files.R for details.
              tags$br(),
              textOutput('upload.message') ),
        conditionalPanel("input.select_dataset=='Gene List'",
                tagList(tags$div(
                   tags$p("You can paste a list of gene symbols and search for over-represented pathways."),
                   tags$p("The list can have one gene per line, or have all genes in one line but separated by comma."),
                   tags$p(tags$i("(Pathway visualization is not available for gene list as logFC data is required for visualization.)")),
                   tags$a(href="Genes_Down_in_ADvsCtrl.txt",  target = "_blank", "Download an example gene list file (624 genes down in ADvsCtrl)"),
                   tags$hr()
                ),
                actionButton("uploadList", "Start Analysis with Gene List")
            ) ),          
       ),
       tabPanel(title="Comparisons", 
         tags$p("Choose cutoff for differentially expressed genes (DEG)."),
         fluidRow(
         column(width=3,numericInput("deg_fccut", label= "Fold Change Threshold",  value = 1.2, min=1, step=0.1)),
         column(width=3,numericInput("deg_pvalcut", label= "P-value Threshold",  value=0.01, min=0, step=0.001)),
         column(width=3,radioButtons("deg_psel", label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),selected="Pval", inline = TRUE))
         ),
         tags$hr(),
         DT::dataTableOutput("deg_counts"),
         tags$p("Click a comparison name to view indivdual gene data below."), tags$hr(),
         selectInput("deg_test", label="Select a Comparison to View Data", choices=NULL),
         radioButtons("subset",label="Display", choices=c("All Genes", "Only DEGs" ), inline = T, selected= "Only DEGs"),
         DT::dataTableOutput("comp_table"))
#       tabPanel(title="Help", htmlOutput('help_input'))
    )
    )
  )
),
##

##########################################################################################################
## Venn Diagram
##########################################################################################################
tabPanel("Venn Diagram",
	fluidRow(
		column(3,
			wellPanel(
			  conditionalPanel("input.venn_combined=='Venn Diagram from Current Project'",
				numericInput("venn_fccut", label= "Choose Fold Change Cutoff", value = 1.5, min=1, step=0.1),
				numericInput("venn_pvalcut", label= "Choose P value Cutoff", value=0.05, min=0, step=0.001),
				radioButtons("venn_psel", label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),selected="Padj", inline = TRUE),
				radioButtons("venn_updown", label= "All, Up or Down?", choices= c("All"="All","Up"="Up","Down"="Down"),inline = TRUE),
				selectInput("venn_test1", label="Select List 1", choices=NULL),
				conditionalPanel("input.venn_tabset=='Venn Diagram'",
				colourInput("col1", "Select colour", "#0000FF",palette = "limited")
				),

				selectInput("venn_test2", label="Select List 2", choices=NULL),
				conditionalPanel("input.venn_tabset=='Venn Diagram'",
				colourInput("col2", "Select colour", "#FF7F00",palette = "limited")
				),

				selectInput("venn_test3", label="Select List 3", choices=NULL),
				conditionalPanel("input.venn_tabset=='Venn Diagram'",
				colourInput("col3", "Select colour", "#00FF00",palette = "limited")
				),

				selectInput("venn_test4", label="Select List 4", choices=NULL),
				conditionalPanel("input.venn_tabset=='Venn Diagram'",
				colourInput("col4", "Select colour", "#FF00FF",palette = "limited")
				),

				selectInput("venn_test5", label="Select List 5", choices=NULL),
				conditionalPanel("input.venn_tabset=='Venn Diagram'",
				colourInput("col5", "Select colour", "#FFFF00", palette = "limited")
				),
				conditionalPanel("input.venn_tabset=='Intersection Output'",
				radioButtons("vennlistname", label= "Label name", choices= c("Gene.Name"="Gene.Name","AC Number"="AC", "UniqueID"="UniqueID"),inline = TRUE, selected = "Gene.Name"))
			  ),
				conditionalPanel("input.venn_combined=='Venn Diagram Across Projects'",
				numericInput("vennP_fccut", label= "Choose Fold Change Cutoff", value = 1.2, min=1, step=0.1),
				numericInput("vennP_pvalcut", label= "Choose P value Cutoff", value=0.01, min=0, step=0.001),
				radioButtons("vennP_psel", label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
				checkboxInput("upperSymbols", "Upper Case Gene Symbols (e.g. mouse vs human)?",  FALSE, width="90%"),
				selectInput("dataset1", "Data set1", choices=NULL),
				selectInput("vennP_test1", label="Select List 1", choices=NULL),
				conditionalPanel("input.vennP_tabset=='Venn Diagram'",
				                 colourInput("vennPcol1", "Select colour", "#0000FF",palette = "limited")
				),
				selectInput("dataset2", "Data set2", choices=NULL),
				selectInput("vennP_test2", label="Select List 2", choices=NULL),
				conditionalPanel("input.vennP_tabset=='Venn Diagram'",
				                 colourInput("vennPcol2", "Select colour", "#FF7F00",palette = "limited")
				),
				selectInput("dataset3", "Data set3", choices=NULL),
				selectInput("vennP_test3", label="Select List 3", choices=NULL),
				conditionalPanel("input.vennP_tabset=='Venn Diagram'",
				                 colourInput("vennPcol3", "Select colour", "#00FF00",palette = "limited")
				),
				selectInput("dataset4", "Data set4", choices=NULL),
				selectInput("vennP_test4", label="Select List 4", choices=NULL),
				conditionalPanel("input.vennP_tabset=='Venn Diagram'",
				                 colourInput("vennPcol4", "Select colour", "#FF00FF",palette = "limited")
				),
				selectInput("dataset5", "Data set5", choices=NULL),
				selectInput("vennP_test5", label="Select List 5", choices=NULL),
				conditionalPanel("input.vennP_tabset=='Venn Diagram'",
				                 colourInput("vennPcol5", "Select colour", "#FFFF00", palette = "limited")
				)
				)
				)
			),
			column(9,
			  tabsetPanel(id="venn_combined",
			    tabPanel(title="Venn Diagram from Current Project",
				tabsetPanel(id="venn_tabset",
					tabPanel(title="Venn Diagram",
						column(9,
							actionButton("vennDiagram", "Save to output"),plotOutput("vennDiagram", height = 800)
						),
						column(3,
							textInput("title", "Title", width = "100%"),
							sliderInput("maincex", "Title Size", min = 0, max = 6, value = 3, width = "100%"),
							sliderInput("alpha", "Opacity", min = 0, max = 1, value = 0.4, width = "100%"),
							sliderInput("lwd", "Line Width", min = 1, max = 4, value = 1, width = "100%"),
							sliderInput("lty", "Line Type", min = 1, max = 6, value = 1, width = "100%"),
							radioButtons("fontface","Number Font face",list("plain", "bold", "italic"),selected = "plain",	inline = TRUE),
							sliderInput("cex", "Font size", min = 1, max = 4, value = 2, width = "100%"),
							radioButtons("catfontface","Label Font face",list("plain", "bold", "italic"),	selected = "plain", inline = TRUE),
							sliderInput("catcex", "Font size", min = 1, max = 2, step=0.1, value = 1.8, width = "100%"),
							sliderInput("margin", "Margin", min = 0, max = 1, step=0.05, value = 0.1, width = "100%")
						)
					),
					tabPanel(title="Venn Diagram (black & white)", plotOutput("SvennDiagram",height = 800, width = 800)),
					tabPanel(title="Intersection Output", htmlOutput("vennHTML")),
					tabPanel(title="DEG Table", actionButton("venn_DEG_data", "Save to output"), DT::dataTableOutput("venn_DEG_Data"))
					#tabPanel(title="Help", htmlOutput('help_venn'))
				)
			    )

				)
			)
		)
),


##########################################################################################################
## footer
##########################################################################################################

footer= HTML(footer_text)


)
)
) #for tagList

server <- function(input, output, session) {
  source("pathwayApp_files/inputdata.R",local = TRUE)
  source("pathwayApp_files/process_uploaded_files.R",local = TRUE)

  source("pathwayApp_files/venn.R",local = TRUE)
  source("pathwayApp_files/genesetmodule_PA.R",local = TRUE)
  insertTab(session=session,  inputId = "menu", target = "Load Data",  position = "after", 
              tabPanel("Gene Set Enrichment", geneset_ui(id = "GS")) )
  geneset_server(id = "GS")
 }

shinyApp(ui, server)