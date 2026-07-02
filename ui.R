##############################################
### shiny explorer.pro app - from REF2 gsva workspaces
### 29-Aug-2025, Anja
### 07-Aug-2025, Anja
### 18-Jul-2025, Anja
### 15-Jul-2025, Anja - first pro version for offline use
##############################################
library(shiny)
library(plotly)
library(DT)
#
### Define UI
##############################################
ui <- fluidPage(
    tags$style(HTML("
    .nav-tabs li a[data-value='Home'],
    .nav-tabs li a[data-value='Gene.Table'],
    .nav-tabs li a[data-value='Gene.Bubbles'],
    .nav-tabs li a[data-value='Gene.BoxPlots'],
    .nav-tabs li a[data-value='Genes.in.Sets'],
    .nav-tabs li a[data-value='X']{
      background-color: #ADD8E6; /* Light Blue */
      color: black;
    }
    .nav-tabs li a[data-value='Project.Home'],
    .nav-tabs li a[data-value='GeneSet.Table'],
    .nav-tabs li a[data-value='GeneSet.Bubbles'],
    .nav-tabs li a[data-value='GeneSet.BoxPlots'],
    .nav-tabs li a[data-value='Scatter']{
      background-color: #90EE90; /* Light Green */
      color: black;
    }
  ")),
  #
  titlePanel("APCRC-Q RNAseq Project Explorer"),
  #
  sidebarLayout(
    sidebarPanel(width=2,
      conditionalPanel(
        condition = "input.tabs == 'Home' || input.tabs == 'Project.Home'",
        p(em("explorer.pro.draft.2025-08-29"), style = "text-align: right;"),
        fileInput("file", "Choose any REF2 gene.set.scoring.Rdata workspace:", accept = c(".Rdata")),
        p(em("... once upload is complete, give it another 10-20 sec for the data to be processed ... the tool is ready when the project description is shown in the Project.Home tab ...")),
        hr()
      ),
      #
      conditionalPanel(
        condition = "input.tabs == 'Gene.Table' || input.tabs == 'Gene.Bubbles' || input.tabs == 'Gene.BoxPlots'",
        textAreaInput("genes", "Enter gene.symbols (one per line):", rows = 10),
        actionButton("plot.genes", "generate gene plots"),br(),
        actionButton("update.gene.table", "subset gene table"),br(),
        actionButton("restore.gene.table", "restore gene table"),br(),br()
      ),  
      conditionalPanel(
        condition = "input.tabs == 'Gene.Bubbles'", 
        actionButton("replot.gene.bubbles", "use uncorrected.pvalues")
      ),
      conditionalPanel(
        condition = "input.tabs == 'Gene.BoxPlots'", 
        actionButton("replot.gene.boxes", "use log2 transformation")
      ),
      #
      conditionalPanel(
        condition = "input.tabs == 'GeneSet.Table' || input.tabs == 'GeneSet.Bubbles' || input.tabs == 'GeneSet.BoxPlots'",
        textAreaInput("sets", "Enter gene.sets (one per line):", rows = 20),
        actionButton("plot.sets", "generate gene.set plots"),br(),
        actionButton("update.set.table", "subset gene.set table"),br(),
        actionButton("restore.set.table", "restore gene.set table"),br(),br()
      ), 
      conditionalPanel(
        condition = "input.tabs == 'GeneSet.Bubbles'", 
        actionButton("replot.set.bubbles", "use uncorrected.pvalues")
      ), 
      #
      conditionalPanel(
        condition = "input.tabs == 'Genes.in.Sets'",
        textAreaInput("oneset", "Enter a gene.set:", rows = 3),
        actionButton("find.genes", "list genes in set")
      ),  
      #
      conditionalPanel(
        condition = "input.tabs == 'Scatter'",
        textAreaInput("entry.a", "X.axis - enter a gene or gene.set:", rows = 3),
        textAreaInput("entry.b", "Y.axis - enter a gene or gene.set:", rows = 1),
        actionButton("plot.scatter", "generate scatter plot")
      )  
    ),
    #
    mainPanel(width=10, style = "padding-left:30px; padding-right:30px;",
       tabsetPanel(
         id = "tabs",
		     tabPanel("Home", value = "Home", br(),
		              p(em("... find our REF2 RNAseq data sets here ... U:/Research/Projects/ihbi/apcrcq/rnaseq_ref2 ...")),
				          p(em("... the gene.set.scoring.Rdata workspaces are in the respective project folders under ... GSVA/R.utilities ...")),
                  HTML('
					<table border="1" style="border-collapse:collapse; width: 100%;">
					  <thead>
						  <tr><th>Data.Set</th><th>Summary</th></tr>
						</thead>
					  <tbody>
						  <tr><td>ADIPOKX</td><td>xxx.</td></tr>
						  <tr><td>ALLO35CRX</td><td>xxx.</td></tr>
						  <tr><td>ALLOLU70X</td><td>xxx.</td></tr>
						  <tr><td>ALU35CRX</td><td>xxx.</td></tr>
						  <tr><td>ALY35CRX</td><td>xxx.</td></tr>
						  <tr><td>ALYMRX</td><td>xxx.</td></tr>
						  <tr><td>APCRCQ</td><td>xxx.</td></tr>
						  <tr><td>ATT1</td><td>xxx.</td></tr>	
						  <tr><td>ATT2</td><td>xxx.</td></tr>
						  <tr><td>ATT3</td><td>xxx.</td></tr>
						  <tr><td>ATT4</td><td>xxx.</td></tr>	
						  <tr><td>ATTX1</td><td>xxx.</td></tr>	
						  <tr><td>ATTX2</td><td>xxx.</td></tr>
						  <tr><td>CBL1078</td><td>xxx.</td></tr>
						  <tr><td>CBLX1</td><td>xxx.</td></tr>	
					  </tbody>
					</table>
				  ')
				  ),
         tabPanel("Project.Home", value = "Project.Home", br(), 
                   p(em("... once the input has been processed, a project overview will show here ...")),
                   verbatimTextOutput("displayDescription"),br(),
                   dataTableOutput("metaTable")),
         #
         tabPanel("Gene.Table", value = "Gene.Table", br(), 
                   p(em("... once the input has been processed, the GENE data table of the project will show here ...")),
                   dataTableOutput("outputTable")),
         tabPanel("Gene.Bubbles", value = "Gene.Bubbles", br(), 
                   p(em("... give it a moment ... this panel will show the differential expression results for the selected GENES of interest ...")),
                   plotlyOutput("bubblePlot")),
         tabPanel("Gene.BoxPlots", value = "Gene.BoxPlots", br(), 
                   p(em("... give it a moment ... this panel will show the expression levels for the selected GENES of interest ...")),
                   plotOutput("boxPlot")),
         #
         tabPanel("GeneSet.Table", value = "GeneSet.Table", br(), 
                   p(em("... once the input has been processed, the GENE.SET data table of the project will show here ...")),
                   dataTableOutput("outputTableSet")),
         tabPanel("GeneSet.Bubbles", value = "GeneSet.Bubbles", br(), 
                   p(em("... give it a moment ... this panel will show the differential expression results for the selected GENE.SETS of interest ...")),
                   plotlyOutput("bubblePlotSet")),
         tabPanel("GeneSet.BoxPlots", value = "GeneSet.BoxPlots", br(), 
                   p(em("... give it a moment ... this panel will show the expression levels for the selected GENE.SETS of interest ...")),
                   plotOutput("boxPlotSet")),
         #
         tabPanel("Genes.in.Sets", value = "Genes.in.Sets", br(), 
                   p(em("... finding the genes in a GENE.SET ...")),
                   verbatimTextOutput("GenesInSets"),
                   plotOutput("dummy1")),
                 #
         tabPanel("Scatter", value = "Scatter", br(), 
                   p(em("... scatter plot with correlation analysis ... works with genes (fpkm expression) and gene.sets (gsva.scores) ...")),
                   plotOutput("scatterPlot")),
                        #
         tabPanel("X", value = "X", br(), 
                   p(em("... you might be interested in exploring these genes and sets ...")),
                  HTML('
					<table border="1" style="border-collapse:collapse; width: 100%;">
					  <thead>
						  <tr><th>Category</th><th>Genes/Sets</th></tr>
						</thead>
					  <tbody>
						  <tr><td>androgen.activated.genes</td><td>KLK3, FKBP5, TMPRSS2,PGC, KLK2, CHRNA2, STEAP4, ORM1, LINC00844, ST6GALNAC1, KLK5, TRPM8, GNMT, SLC45A3, TUBA3D</td></tr>
						  <tr><td>androgen.repressed.genes</td><td>AR, UGT2B17, GRIK1, OPRK1, PLA2G2A, FOLH1, DDC, IFI27L2, CAMK2N1, LRRN1, SI, AMIGO2, DAB1, RP11-1136G4.2, FAM198B</td></tr>
						  <tr><td>androgen.activated.gene.sets</td><td>APCRCQ.REF2-ARPC-UP.221<br>Public-Hieronymus.2006.AR-UP.21<br>C2.CGP_DOANE_RESPONSE_TO_ANDROGEN_UP<br>C2.CGP_WANG_RESPONSE_TO_ANDROGEN_UP<br>C2.CGP_NELSON_RESPONSE_TO_ANDROGEN_UP</td></tr>
						  <tr><td>androgen.repressed.gene.sets</td><td>APCRCQ.REF2-ARPC-DN.199<br>Public-Hieronymus.2006.AR-DN.6<br>C2.CGP_DOANE_RESPONSE_TO_ANDROGEN_DN<br>C2.CGP_MOTAMED_RESPONSE_TO_ANDROGEN_DN<br>C2.CGP_NELSON_RESPONSE_TO_ANDROGEN_DN</td></tr>
						  <tr><td>neuroendocrine.marker.genes</td><td>CHGA,CHGB,ENO2,SYP,NCAM1,ACTL6B,SNAP25,INSM1,ASCL1,CHRNB2,SRRM4,CELF3,PCSK1,SOX2,POU3F2,LMO3,NKX2-1,SCG2,SYT4</td></tr>
						  <tr><td>Labrecque.mCRPC.subtype.gene.sets</td><td>xxx.</td></tr>
						  <tr><td>pan.epithelial.markers</td><td>EPCAM,CDH1,TACSTD2,KRT17,KRT6A,KRT16,KRT6B,KRT15,KRT6C,KRTCAP3,SFN</td></tr>
						  <tr><td>prostate.lineage.markers</td><td>NKX3-1,HOXB13,SOX9,RLN1,AMACR,PSCA,AZGP1,PCGEM1,PCA3</td></tr>
					  </tbody>
					</table>
				  ')
        ))
    )
  )
)
