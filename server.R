##############################################
### shiny explorer.pro app - from REF2 gsva workspaces
### 07-Aug-2025, Anja
##############################################
library(shiny)
library(edgeR)
library(ggplot2)
library(plotly)
library(ggrepel)
library(ggpubr)
library(scales)
library(reshape2)
library(DT)
library(cowplot)
library(GSEABase)
#library(excelR) #NOPE, CAN'T HANDLE LARGE TABLES
#library(gridExtra)
#
#
options(shiny.maxRequestSize = 500 * 1024^2) #increase maximum file upload size 
#
##############################################
##############################################
### define server logic
##############################################
server <- function(input, output, session) {
  rv <- reactiveValues() #for multiple objects
  rvx <- reactiveVal() #one object ... later edgeR.table
  rvy <- reactiveVal() #one object ... later gsva.table
  rvz <- reactiveVal() #one object ... later project.description
  rva <- reactiveVal() #one object ... later gsets
    #
    #
    #######################################################
    ### display project description as soon as input is uploaded 
    #######################################################
    output$displayDescription <- renderText({
      req(rvz)
      project.description<-rvz()
      project.description
    })
    #
    output$metaTable <- renderDT({
          req(rv)  # wait until rv has data
          datatable(rv$meta.data, rownames=F, caption='Meta Data')
    })
    #
    #######################################################
    ### display data table(s) as soon as input is uploaded
    #######################################################
    output$outputTable <- DT::renderDT({
          req(rvx())  # wait until rvx() has data
          edgeR.table<-rvx()
          ncols<-ncol(edgeR.table)
          datatable(edgeR.table,
                    options = list(scrollX = TRUE, pageLength = 10, dom = "Blrtip",
                              buttons = list('colvis', 
                                        list(extend="colvisGroup", text="Show.Names.Only", show=c(0), hide=c(1:(ncols-1))),
                                        list(extend="colvisGroup", text="Show.All.Columns", show=c(0:(ncols-1)), hide=c()))), 
                    filter = 'top',
                    extensions = "Buttons",
                    selection = 'none',  # disable default row selection
                    rownames = FALSE,
                    escape = FALSE)
    })
    #
    #
    output$outputTableSet <- DT::renderDT({
          req(rvy())  # wait until rvy() has data
          gsva.table<-rvy()
          ncols<-ncol(gsva.table)
          datatable(gsva.table,
                    options = list(scrollX = TRUE, pageLength = 10, dom = "Blrtip",
                              buttons = list('colvis', 
                                        list(extend="colvisGroup", text="Show.Names.Only", show=c(0), hide=c(1:(ncols-1))),
                                        list(extend="colvisGroup", text="Show.All.Columns", show=c(0:(ncols-1)), hide=c()))),  
                    filter = 'top', 
                    extensions = "Buttons",
                    selection = 'none',  # disable default row selection
                    rownames = FALSE,
                    escape = FALSE)
    }) 
  
    
#######################################################
### user input - load workspace
#######################################################
observeEvent(input$file, {
    req(input$file)
    load(input$file$datapath)
    #
showNotification("checkpoint:upload.done", duration=5)
    #
    try(edgeR.output<-edgeR.output.H, silent=T) #compatibility with edgeR workspace 
    try(edgeR.output<-output.H, silent=T)
    try(edgeR.results.nc<-edgeR.results.nc.H, silent=T)
    try(edgeR.results.nc<-results.nc.H, silent=T)
    try(edgeR.output<-edgeR.output.M, silent=T) #compatibility with mouse workspace
    try(edgeR.output<-output.M, silent=T)
    try(edgeR.results.nc<-edgeR.results.nc.M, silent=T)
    try(edgeR.results.nc<-results.nc.M, silent=T)
    #
    #
    ### tailoring the edgeR.output to show as data table
    edgeR.output$description<-gsub("\\[.*","",edgeR.output$anno.Description) #shorten the Description
    edgeR.table<-edgeR.output[,grep("FC|pvalue|Name|res|fpkm|description",colnames(edgeR.output)),drop=F]
    #
    index <- grep("^res\\.", colnames(edgeR.table), value = TRUE)
    edgeR.table[index] <- lapply(edgeR.table[index], as.factor) #convert res columns to factors for better table filtering
    #
    index <- grep("^fpkm\\.", colnames(edgeR.table), value = TRUE)
    edgeR.table[index] <- lapply(edgeR.table[index], function(x) round(x, 3)) #round expression values
    #
    #index <- grep("pvalue\\.", colnames(edgeR.table), value = TRUE)
    #edgeR.table[index] <- lapply(edgeR.table[index], scales::label_number_auto()) #meaningful pvalue labels
    #
    colnames(edgeR.table)<-gsub("(?<=FC)\\.|(?<=pvalue)\\.|(?<=fpkm)\\.|(?<=res)\\.|\\.(?=vs)","<br>",colnames(edgeR.table),perl=TRUE) #break headers in DT
    #
    #
    ### tailoring the gsva.output to show as data table
    gsva.table<-gsva.output[,grep("dif|pvalue|res|avg.score",colnames(gsva.output)),drop=F]
    #
    gsva.table$gene.set<-rownames(gsva.table)
    gsva.table<-gsva.table[, c("gene.set", setdiff(names(gsva.table), "gene.set"))]  # Move to first position          
    #
    index <- grep("^res\\.", colnames(gsva.table), value = TRUE)
    gsva.table[index] <- lapply(gsva.table[index], as.factor) #convert res columns to factors for better table filtering
    #
    index <- grep("^avg.score\\.", colnames(gsva.table), value = TRUE)
    gsva.table[index] <- lapply(gsva.table[index], function(x) round(x, 3)) #round scores
    #
    colnames(gsva.table)<-gsub("(?<=dif)\\.|(?<=pvalue)\\.|(?<=avg.score)\\.|(?<=res)\\.|\\.(?=vs)","<br>",colnames(gsva.table),perl=TRUE) #break headers in DT
    #
    #
    ### tailoring the project description
    project.description<-paste0('Project: ', experiment.name,
                                '\nSummary: ', project.description,
                                '\n\nContrasts: ', paste(colnames(contrasts), collapse = ' - '),
                                '\nGene Expression Value: ', expression.value)
    #
    #
    rvx(edgeR.table) #assign to reactive object
    rvy(gsva.table) #assign to reactive object
    rvz(project.description) #assign to reactive object
    rva(gsets)
    #
    rm(list=setdiff(ls(),c('rv','rvx','rvy','rvz','rva','experiment.code','edgeR.output','edgeR.results','edgeR.results.nc','contrasts','level.order','group.colors','meta.data','gsva.output','gsva.results.nc','combined.data'))) #cleanup
    #
    #
    #######################################################      
    ### add all objects into reactive variable
    ####################################################### 
    obj_names <- ls()
    for (obj_name in obj_names) {
      rv[[obj_name]] <- get(obj_name)
    }
    #
#showNotification(paste("stored.objects:", paste(obj_names, collapse = ", ")), duration=10)
#session$flushReact() #errors
}, priority = 10)


####################################################### 
### make plotting objects
#######################################################  
observeEvent(rv$edgeR.output, { 
    req(rv)
    #
    edgeR.output <- rv$edgeR.output
    edgeR.results.nc <- rv$edgeR.results.nc
    contrasts <- rv$contrasts
    gsva.output <- rv$gsva.output
    gsva.results.nc <- rv$gsva.results.nc
    #
showNotification("making.plotting.objects ...", duration=20)
    #
    #######################################################
    ### PREPARE GENE.LEVEL OBJECTS FOR PLOTTING	
    #######################################################
    rownames(edgeR.output)<-paste(edgeR.output$Gene.Name, rownames(edgeR.output), sep="_") #combine symbol_ensemblID
    #
    dif<-edgeR.output[ ,grep("FC.",colnames(edgeR.output)),drop=F]
    dif<-melt(as.matrix(dif))
    colnames(dif)<-c("gene","contrast","fc")
    dif$contrast<-gsub("FC.","",dif$contrast)
    #
    pval<-edgeR.output[ ,grep("adj.pvalue.",colnames(edgeR.output)),drop=F]
    pval<-melt(as.matrix(pval))
    colnames(pval)<-c("gene","contrast","adj.pvalue")
    pval$contrast<-gsub("adj.pvalue.","",pval$contrast)
    #
    res<-edgeR.output[ ,grep("res.",colnames(edgeR.output)),drop=F]
    res<-melt(as.matrix(res))
    colnames(res)<-c("gene","contrast","res")
    res$contrast<-gsub("res.","",res$contrast)
    #
    pvalnc<-edgeR.output[ ,grep("non.corrected.pvalue.",colnames(edgeR.output)),drop=F] #added 19/8/24
    pvalnc<-melt(as.matrix(pvalnc))
    colnames(pvalnc)<-c("gene","contrast","pvalue")
    pvalnc$contrast<-gsub("non.corrected.pvalue.","",pvalnc$contrast)
    #
    rownames(edgeR.results.nc)<-rownames(edgeR.output) #!!! CAREFUL ... converts to symbol_ensemblID
    resnc<-melt(as.matrix(edgeR.results.nc)) #added 07/09/24
    colnames(resnc)<-c("gene","contrast","resnc")
    #
    ## naming it bubble.fx.data to make name unique in case workspace already has a bubble.data object
    bubble.fx.data<-merge(dif,res)
    bubble.fx.data<-merge(bubble.fx.data,pval)
    bubble.fx.data<-merge(bubble.fx.data,resnc)
    bubble.fx.data<-merge(bubble.fx.data,pvalnc)
    rm(dif,pval,res,pvalnc,resnc)
    #
    ## optional, to include overall gene expression level as an extra column (14-Apr-2023)
    bubble.fx.data$exp<-NA
    #
    exp<-edgeR.output[ ,grep("fpkm",colnames(edgeR.output))]
    exp<-data.frame(gene=rownames(exp), contrast="avg.fpkm", fc=0, res=0, adj.pvalue=NA, resnc=0, pvalue=NA, exp=round(rowMeans(exp),1))
    bubble.fx.data<-rbind(bubble.fx.data, exp)
    rm(exp)
    #
    ## amended 20/3/2025 for CBLX2 project to calculate avg per pdx model;
    #exp<-edgeR.output[ ,grep("fpkm.LTL352",colnames(edgeR.output))]
    #exp<-data.frame(gene=rownames(exp), contrast="LTL352.avg.fpkm", fc=0, res=0, adj.pvalue=NA, resnc=0, pvalue=NA, exp=round(rowMeans(exp),1))
    #bubble.fx.data<-rbind(bubble.fx.data, exp)
    #
    #exp<-edgeR.output[ ,grep("fpkm.LTL545",colnames(edgeR.output))]
    #exp<-data.frame(gene=rownames(exp), contrast="LTL545.avg.fpkm", fc=0, res=0, adj.pvalue=NA, resnc=0, pvalue=NA, exp=round(rowMeans(exp),1))
    #bubble.fx.data<-rbind(bubble.fx.data, exp)
    #
    bubble.fx.data$res<-factor(bubble.fx.data$res, levels=c(1,0,-1)) #be careful that legend.labels match factor levels!
    bubble.fx.data$de<-ifelse(bubble.fx.data$res==0,"nonDE","DE")
    bubble.fx.data$dex<-ifelse(bubble.fx.data$res==1,"de.UP",ifelse(bubble.fx.data$res==(-1),"de.DN","non.de")) #added 05/09/2024
    bubble.fx.data$dex<-factor(bubble.fx.data$dex, levels=c("de.UP","non.de","de.DN")) #be careful that legend.labels match factor levels!
    #
    bubble.fx.data$resnc<-factor(bubble.fx.data$resnc, levels=c(1,0,-1)) #be careful that legend.labels match factor levels!
    bubble.fx.data$denc<-ifelse(bubble.fx.data$resnc==0,"nonDE","DE")
    bubble.fx.data$dexnc<-ifelse(bubble.fx.data$resnc==1,"de.UP",ifelse(bubble.fx.data$resnc==(-1),"de.DN","non.de")) #added 05/09/2024
    bubble.fx.data$dexnc<-factor(bubble.fx.data$dexnc, levels=c("de.UP","non.de","de.DN")) #be careful that legend.labels match factor levels!
    #
    #bubble.fx.data$contrast<-factor(bubble.fx.data$contrast, levels=c(colnames(contrasts), "avg.fpkm")) #generic 
    #bubble.fx.data$contrast<-factor(bubble.fx.data$contrast, levels=c(colnames(contrasts), "LTL352.avg.fpkm","LTL545.avg.fpkm"))
    bubble.fx.data$symbol<-gsub("_.*","",bubble.fx.data$gene) #symbol off gene names
    #
    #dim(bubble.fx.data) #364326     13
    #colnames(bubble.fx.data)
    # [1] "gene"       "contrast"   "fc"         "res"        "adj.pvalue" "resnc"      "pvalue"    
    # [8] "exp"        "de"         "dex"        "denc"       "dexnc"      "symbol"   
    #
    rv$bubble.fx.data <- bubble.fx.data #add new object to reactive Var
    #
showNotification("bubble.fx.data.generated", duration=5)
    #
    #
    #######################################################
    ### PREPARE GENE.SET.LEVEL OBJECTS FOR PLOTTING	
    #######################################################
    dif<-gsva.output[ ,grep("dif.",colnames(gsva.output)),drop=F]
    dif<-melt(as.matrix(dif))
    colnames(dif)<-c("gene.set","contrast","dif")
    dif$contrast<-gsub("dif.","",dif$contrast)
    #
    pval<-gsva.output[ ,grep("adj.pvalue.",colnames(gsva.output)),drop=F]
    pval<-melt(as.matrix(pval))
    colnames(pval)<-c("gene.set","contrast","adj.pvalue")
    pval$contrast<-gsub("adj.pvalue.","",pval$contrast)
    #
    res<-gsva.output[ ,grep("res.",colnames(gsva.output)),drop=F]
    res<-melt(as.matrix(res))
    colnames(res)<-c("gene.set","contrast","res")
    res$contrast<-gsub("res.","",res$contrast)
    #
    pvalnc<-gsva.output[ ,grep("^pvalue.",colnames(gsva.output)),drop=F] #added 19/8/24
    pvalnc<-melt(as.matrix(pvalnc))
    colnames(pvalnc)<-c("gene.set","contrast","pvalue")
    pvalnc$contrast<-gsub("pvalue.","",pvalnc$contrast)
    #
    rownames(gsva.results.nc)<-rownames(gsva.output) 
    resnc<-melt(as.matrix(gsva.results.nc)) #added 07/09/24
    colnames(resnc)<-c("gene.set","contrast","resnc")
    #
    ## naming it bubble.fx.data.sets to make name unique in case workspace already has a bubble.data object
    bubble.fx.data.sets<-merge(dif,res)
    bubble.fx.data.sets<-merge(bubble.fx.data.sets,pval)
    bubble.fx.data.sets<-merge(bubble.fx.data.sets,resnc)
    bubble.fx.data.sets<-merge(bubble.fx.data.sets,pvalnc)
    rm(dif,pval,res,pvalnc,resnc)
    #
    bubble.fx.data.sets$res<-factor(bubble.fx.data.sets$res, levels=c(1,0,-1)) #be careful that legend.labels match factor levels!
    bubble.fx.data.sets$de<-ifelse(bubble.fx.data.sets$res==0,"nonDE","DE")
    bubble.fx.data.sets$dex<-ifelse(bubble.fx.data.sets$res==1,"de.UP",ifelse(bubble.fx.data.sets$res==(-1),"de.DN","non.de")) #added 05/09/2024
    bubble.fx.data.sets$dex<-factor(bubble.fx.data.sets$dex, levels=c("de.UP","non.de","de.DN")) #be careful that legend.labels match factor levels!
    #
    bubble.fx.data.sets$resnc<-factor(bubble.fx.data.sets$resnc, levels=c(1,0,-1)) #be careful that legend.labels match factor levels!
    bubble.fx.data.sets$denc<-ifelse(bubble.fx.data.sets$resnc==0,"nonDE","DE")
    bubble.fx.data.sets$dexnc<-ifelse(bubble.fx.data.sets$resnc==1,"de.UP",ifelse(bubble.fx.data.sets$resnc==(-1),"de.DN","non.de")) #added 05/09/2024
    bubble.fx.data.sets$dexnc<-factor(bubble.fx.data.sets$dexnc, levels=c("de.UP","non.de","de.DN")) #be careful that legend.labels match factor levels!
    #
    rv$bubble.fx.data.sets <- bubble.fx.data.sets #add new object to reactive Var
    #
showNotification("bubble.fx.data.sets.generated", duration=5)
    #
}, priority = -10)


#######################################################
### user input - genes of interest
#######################################################
observeEvent(input$plot.genes, {
    req(rv)
    req(input$genes)
    #
    genes <- unlist(strsplit(input$genes, "\n"))
    genes <- trimws(genes)
    genes <- genes[genes != ""]
    #
    bubble.fx.data <- rv$bubble.fx.data  
    experiment.code <- rv$experiment.code
    edgeR.output <- rv$edgeR.output
    contrasts <- rv$contrasts
    meta.data <- rv$meta.data
    #
    ################################################################
    ### bubble.plot for differential gene expression
    ################################################################
    show.contrasts<-c(colnames(contrasts), "avg.fpkm") #more flexible solution
    #show.contrasts<-levels(bubble.fx.data$contrast) #all contrasts
    #show.contrasts<-c("LTL352.CBL.vs.VEH","LTL352.avg.fpkm","LTL545.CBL.vs.VEH","LTL545.avg.fpkm")
    #bubble.fx.data<-bubble.fx.data[bubble.fx.data$contrast %in% show.contrasts,]
    #
    bubble.fx.data.sub<-bubble.fx.data[bubble.fx.data$symbol %in% genes, ] #subset 
    bubble.fx.data.sub$symbol<-factor(bubble.fx.data.sub$symbol, levels=rev(genes)) #match row order
    #
    output$bubblePlot <- renderPlotly({
    #
    bp1<-ggplot(bubble.fx.data.sub, aes(x=contrast, y=symbol, #))+ 
      text=paste(contrast,'\n',gene,'\nFC:',round(fc,2),'\nFDR.pvalue:',round(adj.pvalue,4),'\nDE.result:',dex,'\nuncorrected.pvalue:',round(pvalue,4))))+ #hovertext
      geom_point(aes(color=fc, shape=dex, size=-log10(adj.pvalue+1e-323)))+ 
      geom_text(aes(label=exp),size=3.5)+
      scale_size(range=c(3,10))+ #interactive needs smaller values
      scale_shape_manual(values=c("triangle-up","circle","triangle-down"),drop=F)+ #TRIANGLES
      scale_color_gradient2(low="#12497e",mid="white",high="#79091e", midpoint = 0, trans='pseudo_log')+ 
      scale_x_discrete(limits=show.contrasts)+ #categories and order to plot on x.axis
      theme_bw(base_size=12)+
      theme(axis.title=element_blank(),axis.text.x=element_text(angle=30,hjust=1,face="bold"),axis.text.y=element_text(face="bold"))+ 
      theme(legend.position = "right", legend.box = "horizontal",legend.margin = margin(0,0,0,0))+ #not ggplotly safe
      geom_vline(xintercept=length(colnames(contrasts))+0.5)+
      guides(size='none')+ #remove a particular legend, plotly compatible
      labs(shape=NULL) #remove a particular legend title, plotly compatible
    #print(bp) 
    ggplotly(bp1, tooltip='text', height=length(genes)*25+350)%>%layout(margin=list(t=150),
    title=list(x=0, xref='paper', font=list(size=17), text=paste0(experiment.code,' - differential expression results','<br><em><sub>DE is based on abs.FC>=1.5, FDR.pvalue<=0.05 and count-cutoff<br>shape = differential expression (DE), color = direction of change, size = -log10(pvalue)</em></sub>')))
    })
    #
    #
    ################################################################
    ### box.plots for gene expression
    ################################################################
    output$boxPlot <- renderPlot({
      plot_list <- list()
      #
      for(goi in genes){
        subset<-edgeR.output[edgeR.output$Gene.Name==goi,grep("Gene.Name|fpkm.",colnames(edgeR.output))]
        subset<-melt(subset, variable.name = "Sample", value.name = "fpkm")
        subset$Sample<-gsub('fpkm.','',subset$Sample)
        subset<-plyr::join(subset, meta.data, by='Sample', type='left')
        #
        p<-ggplot(subset,aes(x=Type,y=fpkm,fill=Type))+
          geom_boxplot(lwd=0.4,fatten=2,outlier.shape=NA)+ 
        	geom_point(position=position_jitter(seed=1, width=0.1, height=0.0),size=3,alpha=0.7)+
          scale_x_discrete(limits=rv$level.order)+
          theme_bw(base_size=14)+
          theme(legend.position="none",axis.text.x=element_text(angle=30,hjust=1))+
          geom_text_repel(aes(label=Repeat))+
          scale_fill_manual(values=rv$group.colors)+
          labs(x=NULL, y='fpkm', title=paste(experiment.code, '-', goi))
        plot_list[[goi]] <- p
        }
      #do.call(grid.arrange, c(plot_list, ncol = 1)) #gridExtra
      plot_grid(plotlist = plot_list, ncol = 1, align = 'v') #cowplot
    }, height = length(genes) *450 ) #default base height of plots is 400 px  
})

    
#######################################################
### user input - re-plot gene bubbles of interest
#######################################################
observeEvent(input$replot.gene.bubbles, {
    req(rv)
    req(input$genes)
    #
    genes <- unlist(strsplit(input$genes, "\n"))
    genes <- trimws(genes)
    genes <- genes[genes != ""]
    #
    bubble.fx.data <- rv$bubble.fx.data  
    experiment.code <- rv$experiment.code
    edgeR.output <- rv$edgeR.output
    contrasts <- rv$contrasts
    meta.data <- rv$meta.data
    #
    ################################################################
    ### bubble.plot for differential gene expression
    ################################################################
    show.contrasts<-c(colnames(contrasts), "avg.fpkm") #more flexible solution
    #
    bubble.fx.data.sub<-bubble.fx.data[bubble.fx.data$symbol %in% genes, ] #subset 
    bubble.fx.data.sub$symbol<-factor(bubble.fx.data.sub$symbol, levels=rev(genes)) #match row order
    #
    output$bubblePlot <- renderPlotly({
    #
    bp2<-ggplot(bubble.fx.data.sub, aes(x=contrast, y=symbol, #))+ 
      text=paste(contrast,'\n',gene,'\nFC:',round(fc,2),'\nuncorrected.pvalue:',round(pvalue,4))))+ #hovertext
      geom_point(aes(color=fc, shape=dexnc, size=-log10(pvalue+1e-323)))+ 
      geom_text(aes(label=exp),size=3.5)+
      scale_size(range=c(3,10))+ #interactive needs smaller values
      scale_shape_manual(values=c("triangle-up","circle","triangle-down"),drop=F)+ #TRIANGLES
      scale_color_gradient2(low="#12497e",mid="white",high="#79091e", midpoint = 0, trans='pseudo_log')+ 
      scale_x_discrete(limits=show.contrasts)+ #categories and order to plot on x.axis
      theme_bw(base_size=12)+
      theme(axis.title=element_blank(),axis.text.x=element_text(angle=30,hjust=1,face="bold"),axis.text.y=element_text(face="bold"))+ 
      theme(legend.position = "right", legend.box = "horizontal",legend.margin = margin(0,0,0,0))+ #not ggplotly safe
      geom_vline(xintercept=length(colnames(contrasts))+0.5)+
      guides(size='none')+ #remove a particular legend, plotly compatible
      labs(shape=NULL) #remove a particular legend title, plotly compatible
    #print(bp) 
    ggplotly(bp2, tooltip='text', height=length(genes)*25+350)%>%layout(margin=list(t=150),
    title=list(x=0, xref='paper', font=list(size=17), text=paste0(experiment.code,' - differential expression results','<br><em><sub>DE is based on abs.FC>=1.5, UNCORRECTED.pvalue<=0.05 and count-cutoff<br>shape = differential expression (DE), color = direction of change, size = -log10(pvalue)</em></sub>')))
    })
})
    

#######################################################
### user input - re-plot gene boxes of interest
#######################################################
observeEvent(input$replot.gene.boxes, {
    req(rv)
    req(input$genes)
    #
    genes <- unlist(strsplit(input$genes, "\n"))
    genes <- trimws(genes)
    genes <- genes[genes != ""]
    #
    bubble.fx.data <- rv$bubble.fx.data  
    experiment.code <- rv$experiment.code
    edgeR.output <- rv$edgeR.output
    contrasts <- rv$contrasts
    meta.data <- rv$meta.data
    #
    ################################################################
    ### box.plots for gene expression
    ################################################################
    output$boxPlot <- renderPlot({
      plot_list <- list()
      #
      for(goi in genes){
        subset<-edgeR.output[edgeR.output$Gene.Name==goi,grep("Gene.Name|fpkm.",colnames(edgeR.output))]
        subset<-melt(subset, variable.name = "Sample", value.name = "fpkm")
        subset$Sample<-gsub('fpkm.','',subset$Sample)
        subset<-plyr::join(subset, meta.data, by='Sample', type='left')
        #
        pp<-ggplot(subset,aes(x=Type,y=log2(fpkm+1),fill=Type))+
          geom_boxplot(lwd=0.4,fatten=2,outlier.shape=NA)+ 
        	geom_point(position=position_jitter(seed=1, width=0.1, height=0.0),size=3,alpha=0.7)+
          scale_x_discrete(limits=rv$level.order)+
          theme_bw(base_size=14)+
          theme(legend.position="none",axis.text.x=element_text(angle=30,hjust=1))+
          geom_text_repel(aes(label=Repeat))+
          scale_fill_manual(values=rv$group.colors)+
          labs(x=NULL, y='log2.fpkm', title=paste(experiment.code, '-', goi))
        plot_list[[paste0(goi,".b")]] <- pp
        }
      #do.call(grid.arrange, c(plot_list, ncol = 1)) #gridExtra
      plot_grid(plotlist = plot_list, ncol = 1, align = 'v') #cowplot
    }, height = length(genes) *450 ) #default base height of plots is 400 px 
})
    
    
#######################################################
### user input - subset gene table with selection
#######################################################
observeEvent(input$update.gene.table, {
    req(rvx)
    req(input$genes)
    #
    genes <- unlist(strsplit(input$genes, "\n"))
    genes <- trimws(genes)
    genes <- genes[genes != ""]
    #
    edgeR.table<-rvx()
    ncols<-ncol(edgeR.table)
    #
    output$outputTable <- DT::renderDT({
      datatable(edgeR.table[edgeR.table$Gene.Name %in% genes,],
                    options = list(scrollX = TRUE, pageLength = 10, dom = "Blrtip",
                              buttons = list('colvis', 
                                        list(extend="colvisGroup", text="Show.Names.Only", show=c(0), hide=c(1:(ncols-1))),
                                        list(extend="colvisGroup", text="Show.All.Columns", show=c(0:(ncols-1)), hide=c()))), 
                    filter = 'top',
                    extensions = "Buttons",
                    selection = 'none',  # disable default row selection
                    rownames = FALSE,
                    escape = FALSE)
    })      
})


#######################################################
### user input - restore gene table
#######################################################
  observeEvent(input$restore.gene.table, {
    req(rvx)
    edgeR.table<-rvx()
    ncols<-ncol(edgeR.table)
    #
    output$outputTable <- DT::renderDT({
      datatable(edgeR.table,
                    options = list(scrollX = TRUE, pageLength = 10, dom = "Blrtip",
                              buttons = list('colvis', 
                                        list(extend="colvisGroup", text="Show.Names.Only", show=c(0), hide=c(1:(ncols-1))),
                                        list(extend="colvisGroup", text="Show.All.Columns", show=c(0:(ncols-1)), hide=c()))), 
                    filter = 'top',
                    extensions = "Buttons",
                    selection = 'none',  # disable default row selection
                    rownames = FALSE,
                    escape = FALSE)
    })
})
    
    
#######################################################
### user input - gene.sets of interest
#######################################################
observeEvent(input$plot.sets, {
    req(rv)
    req(input$sets)
    #
    sets <- unlist(strsplit(input$sets, "\n"))
    sets <- trimws(sets)
    sets <- sets[sets != ""]
    #
    bubble.fx.data.sets <- rv$bubble.fx.data.sets  
    experiment.code <- rv$experiment.code
    gsva.output <- rv$gsva.output
    contrasts <- rv$contrasts
    meta.data <- rv$meta.data
    combined.data <- rv$combined.data
    #
    ################################################################
    ### bubble.plot for differential gene.set expression
    ################################################################
    show.contrasts<-c(colnames(contrasts)) #more flexible solution
    #
    bubble.fx.data.sub<-bubble.fx.data.sets[bubble.fx.data.sets$gene.set %in% sets, ] #subset 
    bubble.fx.data.sub$gene.set<-factor(bubble.fx.data.sub$gene.set, levels=rev(sets)) #match row order
    #
    output$bubblePlotSet <- renderPlotly({
    #
    bp1<-ggplot(bubble.fx.data.sub, aes(x=contrast, y=gene.set, #))+ 
      text=paste(contrast,'\n',gene.set,'\nDIF:',round(dif,2),'\nFDR.pvalue:',round(adj.pvalue,4),'\nDE.result:',dex,'\nuncorrected.pvalue:',round(pvalue,4))))+ #hovertext
      geom_point(aes(color=dif, shape=dex, size=-log10(adj.pvalue+1e-323)))+ 
      scale_size(range=c(3,10))+ #interactive needs smaller values
      scale_shape_manual(values=c("triangle-up","circle","triangle-down"),drop=F)+ #TRIANGLES
      scale_color_gradient2(low="#12497e",mid="white",high="#79091e", midpoint = 0, limits = c(min(bubble.fx.data.sets$dif),max(bubble.fx.data.sets$dif)), oob = scales::squish)+ 
      scale_x_discrete(limits=show.contrasts)+ #categories and order to plot on x.axis
      theme_bw(base_size=12)+
      theme(axis.title=element_blank(),axis.text.x=element_text(angle=30,hjust=1,face="bold"),axis.text.y=element_text(face="bold"))+ 
      theme(legend.position = "right", legend.box = "horizontal",legend.margin = margin(0,0,0,0))+ #not ggplotly safe
      #geom_vline(xintercept=length(colnames(contrasts))+0.5)+
      guides(size='none')+ #remove a particular legend, plotly compatible
      labs(shape=NULL) #remove a particular legend title, plotly compatible
    #print(bp) 
    ggplotly(bp1, tooltip='text', height=length(sets)*25+350)%>%layout(margin=list(t=150),
    title=list(x=0, xref='paper', font=list(size=17), text=paste0(experiment.code,' - differential expression results','<br><em><sub>DE is based on FDR.pvalue<=0.05 only<br>shape = differential expression (DE), color = direction of change, size = -log10(pvalue)</em></sub>')))
    })
    #
    #
    ################################################################
    ### box.plots for gene.set expression
    ################################################################
    output$boxPlotSet <- renderPlot({
      plot_list <- list()
      #
      for(gsoi in sets){
        subset<-combined.data[,c(gsoi,"Type","Repeat"),drop=F]
        #
        p<-ggplot(subset,aes(x=Type,y=!!ensym(gsoi),fill=Type))+
          geom_boxplot(lwd=0.4,fatten=2,outlier.shape=NA)+ 
        	geom_point(position=position_jitter(seed=1, width=0.1, height=0.0),size=3,alpha=0.7)+
          scale_x_discrete(limits=rv$level.order)+
          theme_bw(base_size=14)+
          theme(legend.position="none",axis.text.x=element_text(angle=30,hjust=1))+
          geom_text_repel(aes(label=Repeat))+
          scale_fill_manual(values=rv$group.colors)+
          labs(x=NULL, y='gsva.score', title=paste(experiment.code, '-', gsoi))
        #
        plot_list[[gsoi]] <- p
        }
      #do.call(grid.arrange, c(plot_list, ncol = 1)) #gridExtra
      plot_grid(plotlist = plot_list, ncol = 1, align = 'v') #cowplot
    }, height = length(sets) *450 ) #default base height of plots is 400 px  
})

    
#######################################################
### user input - re-plot gene.set bubbles of interest
#######################################################
observeEvent(input$replot.set.bubbles, {
    req(rv)
    req(input$sets)
    #
    sets <- unlist(strsplit(input$sets, "\n"))
    sets <- trimws(sets)
    sets <- sets[sets != ""]
    #
    bubble.fx.data.sets <- rv$bubble.fx.data.sets  
    experiment.code <- rv$experiment.code
    gsva.output <- rv$gsva.output
    contrasts <- rv$contrasts
    meta.data <- rv$meta.data
    #
    ################################################################
    ### bubble.plot for differential gene.set expression
    ################################################################
    show.contrasts<-c(colnames(contrasts)) #more flexible solution
    #
    bubble.fx.data.sub<-bubble.fx.data.sets[bubble.fx.data.sets$gene.set %in% sets, ] #subset 
    bubble.fx.data.sub$gene.set<-factor(bubble.fx.data.sub$gene.set, levels=rev(sets)) #match row order
    #
    output$bubblePlotSet <- renderPlotly({
    #
    bp2<-ggplot(bubble.fx.data.sub, aes(x=contrast, y=gene.set, #))+ 
      text=paste(contrast,'\n',gene.set,'\nDIF:',round(dif,2),'\nuncorrected.pvalue:',round(pvalue,4))))+ #hovertext
      geom_point(aes(color=dif, shape=dexnc, size=-log10(pvalue+1e-323)))+ 
      scale_size(range=c(3,10))+ #interactive needs smaller values
      scale_shape_manual(values=c("triangle-up","circle","triangle-down"),drop=F)+ #TRIANGLES
      scale_color_gradient2(low="#12497e",mid="white",high="#79091e", midpoint = 0, limits = c(min(bubble.fx.data.sets$dif),max(bubble.fx.data.sets$dif)), oob = scales::squish)+ 
      scale_x_discrete(limits=show.contrasts)+ #categories and order to plot on x.axis
      theme_bw(base_size=12)+
      theme(axis.title=element_blank(),axis.text.x=element_text(angle=30,hjust=1,face="bold"),axis.text.y=element_text(face="bold"))+ 
      theme(legend.position = "right", legend.box = "horizontal",legend.margin = margin(0,0,0,0))+ #not ggplotly safe
      #geom_vline(xintercept=length(colnames(contrasts))+0.5)+
      guides(size='none')+ #remove a particular legend, plotly compatible
      labs(shape=NULL) #remove a particular legend title, plotly compatible
    #print(bp) 
    ggplotly(bp2, tooltip='text', height=length(sets)*25+350)%>%layout(margin=list(t=150),
    title=list(x=0, xref='paper', font=list(size=17), text=paste0(experiment.code,' - differential expression results','<br><em><sub>DE is based on UNCORRECTED.pvalue<=0.05<br>shape = differential expression (DE), color = direction of change, size = -log10(pvalue)</em></sub>')))
    })
})

    
#######################################################
### user input - subset sets table with selection
#######################################################
observeEvent(input$update.set.table, {
    req(rvy)
    req(input$sets)
    #
    sets <- unlist(strsplit(input$sets, "\n"))
    sets <- trimws(sets)
    sets <- sets[sets != ""]
    #
    gsva.table<-rvy()
    ncols<-ncol(gsva.table)
    #
    output$outputTable <- DT::renderDT({
      datatable(gsva.table[gsva.table$gene.set %in% sets,],
                    options = list(scrollX = TRUE, pageLength = 10, dom = "Blrtip",
                              buttons = list('colvis', 
                                        list(extend="colvisGroup", text="Show.Names.Only", show=c(0), hide=c(1:(ncols-1))),
                                        list(extend="colvisGroup", text="Show.All.Columns", show=c(0:(ncols-1)), hide=c()))), 
                    filter = 'top',
                    extensions = "Buttons",
                    selection = 'none',  # disable default row selection
                    rownames = FALSE,
                    escape = FALSE)
    })      
})


#######################################################
### user input - restore sets table
#######################################################
observeEvent(input$restore.set.table, {
    req(rvy)
    gsva.table<-rvy()
    ncols<-ncol(gsva.table)
    #
    output$outputTable <- DT::renderDT({
      datatable(gsva.table,
                    options = list(scrollX = TRUE, pageLength = 10, dom = "Blrtip",
                              buttons = list('colvis', 
                                        list(extend="colvisGroup", text="Show.Names.Only", show=c(0), hide=c(1:(ncols-1))),
                                        list(extend="colvisGroup", text="Show.All.Columns", show=c(0:(ncols-1)), hide=c()))), 
                    filter = 'top',
                    extensions = "Buttons",
                    selection = 'none',  # disable default row selection
                    rownames = FALSE,
                    escape = FALSE)
    })
})

    
#######################################################
### user input - genes in set
#######################################################  
observeEvent(input$find.genes, {
    req(rva)
    gsets<-rva()
    #
    req(input$oneset)
    oneset<- unlist(strsplit(input$oneset, "\n"))
    oneset<- trimws(oneset)
    oneset<- oneset[oneset!= ""]
    #
    output$GenesInSets <- renderText({
      setgenes<-GSEABase::geneIds(gsets[[oneset]]) 
      paste0(oneset,' contains ',length(setgenes),' genes:\n\n',paste0(sort(setgenes),collapse='\n'))
    })
})
    
    
#######################################################
### user input - scatter
#######################################################
observeEvent(input$plot.scatter, {
    req(rv)
    req(input$entry.a)
    req(input$entry.b)
    #
    entry.a <- unlist(strsplit(input$entry.a, "\n"))
    entry.a <- trimws(entry.a)
    entry.a <- entry.a[entry.a != ""]
    #
    entry.b <- unlist(strsplit(input$entry.b, "\n"))
    entry.b <- trimws(entry.b)
    entry.b <- entry.b[entry.b != ""]
    #
    combined.data <- rv$combined.data  
        #
    ################################################################
    ### scatter plots
    ################################################################
    output$scatterPlot <- renderPlot({
    #
    validate(
    need(entry.a %in% colnames(combined.data), "Error: non valid X entry"),
    need(entry.b %in% colnames(combined.data), "Error: non valid Y entry"),
    )
    #
    sp<-ggplot(combined.data, aes(x=!!ensym(entry.a), y=!!ensym(entry.b)))+ 
      ggpubr::stat_cor(method='spearman', color='grey70', size=5, aes(label = paste("Spearman:",..r.label..,..rr.label..,..p.label..,sep = "~` `~")))+ 
      geom_smooth(method='lm', se=FALSE, color="grey70", linewidth=1, formula = y ~ x)+
      geom_point(aes(color=Type), size=10, alpha=0.5)+ 
      scale_color_manual(values=rv$group.colors)+ 
      geom_text(aes(label=Repeat),size=4)+ 
      theme_bw(base_size=12)
    print(sp)
    }, width=800, height=600)
})
        
}
##############################################
##############################################
