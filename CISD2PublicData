## Packages for this study.
    library(tidyverse)              # Data modulation and plotting.
    library(data.table)             # Data convenience features.
    require(readxl)                 # Reading and writing from excel.
    require(Rmisc)                  # Calculation SE for ploting SE bar.
    require(ggsignif)               # Significance of ploting using ggplot2.
    require(scales)                 # Scale Functions for Visualization.
    color = c("#A1A9D0", "#F0988C", "#96CCCB", "#F5DC75", '#579AC3',
              "#C4A5DE", "#A5CCC7", "#CFEAF1", "#9E9E9E", "#7E7DB4")
## Figure01 CISD2 in cellaular senescence
#### Seurat analysis GSE227136
    options(future.globals.maxSize = 4096*1024^2 )
    set.seed(2811)
    setwd("~/Studio/01IPFfibroblast/")
    library(data.table)
    library(googlesheets4)
    library(RCurl)
    library(qdap)
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(ade4)
    library(Matrix)
    library(harmony)
    source("utilities.R")
    source("lung_celltype_markers.R")
    source("lung_celltype_colors.R")
    ild_all <- readRDS("~/Studio/01IPFfibroblast/GSE227136_ILD_all_celltypes_Seurat.rds")
    umap <- as.data.frame(ild_all@reductions$umap@cell.embeddings)
    celltype <- as.data.frame(ild_all@meta.data  %>%  select(c('lineage', 'manual_annotation_1', 'Diagnosis')))
    meta <- cbind(celltype, umap)
    plot1 = ggplot(data.frame(x=meta$UMAP_1, y=meta$UMAP_2)) +
        scattermore::geom_scattermore(aes(x, y, color = meta$'lineage'),
                   pointsize = 3, alpha = 0.3, pixels = c(1000,1000), interpolate=TRUE) + 
            scale_color_manual(values = color) +
            labs(x= "UMAP_1", y = "UMAP_2", title = 'GSE227136') + 
            theme_classic() + 
            theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    plot2 = ggplot(data.frame(x=meta$UMAP_1, y=meta$UMAP_2)) +
        scattermore::geom_scattermore(aes(x, y, color = meta$'Diagnosis'),
                   pointsize = 3, alpha = 1, pixels = c(1000,1000), interpolate = TRUE) + 
            scale_color_manual(values = color) +
            labs(x= "UMAP_1", y = "UMAP_2", title = 'GSE227136') + 
            theme_classic() + 
            theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    FeaturePlot(object = ild_all, features = "CISD2") 
    mydata = FetchData(ild_all,vars = c('CISD2', 'UMAP_1','UMAP_2'))
    plot3 = ggplot(data.frame(x = mydata$UMAP_1, y = mydata$UMAP_2)) +
        scattermore::geom_scattermore(aes(x, y, color = mydata$CISD2),
                   pointsize = 3, alpha = 1, pixels = c(1000,1000), interpolate=TRUE) + 
            scale_color_gradient(low = "white", high = "#F0988C",) +
            labs(x= "UMAP_1", y = "UMAP_2", title = 'GSE227136', color = 'CISD2') + 
            theme_classic() + 
            theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 12),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    cell1 = celltype %>% 
        na.omit() %>%  
        rownames_to_column('x') %>% 
        inner_join(mydata, by = 'x') %>% 
        # filter(manual_annotation_1 == 'AT2') %>% 
        mutate(Diagnosis = ifelse(Diagnosis != 'Control', 'IPF', 'Control')) %>% 
        dplyr::rename('celltype' = 'manual_annotation_1')
    x = data.frame()
    for (i in cell1$celltype %>% unique()) {
        x1 = t.test(cell1$CISD2 ~ cell1$Diagnosis, subset = cell1$celltype == i)
        x2 = data.frame(t = x1$statistic, pvalue = x1$p.value, celltype = i)
        x = rbind(x2, x)
    }    
    require(Rmisc)
    y = x %>% 
        filter(pvalue < 0.05) %>% 
        inner_join(cell1, by = "celltype") %>% 
        mutate(stats = ifelse(pvalue >= 0.01, '*', ifelse(pvalue >= 0.001, "**", "***"))) %>% 
        arrange(celltype)
    y$celltype = factor(y$celltype, levels = unique(y$celltype))
    plot4 = y %>% 
        # summarySE(., measurevar = "AT2", groupvars = "group") %>%
        ggplot(aes(x = Diagnosis, y = CISD2, fill = Diagnosis)) + 
            geom_violin(alpha = 0.8, linewidth = 0.2)  +
            theme_classic() +
            # ggsignif::geom_signif(aes(x = Diagnosis, y = CISD2, group = Diagnosis), map_signif_level = T, comparisons = list(c('Control', 'IPF')), size = 0.5, textsize = 4.5, vjust = 0.1, tip_length = 0.01, test = 't.test') +
            scale_fill_manual(values = color) +
            geom_point(size = 2, shape = 21, stroke = 0, color = 'black', alpha = 0.8) +
            labs(x= "", y = "CISD2 expression", title = 'GSE227136') + 
            theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'none',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 0),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    t.test(cell1$CISD2 ~ cell1$Diagnosis, data = y)
    xstat = x %>% 
        inner_join(y %>% select(c('stats', 'celltype')), by = "celltype") %>% 
        distinct(celltype, .keep_all = T) %>% 
        arrange((celltype))
    xstat$celltype = factor(xstat$celltype, levels = unique(xstat$celltype))
    ploty = y %>% 
        ggplot(aes(x = celltype, y = CISD2, fill = Diagnosis)) +
            geom_violin(alpha = 0.8, linewidth = 0.3) + 
            ggsignif::geom_signif(annotations = xstat$stats, xmin = c(1:20), xmax = c(1:20), y_position = c(2.25), size = 0, textsize = 4.5, vjust = 0.1, tip_length = 0.00) +
            scale_fill_manual(values = color) +
            theme_classic() + 
            coord_flip() +
            labs(x= "", y = "CISD2 expression", title = 'GSE227136') +
            theme(axis.text.y = element_text(size = 14, angle = 0, hjust = 1),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Studio/01IPFfibroblast/F016.jpg', width = 7, height = 6.5)
    save.image()
    require(patchwork)
    design = '#AAAAAAAAAAAAAAAA'
    plot1 / plot3
    ggsave('~/Studio/01IPFfibroblast/F015.jpg', width = 5, height = 6.5)
    head(y)
    ploty + plot_layout(design = design) 
    ggsave('~/Studio/01IPFfibroblast/F014.jpg', width = 12, height = 5)
    ydata = FetchData(ild_all, vars = c('CISD2', 'BCL2', 'BECN1', 'GIMAP5', 'SDHD', 'CANX', 'COX4I1', 'OSBP',  'UMAP_1','UMAP_2')) %>% 
        rownames_to_column('x')
    yydata = ydata %>% 
        na.omit() %>% 
        inner_join(y, by = 'x') %>% 
        filter(celltype == 'AT2') %>% 
        filter(Diagnosis == 'IPF')
    cor.test(yydata$CISD2.x, yydata$BCL2)
#### CISD2 in AT2
    AT2 = c(0.262972376417385, 0.274877822229537, 0.0929553631671417, 0.485391326689285, 0.333813838538755, 
        0.321787346088448, 0.354475887245615, 0.324140210264222, 0.259074842652597, 0.353314515201445, 0.601525268386735, 0.276136657916767, 0.203078173673341, 0.336702982477433, 0.314828231015725, 0.227587471962018, 0.200021822289487, 0.414602181526345, 0.139816739239797, 0.154443604146954, 0.148342811679714, 0.191083089920226, 0.279006852998741, 0.122408220707511, 0.104161857842754, 0.305982058939666, 0.297185068686748, 0.275496227295779, 0.218946225749429, 0.19333745377012, 0, 0.231877465172906, 0.0440990508153903, 0.223402743549885, 0.137367189005107, 0.218869454905241, 0.404249350925755, 0.309692945595813, 0.19330083033286, 0.140331197269569, 0, 0.171711325683941, 0.436779224603534, 0.613747190389058, 0.531700946197026, 0.56820112899566, 0.587314599955271, 0.358454508380118, 0.335088872550904, 0.27965008567785, 0.308640993296377, 0.16709327335554, 0.445087993777397, 0.452292992620164, 0.377764199961168, 0.3583339739785, 0.322405391101069, 0.42152961339776, 0.285303495250326, 0.299185801845793, 0.34799918938656, 0.251191203957009, 0.407354739948511, 0.314596677832416, 0.109386039923885, 0.528521918422999, 0.275009865395547, 0.320349837927847, 0.356028001599325, 0.416743491841915, 0.230792709230611, 0.287362185648153, 0.156097691549803, 0.296263815446949, 0.141904910387201, 0.397094598196361, 0.462259436244865, 0.290977289221756, 0.257349149579741, 0.360080997440221, 0.408852128779994, 0.27858544609953, 0.335154682550706, 0.386068562972179, 0.157208226187019, 0.103985041450719, 0.260898838005722, 0.39919483737885, 0.26292389793547, 0.368417404490614, 0.33160995031778, 0.265300502484885, 0.310449796859499, 0.290287872972824, 0.22337620369201, 0.243534752834763, 0.281228614998641, 0.314616854001022, 0.234025390510249, 0.303828002612311, 0.21932155166429, 0.458143164268815, 0.0695370757643083, 0.231105728537657, 0.0834526463003016, 0.338927548776011, 0.257086179675496, 0.349663969297882, 0.393979271415836)
    patients = c('THD0001', 'THD0002', 'THD0005', 'THD0006', 'THD0007', 'THD0008', 'THD0009', 'THD0010', 
        'THD0012', 'THD0014', 'THD0015', 'THD0016', 'THD0017', 'THD0019', 'THD0021', 'THD0022', 'THD0023', 'THD0026', 'TILD006', 'TILD010', 'TILD015', 'TILD019', 'TILD028', 'TILD030', 'TILD037', 'TILD039', 'TILD041', 'TILD049', 'TILD051', 'TILD055', 'TILD059', 'TILD062', 'TILD084', 'TILD102', 'TILD103', 'TILD109', 'TILD111', 'TILD113', 'TILD116', 'TILD123', 'TILD126', 'TILD136', 'VUHD071', 'VUHD072', 'VUHD073', 'VUHD078', 'VUHD080', 'VUHD101', 'VUHD103', 'VUHD104', 'VUHD106', 'VUHD107', 'VUHD108', 'VUHD109', 'VUHD111', 'VUHD112', 'VUHD113', 'VUHD114', 'VUHD115', 'VUHD117', 'VUHD65', 'VUHD66', 'VUHD67', 'VUHD68', 'VUHD69', 'VUHD70', 'VUHD84', 'VUHD85', 'VUHD92', 'VUHD94', 'VUHD95', 'VUHD98', 'VUILD100', 'VUILD102', 'VUILD104', 'VUILD108', 'VUILD53', 'VUILD54', 'VUILD55', 'VUILD57', 'VUILD58', 'VUILD59', 'VUILD60', 'VUILD61', 'VUILD62', 'VUILD63', 'VUILD64', 'VUILD67', 'VUILD69', 'VUILD73', 'VUILD77', 'VUILD79', 'VUILD82', 'VUILD83', 'VUILD84', 'VUILD85', 'VUILD86', 'VUILD87', 'VUILD88', 'VUILD89', 'VUILD90', 'VUILD91', 'VUILD92', 'VUILD93', 'VUILD95', 'VUILD96', 'VUILD97', 'VUILD98', 'VUILD99')    
    GSE227136 = data.frame(patients, AT2) %>% 
        mutate(group = ifelse(grepl('HD', patients), 'HD', 'ILD'))
    plotAT2 = GSE227136 %>% 
        # summarySE(., measurevar = "AT2", groupvars = "group") %>%
        ggplot(aes(x = group, y = AT2, fill = group)) + 
            geom_violin(alpha = 0.8, linewidth = 0)  +
            theme_classic() +
            ggsignif::geom_signif(aes(x = group, y = AT2, group = group), map_signif_level = T, comparisons = list(c('HD', 'ILD')), size = 0.5, textsize = 4.5, vjust = 0.1, tip_length = 0.01, test = 'wilcox.test') +
            scale_fill_manual(values = color) +
            geom_point(size = 2, shape = 21, stroke = 0, color = 'black', alpha = 0.8) +
            scale_x_discrete(labels = c('HD' = 'Control', 'ILD' = 'IPF')) +
            labs(x= "AT2", y = "CISD2 expression", title = 'GSE227136') + 
            theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'none',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 0),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    require(patchwork)
    design = 'AAAAABBBBBCCC'
    plot1 + plot3 + plotAT2 + plot_layout(design = design)
    ggsave('~/Studio/01IPFfibroblast/F010.jpg', width = 12, height = 4)
    head(ild_all@meta.data)
#### Correction plot
    g132607 = getGEO(filename = ("~/Studio/01IPFfibroblast/GSE132607_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F)
    anno = data.table::fread("~/Studio/01IPFfibroblast/GPL15207-17536.txt") %>% 
        select(c('ID', 'Gene Symbol'))
    g132607c = pData(g132607)
    g132607c1 = g132607c %>% 
        filter(characteristics_ch1.5 == 'diagnosis: IPF') %>% 
        select(geo_accession)
    g132607e = exprs(g132607) %>% 
        as.data.frame() %>% 
        rownames_to_column('ID') %>% 
        mutate_at('ID', as.factor) %>% 
        inner_join(anno, by = 'ID') %>% 
        dplyr::rename('symbol' = 'Gene Symbol') %>% 
        select(- c('ID')) %>% 
        select('symbol', everything()) %>% 
        mutate(rowMean = rowMeans(.[grep('GSM', names(.))])) %>%
        dplyr::arrange(desc(rowMean)) %>%
        dplyr::distinct(symbol, .keep_all = T) %>%
        dplyr::select(-c('rowMean')) %>% 
        select(g132607c1$geo_accession) %>% 
        filter(symbol == 'CISD2' | symbol == 'TIMM10' | symbol == 'TIMM10B' | symbol == 'CISD3' | symbol == 'TIMM9' | symbol == 'NDUFA8' | symbol == 'WFS1' | symbol == 'CANX' | symbol == 'BECN1' | symbol == 'BCL2' | symbol == 'GIMAP5') %>% 
        t()
    g132607c1 = g132607c %>% 
        filter(characteristics_ch1.5 == 'diagnosis: IPF') %>% 
        select(geo_accession)
    colnames(g132607e) = g132607e[1, ]
    g132607e = g132607e[-1,] %>% as.data.frame() %>% mutate_all(as.numeric)
    require(reshape)
    cormat <- round(cor(g132607e),2)
    head(cormat)
    get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)}
    lower_tri <- get_lower_tri(cormat)
    melted_cormat <- melt(lower_tri, na.rm = TRUE)
    ggplot(data = melted_cormat, aes(X2, X1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(high = "#F0988C", low = "#A1A9D0", mid = "white", 
    midpoint = 0, limit = c(-1,1), space = "Lab", 
    name="Pearson\nCorrelation") +
    theme_minimal()+ 
    geom_text(aes(X2, X1, label = value), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.8, 0.1),
      legend.direction = "horizontal")+
      guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                title.position = "top", title.hjust = 0.5))
      theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                size = 12, hjust = 1))+
      coord_fixed()
    head(g132607e)
    cor.test(g132607e$BECN1, g132607e$CISD2, method = c("pearson"))
    ggplot(data = g132607e, aes(x = CISD2, y = BECN1)) + 
        geom_point(size = 3, shape = 21) + 
        geom_smooth(method = "lm", color = "black") + 
        labs(title = "GSE132607", y = 'BECN1 expression', x = 'CISD2 expression') + 
        theme_classic() +
        theme(legend.position = 'none',
        legend.title = element_blank(),
        axis.text = element_text(size = 14, vjust = 0.5),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16))
