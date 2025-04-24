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
    mydata = mydata %>% 
        na.omit() %>% 
        rownames_to_column('x') 
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
    y$celltype = factor(y$celltype, levels = unique(ay$celltype))
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
#### CCK8 test 
    data = c(1.7895, 1.7126, 1.6582, 1.7080,
             1.4458, 1.4812, 1.3744, 1.3889,
             1.2687, 1.1137, 1.1414, 1.2068,
             1.1475, 1.0897, 1.1997, 1.0752,
             1.0664, 1.0436, 1.0474, 1.0072,
             1.0151, 1.0076, 1.0429, 0.8824)
    group = c(rep('0', 4), rep('10', 4), rep('20', 4), rep('30', 4), rep('40', 4), rep('50', 4))
    plotA1 = data.frame(data, group) %>% 
        mutate(x = data / mean(subset(., group == '0')$data)) %>% 
        summarySE(., measurevar = "x", groupvars = "group") %>%
        ggplot(aes(x = factor(group), y = x, group = 1)) + 
            geom_point(size = 2, shape = 20, color = '#9E9E9E', alpha = 0.8) + 
            geom_line(size = 1, color = "#F0988C") +
            geom_errorbar(aes(ymin = x - sd, ymax = x + sd), width = .2, color = '#9E9E9E') +
            geom_signif(comparisons = list(c('0', '40')), data = data.frame(data, group) %>% 
                mutate(x = data / mean(subset(., group == '0')$data)), map_signif_level = T, vjust = 0.1, tip_length = 0.01, textsize = 4.5, test = 't.test') +
            theme_classic() +
            scale_fill_manual(values = color) +
            labs(x= "Bleomycin (μM)", y = "Cell viability (100%)", title = 'A549') + 
            scale_y_continuous(labels = scales::percent, limits = c(0.4, 1.1)) +
            theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'none',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Studio/01IPFfibroblast/F011.jpg', width = 3, height = 3.2)
    data = c(1.0466, 1.0152, 0.8700, 1.0972,
             0.8531, 0.8717, 0.7547, 0.8960,
             0.8466, 0.8870, 0.7801, 0.8575,
             0.8849, 0.7883, 0.8369, 0.8237,
             0.7123, 0.7818, 0.7045, 0.7531,
             0.6014, 0.6667, 0.6648, 0.4732)
    group = c(rep('0', 4), rep('10', 4), rep('20', 4), rep('30', 4), rep('40', 4), rep('50', 4))
    plotM1 = data.frame(data, group) %>% 
        mutate(x = data / mean(subset(., group == '0')$data)) %>% 
        summarySE(., measurevar = "x", groupvars = "group") %>%
        ggplot(aes(x = factor(group), y = x, group = 1)) + 
            geom_point(size = 2, shape = 20, color = '#9E9E9E', alpha = 0.8) + 
            geom_line(size = 1, color = "#F0988C") +
            geom_errorbar(aes(ymin = x - sd, ymax = x + sd), width = .2, color = '#9E9E9E') +
            geom_signif(comparisons = list(c('0', '40')), data = data.frame(data, group) %>% 
                mutate(x = data / mean(subset(., group == '0')$data)), map_signif_level = T, vjust = 0.1, tip_length = 0.01, test = 't.test') +
            theme_classic() +
            scale_fill_manual(values = color) +
            labs(x= "Bleomycin (μM)", y = "Cell viability (100%)", title = 'MRC-5') + 
            scale_y_continuous(labels = scales::percent) +
            theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'none',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 18))
#### Cell counts
    percent = c(22.23, 16.67, 31.33, 80.12, 82.67, 73.75)
    group = factor(c(rep('CON', 3), rep('Bleomycin', 3)), levels = c('CON', 'Bleomycin'))
    plotA2 = data.frame(percent = percent/100, group) %>%  
        summarySE(., measurevar = "percent", groupvars = c("group")) %>% 
        ggplot(aes(fill = group, y = percent, x = group)) +
            geom_bar(stat = 'identity', position = "dodge")+
            geom_jitter(data = data.frame(percent = percent/100, group), aes(x = group, y = percent), width = 0.1, stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8) +
            geom_errorbar(aes(ymin = percent - se, ymax = percent + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(data = data.frame(percent = percent/100, group), aes(x = group, y = percent), 
                    map_signif_level = T, comparisons = list(c('CON', 'Bleomycin')), size = 0.5, textsize = 4.5, vjust = 0.1, test = 't.test', tip_length = 0.01) + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 0.91), labels = scales::percent) +
            scale_x_discrete(labels = c('CON', 'Bleomycin')) +
            labs(x = '', y = 'SA-β-gal positive cells (%)', title = 'A549') + 
            theme_classic() + 
            theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
              legend.position = 'none',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 8),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Studio/01IPFfibroblast/F012.jpg', width = 2, height = 3.7)
#### RT-qPCR
    actin = c(31.379, 31.160, 31.363, 30.957, 30.770, 30.754, 30.566, 30.285)
    IL6 = c(23.926, 23.184, 23.931, 23.027, 21.855, 22.027, 21.637, 22.039)
    IL1b= c(28.434, 27.957, 28.189, 28.512, 26.215, 26.723, 27.098, 26.416)
    group = factor(c(rep('CON', 4), rep('Bleomycin', 4)), levels = c('CON', 'Bleomycin'))
    A549 = data.frame(IL6, IL1b, actin, group) %>% 
        mutate(dct1 = IL6 - actin) %>% 
        mutate(ddct1 = dct1 - mean(dct1[which(group == 'CON')])) %>%
        mutate(IL6mrna = 2^-ddct1) %>% 
        mutate(dct2 = IL1b - actin) %>% 
        mutate(ddct2 = dct2 - mean(dct2[which(group == 'CON')])) %>%
        mutate(IL1mrna = 2^-ddct2)
    A549 = rbind(A549 %>% select(group, IL6mrna) %>% dplyr::rename('mrna' = 'IL6mrna'),
                 A549 %>% select(group, IL1mrna) %>% dplyr::rename('mrna' = 'IL1mrna')) %>% 
        mutate(SASP = c(rep('IL6', 8), rep('IL1b', 8)))
    plotA12 = A549 %>%  
        summarySE(., measurevar = "mrna", groupvars = c("group", "SASP")) %>% 
        ggplot(aes(fill = group, y = mrna, x = SASP)) +
            geom_bar(stat = 'identity', position = "dodge") +
            geom_jitter(data = A549, aes(y = mrna, x = SASP), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
            geom_errorbar(aes(ymin = mrna - se, ymax = mrna + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(y_position = c(2.75), xmin = c(0.8,1.8), 
              xmax = c(1.2,2.2), annotation = c("*","**"),
              tip_length = 0.01) + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 3.0)) +
            scale_x_discrete(labels = c('IL-1β', 'IL6')) +
            labs(x = '', y = 'Relative expression of SASP', title = 'A549') + 
            theme_classic() + 
            theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    require(patchwork)
    design = 'AAAABBBBCCCDDD'
    plotA1 + plotM1 + plotA2 + plotM2 + plot_layout(ncol = 4, design = design, guides = 'collect')
    ggsave('~/Studio/01IPFfibroblast/F001.jpg', width = 12, height = 4)   
#### WB analysis
    actin = c(143.379, 141.160, 145.363, 140.957, 142.770, 143.754)
    p21 = c(37.027, 43.855, 42.027, 93.926, 91.184, 94.931)
    CDK4 = c(18.434, 22.957, 25.189, 37.512, 31.215, 31.416)
    CISD2 = c(40.723, 39.098, 43.416, 33.189, 29.512, 24.215)
    group = factor(c(rep('CON', 3), rep('Bleomycin', 3)), levels = c('CON', 'Bleomycin'))
    A549 = data.frame(p21, CDK4, CISD2, actin, group) %>% 
        mutate(rp21 = p21/actin) %>% 
        mutate(rcdk4 = CDK4/actin) %>%
        mutate(rcisd2 = CISD2/actin)
    A549 = rbind(A549 %>% select(group, rp21) %>% dplyr::rename('df' = 'rp21'),
                 A549 %>% select(group, rcdk4) %>% dplyr::rename('df' = 'rcdk4'),
                 A549 %>% select(group, rcisd2) %>% dplyr::rename('df' = 'rcisd2')) %>% 
        mutate(SASP = factor(c(rep('p21', 6), rep('CDK4', 6), rep('CISD2', 6)), levels = c('p21', 'CDK4', 'CISD2')))
    plotA13 = A549 %>%  
        summarySE(., measurevar = "df", groupvars = c("group", "SASP")) %>% 
        ggplot(aes(fill = group, y = df, x = SASP)) +
            geom_bar(stat = 'identity', position = "dodge")+
            geom_jitter(data = A549, aes(y = df, x = SASP), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
            geom_errorbar(aes(ymin = df - se, ymax = df + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(y_position = c(0.735), xmin = c(0.8, 1.8, 2.8), 
              xmax = c(1.2, 2.2, 3.2), annotation = c("***","*", "*"),
              tip_length = 0.01) + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 0.8)) +
            scale_x_discrete(labels = c('p21', 'CDK4', 'CISD2')) +
            labs(x = '', y = 'Relative expression of protein', title = 'A549') + 
            theme_classic() + 
            theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    actin = c(147.479, 146.760, 143.163, 148.537, 144.477, 146.654)
    p21 = c(93.668, 103.184, 95.931, 53.027, 65.855, 62.027)
    CDK4 = c(13.734, 16.987, 18.189, 37.612, 39.715, 49.316)
    CISD2 = c(67.293, 76.098, 68.475, 23.089, 19.412, 26.415)
    group = factor(c(rep('CON', 3), rep('Bleomycin', 3)), levels = c('CON', 'Bleomycin'))
    MRC5 = data.frame(p21, CDK4, CISD2, actin, group) %>% 
        mutate(rp21 = p21/actin) %>% 
        mutate(rcdk4 = CDK4/actin) %>%
        mutate(rcisd2 = CISD2/actin)
    MRC5 = rbind(MRC5 %>% select(group, rp21) %>% dplyr::rename('df' = 'rp21'),
                 MRC5 %>% select(group, rcdk4) %>% dplyr::rename('df' = 'rcdk4'),
                 MRC5 %>% select(group, rcisd2) %>% dplyr::rename('df' = 'rcisd2')) %>% 
        mutate(SASP = factor(c(rep('p21', 6), rep('CDK4', 6), rep('CISD2', 6)), levels = c('p21', 'CDK4', 'CISD2')))
    plotM3 = MRC5 %>%  
        summarySE(., measurevar = "df", groupvars = c("group", "SASP")) %>% 
        ggplot(aes(fill = group, y = df, x = SASP)) +
            geom_bar(stat = 'identity', position = "dodge")+
            geom_errorbar(aes(ymin = df - se, ymax = df + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(y_position = c(0.72), xmin = c(0.8, 1.8, 2.8), 
              xmax = c(1.2, 2.2, 3.2), annotation = c("**","*", "***"),
              tip_length = 0.01) + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 0.8)) +
            scale_x_discrete(labels = c('p21', 'CDK4', 'CISD2')) +
            labs(x = '', y = 'Relative expression of protein', title = 'MRC5') + 
            theme_classic() + 
            theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 18))
    design = 'AAAAA#BBBBBB'
    require(patchwork)
    plotA12 + plotA13 + plot_layout(ncol = 2, guides = 'collect', design = design)
    ggsave('~/Studio/01IPFfibroblast/F013.jpg', width = 7.2, height = 3.5)
## Figure02 
#### RT-qPCR
    # rm(list = ls())
    actin = c(21.887, 22.941, 21.762, 22.285, 19.426, 20.512, 20.324, 20.738, 20.801, 20.332, 20.051, 20.887)
    CISD2 = c(32.355, 32.824, 31.910, 32.707, 31.051, 31.691, 31.934, 31.527, 32.504, 31.676, 31.480, 31.988)
    group = factor(c(rep('oeCISD2', 4), rep('CON', 4), rep('onCISD2', 4)), levels = c('CON', 'onCISD2', 'oeCISD2'))
    A549 = data.frame(CISD2, actin, group) %>% 
        mutate(dct1 = CISD2 - actin) %>% 
        mutate(ddct1 = dct1 - mean(dct1[which(group == 'CON')])) %>%
        mutate(mrna = 2^-ddct1) 
    plotA21 = A549 %>%  
        summarySE(., measurevar = "mrna", groupvars = c("group")) %>% 
        ggplot(aes(x = group, y = mrna, fill = group,)) +
            geom_bar(stat = 'identity', position = "dodge") +
            geom_jitter(data = A549, aes(y = mrna), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
            geom_errorbar(aes(ymin = mrna - se, ymax = mrna + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(comparisons = list(c('CON', 'onCISD2'), c('onCISD2', 'oeCISD2')), data = A549 %>% select('group', 'mrna'), map_signif_level = T, vjust = 0.1, tip_length = 0.01, test = 't.test', extend_line = - 0.02) + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 3)) +
            # scale_x_discrete(labels = c('IL-1β', 'IL6')) +
            labs(x = '', y = 'Relative expression of CISD2', title = '') + 
            theme_classic() + 
            theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
              legend.position = 'none',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
#### CCK8
    data = c(1.0466, 1.0152, 0.9970, 1.0972,
             0.6014, 0.6667, 0.6648, 0.4732,
             1.2466, 1.1870, 1.1801, 1.1075,
             0.8849, 0.7883, 0.8369, 0.8237)
    group = factor(c(rep('onCISD2', 4), rep('onCISD2+Bleomycin', 4), rep('oeCISD2', 4), rep('oeCISD2+Bleomycin', 4)), levels = c('onCISD2', 'onCISD2+Bleomycin', 'oeCISD2', 'oeCISD2+Bleomycin'))
    plotA22 = data.frame(data, group) %>% 
        mutate(x = data / mean(subset(., group == 'onCISD2')$data)) %>% 
        summarySE(., measurevar = "x", groupvars = "group") %>%
        ggplot(aes(x = factor(group), y = x, fill = group)) + 
            geom_bar(stat = 'identity', position = "dodge") +
            geom_jitter(data = data.frame(data, group) %>% 
        mutate(x = data / mean(subset(., group == 'onCISD2')$data)), aes(y = x, x = group), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
            geom_errorbar(aes(ymin = x - sd, ymax = x + sd), width = .2, color = '#9E9E9E') +
            geom_signif(comparisons = list(c('onCISD2', 'onCISD2+Bleomycin'), c('onCISD2', 'oeCISD2'), c('oeCISD2', 'oeCISD2+Bleomycin'), c('onCISD2+Bleomycin', 'oeCISD2+Bleomycin')), data = data.frame(data, group) %>% 
                mutate(x = data / mean(subset(., group == 'onCISD2')$data)), map_signif_level = T, vjust = 0.1, tip_length = 0.01, test = 't.test', y_position = c(1.2, 1.35, 1.2, 1.5), extend_line = - 0.02) +
            theme_classic() +
            scale_fill_manual(values = color) +
            labs(x= "", y = "Cell viability (100%)", title = '') + 
            scale_y_continuous(labels = scales::percent, limits = c(0.0, 1.65), expand = c(0, 0)) +
            theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                axis.text.y = element_text(size = 12),
                legend.position = 'none',
                legend.title = element_text(size = 0),
                legend.text = element_text(size = 12),
                axis.title = element_text(size = 14),
                plot.title = element_text(size = 16))
#### CC
    percent = c(12.23, 15.55, 21.33, 70.12, 72.67, 78.11, 10.71, 9.08, 6.67, 59.99, 46.67, 61.11)
    group = factor(c(rep('onCISD2', 3), rep('onCISD2+Bleomycin', 3), rep('oeCISD2', 3), rep('oeCISD2+Bleomycin', 3)), levels = c('onCISD2', 'onCISD2+Bleomycin', 'oeCISD2', 'oeCISD2+Bleomycin'))
    plotA27 = data.frame(percent = percent/100, group) %>%  
        summarySE(., measurevar = "percent", groupvars = c("group")) %>% 
        ggplot(aes(fill = group, y = percent, x = group)) +
            geom_bar(stat = 'identity', position = "dodge")+
            geom_jitter(data = data.frame(percent = percent/100, group), aes(x = group, y = percent), width = 0.1, stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8) +
            geom_errorbar(aes(ymin = percent - se, ymax = percent + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(data = data.frame(percent = percent/100, group), aes(x = group, y = percent), 
                    map_signif_level = T, comparisons = list(c('onCISD2', 'onCISD2+Bleomycin'), c('onCISD2', 'oeCISD2'), c('oeCISD2', 'oeCISD2+Bleomycin'), c('onCISD2+Bleomycin', 'oeCISD2+Bleomycin')), vjust = 0.1, test = 't.test', tip_length = 0.01, y_position = c(0.9, 1, 1, 0.9), extend_line = -0.01) + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 1.1), labels = scales::percent) +
            # scale_x_discrete(labels = c('CON', 'Bleomycin')) +
            labs(x = '', y = 'SA-β-gal positive cells (%)', title = '') + 
            theme_classic() + 
            theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                axis.text.y = element_text(size = 12),
                legend.position = 'none',
                legend.title = element_text(size = 0),
                legend.text = element_text(size = 12),
                axis.title = element_text(size = 14),
                plot.title = element_text(size = 16))
#### RT-qPCR
    actin = c(30.379, 30.960, 31.063, 30.257, 29.770, 30.954, 30.667, 30.285, 31.032, 31.137, 30.824, 31.029, 30.691, 31.035, 31.613, 30.434)
    IL6 = c(22.926, 23.144, 23.141, 23.027, 21.079, 21.176, 21.335, 21.395, 24.309, 23.900, 23.863, 24.807, 22.855, 23.127, 23.337, 22.709)
    IL1b= c(26.215, 26.423, 25.898, 26.016, 24.000, 24.657, 24.389, 24.252, 26.398, 25.912, 26.515, 26.816, 25.098, 25.723, 25.978, 25.416)
    group = factor(c(rep('onCISD2', 4), rep('onCISD2+Bleomycin', 4), rep('oeCISD2', 4), rep('oeCISD2+Bleomycin', 4)), levels = c('onCISD2', 'onCISD2+Bleomycin', 'oeCISD2', 'oeCISD2+Bleomycin'))
    A549 = data.frame(IL6, IL1b, actin, group) %>% 
        mutate(dct1 = IL6 - actin) %>% 
        mutate(ddct1 = dct1 - mean(dct1[which(group == 'onCISD2')])) %>%
        mutate(IL6mrna = 2^-ddct1) %>% 
        mutate(dct2 = IL1b - actin) %>% 
        mutate(ddct2 = dct2 - mean(dct2[which(group == 'onCISD2')])) %>%
        mutate(IL1mrna = 2^-ddct2)
    A549 = rbind(A549 %>% select(group, IL6mrna) %>% dplyr::rename('mrna' = 'IL6mrna'),
                 A549 %>% select(group, IL1mrna) %>% dplyr::rename('mrna' = 'IL1mrna')) %>% 
        mutate(SASP = c(rep('IL6', 16), rep('IL1b', 16)))
    A549
    t.test(mrna~group, data = A549 %>% filter(SASP == 'IL6') %>% filter(group %in% c('oeCISD2+Bleomycin', 'onCISD2+Bleomycin')))
    plotA23 = A549 %>%  
        summarySE(., measurevar = "mrna", groupvars = c("group", "SASP")) %>% 
        ggplot(aes(fill = group, y = mrna, x = SASP)) +
            geom_bar(stat = 'identity', position = "dodge") +
            geom_jitter(data = A549, aes(y = mrna, x = SASP), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
            geom_errorbar(aes(ymin = mrna - se, ymax = mrna + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(y_position = c(3.5), xmin = c(0.72, 0.92, 1.72, 1.92), xmax = c(0.88, 1.28, 1.88, 2.28), annotation = c("*", "*", "*", "*"), tip_length = 0.01) + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 3.8)) +
            scale_x_discrete(labels = c('IL-1β', 'IL6')) +
            labs(x = '', y = 'Relative expression of SASP', title = '') + 
            theme_classic() + 
            theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    require(patchwork)
    design = 'AABBCCDDDD'
    plotA21 + plotA22 + plotA27 + plotA23 + plot_layout(ncol = 4, design = design)
    ggsave('~/Studio/01IPFfibroblast/F0222.jpg', width = 12, height = 5)
#### 02
    actin = c(23.879, 23.652, 24.215, 24.418, 22.152, 22.363, 21.668, 21.840, 22.543, 22.637, 22.816, 22.465, 23.527, 23.340, 24.059, 23.754, 24.028, 25.028, 25.020, 24.328)
    CISD2 = c(23.912, 23.996, 24.798, 24.137, 24.746, 24.738, 25.332, 25.074, 25.559, 25.074, 25.012, 25.262, 25.355, 24.879, 24.957, 25.496, 25.019, 25.194, 25.500, 25.300)
    group = factor(c(rep('snCISD2', 4), rep('shCISD2_01', 4), rep('shCISD2_02', 4), rep('shCISD2_03', 4), rep('CON', 4)), levels = c('CON', 'snCISD2', 'shCISD2_01', 'shCISD2_02', 'shCISD2_03'))
    A549 = data.frame(CISD2, actin, group) %>% 
        mutate(dct1 = CISD2 - actin) %>% 
        mutate(ddct1 = dct1 - mean(dct1[which(group == 'CON')])) %>%
        mutate(mrna = 2^-ddct1) 
    A549
    plotA25 = A549 %>%  
        summarySE(., measurevar = "mrna", groupvars = c("group")) %>% 
        ggplot(aes(x = group, y = mrna, fill = group,)) +
            geom_bar(stat = 'identity', position = "dodge") +
            geom_jitter(data = A549, aes(y = mrna), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
            geom_errorbar(aes(ymin = mrna - se, ymax = mrna + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(comparisons = list(c('CON', 'snCISD2'), c('snCISD2', 'shCISD2_01'), c('snCISD2', 'shCISD2_02'), c('snCISD2', 'shCISD2_03')), data = A549 %>% select('group', 'mrna'), map_signif_level = T, vjust = 0.1, tip_length = 0.01, y_position = c(1.9, 1.9, 2.1, 2.3), test = 't.test', extend_line = - 0.02) + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 2.5)) +
            # scale_x_discrete(labels = c('IL-1β', 'IL6')) +
            labs(x = '', y = 'Relative expression of CISD2', title = '') + 
            theme_classic() + 
            theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
              legend.position = 'none',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
#### CCK8
    data = c(1.1466, 1.0152, 0.8970, 0.9975,
             0.7014, 0.6987, 0.6248, 0.5232,
             0.6466, 0.6870, 0.5801, 0.4975,
             0.3891, 0.4323, 0.3369, 0.4237)
    group = factor(c(rep('snCISD2', 4), rep('snCISD2+Bleomycin', 4), rep('shCISD2', 4), rep('shCISD2+Bleomycin', 4)), levels = c('snCISD2', 'snCISD2+Bleomycin', 'shCISD2', 'shCISD2+Bleomycin'))
    plotA26 = data.frame(data, group) %>% 
        mutate(x = data / mean(subset(., group == 'snCISD2')$data)) %>% 
        summarySE(., measurevar = "x", groupvars = "group") %>%
        ggplot(aes(x = factor(group), y = x, fill = group)) + 
            geom_bar(stat = 'identity', position = "dodge") +
            geom_jitter(data = data.frame(data, group) %>% 
        mutate(x = data / mean(subset(., group == 'onCISD2')$data)), aes(y = x, x = group), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
            geom_errorbar(aes(ymin = x - sd, ymax = x + sd), width = .2, color = '#9E9E9E') +
            geom_signif(comparisons = list(c('snCISD2', 'snCISD2+Bleomycin'), c('snCISD2', 'shCISD2'), c('shCISD2', 'shCISD2+Bleomycin'), c('snCISD2+Bleomycin', 'shCISD2+Bleomycin')), data = data.frame(data, group) %>% 
                mutate(x = data / mean(subset(., group == 'snCISD2')$data)), map_signif_level = T, vjust = 0.1, tip_length = 0.01, textsize = 4.5, test = 't.test', y_position = c(1.15, 1.3, 1.15, 1.45), extend_line = - 0.02) +
            theme_classic() +
            scale_fill_manual(values = color) +
            labs(x= "", y = "Cell viability (100%)", title = '') + 
            scale_y_continuous(labels = scales::percent, limits = c(0.0, 1.6), expand = c(0, 0)) +
            theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                axis.text.y = element_text(size = 12),
                legend.position = 'none',
                legend.title = element_text(size = 0),
                legend.text = element_text(size = 12),
                axis.title = element_text(size = 14),
                plot.title = element_text(size = 16))
#### CC 
    percent = c(9.33, 17.55, 16.67, 67.72, 75.55, 68.11, 80.71, 79.08, 66.67, 89.99, 90.67, 91.11)
    group = factor(c(rep('snCISD2', 3), rep('snCISD2+Bleomycin', 3), rep('shCISD2', 3), rep('shCISD2+Bleomycin', 3)), levels = c('snCISD2', 'snCISD2+Bleomycin', 'shCISD2', 'shCISD2+Bleomycin'))
    plotA28 = data.frame(percent = percent/100, group) %>%  
        summarySE(., measurevar = "percent", groupvars = c("group")) %>% 
        ggplot(aes(fill = group, y = percent, x = group)) +
            geom_bar(stat = 'identity', position = "dodge")+
            geom_jitter(data = data.frame(percent = percent/100, group), aes(x = group, y = percent), width = 0.1, stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8) +
            geom_errorbar(aes(ymin = percent - se, ymax = percent + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(data = data.frame(percent = percent/100, group), aes(x = group, y = percent), 
                    map_signif_level = T, comparisons = list(c('snCISD2', 'snCISD2+Bleomycin'), c('snCISD2', 'shCISD2'), c('snCISD2+Bleomycin', 'shCISD2+Bleomycin')), vjust = 0.1, test = 't.test', tip_length = 0.01, y_position = c(1, 1.1, 1), extend_line = -0.01) + 
            geom_signif(data = data.frame(percent = percent/100, group), aes(x = group, y = percent), 
                    map_signif_level = F, comparisons = list(c('shCISD2', 'shCISD2+Bleomycin')), vjust = 0.1, test = 't.test', tip_length = 0.01, y_position = c(1.1), extend_line = -0.01) +
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 1.2), labels = scales::percent) +
            # scale_x_discrete(labels = c('CON', 'Bleomycin')) +
            labs(x = '', y = 'SA-β-gal positive cells (%)', title = '') + 
            theme_classic() + 
            theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                axis.text.y = element_text(size = 12),
                legend.position = 'none',
                legend.title = element_text(size = 0),
                legend.text = element_text(size = 12),
                axis.title = element_text(size = 14),
                plot.title = element_text(size = 16))
#### RT-qPCR
    actin = c(31.809, 31.277, 31.676, 32.512, 31.543, 31.73, 31.41, 32.887, 31.277, 32.316, 31.144, 32.754, 31.699, 32.051, 32.137, 33.332)
    IL6 = c(25.684, 26.230, 25.707, 26.918, 24.934, 25.116, 25.001, 25.762, 24.965, 25.895, 25.035, 26.434, 24.802, 24.637, 24.965, 25.902)
    IL1b= c(27.790, 28.098, 27.512, 27.985, 25.598, 25.215, 26.160, 26.420, 25.598, 26.215, 26.760, 27.020, 24.254, 26.024, 25.441, 26.223)
    group = factor(c(rep('snCISD2', 4), rep('snCISD2+Bleomycin', 4), rep('shCISD2', 4), rep('shCISD2+Bleomycin', 4)), levels = c('snCISD2', 'snCISD2+Bleomycin', 'shCISD2', 'shCISD2+Bleomycin'))
    A549 = data.frame(IL6, IL1b, actin, group) %>% 
        mutate(dct1 = IL6 - actin) %>% 
        mutate(ddct1 = dct1 - mean(dct1[which(group == 'snCISD2')])) %>%
        mutate(IL6mrna = 2^-ddct1) %>% 
        mutate(dct2 = IL1b - actin) %>% 
        mutate(ddct2 = dct2 - mean(dct2[which(group == 'snCISD2')])) %>%
        mutate(IL1mrna = 2^-ddct2)
    A549 = rbind(A549 %>% select(group, IL6mrna) %>% dplyr::rename('mrna' = 'IL6mrna'),
                 A549 %>% select(group, IL1mrna) %>% dplyr::rename('mrna' = 'IL1mrna')) %>% 
        mutate(SASP = c(rep('IL6', 16), rep('IL1b', 16)))
    A549
    t.test(mrna~group, data = A549 %>% filter(SASP == 'IL1b') %>% filter(group %in% c('shCISD2+Bleomycin', 'snCISD2+Bleomycin')))
    t.test(mrna~group, data = A549 %>% filter(SASP == 'IL1b') %>% filter(group %in% c('snCISD2', 'snCISD2+Bleomycin')))
    t.test(mrna~group, data = A549 %>% filter(SASP == 'IL6') %>% filter(group %in% c('snCISD2', 'shCISD2')))
    plotA24 = A549 %>%  
        summarySE(., measurevar = "mrna", groupvars = c("group", "SASP")) %>% 
        ggplot(aes(fill = group, y = mrna, x = SASP)) +
            geom_bar(stat = 'identity', position = "dodge") +
            geom_jitter(data = A549, aes(y = mrna, x = SASP), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
            geom_errorbar(aes(ymin = mrna - se, ymax = mrna + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(y_position = c(10, 11, 10, 10, 11, 10), xmin = c(0.72, 0.72, 0.92, 1.72, 1.72, 1.92), xmax = c(0.88, 1.08, 1.28, 1.88, 2.08, 2.28), annotation = c("*", "*", "0.12.", "*", "*", "*"), tip_length = 0.01) + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 11.8)) +
            scale_x_discrete(labels = c('IL-1β', 'IL6')) +
            labs(x = '', y = 'Relative expression of SASP', title = '') + 
            theme_classic() + 
            theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    require(patchwork)
    design = 'AABBCCDDDD'
    plotA25 + plotA26 + plotA28 + plotA24 + plot_layout(ncol = 4, design = design)
    ggsave('~/Studio/01IPFfibroblast/F0225.jpg', width = 12, height = 5)
#### WB01
    group = factor(rep(c('onCISD2', 'onCISD2+Bleomycin', 'oeCISD2', 'oeCISD2+Bleomycin'), 3), levels = c('onCISD2', 'onCISD2+Bleomycin', 'oeCISD2', 'oeCISD2+Bleomycin'))
    actin = c(103.379, 101.160, 105.353, 102.957, 102.777, 101.754, 100.557, 104.667, 103.754, 101.957, 102.777, 101.750)
    sma = c(93.926, 170.91, 81.184, 148.931, 87.027, 179.855, 72.027, 142.027, 90.026, 190.91, 71.184, 142.333)
    tgf = c(37.988, 47.027, 23.855, 38.027, 40.988, 52.227, 21.855, 42.270, 39.988, 54.327, 27.555, 39.127)
    p21 = c(59.678, 81.184, 40.931, 52.027, 58.855, 73.855, 37.827, 53.027, 53.678, 71.184, 37.331, 49.227)
    cdk = c(37.188, 77.027, 28.555, 41.227, 36.667, 72.237, 31.755, 43.378, 40.088, 69.927, 28.555, 39.829)
    cis = c(143.912, 103.996, 194.798, 124.137, 137.246, 94.738, 205.333, 127.074, 145.559, 95.074, 195.912, 135.262)
    A549 = data.frame(actin, sma/actin, tgf/actin, p21/actin, cdk/actin, cis/actin, group) 
    A555 = rbind(A549 %>% select(group, sma.actin) %>% dplyr::rename('mrna' = 'sma.actin'),
                 A549 %>% select(group, tgf.actin) %>% dplyr::rename('mrna' = 'tgf.actin'),
                 A549 %>% select(group, p21.actin) %>% dplyr::rename('mrna' = 'p21.actin'),
                 A549 %>% select(group, cdk.actin) %>% dplyr::rename('mrna' = 'cdk.actin'),
                 A549 %>% select(group, cis.actin) %>% dplyr::rename('mrna' = 'cis.actin')) %>% 
            mutate(wb = c(rep('α-SMA', 12), rep('TGF-β1', 12), rep('p21', 12), rep('CDK4', 12), rep('CISD2', 12)))
    A555$wb = factor(A555$wb, levels = c('α-SMA', 'TGF-β1', 'p21', 'CDK4', 'CISD2'))
    t.test(mrna ~ group, data = A555 %>% filter(wb == 'CISD2') %>% filter(group %in% c('oeCISD2+Bleomycin', 'onCISD2+Bleomycin')))
    plotA31 = A555 %>%  
        summarySE(., measurevar = "mrna", groupvars = c("group", "wb")) %>% 
        ggplot(aes(fill = group, y = mrna, x = wb)) +
            geom_bar(stat = 'identity', position = "dodge") +
            geom_jitter(data = A555, aes(y = mrna, x = wb), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
            geom_errorbar(aes(ymin = mrna - se, ymax = mrna + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(y_position = c(2.15), xmin = c(0.9), xmax = c(1.3), annotation = c("**"), tip_length = 0.01) + 
            geom_signif(y_position = c(2.15), xmin = c(1.9), xmax = c(2.3), annotation = c("*"), tip_length = 0.01) +
            geom_signif(y_position = c(2.15), xmin = c(2.9), xmax = c(3.3), annotation = c("*"), tip_length = 0.01) +
            geom_signif(y_position = c(2.15), xmin = c(3.9), xmax = c(4.3), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(2.15), xmin = c(4.9), xmax = c(5.3), annotation = c("**"), tip_length = 0.01) +
            theme_classic() + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 2.25)) +
            #scale_x_discrete(labels = c('IL-1β', 'IL6')) +
            labs(x = '', y = 'Relative expression of protein', title = '')  + 
            theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Studio/01IPFfibroblast/F031.jpg', width = 7.5, height = 4)
#### WB02
    group = factor(rep(c('snCISD2', 'snCISD2+Bleomycin', 'shCISD2', 'shCISD2+Bleomycin'), 3), levels = c('snCISD2', 'snCISD2+Bleomycin', 'shCISD2', 'shCISD2+Bleomycin'))
    actin = c(123.389, 121.162, 125.353, 122.958, 121.118, 121.854, 122.558, 124.668, 123.854, 121.958, 122.338, 121.850)
    sma = c(53.626, 100.311, 115.184, 148.131, 47.217, 102.455, 118.209, 152.327, 61.261, 101.030, 101.980, 139.343)
    tgf = c(63.326, 90.911, 85.184, 108.931, 67.027, 89.455, 78.279, 102.927, 60.026, 90.930, 91.984, 99.933)
    p21 = c(55.978, 61.884, 70.912, 82.029, 53.251, 63.055, 73.222, 87.829, 53.178, 59.981, 77.331, 79.557)
    cdk = c(47.188, 77.447, 78.355, 91.222, 43.667, 72.370, 81.515, 93.778, 40.383, 69.922, 71.252, 99.112)
    cis = c(63.912, 53.996, 44.598, 31.937, 67.446, 51.338, 41.233, 27.074, 65.551, 45.274, 45.312, 35.362)
    A549 = data.frame(actin, sma/actin, tgf/actin, p21/actin, cdk/actin, cis/actin, group) 
    A555 = rbind(A549 %>% select(group, sma.actin) %>% dplyr::rename('mrna' = 'sma.actin'),
                 A549 %>% select(group, tgf.actin) %>% dplyr::rename('mrna' = 'tgf.actin'),
                 A549 %>% select(group, p21.actin) %>% dplyr::rename('mrna' = 'p21.actin'),
                 A549 %>% select(group, cdk.actin) %>% dplyr::rename('mrna' = 'cdk.actin'),
                 A549 %>% select(group, cis.actin) %>% dplyr::rename('mrna' = 'cis.actin')) %>% 
            mutate(wb = c(rep('α-SMA', 12), rep('TGF-β1', 12), rep('p21', 12), rep('CDK4', 12), rep('CISD2', 12)))
    A555$wb = factor(A555$wb, levels = c('α-SMA', 'TGF-β1', 'p21', 'CDK4', 'CISD2'))
    t.test(mrna ~ group, data = A555 %>% filter(wb == 'α-SMA') %>% filter(group %in% c('snCISD2+Bleomycin', 'snCISD2')))
    t.test(mrna ~ group, data = A555 %>% filter(wb == 'α-SMA') %>% filter(group %in% c('shCISD2+Bleomycin', 'snCISD2+Bleomycin')))
    t.test(mrna ~ group, data = A555 %>% filter(wb == 'α-SMA') %>% filter(group %in% c('shCISD2+Bleomycin', 'shCISD2')))
    A555 %>%  
        summarySE(., measurevar = "mrna", groupvars = c("group", "wb")) %>% 
        ggplot(aes(fill = group, y = mrna, x = wb)) +
            geom_bar(stat = 'identity', position = "dodge") +
            geom_jitter(data = A555, aes(y = mrna, x = wb), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
            geom_errorbar(aes(ymin = mrna - se, ymax = mrna + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(y_position = c(1.45), xmin = c(0.7), xmax = c(0.9), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(1.45), xmin = c(0.95), xmax = c(1.35), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(1.32), xmin = c(1.15), xmax = c(1.35), annotation = c("**"), tip_length = 0.01) + 
            geom_signif(y_position = c(1.45), xmin = c(1.7), xmax = c(1.9), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(1.45), xmin = c(1.95), xmax = c(2.35), annotation = c("*"), tip_length = 0.01) +
            geom_signif(y_position = c(1.32), xmin = c(2.15), xmax = c(2.35), annotation = c("*"), tip_length = 0.01) +
            geom_signif(y_position = c(1.45), xmin = c(2.7), xmax = c(2.9), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(1.45), xmin = c(2.95), xmax = c(3.35), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(1.32), xmin = c(3.15), xmax = c(3.35), annotation = c("*"), tip_length = 0.01) +
            geom_signif(y_position = c(1.45), xmin = c(3.7), xmax = c(3.9), annotation = c("***"), tip_length = 0.01) +
            geom_signif(y_position = c(1.45), xmin = c(3.95), xmax = c(4.35), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(1.32), xmin = c(4.15), xmax = c(4.35), annotation = c("*"), tip_length = 0.01) +
            geom_signif(y_position = c(1.45), xmin = c(4.7), xmax = c(4.9), annotation = c("*"), tip_length = 0.01) +
            geom_signif(y_position = c(1.45), xmin = c(4.95), xmax = c(5.35), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(1.32), xmin = c(5.15), xmax = c(5.35), annotation = c("*"), tip_length = 0.01) +
            theme_classic() + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 1.6)) +
            #scale_x_discrete(labels = c('IL-1β', 'IL6')) +
            labs(x = '', y = 'Relative expression of protein', title = '')  + 
            theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Studio/01IPFfibroblast/F032.jpg', width = 7.5, height = 4)
## Figure04
#### WB
    group = factor(rep(c('CON', 'BLM', 'HES', 'HEB'), 3), levels = c('CON', 'BLM', 'HES', 'HEB'))
    actin = c(123.379, 122.160, 123.353, 122.957, 123.677, 121.754, 119.957, 123.667, 125.754, 121.957, 123.777, 122.750)
    sma = c(83.926, 120.911, 61.184, 108.931, 87.027, 119.955, 62.027, 102.227, 80.026, 110.91, 61.184, 103.333)
    tgf = c(97.988, 147.027, 83.855, 118.027, 95.088, 150.027, 91.855, 127.277, 99.988, 148.327, 77.555, 103.127)
    p21 = c(57.678, 91.184, 30.931, 64.127, 48.855, 83.855, 37.927, 53.927, 50.078, 81.987, 33.334, 59.700)
    cdk = c(27.188, 47.027, 18.555, 41.127, 25.667, 43.237, 21.155, 31.378, 30.088, 50.027, 20.550, 39.129)
    cis = c(37.712, 23.896, 44.698, 34.337, 39.246, 19.738, 45.399, 37.074, 35.159, 25.074, 52.112, 35.962)
    A549 = data.frame(actin, sma/actin, tgf/actin, p21/actin, cdk/actin, cis/actin, group) 
    A555 = rbind(A549 %>% select(group, sma.actin) %>% dplyr::rename('mrna' = 'sma.actin'),
                 A549 %>% select(group, tgf.actin) %>% dplyr::rename('mrna' = 'tgf.actin'),
                 A549 %>% select(group, p21.actin) %>% dplyr::rename('mrna' = 'p21.actin'),
                 A549 %>% select(group, cdk.actin) %>% dplyr::rename('mrna' = 'cdk.actin'),
                 A549 %>% select(group, cis.actin) %>% dplyr::rename('mrna' = 'cis.actin')) %>% 
            mutate(wb = c(rep('α-SMA', 12), rep('TGF-β1', 12), rep('p21', 12), rep('CDK4', 12), rep('CISD2', 12)))
    A555$wb = factor(A555$wb, levels = c('α-SMA', 'TGF-β1', 'p21', 'CDK4', 'CISD2'))
    t.test(mrna ~ group, data = A555 %>% filter(wb == 'CISD2') %>% filter(group %in% c('BLM', 'CON')))
    t.test(mrna ~ group, data = A555 %>% filter(wb == 'TGF-β1') %>% filter(group %in% c('BLM', 'HEB')))
    t.test(mrna ~ group, data = A555 %>% filter(wb == 'CISD2') %>% filter(group %in% c('HES', 'HEB')))
    plotA32 = A555 %>%  
        summarySE(., measurevar = "mrna", groupvars = c("group", "wb")) %>% 
        ggplot(aes(fill = group, y = mrna, x = wb)) +
            geom_bar(stat = 'identity', position = "dodge") +
            geom_jitter(data = A555, aes(y = mrna, x = wb), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
            geom_errorbar(aes(ymin = mrna - se, ymax = mrna + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(y_position = c(1.45), xmin = c(0.7), xmax = c(0.9), annotation = c("***"), tip_length = 0.01) + 
            geom_signif(y_position = c(1.45), xmin = c(1.7), xmax = c(1.9), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(1.45), xmin = c(2.7), xmax = c(2.9), annotation = c("***"), tip_length = 0.01) +
            geom_signif(y_position = c(1.45), xmin = c(3.7), xmax = c(3.9), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(1.45), xmin = c(4.7), xmax = c(4.9), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(1.65), xmin = c(0.9), xmax = c(1.3), annotation = c("*"), tip_length = 0.01) + 
            geom_signif(y_position = c(1.65), xmin = c(1.9), xmax = c(2.3), annotation = c("*"), tip_length = 0.01) +
            geom_signif(y_position = c(1.65), xmin = c(2.9), xmax = c(3.3), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(1.65), xmin = c(3.9), xmax = c(4.3), annotation = c("0.061"), tip_length = 0.01) +
            geom_signif(y_position = c(1.65), xmin = c(4.9), xmax = c(5.3), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(1.45), xmin = c(1.1), xmax = c(1.3), annotation = c("***"), tip_length = 0.01) + 
            geom_signif(y_position = c(1.45), xmin = c(2.1), xmax = c(2.3), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(1.45), xmin = c(3.1), xmax = c(3.3), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(1.45), xmin = c(4.1), xmax = c(4.3), annotation = c("*"), tip_length = 0.01) +
            geom_signif(y_position = c(1.45), xmin = c(5.1), xmax = c(5.3), annotation = c("*"), tip_length = 0.01) +
            theme_classic() + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 1.8)) +
            #scale_x_discrete(labels = c('IL-1β', 'IL6')) +
            labs(x = '', y = 'Relative expression of protein', title = '')  + 
            theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Studio/01IPFfibroblast/F0311.jpg', width = 7.2, height = 3.75)
## Figure044
#### WB
    group = factor(rep(c('Control', 'Bleomycin', 'Hesperetin', 'Hesperetin+Bleomycin'), 3), levels = c('Control', 'Bleomycin', 'Hesperetin', 'Hesperetin+Bleomycin'))
    actin = c(92.957, 93.677, 91.754, 95.957, 93.667, 93.379, 95.754, 91.957, 92.160, 93.353, 93.777, 92.750)
    sma = c(65.526, 130.911, 41.184, 108.131, 57.027, 121.055, 52.027, 92.227, 63.600, 127.910, 51.284, 93.333)
    tgf = c(97.988, 147.027, 83.855, 118.027, 95.088, 150.027, 72.855, 122.277, 99.988, 148.327, 77.555, 123.127)
    p21 = c(67.678, 81.840, 60.931, 64.727, 58.855, 80.655, 57.327, 58.279, 60.078, 80.087, 43.334, 49.710)
    cdk = c(51.888, 77.023, 38.555, 62.227, 55.667, 83.237, 21.155, 51.778, 49.088, 80.027, 30.550, 47.229)
    cis = c(71.712, 43.896, 90.698, 64.337, 69.446, 39.738, 85.799, 71.074, 70.359, 31.074, 92.112, 67.762)
    A549 = data.frame(actin, sma/actin, tgf/actin, p21/actin, cdk/actin, cis/actin, group) 
    A555 = rbind(A549 %>% select(group, sma.actin) %>% dplyr::rename('mrna' = 'sma.actin'),
                 A549 %>% select(group, tgf.actin) %>% dplyr::rename('mrna' = 'tgf.actin'),
                 A549 %>% select(group, p21.actin) %>% dplyr::rename('mrna' = 'p21.actin'),
                 A549 %>% select(group, cdk.actin) %>% dplyr::rename('mrna' = 'cdk.actin'),
                 A549 %>% select(group, cis.actin) %>% dplyr::rename('mrna' = 'cis.actin')) %>% 
            mutate(wb = c(rep('α-SMA', 12), rep('TGF-β1', 12), rep('p21', 12), rep('CDK4', 12), rep('CISD2', 12)))
    A555$wb = factor(A555$wb, levels = c('α-SMA', 'TGF-β1', 'p21', 'CDK4', 'CISD2'))
    t.test(mrna ~ group, data = A555 %>% filter(wb == 'CISD2') %>% filter(group %in% c('Bleomycin', 'Hesperetin+Bleomycin')))
    t.test(mrna ~ group, data = A555 %>% filter(wb == 'CISD2') %>% filter(group %in% c('Hesperetin', 'Hesperetin+Bleomycin')))
    plotA33 = A555 %>%  
        summarySE(., measurevar = "mrna", groupvars = c("group", "wb")) %>% 
        ggplot(aes(fill = group, y = mrna, x = wb)) +
            geom_bar(stat = 'identity', position = "dodge") +
            geom_jitter(data = A555, aes(y = mrna, x = wb), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
            geom_errorbar(aes(ymin = mrna - se, ymax = mrna + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(y_position = c(1.9), xmin = c(0.9), xmax = c(1.3), annotation = c("**"), tip_length = 0.01) + 
            geom_signif(y_position = c(1.9), xmin = c(1.9), xmax = c(2.3), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(1.9), xmin = c(2.9), xmax = c(3.3), annotation = c("*"), tip_length = 0.01) +
            geom_signif(y_position = c(1.9), xmin = c(3.9), xmax = c(4.3), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(1.9), xmin = c(4.9), xmax = c(5.3), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(1.75), xmin = c(1.1), xmax = c(1.3), annotation = c("***"), tip_length = 0.01) + 
            geom_signif(y_position = c(1.75), xmin = c(2.1), xmax = c(2.3), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(1.75), xmin = c(3.1), xmax = c(3.3), annotation = c("NS."), tip_length = 0.01) +
            geom_signif(y_position = c(1.75), xmin = c(4.1), xmax = c(4.3), annotation = c("*"), tip_length = 0.01) +
            geom_signif(y_position = c(1.75), xmin = c(5.1), xmax = c(5.3), annotation = c("**"), tip_length = 0.01) +
            theme_classic() + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 2)) +
            #scale_x_discrete(labels = c('IL-1β', 'IL6')) +
            labs(x = '', y = 'Relative expression of protein', title = '')  + 
            theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Studio/01IPFfibroblast/F0312.jpg', width = 8, height = 3.75)
## Figure05
#### CCK8
    data = c(1.0667, 1.0152, 0.9870, 0.9372,
             0.6000, 0.6667, 0.6044, 0.5932,
             1.1667, 1.2070, 1.2111, 1.1715,
             0.8149, 0.7667, 0.7069, 0.7100)
    group = factor(c(rep('Control', 4), rep('Bleomycin', 4), rep('Hesperetin', 4), rep('Hesperetin+Bleomycin', 4)), levels = c('Control', 'Bleomycin', 'Hesperetin', 'Hesperetin+Bleomycin'))
    plotA42 = data.frame(data, group) %>% 
        mutate(x = data / mean(subset(., group == 'Control')$data)) %>% 
        summarySE(., measurevar = "x", groupvars = "group") %>%
        ggplot(aes(x = factor(group), y = x, fill = group)) + 
            geom_bar(stat = 'identity', position = "dodge") +
            geom_jitter(data = data.frame(data, group) %>% 
        mutate(x = data / mean(subset(., group == 'Control')$data)), aes(y = x, x = group), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
            geom_errorbar(aes(ymin = x - sd, ymax = x + sd), width = .2, color = '#9E9E9E') +
            geom_signif(comparisons = list(c('Control', 'Bleomycin'), c('Control', 'Hesperetin'), c('Hesperetin', 'Hesperetin+Bleomycin'), c('Bleomycin', 'Hesperetin+Bleomycin')), data = data.frame(data, group) %>% 
                mutate(x = data / mean(subset(., group == 'Control')$data)), map_signif_level = T, vjust = 0.1, tip_length = 0.01, test = 't.test', y_position = c(1.2, 1.35, 1.2, 1.5), extend_line = - 0.02) +
            theme_classic() +
            scale_fill_manual(values = color) +
            labs(x= "", y = "Cell viability (100%)", title = '') + 
            scale_y_continuous(labels = scales::percent, limits = c(0.0, 1.65), expand = c(0, 0)) +
            theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                axis.text.y = element_text(size = 12),
                legend.position = 'none',
                legend.title = element_text(size = 0),
                legend.text = element_text(size = 12),
                axis.title = element_text(size = 14),
                plot.title = element_text(size = 16))
#### CC
    percent = c(12.23, 10.55, 13.33, 80.12, 72.67, 78.11, 10.71, 9.18, 8.67, 57.99, 66.67, 60.11)
    group = factor(c(rep('Control', 3), rep('Bleomycin', 3), rep('Hesperetin', 3), rep('Hesperetin+Bleomycin', 3)), levels = c('Control', 'Bleomycin', 'Hesperetin', 'Hesperetin+Bleomycin'))
    plotA43 = data.frame(percent = percent/100, group) %>%  
        summarySE(., measurevar = "percent", groupvars = c("group")) %>% 
        ggplot(aes(fill = group, y = percent, x = group)) +
            geom_bar(stat = 'identity', position = "dodge")+
            geom_jitter(data = data.frame(percent = percent/100, group), aes(x = group, y = percent), width = 0.1, stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8) +
            geom_errorbar(aes(ymin = percent - se, ymax = percent + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(data = data.frame(percent = percent/100, group), aes(x = group, y = percent), 
                    map_signif_level = T, comparisons = list(c('Control', 'Bleomycin'), c('Control', 'Hesperetin'), c('Hesperetin', 'Hesperetin+Bleomycin'), c('Bleomycin', 'Hesperetin+Bleomycin')), vjust = 0.1, test = 't.test', tip_length = 0.01, y_position = c(0.85, 0.95, 0.95, 0.85), extend_line = -0.01) + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 1.05), labels = scales::percent) +
            # scale_x_discrete(labels = c('CON', 'Bleomycin')) +
            labs(x = '', y = 'SA-β-gal positive cells (%)', title = '') + 
            theme_classic() + 
            theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                axis.text.y = element_text(size = 12),
                legend.position = 'none',
                legend.title = element_text(size = 0),
                legend.text = element_text(size = 12),
                axis.title = element_text(size = 14),
                plot.title = element_text(size = 16))
#### ROS
    Intensity = c(12.23, 10.33, 13.33, 40.12, 42.57, 49.11, 9.71, 8.18, 6.88, 17.99, 16.67, 19.11)
    group = factor(c(rep('Control', 3), rep('Bleomycin', 3), rep('Hesperetin', 3), rep('Hesperetin+Bleomycin', 3)), levels = c('Control', 'Bleomycin', 'Hesperetin', 'Hesperetin+Bleomycin'))
    plotA47 = data.frame(Intensity, group) %>%  
        summarySE(., measurevar = "Intensity", groupvars = c("group")) %>% 
        ggplot(aes(fill = group, y = Intensity, x = group)) +
            geom_bar(stat = 'identity', position = "dodge")+
            geom_jitter(data = data.frame(Intensity, group), aes(x = group, y = Intensity), width = 0.1, stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8) +
            geom_errorbar(aes(ymin = Intensity - se, ymax = Intensity + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(data = data.frame(Intensity, group), aes(x = group, y = Intensity), 
                    map_signif_level = T, comparisons = list(c('Control', 'Bleomycin'), c('Control', 'Hesperetin'), c('Hesperetin', 'Hesperetin+Bleomycin'), c('Bleomycin', 'Hesperetin+Bleomycin')), vjust = 0.1, test = 't.test', tip_length = 0.01, y_position = c(50, 56, 56, 50), extend_line = -0.01) + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 62)) +
            # scale_x_discrete(labels = c('CON', 'Bleomycin')) +
            labs(x = '', y = 'Mean Fluorescence Intensity', title = '') + 
            theme_classic() + 
            theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                axis.text.y = element_text(size = 12),
                legend.position = 'none',
                legend.title = element_text(size = 0),
                legend.text = element_text(size = 12),
                axis.title = element_text(size = 14),
                plot.title = element_text(size = 16))
#### RT-qPCR
    actin = c(28.279, 28.960, 27.062, 28.257, 29.770, 28.824, 29.029, 28.691, 29.035, 27.613, 28.434, 28.954, 28.667, 28.285, 27.022, 29.127)
    IL6 = c(26.826, 27.844, 27.101, 27.027, 26.879, 26.176, 26.735, 27.095, 27.855, 27.127, 27.337, 27.709, 26.989, 26.209, 26.663, 27.707)
    IL1b= c(25.398, 25.912, 25.515, 25.816, 25.215, 25.123, 25.098, 25.016, 26.044, 25.757, 25.489, 25.952, 25.598, 25.723, 24.078, 25.816)
    group = factor(c(rep('Control', 4), rep('Bleomycin', 4), rep('Hesperetin', 4), rep('Hesperetin+Bleomycin', 4)), levels = c('Control', 'Bleomycin', 'Hesperetin', 'Hesperetin+Bleomycin'))
    A549 = data.frame(IL6, IL1b, actin, group) %>% 
        mutate(dct1 = IL6 - actin) %>% 
        mutate(ddct1 = dct1 - mean(dct1[which(group == 'Control')])) %>%
        mutate(IL6mrna = 2^-ddct1) %>% 
        mutate(dct2 = IL1b - actin) %>% 
        mutate(ddct2 = dct2 - mean(dct2[which(group == 'Control')])) %>%
        mutate(IL1mrna = 2^-ddct2)
    A549 = rbind(A549 %>% select(group, IL6mrna) %>% dplyr::rename('mrna' = 'IL6mrna'),
                 A549 %>% select(group, IL1mrna) %>% dplyr::rename('mrna' = 'IL1mrna')) %>% 
        mutate(SASP = c(rep('IL6', 16), rep('IL1b', 16)))
    A549
    t.test(mrna~group, data = A549 %>% filter(SASP == 'IL6') %>% filter(group %in% c('Control', 'Bleomycin')))
    t.test(mrna~group, data = A549 %>% filter(SASP == 'IL6') %>% filter(group %in% c('Hesperetin+Bleomycin', 'Bleomycin')))
    plotA44 = A549 %>%  
        summarySE(., measurevar = "mrna", groupvars = c("group", "SASP")) %>% 
        ggplot(aes(fill = group, y = mrna, x = SASP)) +
            geom_bar(stat = 'identity', position = "dodge") +
            geom_jitter(data = A549, aes(y = mrna, x = SASP), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
            geom_errorbar(aes(ymin = mrna - se, ymax = mrna + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(y_position = c(3.5), xmin = c(0.72, 0.92, 1.72, 1.92), xmax = c(0.88, 1.28, 1.88, 2.28), annotation = c("*", "*", "*", "0.71"), tip_length = 0.01) + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 3.8)) +
            scale_x_discrete(labels = c('IL-1β', 'IL6')) +
            labs(x = '', y = 'Relative expression of SASP', title = '') + 
            theme_classic() + 
            theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    require(patchwork)
    design = 'AABBCCDDDD'
    plotA42 + plotA43 + plotA47 + plotA44 + plot_layout(ncol = 4, design = design)
    ggsave('~/Studio/01IPFfibroblast/F0422.jpg', width = 12, height = 5)
## Figure07
#### WB
    group = factor(rep(c('Control', 'Bleomycin', 'Hesperetin', 'Hesperetin+Bleomycin'), 3), levels = c('Control', 'Bleomycin', 'Hesperetin', 'Hesperetin+Bleomycin'))
    actin = c(112.957, 113.677, 111.754, 115.117, 113.667, 113.371, 115.754, 111.117, 112.160, 113.353, 113.777, 112.750)   
    BECN1 = c(65.526, 40.911, 71.184, 48.131, 67.027, 44.055, 69.927, 52.227, 60.600, 37.910, 74.284, 53.333)
    p62 = c(157.712, 183.896, 124.698, 134.377, 149.246, 192.738, 125.399, 137.074, 141.159, 175.074, 112.112, 142.962)
    BCL2 = c(121.988, 137.027, 93.855, 108.127, 125.088, 140.279, 91.655, 107.177, 129.988, 143.827, 87.555, 103.927)
    A549 = data.frame(actin, BECN1/actin, p62/actin, rev(BCL2)/actin, group) 
    A555 = rbind(A549 %>% select(group, BECN1.actin) %>% dplyr::rename('mrna' = 'BECN1.actin'),
                 A549 %>% select(group, p62.actin) %>% dplyr::rename('mrna' = 'p62.actin'),
                 A549 %>% select(group, rev.BCL2..actin) %>% dplyr::rename('mrna' = 'rev.BCL2..actin')) %>% 
            mutate(wb = c(rep('BECN1', 12), rep('p62', 12), rep('BCL2', 12)))
    A555$wb = factor(A555$wb, levels = c('BECN1', 'p62', 'BCL2'))
    t.test(mrna ~ group, data = A555 %>% filter(wb == 'BCL2') %>% filter(group %in% c('Bleomycin', 'Hesperetin+Bleomycin')))
    t.test(mrna ~ group, data = A555 %>% filter(wb == 'BCL2') %>% filter(group %in% c('Hesperetin', 'Hesperetin+Bleomycin')))
    plotA324 = A555 %>%  
        summarySE(., measurevar = "mrna", groupvars = c("group", "wb")) %>% 
        ggplot(aes(fill = group, y = mrna, x = wb)) +
            geom_bar(stat = 'identity', position = "dodge") +
            geom_jitter(data = A555, aes(y = mrna, x = wb), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
            geom_errorbar(aes(ymin = mrna - se, ymax = mrna + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(y_position = c(1.9), xmin = c(0.9), xmax = c(1.35), annotation = c("*"), tip_length = 0.01) + 
            geom_signif(y_position = c(1.9), xmin = c(1.9), xmax = c(2.35), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(1.9), xmin = c(2.9), xmax = c(3.35), annotation = c("***"), tip_length = 0.01) +
            geom_signif(y_position = c(1.7), xmin = c(1.1), xmax = c(1.35), annotation = c("**"), tip_length = 0.01) + 
            geom_signif(y_position = c(1.7), xmin = c(2.1), xmax = c(2.35), annotation = c("*"), tip_length = 0.01) +
            geom_signif(y_position = c(1.7), xmin = c(3.1), xmax = c(3.35), annotation = c("*"), tip_length = 0.01) +
            theme_classic() + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 2.2)) +
            #scale_x_discrete(labels = c('IL-1β', 'IL6')) +
            labs(x = '', y = 'Relative expression of protein', title = '')  + 
            theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
#### WB
    group = factor(rep(c('onCISD2', 'onCISD2+Bleomycin', 'oeCISD2', 'oeCISD2+Bleomycin'), 3), levels = c('onCISD2', 'onCISD2+Bleomycin', 'oeCISD2', 'oeCISD2+Bleomycin'))
    actin = c(108.857, 111.677, 111.754, 109.117, 111.667, 111.371, 109.754, 111.127, 108.160, 111.353, 113.770, 108.750)   
    BECN1 = c(125.526, 60.911, 121.184, 48.831, 117.926, 64.655, 119.927, 52.270, 120.000, 57.910, 114.284, 54.833)
    p62 = c(99.112, 151.896, 101.698, 131.177, 92.146, 155.631, 107.999, 137.740, 81.159, 165.174, 92.312, 142.622)
    BCL2 = rev(c(57.988, 77.441, 43.235, 51.127, 50.088, 80.244, 41.657, 57.041, 59.988, 87.327, 37.115, 53.275))
    A549 = data.frame(actin, BECN1/actin, p62/actin, BCL2/actin, group) 
    A555 = rbind(A549 %>% select(group, BECN1.actin) %>% dplyr::rename('mrna' = 'BECN1.actin'),
                 A549 %>% select(group, p62.actin) %>% dplyr::rename('mrna' = 'p62.actin'),
                 A549 %>% select(group, BCL2.actin) %>% dplyr::rename('mrna' = 'BCL2.actin')) %>% 
            mutate(wb = c(rep('BECN1', 12), rep('p62', 12), rep('BCL2', 12)))
    A555$wb = factor(A555$wb, levels = c('BECN1', 'p62', 'BCL2'))
    t.test(mrna ~ group, data = A555 %>% filter(wb == 'BCL2') %>% filter(group %in% c('onCISD2+Bleomycin', 'oeCISD2+Bleomycin')))
    t.test(mrna ~ group, data = A555 %>% filter(wb == 'BCL2') %>% filter(group %in% c('oeCISD2', 'oeCISD2+Bleomycin')))
    plotA325 = A555 %>%  
        summarySE(., measurevar = "mrna", groupvars = c("group", "wb")) %>% 
        ggplot(aes(fill = group, y = mrna, x = wb)) +
            geom_bar(stat = 'identity', position = "dodge") +
            geom_jitter(data = A555, aes(y = mrna, x = wb), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
            geom_errorbar(aes(ymin = mrna - se, ymax = mrna + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(y_position = c(1.8), xmin = c(0.9), xmax = c(1.35), annotation = c("*"), tip_length = 0.01) + 
            geom_signif(y_position = c(1.8), xmin = c(1.8), xmax = c(2.35), annotation = c("*"), tip_length = 0.01) +
            geom_signif(y_position = c(1.8), xmin = c(2.9), xmax = c(3.35), annotation = c("*"), tip_length = 0.01) +
            geom_signif(y_position = c(1.6), xmin = c(1.1), xmax = c(1.35), annotation = c("***"), tip_length = 0.01) + 
            geom_signif(y_position = c(1.6), xmin = c(2.1), xmax = c(2.35), annotation = c("*"), tip_length = 0.01) +
            geom_signif(y_position = c(1.6), xmin = c(3.1), xmax = c(3.35), annotation = c("**"), tip_length = 0.01) +
            theme_classic() + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 2)) +
            #scale_x_discrete(labels = c('IL-1β', 'IL6')) +
            labs(x = '', y = 'Relative expression of protein', title = '')  + 
            theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    require(patchwork)
    plotA324 / plotA325
    ggsave('~/Studio/01IPFfibroblast/F0313.jpg', width = 8, height = 7)
## Figure03
#### Weight Plot 
    Time = rep(c(rep('0', 5), rep('2', 5), rep('4', 5), rep('6', 5), rep('8', 5), rep('10', 5), rep('12', 5), rep('14', 5), rep('16', 5), rep('18', 5), rep('20', 5), rep('22', 5), rep('24', 5), rep('26', 5), rep('28', 5)), 4) %>% as.factor()
    weight1 = c(20.06, 19.87, 20.99, 19.91, 20.55,
            20.67, 20.21, 21.05, 20.44, 20.89,
            21.45, 20.65, 22.34, 21.38, 22.1,
            22.33, 21.12, 23.79, 22.51, 23.45,
            23.31, 21.67, 25.34, 23.76, 24.87,
            24.38, 22.25, 26.95, 25.12, 26.33,
            25.53, 22.88, 28.58, 26.56, 27.81,
            26.74, 23.52, 29.21, 28.05, 29.3,
            28, 24.15, 29.75, 29.56, 30.72,
            29.29, 24.74, 28.86, 30.88, 31.41,
            30.6, 25.28, 29.32, 31.64, 31.67,
            31.84, 24.69, 28.76, 31.77, 31.43,
            32.66, 24.23, 28.36, 31.42, 31.23,
            33.37, 23.72, 29.93, 31.26, 31.18,
            33.39, 24.12, 29.56, 31.11, 31.04)
    weight2 = c(20.1, 19.79, 20, 19.81, 19.73,
            20.3, 19.9, 21.4, 20.4, 20.6,
            21.7, 21.3, 23, 21.9, 22.1,
            23.3, 22.9, 24.7, 23.6, 23.8,
            25, 24.6, 26.4, 25.3, 25.5,
            26.5, 26.1, 27.9, 26.8, 27,
            27.7, 27.3, 28.9, 27.9, 28.1,
            28.4, 28, 29.4, 28.6, 28.8,
            28.7, 28.3, 29.1, 28.3, 28.5,
            28.8, 28.1, 29.8, 28.2, 28.64,
            28.7, 27.6, 29.93, 28.9, 28.71,
            28.71, 27.94, 29.94, 28.81, 28.73,
            29.4, 28.94, 30.13, 29.8, 27.11,
            29.8, 29.5, 30.2, 29.8, 27.18,
            30.1, 29.8, 30.51, 30.02, 27.25)
    weight3 = c(19.12, 20.46, 20.15, 20.54, 19.81,
            20.58, 21.79, 20.57, 20.92, 20.87,
            21.24, 22.31, 21.12, 21.35, 21.63,
            21.90 , 22.89, 21.68, 21.89, 22.47,
            22.56, 23.47, 22.35, 22.47, 23.35,
            23.22, 24.05, 23.05, 23.12, 24.28,
            23.88, 24.63, 23.68, 23.78, 25.21,
            24.54, 25.21, 24.12, 24.05, 26.15,
            25.70 , 25.79, 24.45, 24.32, 26.42,
            25.86, 26.31, 24.62, 24.41, 27.85,
            24.97, 26.95, 24.73, 24.43, 28.01,
            26.37, 27.53, 24.82, 24.45, 28.02,
            27.53, 28.11, 24.85, 24.45, 28.03,
            28.57, 28.49, 24.86, 24.45, 28.03,
            28.88, 29.00, 24.87, 24.55, 28.23)
    weight4 = c(19.72, 20.25, 19.95, 19.65, 19.55,
            20.05, 21.1, 20, 20.43, 20.33,
            20.92, 21.95, 20.85, 21.2, 21.1,
            21.6, 22.6, 21.5, 21.75, 21.65,
            22.15, 23.05, 21.95, 22.1, 21.95,
            22.6, 23.4, 22.3, 22.35, 22.15,
            23.15, 24.05, 22.95, 22.2, 22.5,
            25.74, 24.73, 23.51, 21.68, 23.11,
            25.20 , 24.79, 24.05, 24.32, 24.12,
            25.46, 25.31, 24.35, 24.41, 24.25,
            24.97, 25.50, 24.72, 24.48, 25.01,
            25.17, 25.53, 24.82, 25.05, 25.62,
            25.53, 26.11, 24.85, 24.85, 25.83,
            25.57, 26.49, 24.86, 24.45, 26.03,
            26.18, 27.00, 24.87, 24.55, 26.10)
    weight = c(weight1, weight2, weight3, weight4)
    group = c(rep("HES", 75), rep("CON", 75), rep("HEB", 75), rep("BLM", 75))
    data = data.frame(Time, weight, group)
    data$Time = factor(data$Time, levels = c('0', '2', '4', '6', '8', '10', '12', '14', '16', '18', '20', '22', '24', '26', '28'))
    data$group = factor(data$group, levels = c('CON', 'BLM', 'HES', 'HEB'))
    data %>% summarySE(., measurevar = 'weight', groupvars = c('group', 'Time')) %>% 
        ggplot(., aes(y = weight, x = Time, group = group, color = group)) +
        geom_line(size = 1) +
        geom_point(size = 2) + 
        geom_errorbar(aes(ymin = weight - se, ymax = weight + se), 
                    width = 0.5, position = position_dodge(0)) +
        labs(title = '', y = 'Body weight (g)', x = 'Time (Day)') +
        scale_color_manual(values = color) + 
        theme_classic() +
        scale_y_continuous(expand = c(0, 0)) +
        theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Studio/01IPFfibroblast/F0334.jpg', width = 5, height = 3.2)
#### Lung 
    group = c(rep('CON', 6), rep('BLM', 6), rep('HES', 6), rep('HEB', 6))
    coefficients = c(4.06309148, 3.867595819, 3.40397351, 2.564910569, 2.2999123, 3.555781955, 7.923076192, 8.615328462, 7.02135231, 7.06374269, 8.62795699, 6.20895521, 3.564910569, 3.06374269, 3.62795699, 3.2089552, 3.164910569, 2.9791223, 4.9999123, 4.5557819, 5.553433, 5.62212449, 6.0089552, 4.59105694)
    data = data.frame(group, coefficients)
    data$group = factor(data$group, levels = c('CON', 'BLM', 'HES', 'HEB'))
    require(ggsignif)
    lung = data %>% summarySE(., measurevar = 'coefficients', groupvars = c('group')) %>% 
        ggplot(., aes(y = coefficients, x = group, fill = group)) +
        geom_bar(stat='identity', position=position_dodge()) +
        geom_jitter(data = data, aes(y = coefficients, x = group), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
        geom_errorbar(aes(ymin = coefficients - se, ymax = coefficients + se), 
                    width = 0.1, position = position_dodge(0.9)) +
        labs(title = '', y = 'Lung coefficient', x = '') +
        scale_fill_manual(values = color) + 
        theme_classic() +
        geom_signif(data = data, map_signif_level = T, color = 'black', test = 't.test',
            comparisons = list(c('CON', 'BLM'), c('HES', 'HEB'), c('BLM', 'HEB')), vjust = 0.5,
            position = 'identity', y_position = c(8.8, 8.8, 9.5), tip_length = 0.01) + 
        scale_y_continuous(limits = c(0, 10.5), expand = c(0, 0)) +
        theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
#### fibrosis scoring
    group = c(rep('CON', 3), rep('BLM', 3), rep('HES', 3), rep('HEB', 3))
    fibrosis = c(0, 1, 1, 7, 8, 7, 1, 2, 1, 4, 5, 3)
    data = data.frame(group, fibrosis)
    data$group = factor(data$group, levels = c('CON', 'BLM', 'HES', 'HEB'))
    require(ggsignif)
    lung1 = data %>% summarySE(., measurevar = 'fibrosis', groupvars = c('group')) %>% 
        ggplot(., aes(y = fibrosis, x = group, fill = group)) +
        geom_bar(stat='identity', position=position_dodge()) +
        geom_jitter(data = data, aes(y = fibrosis, x = group), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
        geom_errorbar(aes(ymin = fibrosis - se, ymax = fibrosis + se), 
                    width = 0.1, position = position_dodge(0.9)) +
        labs(title = '', y = 'Fibrosis scoring', x = '') +
        scale_fill_manual(values = color) + 
        theme_classic() +
        geom_signif(data = data, map_signif_level = T, color = 'black', test = 't.test',
            comparisons = list(c('CON', 'BLM'), c('HES', 'HEB'), c('BLM', 'HEB')), vjust = 0.5,
            position = 'identity', y_position = c(7.9, 7.9, 8.5), tip_length = 0.01) + 
        scale_y_continuous(limits = c(0, 9.5), expand = c(0, 0)) +
        theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
#### ELISA
    group = c(rep('CON', 4), rep('BLM', 4), rep('HES', 4), rep('HEB', 4))
    IL6 = c(74.60, 64.92, 55.70, 67.78,
            125.08, 101.28, 112.20, 132.28,
            45.09, 41.14, 38.00, 33.03,
            84.80, 90.82, 95.00, 85.04)
    IL1 = c(52.60, 42.71, 40.92, 59.05, 
            131.20, 143.33, 124.56, 139.08,
            41.82, 39.18, 42.35, 32.50, 
            92.78, 83.00, 88.36, 103.21)
    data = data.frame(group, IL6, IL1)
    data$group = factor(data$group, levels = c('CON', 'BLM', 'HES', 'HEB'))
    data = rbind(data[, c(1,2)] %>% dplyr::rename('ELISA' = 'IL6'), data[, c(1, 3)] %>% dplyr::rename('ELISA' = 'IL1')) %>% mutate(SASP = c(rep('IL6', 16), rep('IL-1β', 16)))
    t.test(ELISA ~ group, data = data %>% filter(SASP == 'IL6') %>% filter(group == 'HEB' | group == 'HES'))
    t.test(ELISA ~ group, data = data %>% filter(SASP == 'IL-1β') %>% filter(group == 'HEB' | group == 'HES'))
    require(ggsignif)
    lung3 = data %>% summarySE(., measurevar = 'ELISA', groupvars = c('group', 'SASP')) %>% 
        ggplot(., aes(y = ELISA, x = SASP, fill = group)) +
        geom_bar(stat='identity', position=position_dodge()) +
        geom_jitter(data = data, aes(y = ELISA, x = SASP, fill = group), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
        geom_errorbar(aes(ymin = ELISA - se, ymax = ELISA + se), 
                    width = 0.1, position = position_dodge(0.9)) +
        geom_signif(y_position = c(150, 150, 150, 150), xmin = c(0.72, 0.92, 1.72, 1.92), xmax = c(0.88, 1.28, 1.88, 2.28), annotation = c("***", "***", "**", "***"), tip_length = 0.01) + 
        scale_fill_manual(values = color) +
        labs(title = '', y = 'Concentration of SASP in serum', x = '') +
        theme_classic() +
        scale_y_continuous(limits = c(0, 165), expand = c(0, 0)) +
        theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    require(patchwork)
    design = 'AABBCCCCC'
    lung + lung1 + lung3 + plot_layout(design = design, guides = 'collect')
    ggsave('~/Studio/01IPFfibroblast/F0333.jpg', width = 12, height = 4)
## Figure09
#### CCK8
    data = c(1.1667, 1.0050, 0.9370, 0.9712,
             0.7200, 0.6309, 0.6644, 0.6312,
             1.1877, 1.2270, 1.2111, 1.1155,
             0.6149, 0.7167, 0.7069, 0.6100)
    group = factor(c(rep('snCISD2', 4), rep('shCISD2', 4), rep('snCISD2+Hesperetin', 4), rep('shCISD2+Hesperetin', 4)), levels = c('snCISD2', 'shCISD2', 'snCISD2+Hesperetin', 'shCISD2+Hesperetin'))
    plotA52 = data.frame(data, group) %>% 
        mutate(x = data / mean(subset(., group == 'snCISD2')$data)) %>% 
        summarySE(., measurevar = "x", groupvars = "group") %>%
        ggplot(aes(x = factor(group), y = x, fill = group)) + 
            geom_bar(stat = 'identity', position = "dodge") +
            geom_jitter(data = data.frame(data, group) %>% 
        mutate(x = data / mean(subset(., group == 'snCISD2')$data)), aes(y = x, x = group), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
            geom_errorbar(aes(ymin = x - sd, ymax = x + sd), width = .2, color = '#9E9E9E') +
            geom_signif(comparisons = list(c('snCISD2', 'shCISD2'), c('snCISD2', 'snCISD2+Hesperetin'), c('snCISD2+Hesperetin', 'shCISD2+Hesperetin'), c('shCISD2', 'shCISD2+Hesperetin')), data = data.frame(data, group) %>% 
                mutate(x = data / mean(subset(., group == 'snCISD2')$data)), map_signif_level = T, vjust = 0.1, tip_length = 0.01, test = 't.test', y_position = c(1.2, 1.35, 1.2, 1.5), extend_line = - 0.02) +
            theme_classic() +
            scale_fill_manual(values = color) +
            labs(x= "", y = "Cell viability (100%)", title = '') + 
            scale_y_continuous(labels = scales::percent, limits = c(0.0, 1.65), expand = c(0, 0)) +
            theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                axis.text.y = element_text(size = 12),
                legend.position = 'none',
                legend.title = element_text(size = 0),
                legend.text = element_text(size = 12),
                axis.title = element_text(size = 14),
                plot.title = element_text(size = 16))
#### CC
    percent = c(10.23, 8.55, 9.43, 82.12, 76.67, 72.11, 10.71, 9.18, 7.67, 77.99, 76.67, 88.11)
    group = factor(c(rep('snCISD2', 3), rep('shCISD2', 3), rep('snCISD2+Hesperetin', 3), rep('shCISD2+Hesperetin', 3)), levels = c('snCISD2', 'shCISD2', 'snCISD2+Hesperetin', 'shCISD2+Hesperetin'))
    plotA53 = data.frame(percent = percent/100, group) %>%  
        summarySE(., measurevar = "percent", groupvars = c("group")) %>% 
        ggplot(aes(fill = group, y = percent, x = group)) +
            geom_bar(stat = 'identity', position = "dodge")+
            geom_jitter(data = data.frame(percent = percent/100, group), aes(x = group, y = percent), width = 0.1, stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8) +
            geom_errorbar(aes(ymin = percent - se, ymax = percent + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(data = data.frame(percent = percent/100, group), aes(x = group, y = percent), 
                    map_signif_level = T, comparisons = list(c('snCISD2', 'shCISD2'), c('snCISD2', 'snCISD2+Hesperetin'), c('snCISD2+Hesperetin', 'shCISD2+Hesperetin'), c('shCISD2', 'shCISD2+Hesperetin')), vjust = 0.1, test = 't.test', tip_length = 0.01, y_position = c(0.85, 0.95, 0.95, 0.85), extend_line = -0.01) + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 1.05), labels = scales::percent) +
            # scale_x_discrete(labels = c('CON', 'Bleomycin')) +
            labs(x = '', y = 'SA-β-gal positive cells (%)', title = '') + 
            theme_classic() + 
            theme(axis.text.x = element_text(size = 0, angle = 45, hjust = 1),
                axis.text.y = element_text(size = 12),
                legend.position = 'none',
                legend.title = element_text(size = 0),
                legend.text = element_text(size = 12),
                axis.title = element_text(size = 14),
                plot.title = element_text(size = 16))
#### ROS
    Intensity = c(14.23, 12.33, 10.33, 30.12, 36.67, 39.11, 9.91, 8.81, 11.88, 27.99, 36.67, 39.72)
    group = factor(c(rep('snCISD2', 3), rep('shCISD2', 3), rep('snCISD2+Hesperetin', 3), rep('shCISD2+Hesperetin', 3)), levels = c('snCISD2', 'shCISD2', 'snCISD2+Hesperetin', 'shCISD2+Hesperetin'))
    plotA57 = data.frame(Intensity, group) %>%  
        summarySE(., measurevar = "Intensity", groupvars = c("group")) %>% 
        ggplot(aes(fill = group, y = Intensity, x = group)) +
            geom_bar(stat = 'identity', position = "dodge")+
            geom_jitter(data = data.frame(Intensity, group), aes(x = group, y = Intensity), width = 0.1, stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8) +
            geom_errorbar(aes(ymin = Intensity - se, ymax = Intensity + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(data = data.frame(Intensity, group), aes(x = group, y = Intensity), 
                    map_signif_level = T, comparisons = list(c('snCISD2', 'shCISD2'), c('snCISD2', 'snCISD2+Hesperetin'), c('snCISD2+Hesperetin', 'shCISD2+Hesperetin'), c('shCISD2', 'shCISD2+Hesperetin')), vjust = 0.1, test = 't.test', tip_length = 0.01, y_position = c(41, 46, 46, 41), extend_line = -0.01) + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 50)) +
            # scale_x_discrete(labels = c('CON', 'Bleomycin')) +
            labs(x = '', y = 'Mean Fluorescence Intensity', title = '') + 
            theme_classic() + 
            theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                axis.text.y = element_text(size = 12),
                legend.position = 'none',
                legend.title = element_text(size = 0),
                legend.text = element_text(size = 12),
                axis.title = element_text(size = 14),
                plot.title = element_text(size = 16))
#### RT-qPCR
    actin = c(19.035, 18.613, 18.434, 18.954, 18.667, 18.185, 19.011, 19.117, 18.179, 18.960, 19.061, 18.157, 19.770, 18.814, 19.019, 18.691)
    IL6 = c(17.816, 16.935, 17.095, 17.855, 16.117, 16.344, 16.101, 16.017, 17.279, 17.976, 17.937, 17.709, 16.989, 16.109, 16.663, 17.017)
    IL1b= c(16.898, 16.911, 16.515, 16.816, 16.116, 15.798, 16.016, 16.044, 16.865, 17.157, 17.423, 16.951, 16.598, 16.716, 16.078, 16.816)
    group = factor(c(rep('snCISD2', 4), rep('shCISD2', 4), rep('snCISD2+Hesperetin', 4), rep('shCISD2+Hesperetin', 4)), levels = c('snCISD2', 'shCISD2', 'snCISD2+Hesperetin', 'shCISD2+Hesperetin'))
    A549 = data.frame(IL6, IL1b, actin, group) %>% 
        mutate(dct1 = IL6 - actin) %>% 
        mutate(ddct1 = dct1 - mean(dct1[which(group == 'snCISD2')])) %>%
        mutate(IL6mrna = 2^-ddct1) %>% 
        mutate(dct2 = IL1b - actin) %>% 
        mutate(ddct2 = dct2 - mean(dct2[which(group == 'snCISD2')])) %>%
        mutate(IL1mrna = 2^-ddct2)
    A549 = rbind(A549 %>% select(group, IL6mrna) %>% dplyr::rename('mrna' = 'IL6mrna'),
                 A549 %>% select(group, IL1mrna) %>% dplyr::rename('mrna' = 'IL1mrna')) %>% 
        mutate(SASP = c(rep('IL6', 16), rep('IL1b', 16)))
    A549
    t.test(mrna~group, data = A549 %>% filter(SASP == 'IL1b') %>% filter(group %in% c('snCISD2', 'shCISD2')))
    t.test(mrna~group, data = A549 %>% filter(SASP == 'IL1b') %>% filter(group %in% c('shCISD2', 'shCISD2+Hesperetin')))
    plotA54 = A549 %>%  
        summarySE(., measurevar = "mrna", groupvars = c("group", "SASP")) %>% 
        ggplot(aes(fill = group, y = mrna, x = SASP)) +
            geom_bar(stat = 'identity', position = "dodge") +
            geom_jitter(data = A549, aes(y = mrna, x = SASP), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
            geom_errorbar(aes(ymin = mrna - se, ymax = mrna + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(y_position = c(3.5), xmin = c(0.72, 0.92, 1.72, 1.92), xmax = c(0.88, 1.28, 1.88, 2.28), annotation = c("*", "NS.", "*", "NS."), tip_length = 0.01) + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 3.8)) +
            scale_x_discrete(labels = c('IL-1β', 'IL6')) +
            labs(x = '', y = 'Relative expression of SASP', title = '') + 
            theme_classic() + 
            theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    require(patchwork)
    design = 'AABBCCDDDD'
    design = 'AABBBB'
    plotA53 + plotA54 + plot_layout(ncol = 2, design = design)
    plotA52 + plotA53 + plotA57 + plotA54 + plot_layout(ncol = 4, design = design)
    ggsave('~/Studio/01IPFfibroblast/F0522.jpg', width = 12, height = 5)
    ggsave('~/Studio/01IPFfibroblast/F0523.jpg', width = 12, height = 4)
#### WB
    group = factor(rep(c('snCISD2', 'shCISD2', 'snCISD2+Hesperetin', 'shCISD2+Hesperetin'), 3), levels = c('snCISD2', 'shCISD2', 'snCISD2+Hesperetin', 'shCISD2+Hesperetin'))
    actin = c(148.857, 151.077, 149.754, 149.117, 150.667, 149.371, 149.454, 150.927, 148.160, 150.353, 149.778, 148.175)
    sma = c(73.926, 100.931, 63.184, 98.931, 67.217, 99.845, 70.027, 92.273, 80.026, 100.91, 71.841, 92.333)
    tgf = c(66.088, 87.027, 63.855, 78.027, 60.188, 82.227, 69.155, 79.170, 59.908, 94.327, 67.555, 89.022)
    p21 = c(84.078, 111.084, 70.931, 92.027, 78.055, 93.852, 67.729, 92.275, 73.679, 101.184, 87.031, 99.327)
    cdk = c(27.188, 41.025, 18.955, 37.226, 26.667, 42.239, 29.055, 38.378, 24.088, 39.274, 22.450, 39.329)   
    BECN1 = c(19.776, 50.121, 21.084, 48.231, 17.826, 44.555, 19.234, 52.070, 20.000, 47.910, 18.841, 50.333)
    p62 = c(83.212, 101.960, 78.098, 101.177, 80.146, 95.931, 87.095, 107.011, 81.590, 105.744, 72.912, 100.022)
    BCL2 = c(97.088, 77.441, 93.735, 61.927, 90.078, 70.345, 98.658, 67.941, 99.088, 72.374, 97.015, 73.475)
    cis = c(103.012, 53.096, 104.980, 54.317, 107.949, 50.238, 105.233, 47.274, 105.659, 45.104, 111.320, 52.263)
    A549 = data.frame(actin, sma/actin, tgf/actin, p21/actin, cdk/ actin, BECN1/actin, p62/actin, BCL2/actin, cis/actin,group) 
    A555 = rbind(A549 %>% select(group, sma.actin) %>% dplyr::rename('mrna' = 'sma.actin'),
                 A549 %>% select(group, tgf.actin) %>% dplyr::rename('mrna' = 'tgf.actin'),
                 A549 %>% select(group, p21.actin) %>% dplyr::rename('mrna' = 'p21.actin'),
                 A549 %>% select(group, cdk.actin) %>% dplyr::rename('mrna' = 'cdk.actin'),
                 A549 %>% select(group, BECN1.actin) %>% dplyr::rename('mrna' = 'BECN1.actin'),
                 A549 %>% select(group, p62.actin) %>% dplyr::rename('mrna' = 'p62.actin'),
                 A549 %>% select(group, BCL2.actin) %>% dplyr::rename('mrna' = 'BCL2.actin'),
                 A549 %>% select(group, cis.actin) %>% dplyr::rename('mrna' = 'cis.actin')) %>% 
            mutate(wb = c(rep('α-SMA', 12), rep('TGF-β1', 12), rep('p21', 12), rep('CDK4', 12), rep('BECN1', 12), rep('p62', 12), rep('BCL2', 12), rep('CISD2', 12)))
    A555$wb = factor(A555$wb, levels = c('α-SMA', 'TGF-β1', 'p21', 'CDK4', 'BECN1', 'p62', 'BCL2', 'CISD2'))
    t.test(mrna ~ group, data = A555 %>% filter(wb == 'p21') %>% filter(group %in% c('snCISD2', 'shCISD2')))
    t.test(mrna ~ group, data = A555 %>% filter(wb == 'p21') %>% filter(group %in% c('shCISD2', 'shCISD2+Hesperetin')))
    p1 = A555[c(1:48), ] %>%  
        summarySE(., measurevar = "mrna", groupvars = c("group", "wb")) %>% 
        ggplot(aes(fill = group, y = mrna, x = wb)) +
            geom_bar(stat = 'identity', position = "dodge") +
            geom_jitter(data = A555[c(1:48), ], aes(y = mrna, x = wb), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
            geom_errorbar(aes(ymin = mrna - se, ymax = mrna + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(y_position = c(0.9), xmin = c(0.7), xmax = c(0.85), annotation = c("***"), tip_length = 0.01) + 
            geom_signif(y_position = c(0.9), xmin = c(0.9), xmax = c(1.35), annotation = c("NS."), tip_length = 0.01) +
            geom_signif(y_position = c(0.9), xmin = c(1.7), xmax = c(1.85), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(0.9), xmin = c(1.9), xmax = c(2.35), annotation = c("NS."), tip_length = 0.01) + 
            geom_signif(y_position = c(0.9), xmin = c(2.7), xmax = c(2.85), annotation = c("*"), tip_length = 0.01) +
            geom_signif(y_position = c(0.9), xmin = c(2.9), xmax = c(3.35), annotation = c("NS."), tip_length = 0.01) +
            geom_signif(y_position = c(0.9), xmin = c(3.7), xmax = c(3.85), annotation = c("***"), tip_length = 0.01) +
            geom_signif(y_position = c(0.9), xmin = c(3.9), xmax = c(4.35), annotation = c("NS."), tip_length = 0.01) +
            theme_classic() + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
            #scale_x_discrete(labels = c('IL-1β', 'IL6')) +
            labs(x = '', y = 'Relative expression of protein', title = '')  + 
            theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 12),
              plot.title = element_text(size = 16))
    p2 = A555[c(49:96), ] %>%  
        summarySE(., measurevar = "mrna", groupvars = c("group", "wb")) %>% 
        ggplot(aes(fill = group, y = mrna, x = wb)) +
            geom_bar(stat = 'identity', position = "dodge") +
            geom_jitter(data = A555[c(49:96), ], aes(y = mrna, x = wb), stroke = 0.2, size = 2, shape = 21, color = 'black', alpha = 0.8, position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2)) +
            geom_errorbar(aes(ymin = mrna - se, ymax = mrna + se), 
                    width = 0.1, position = position_dodge(0.9), color = '#9E9E9E') + 
            geom_signif(y_position = c(0.8), xmin = c(0.7), xmax = c(0.85), annotation = c("*"), tip_length = 0.01) + 
            geom_signif(y_position = c(0.8), xmin = c(0.9), xmax = c(1.35), annotation = c("NS."), tip_length = 0.01) +
            geom_signif(y_position = c(0.8), xmin = c(1.7), xmax = c(1.85), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(0.8), xmin = c(1.9), xmax = c(2.35), annotation = c("NS."), tip_length = 0.01) + 
            geom_signif(y_position = c(0.8), xmin = c(2.7), xmax = c(2.85), annotation = c("**"), tip_length = 0.01) +
            geom_signif(y_position = c(0.8), xmin = c(2.9), xmax = c(3.35), annotation = c("NS."), tip_length = 0.01) +
            geom_signif(y_position = c(0.8), xmin = c(3.7), xmax = c(3.85), annotation = c("***"), tip_length = 0.01) +
            geom_signif(y_position = c(0.8), xmin = c(3.9), xmax = c(4.35), annotation = c("NS."), tip_length = 0.01) +
            theme_classic() + 
            scale_fill_manual(values = color) +
            scale_y_continuous(expand = c(0,0), limits = c(0, 0.9)) +
            #scale_x_discrete(labels = c('IL-1β', 'IL6')) +
            labs(x = '', y = 'Relative expression of protein', title = '')  + 
            theme(axis.text = element_text(size = 12, vjust = 0.5),
              legend.position = 'right',
              legend.title = element_text(size = 0),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 12),
              plot.title = element_text(size = 16))
    require(patchwork)
    p1/p2 + plot_layout(guides = 'collect')
    ggsave('~/Studio/01IPFfibroblast/F0323.jpg', width = 7.5, height = 5.5)
