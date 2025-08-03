## Packages for this study.
    library(tidyverse)              # Data modulation and plotting.
    library(data.table)             # Data convenience features.
    require(readxl)                 # Reading and writing from excel.
    require(Rmisc)                  # Calculation SE for ploting SE bar.
    require(ggsignif)               # Significance of ploting using ggplot2.
    require(scales)                 # Scale Functions for Visualization.
    require(ggplot2)
    require(ggpubr)
    color = c("#D4E5F4", "#FBB9BA", "#C2E3D0", "#F6E4D0", "#A1A9D0", "#DBB8B2", 
              "#E7A9C5", "#E3BBED", "#A9C287", "#B7B7EB", "#FFEBAD", '#849184', 
              "#579AC3", "#A5CCC7", "#C4A5DE", '#E64B35', "#F5EBF4", '#00A087', 
              "#D6E2E2", '#4DBBD5', "#D4E6BC", "#AA8984", "#BFB1D0", "#85BF67", 
              "#CCCC99", "#DCD7C1", "#A8CBDF", "#F0988C", "#96CCCB", "#F5DC75",
              '#3C5488', '#F39B7F', '#91D1C2', '#DC0000', '#7E6148', "#9E9E9E", 
              '#B09C85', "#EAB883", "#FFE6B4", "#A5AAF9")
    #save.image(file = '~/OneDrive/Datasets/01MetaRA1.RData')
    #load(file = '~/OneDrive/Datasets/01MetaRA1.RData')
## Figure1 GSVA analysis
#### Metabolic Pathway enrolled (Wu et al.)
    gs_th = readxl::read_excel('~/Studio/01MetaRA/MetabolicPathways.xlsx', sheet = 1)
    gs_path = split(gs_th$symbol, list(gs_th$Pathway))
#### GMT building
    name <- unique(gs_th$Pathway)
    description <- rep(NA,length(name))
    names(description) <- name
    genes <- lapply(name, function(name){
        as.vector(gs_th[gs_th$Pathway == name, "symbol"])})
    names(genes) <- name
    gmtinput <- list(name=name,description=description,genes=genes)
    get_gmt <- function(gmtinput,filename){
        output <- file(filename, open="wt")
        lapply(gmtinput[["name"]],function(name){
        outlines = paste0(c(name, gmtinput[["description"]][[name]],
                        gmtinput[["genes"]][[name]]),collapse='\t')
        writeLines(outlines, con=output)})
    close(output)}
    get_gmt(gmtinput = gmtinput,filename = "~/Studio/01MetaRA/catabolism.gmt")
#### GSE89408 Adjustment
    require(GEOquery)
    require(GSVA)
    g89408c = getGEO(filename = ("~/Studio/01MetaRA/GSE89408_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        pData(.)
    G89408 = data.table::fread('~/Studio/01MetaRA/GSE89408_GEO_count_matrix_rename.txt', header = T)[, c(1, 2:29, 68:219)] %>% 
    tibble::column_to_rownames('V1')
    G89408 = log2 (edgeR::cpm(G89408) +1)
    gg89408 = data.table::fread('~/Studio/01MetaRA/GSE89408_GEO_count_matrix_rename.txt', header = T) %>% 
    tibble::column_to_rownames('V1')
    log2 (edgeR::cpm(gg89408) +1) %>% write.table(., '~/Studio/01MetaRA/GSE89408_total.txt', sep = '\t')
#### Heatmap of GSVA GSE89408
    require(limma)
    require(GSVA)
    g89408c$'disease:ch1'[c(1:28, 67:218)]
    g = factor(c(rep('Normal', 28), rep('Rheumatoid arthritis', 152)), levels = c('Normal', 'Rheumatoid arthritis'))
    design = model.matrix( ~ 1 + g)
    gsva89408 = gsva(gsvaParam(G89408, gs_path, maxDiff = TRUE)) %>% 
        lmFit(., design) %>% 
        eBayes(.) %>% 
        topTable(., coef = 2, number = Inf) %>% 
        mutate(sign = ifelse(adj.P.Val <= 0.001, '***', ifelse(adj.P.Val <= 0.01, '**', ifelse(adj.P.Val <= 0.05, '*', 'ns.')))) %>% 
        tibble::rownames_to_column('Metabolism_pathway') 
    fit1 = gsva(gsvaParam(as.matrix(G89408), gs_path, maxDiff = TRUE)) %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('geo_accession') %>% 
        pivot_longer(!geo_accession,
               names_to = "Metabolism_pathway",
               values_to = 'value') %>% 
        inner_join(gsva89408, by = 'Metabolism_pathway') %>% 
        mutate(names = paste0(.$Metabolism_pathway, '(', .$sign, ')')) %>% 
        na.omit() 
    fit1$names = factor(fit1$names, levels = rev(unique(fit1$names)))
    plot1 = fit1 %>% 
        ggplot(aes(x = geo_accession, y = names, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_gradient(low = "#BEDADA33", high = "#E89B7E") +
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        scale_y_discrete(position = "left") +
        theme_minimal() +
        labs(x = "", y = "") +
        guides(fill = guide_legend("Value")) + 
        theme(legend.position = "right",
              legend.text = element_text(size = 14),
              panel.grid = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 14),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title = element_text(face = 'bold', size = 16))
    fit0 = gsva(gsvaParam(G89408, gs_path, maxDiff = TRUE)) %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('geo_accession') %>% 
        pivot_longer(!geo_accession,
               names_to = "Metabolism_pathway",
               values_to = 'value') %>% 
        inner_join(gsva89408, by = 'Metabolism_pathway') %>% 
        mutate(x = c(1:dim(.)[1])) %>% 
        mutate(y = c("group")) %>% 
        mutate(state = c(rep('Healthy control', 28*7), rep('Rheumatoid arthritis', 152*7))) %>% 
        na.omit()
    fit0$geo_accession = factor(fit0$geo_accession, levels = unique(fit0$geo_accession))
    plot0 = fit0 %>% 
        ggplot(aes(x = geo_accession, y = y, fill = state)) +
        geom_tile(aes(fill = state)) +
        scale_fill_manual(values = color) +
        scale_y_discrete(position = "left") +
        theme_minimal() +
        xlab(NULL) + ylab(NULL) +
        labs(fill = 'state', title = 'GSE89408') + 
        guides(fill = guide_legend("Disease")) + 
        theme(panel.border = element_blank(),
              panel.grid = element_blank(),
              legend.text = element_text(size = 14),
              axis.ticks = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 14),
              axis.text.x = element_blank(),
              plot.title = element_text(face = 'bold', size = 20),
              legend.position = "right")
    require(aplot)
    plot0 %>% insert_bottom(plot1, height = 7)
    ggsave('~/Studio/01MetaRA/F001.jpg', width = 15, height = 3.6) 
#### Woetzel study extracted
    require(GEOquery)
    if (length(getGEO(filename = ("~/Studio/01MetaRA/Woetzel/GSE55584_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F)) > 1) {
        idx = grep("GPL96", attr(getGEO(filename = ("~/Studio/01MetaRA/Woetzel/GSE55584_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F), "names"))} else {
        idx = 1}
    x1 = getGEO(filename = ("~/Studio/01MetaRA/Woetzel/GSE55584_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        exprs(.)
    qx = as.numeric(quantile(x1, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    LogC = (qx[5] > 100) ||
        (qx[6] - qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) {
        x1[which(x1 <= 0)] = NaN
        x1 = log2(x1 + 1)
        print("log2 transform finished")} else {
        print("log2 transform not needed")}
    if (length(getGEO(filename = ("~/Studio/01MetaRA/Woetzel/GSE55457_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F)) > 1) {
        idx = grep("GPL96", attr(getGEO(filename = ("~/Studio/01MetaRA/Woetzel/GSE55457_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F), "names"))} else {
        idx = 1}
    x2 = getGEO(filename = ("~/Studio/01MetaRA/Woetzel/GSE55457_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        exprs(.)
    qx = as.numeric(quantile(x2, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    LogC = (qx[5] > 100) ||
        (qx[6] - qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) {
        x2[which(x2 <= 0)] = NaN
        x2 = log2(x2 + 1)
        print("log2 transform finished")} else {
        print("log2 transform not needed")}
    x3 = getGEO(filename = ("~/Studio/01MetaRA/Woetzel/GSE55235_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        exprs(.)
    qx = as.numeric(quantile(x3, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    LogC = (qx[5] > 100) ||
        (qx[6] - qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) {
        x3[which(x3 <= 0)] = NaN
        x3 = log2(x3 + 1)
        print("log2 transform finished")} else {
        print("log2 transform not needed")}
    x0 = merge((as.data.frame(x1) %>% tibble::rownames_to_column('probe_id')), (as.data.frame(x2) %>% 
        tibble::rownames_to_column('probe_id')), by = 'probe_id') %>% 
        merge(., (as.data.frame(x3) %>% tibble::rownames_to_column('probe_id')), by = 'probe_id') %>% 
        inner_join(AnnoProbe::idmap('GPL96'), by = 'probe_id') %>% 
        select('symbol', everything()) %>% 
        mutate(rowMean = rowMeans(.[grep('GSM', names(.))])) %>%
        dplyr::arrange(desc(rowMean)) %>%
        dplyr::rename('symbol' = 'symbol') %>% 
        dplyr::distinct(symbol, .keep_all = T) %>%
        dplyr::select(-c('rowMean', 'probe_id'))
    Woetzelc = rbind(read.csv('~/Studio/01MetaRA/Woetzel/GSE55457RABC11_PhenotypeFile.csv'), read.csv('~/Studio/01MetaRA/Woetzel/GSE55235RABC10_PhenotypeFile.csv'), read.csv('~/Studio/01MetaRA/Woetzel/GSE55584RABC51_PhenotypeFile.csv') %>% filter(status == '1')) %>% mutate(group = c(rep('GSE55457', 23), rep('GSE55235', 20), rep('GSE55584', 10)))
    Woetzel = x0 %>% 
        tibble::column_to_rownames('symbol') %>% 
        select(Woetzelc$sample)
    require(limma)
    batch = data.frame(sample = names(Woetzel)) %>% 
        inner_join(Woetzelc, by = 'sample')
    eWoetzel = removeBatchEffect(Woetzel, batch$group) %>% 
        na.omit()
    plot2 = eWoetzel %>% 
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('gsm') %>% 
        pivot_longer(!gsm, names_to = 'symbol', values_to = "value") %>% 
        ggplot(aes(x = gsm, y = value)) + 
            geom_boxplot(width = 0.5, position = position_dodge(0.75), fill = '#FBB9BA') +
            theme_classic2() +
            labs(title = "Woetzel\'s study", y = "Gene expression", x = "") +
            scale_fill_manual(values = color) + 
            theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 8),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 18, face = 'bold'))
    ggsave('~/Studio/01MetaRA/S002.jpg', width = 11.5, height = 4)
    boxplot(eWoetzel)
    head(eWoetzel)
    write.table(as.data.frame(eWoetzel), '~/Studio/01MetaRA/Woetzel/Totalexpression.txt', row.names = TRUE, sep = '\t')
    write.table(Woetzelc, '~/Studio/01MetaRA/Woetzel/Totalclinic.txt', row.names = TRUE)
#### GSVA analysis for Woetzel study
    require(GSVA)
    require(aplot)
    eWoetzel = read.table('~/Studio/01MetaRA/Woetzel/Totalexpression.txt', header =  T)
    cWoetzel = read.table('~/Studio/01MetaRA/Woetzel/Totalclinic.txt', header =  T)
    g = factor(cWoetzel$status, levels = c('0', '1'))
    design = model.matrix( ~ 1 + g)
    gsvaWoetzel = gsva(gsvaParam(as.matrix(eWoetzel), gs_path, maxDiff = TRUE)) %>% 
        lmFit(., design) %>% 
        eBayes(.) %>% 
        topTable(., coef = 2, number = Inf) %>% 
        mutate(sign = ifelse(adj.P.Val <= 0.001, '***', ifelse(adj.P.Val <= 0.01, '**', ifelse(adj.P.Val <= 0.05, '*', 'ns.')))) %>% 
        tibble::rownames_to_column('Metabolism_pathway') 
#### Heatmap of GSVA Woetzelff
    fit1 = gsva(gsvaParam(as.matrix(eWoetzel), gs_path, maxDiff = TRUE)) %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('geo_accession') %>% 
        pivot_longer(!geo_accession,
               names_to = "Metabolism_pathway",
               values_to = 'value') %>% 
        inner_join(gsvaWoetzel, by = 'Metabolism_pathway') %>% 
        mutate(names = paste0(.$Metabolism_pathway, '(', .$sign, ')')) %>% 
        dplyr::rename('sample' = 'geo_accession') %>% 
        inner_join(cWoetzel, by = 'sample') %>% 
        arrange(status) %>% 
        na.omit() 
    fit1$sample = factor(fit1$sample, levels = unique(fit1$sample))
    fit1$names = factor(fit1$names, levels = rev(unique(fit1$names)))
    plot1 = fit1 %>% 
        ggplot(aes(x = sample, y = names, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_gradient(low = "#FBF8B422", high = "#FBB463") +
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        scale_y_discrete(position = "left") +
        theme_minimal() +
        labs(x = "", y = "") +
        guides(fill = guide_legend("Value")) + 
        theme(legend.position = "right",
              legend.text = element_text(size = 14),
              panel.grid = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 14),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title = element_text(face = 'bold', size = 16))
    plot0 = fit1 %>% 
        mutate(y = c("group")) %>% 
        mutate(state = ifelse(status == '0', 'Healthy control', 'Rheumatoid arthritis')) %>% 
        ggplot(aes(x = sample, y = y, fill = state)) +
        geom_tile(aes(fill = state)) +
        scale_fill_manual(values = color) +
        scale_y_discrete(position = "left") +
        theme_minimal() +
        xlab(NULL) + ylab(NULL) +
        labs(fill = 'state', title = 'Woetzel\'s cohort') + 
        guides(fill = guide_legend("Disease")) + 
        theme(panel.border = element_blank(),
              panel.grid = element_blank(),
              legend.text = element_text(size = 14),
              axis.ticks = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 14),
              axis.text.x = element_blank(),
              plot.title = element_text(face = 'bold', size = 20),
              legend.position = "right")
    plot0 %>% insert_bottom(plot1, height = 7)
    ggsave('~/Studio/01MetaRA/F002.jpg', width = 15, height = 3.6) 
#### GSE93272 selected 
    require(GEOquery)
    require(limma)
    require(GSVA)
    g93272c = getGEO(filename = ("~/Studio/01MetaRA/GSE93272_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        pData(.) %>% 
        mutate(patients = substr(.$characteristics_ch1, start = 16, nchar(.$characteristics_ch1))) 
    g93272c$title = gsub('Whole blood from rheumatoid arthritis', '', g93272c$title)
    g93272c$title = gsub('Whole blood from healthy control', '', g93272c$title)
    drugnaive = ((read.table("~/Studio/01MetaRA/drug.naive_metadata.txt", header = T) %>% dplyr::rename('patients' = 'INDIVIDUAL.ID'))) %>% mutate(title = paste0('(', SAMPLE.ID, ')')) %>% 
        inner_join(g93272c[, c('title', 'geo_accession')], by = 'title')
    e93272 = getGEO(filename = ("~/Studio/01MetaRA/GSE93272_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        exprs(.) %>% 
        as.data.frame() %>% 
        select(drugnaive$geo_accession) %>% 
        tibble::rownames_to_column('probe_id') %>% 
        inner_join(AnnoProbe::idmap('GPL570'), by = 'probe_id') %>% 
        select('symbol', everything()) %>% 
        mutate(rowMean = rowMeans(.[grep('GSM', names(.))])) %>%
        dplyr::arrange(desc(rowMean)) %>%
        dplyr::rename('symbol' = 'symbol') %>% 
        dplyr::distinct(symbol, .keep_all = T) %>%
        dplyr::select(-c('rowMean', 'probe_id')) %>% 
        tibble::column_to_rownames('symbol')
    e93272 %>% tibble::rownames_to_column('symbol') %>% write.table(., '~/Studio/01MetaRA/GSE93272_series.txt', sep = '\t', row.names = F)
#### limma GSE93272
    status = factor(drugnaive$DISEASE, levels = c('HC', 'RA'))
    design = model.matrix(~ 0 + status)
    colnames(design) = c("HC", "RA")
    contrast.matrix = makeContrasts(RA - HC, levels = colnames(design))
    fit = e93272 %>% 
        lmFit(., design) %>% 
        contrasts.fit(., contrast.matrix) %>% 
        eBayes(.)
    sig = topTable(fit, n = Inf, adjust = "fdr") %>% 
        na.omit() %>% 
        dplyr::rename('log2FC'="logFC") %>% 
        mutate(group = ifelse(adj.P.Val <= 0.05 & log2FC > 0.5, 'Up', ifelse(adj.P.Val <= 0.05 & log2FC < - 0.5, 'Down', 'Not'))) %>% 
        tibble::rownames_to_column('symbol') %>% 
        mutate(log10P = -log10(adj.P.Val)) %>% 
        arrange(group)
    x1 = sig %>% filter(group == 'Up' | group == 'Down')
    intersect(x1$symbol, key)
    write.csv(x1, '~/Studio/01MetaRA/x1.csv')
    require(ggpubr)
    require(ggsci)
    plot1 = ggscatter(sig, "log2FC", "log10P",
        combine = F, merge = T, color = "group", shape = 20, size = 1,
        point = TRUE, font.label = 20, palette = c("#D4E5F4", '#9E9E9E', "#FBB9BA")) +
        geom_vline(xintercept = c(-1, 1), colour = "black", linetype = "dashed") +
        geom_hline(yintercept = -log10(0.05), colour = "black", linetype = "dashed") +
        labs(title = "GSE89408", x = 'log2(FoldChange)', y = 'log10(adj.p.value)') + 
        theme_classic2() +
        guides(color = guide_legend(title = "")) +
        theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 1),
                  axis.text.y = element_text(size = 12),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 10),
                  legend.position = 'right',
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 18, face = 'bold'))
    ggsave('~/Studio/01MetaRA/F004a.jpg', width = 3.6, height = 3.2)
#### Heatmap of GSVA GSE93272
    g = factor(drugnaive$DISEASE, levels = c('HC', 'RA'))
    design = model.matrix( ~ 1 + g)
    gsva93272 = gsva(gsvaParam(as.matrix(e93272), gs_path, maxDiff = TRUE)) %>% 
        lmFit(., design) %>% 
        eBayes(.) %>% 
        topTable(., coef = 2, number = Inf) %>% 
        mutate(sign = ifelse(adj.P.Val <= 0.001, '***', ifelse(adj.P.Val <= 0.01, '**', ifelse(adj.P.Val <= 0.05, '*', 'ns.')))) %>% 
        tibble::rownames_to_column('Metabolism_pathway') 
    fit1 = gsva(gsvaParam(as.matrix(e93272), gs_path, maxDiff = TRUE)) %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('geo_accession') %>% 
        pivot_longer(!geo_accession,
               names_to = "Metabolism_pathway",
               values_to = 'value') %>% 
        inner_join(gsva93272, by = 'Metabolism_pathway') %>% 
        mutate(names = paste0(.$Metabolism_pathway, '(', .$sign, ')')) %>% 
        na.omit() 
    fit1$names = factor(fit1$names, levels = rev(unique(fit1$names)))
    plot1 = fit1 %>% 
        ggplot(aes(x = geo_accession, y = names, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_gradient(low = "#FBF8B422", high = "#FBB463") +
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        scale_y_discrete(position = "left") +
        theme_minimal() +
        labs(x = "", y = "") +
        guides(fill = guide_legend("Value")) + 
        theme(legend.position = "right",
              legend.text = element_text(size = 14),
              panel.grid = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 14),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title = element_text(face = 'bold', size = 16))
    fit0 = gsva(gsvaParam(as.matrix(e93272), gs_path, maxDiff = TRUE)) %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('geo_accession') %>% 
        pivot_longer(!geo_accession,
               names_to = "Metabolism_pathway",
               values_to = 'value') %>% 
        inner_join(gsva93272, by = 'Metabolism_pathway') %>% 
        mutate(x = c(1:dim(.)[1])) %>% 
        mutate(y = c("group")) %>% 
        inner_join(drugnaive[, c('geo_accession', 'DISEASE')], by = 'geo_accession') %>% 
        na.omit() %>% 
        arrange(DISEASE)
    fit0$geo_accession = factor(fit0$geo_accession, levels = unique(fit0$geo_accession))
    fit0$dis = ifelse(fit0$DISEASE == 'HC', 'Healthy control', 'Rheumatoid arthritis')
    plot0 = fit0 %>% 
        ggplot(aes(x = geo_accession, y = y, fill = dis)) +
        geom_tile(aes(fill = dis)) +
        scale_fill_manual(values = color) +
        scale_y_discrete(position = "left") +
        theme_minimal() +
        xlab(NULL) + ylab(NULL) +
        labs(fill = 'dis', title = 'GSE93272') + 
        guides(fill = guide_legend("Disease")) + 
        theme(panel.border = element_blank(),
              panel.grid = element_blank(),
              legend.text = element_text(size = 14),
              axis.ticks = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 14),
              axis.text.x = element_blank(),
              plot.title = element_text(face = 'bold', size = 20),
              legend.position = "right")
    require(aplot)
    plot0 %>% insert_bottom(plot1, height = 7)
    ggsave('~/Studio/01MetaRA/F002.jpg', width = 15, height = 3.6) 
#### ROC plot
    require(plotROC) 
    require(pROC)
    x89408 = gsva(gsvaParam(G89408, gs_path, maxDiff = TRUE)) %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('geo_accession') %>% 
        mutate(status = c(rep('0', 28), rep('1', 152)))
    x0 = data.frame()
    for(i in 2:8){
        x89408$status = factor(x89408$status, levels = c('1', '0'))
        x1 = roc(response = x89408$status, predictor = x89408[, i])
        x2 = round(auc(response = x89408$status, predictor = x89408[, i]), 4)
        x3 = data.frame(symbol = names(x89408)[i], value = x2)
        x0 = rbind(x3, x0)}
    auc89408 = x0
    y89408 = x89408[, - 9] %>% 
        pivot_longer(!geo_accession,
            names_to = "Meta",
            values_to = 'value')  %>% 
        inner_join(x89408 %>% select(c('geo_accession', 'status')), by = 'geo_accession') %>% 
        group_split(Meta)
    str(y89408)
    map_function <- function(y89408){
        pROC::roc(
        data = y89408,
        response = status,
        predictor = value)}
    roc.list = map(y89408, map_function)
    names(roc.list) = colnames(x89408)[2:8]
    plot1 = pROC::ggroc(roc.list, alpha = 1, size = 0.75, legacy = T) +
        geom_abline(linetype = 3, intercept = 0) +
        labs(title = "GSE89408", color = 'Metabolism') +
        theme_classic() +
        scale_color_manual(values = color) +
        theme(axis.text = element_text(size = 12, vjust = 0.5),
              # legend.position = 'left',
              legend.title = element_text(size = 16, face = 'bold'),
              legend.text = element_text(size = 12),
              axis.title = element_text(size = 14, face = 'bold'),
              plot.title = element_text(size = 16, face = 'bold'))    
    x93272 = gsva(gsvaParam(as.matrix(e93272), gs_path, maxDiff = TRUE)) %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('geo_accession') %>% 
        inner_join(drugnaive[, c('geo_accession', 'DISEASE')], by = 'geo_accession')
    x0 = data.frame()
    for(i in 2:8){
        x93272$DISEASE = factor(x93272$DISEASE, levels = c('HC', 'RA'))
        x1 = roc(response = x93272$DISEASE, predictor = x93272[, i])
        x2 = round(auc(response = x93272$DISEASE, predictor = x93272[, i]), 4)
        x3 = data.frame(symbol = names(x93272)[i], value = x2)
        x0 = rbind(x3, x0)}
    auc93272 = x0
    y93272 = x93272[, - 9] %>% 
        pivot_longer(!geo_accession,
            names_to = "Meta",
            values_to = 'value')  %>% 
        inner_join(drugnaive[, c('geo_accession', 'DISEASE')], by = 'geo_accession') %>% 
        group_split(Meta)
    str(y93272)
    map_function <- function(y93272){
        pROC::roc(
        data = y93272,
        response = DISEASE,
        predictor = value)}
    roc.list = map(y93272, map_function)
    names(roc.list) = colnames(x93272)[2:8]
    plot2 = pROC::ggroc(roc.list, alpha = 1, size = 0.75, legacy = T) +
        geom_abline(linetype = 3, intercept = 0) +
        labs(title = "GSE93272", color = 'Metabolism') +
        theme_classic() +
        scale_color_manual(values = color) +
        theme(axis.text = element_text(size = 12, vjust = 0.5),
              # legend.position = 'left',
              legend.title = element_text(size = 16, face = 'bold'),
              legend.text = element_text(size = 12),
              axis.title = element_text(size = 14, face = 'bold'),
              plot.title = element_text(size = 16, face = 'bold'))    
    require(patchwork)
    plot1 + plot2 + plot_layout(guides = 'collect')
    ggsave('~/Studio/01MetaRA/F003.jpg', width = 8, height = 3.6) 
    merge(auc89408 %>% dplyr::rename('GSE89408' = 'value'), auc93272 %>% dplyr::rename('GSE93272' = 'value'), by = 'symbol') %>% write.table(., '~/Studio/01MetaRA/newAUC.txt')
#### Venn plot
    require(ggvenn)
    a = list(GSE89048 = (auc89408 %>% filter(value > 0.7))$symbol, GSE93272 = (auc93272 %>% filter(value > 0.7))$symbol)
    ggvenn(a, stroke_size = 0, set_name_color = "black", set_name_size = 6, show_outside = c("none"),
          fill_color = color, text_size = 6, auto_scale = FALSE)
    ggsave('~/Studio/01MetaRA/F003.jpg', width = 6, height = 4) 
    intersect((auc89408 %>% filter(value > 0.7))$symbol, (auc93272 %>% filter(value > 0.7))$symbol)
#### ROC Woetzel
    xWoetzel = gsva(gsvaParam(as.matrix(eWoetzel), gs_path, maxDiff = TRUE)) %>% 
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('sample') %>% 
        inner_join(cWoetzel, by = 'sample')
    x0 = data.frame()
    for(i in 2:8){
        xWoetzel$status = factor(xWoetzel$status, levels = c('1', '0'))
        x1 = roc(response = xWoetzel$status, predictor = xWoetzel[, i])
        x2 = round(auc(response = xWoetzel$status, predictor = xWoetzel[, i]), 4)
        x3 = data.frame(symbol = names(xWoetzel)[i], value = x2)
        x0 = rbind(x3, x0)}
    aucWoetzel = x0
    yWoetzel = xWoetzel[, - c(9, 10)] %>% 
        pivot_longer(!sample,
            names_to = "Meta",
            values_to = 'value')  %>% 
        inner_join(cWoetzel, by = 'sample') %>% 
        group_split(Meta)
    map_function <- function(yWoetzel){
        pROC::roc(
        data = yWoetzel,
        response = status,
        predictor = value)}
    roc.list = map(yWoetzel, map_function)
    names(roc.list) = colnames(xWoetzel)[2:8]
    plot2 = pROC::ggroc(roc.list, alpha = 1, size = 0.75, legacy = T) +
        geom_abline(linetype = 3, intercept = 0) +
        labs(title = "Woetzel\'s cohort", color = 'Metabolism') +
        theme_classic() +
        scale_color_manual(values = color) +
        theme(axis.text = element_text(size = 12, vjust = 0.5),
              # legend.position = 'left',
              legend.title = element_text(size = 16, face = 'bold'),
              legend.text = element_text(size = 12),
              axis.title = element_text(size = 14, face = 'bold'),
              plot.title = element_text(size = 16, face = 'bold'))    
    merge(auc89408 %>% dplyr::rename('GSE89408' = 'value'), aucWoetzel %>% dplyr::rename('Woetzel' = 'value'), by = 'symbol') %>% dplyr::rename('Metabolism' = 'symbol') %>% write.csv(., '~/Studio/01MetaRA/ROCMetabolism7.csv')
    require(patchwork)
    plot1 + plot2 + plot_layout(guide = 'collect')    
    ggsave('~/Studio/01MetaRA/F003.jpg', width = 8, height = 4) 
#### Venn plot
    require(ggvenn)
    a = list('GSE89408' = (auc89408 %>% filter(value > 0.7))$symbol, 
             'Woetzel\'s cohort' = (aucWoetzel %>% filter(value > 0.7))$symbol)
    ggvenn(a, stroke_size = 1, set_name_color = "black", set_name_size = 10, show_outside = c("none"),
          fill_color = color, text_size = 10, auto_scale = FALSE)     
    ggsave('~/Studio/01MetaRA/F003a.jpg', width = 10, height = 12) 
    intersect((auc89408 %>% filter(value > 0.7))$symbol, (aucWoetzel %>% filter(value > 0.7))$symbol)
#### GSE12021 combined
    require(GEOquery)
    require(AnnoProbe)
    if (length(getGEO(filename = ("~/Studio/01MetaRA/GSE12021-GPL96_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F)) > 1) {
        idx = grep("GPL96", attr(getGEO(filename = ("~/Studio/01MetaRA/GSE12021-GPL96_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F), "names"))} else {
        idx = 1}
    y1 = getGEO(filename = ("~/Studio/01MetaRA/GSE12021-GPL96_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        exprs(.)
    qx = as.numeric(quantile(y1, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    LogC = (qx[5] > 100) ||
        (qx[6] - qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) {
        y1[which(y1 <= 0)] = NaN
        y1 = log2(y1 + 1)
        print("log2 transform finished")} else {
        print("log2 transform not needed")}
    if (length(getGEO(filename = ("~/Studio/01MetaRA/GSE12021-GPL97_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F)) > 1) {
        idx = grep("GPL97", attr(getGEO(filename = ("~/Studio/01MetaRA/GSE12021-GPL97_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F), "names"))} else {
        idx = 1}
    y2 = getGEO(filename = ("~/Studio/01MetaRA/GSE12021-GPL97_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        exprs(.)
    qx = as.numeric(quantile(y2, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    LogC = (qx[5] > 100) ||
        (qx[6] - qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) {
        y2[which(y2 <= 0)] = NaN
        y2 = log2(y2 + 1)
        print("log2 transform finished")} else {
        print("log2 transform not needed")}
    y0 = merge((as.data.frame(y1) %>% tibble::rownames_to_column('probe_id')), (as.data.frame(y2) %>%   
        tibble::rownames_to_column('probe_id')), by = 'probe_id') %>% 
        inner_join(AnnoProbe::idmap('GPL96'), by = 'probe_id') %>% 
        select('symbol', everything()) %>% 
        mutate(rowMean = rowMeans(.[grep('GSM', names(.))])) %>%
        dplyr::arrange(desc(rowMean)) %>%
        dplyr::rename('symbol' = 'symbol') %>% 
        dplyr::distinct(symbol, .keep_all = T) %>%
        dplyr::select(-c('rowMean', 'probe_id'))
    boxplot(y0[, -1])
    G12021c1 = getGEO(filename = ("~/Studio/01MetaRA/GSE12021-GPL96_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        pData(.) %>% 
        select(c('characteristics_ch1.2', 'geo_accession')) %>% 
        dplyr::rename('disease' = 'characteristics_ch1.2')  %>% 
        mutate(gse = 'GPL96')
    G12021c2 = getGEO(filename = ("~/Studio/01MetaRA/GSE12021-GPL97_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        pData(.) %>% 
        select(c('disease:ch1', 'geo_accession')) %>% 
        dplyr::rename('disease' = 'disease:ch1')  %>% 
        mutate(gse = 'GPL97')
    G12021c = rbind(G12021c1, G12021c2)
    G12021c$disease = gsub('disease: ', '', G12021c$disease)
    G12021c$disease = gsub('diesease: ', '', G12021c$disease)
    G12021c$disease = gsub('normal control', 'normal control', G12021c$disease)
    G12021c$group = c(G12021c$disease[1:53], rep('normal control', 4))
    G12021cselected = G12021c %>% 
        filter(group != 'osteoarthritis') %>% 
        mutate(dis = ifelse(group == 'normal control', 'Normal control', 'Rheumatoid arthritis'))
    head(G12021cselected)
    y12021 = y0 %>% 
        tibble::column_to_rownames('symbol') %>% 
        select(G12021cselected$geo_accession)
    boxplot(y12021)
    require(limma)
    batch = G12021cselected$gse
    e12021 = removeBatchEffect(y12021, batch)
    write.csv(e12021, '~/Studio/01MetaRA/GSE12021expression.csv', row.names = TRUE)
    write.csv(G12021cselected, '~/Studio/01MetaRA/GSE12021clinic.csv', row.names = TRUE)
    plot1 = e12021 %>% 
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('geo') %>% 
        pivot_longer(!geo, names_to = 'symbol', values_to = 'value') %>% 
        ggplot(aes(x = geo, y = value)) + 
            geom_boxplot(width = 0.5, position = position_dodge(0.75), fill = '#FBB9BA') +
            theme_classic2() +
            # stat_compare_means(label = "p.signif", method = "wilcox.test") +
            # facet_wrap(~ symbol, nrow = 1) +
            # ggsignif::geom_signif(aes(x = symbol, y = value, group = group), map_signif_level = F, parse = T, y_position = 9.5, comparisons = list(c('Rheumatoid arthritis', 'Healthy control')), size = 0.5, textsize = 2.5, vjust = 0.1, tip_length = 0.01) +
            labs(title = "GSE12021", y = "Gene expression", x = "") +
            scale_fill_manual(values = color) + 
            theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 8),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 18, face = 'bold'))
    ggsave('~/Studio/01MetaRA/S001.jpg', width = 11.5, height = 4)
#### GSVA for GSE12021
    require(limma)
    require(GSVA)
    e12021 = read.csv('~/Studio/01MetaRA/GSE12021expression.csv', header = T) %>% 
        tibble::column_to_rownames('X')
    c12021 = read.csv('~/Studio/01MetaRA/GSE12021clinic.csv', header = T)
    g = factor(c12021$dis, levels = c('Rheumatoid arthritis', 'Normal control'))
    design = model.matrix( ~ 1 + g)
    gsva12021 = gsva(gsvaParam(as.matrix(e12021), gs_path, maxDiff = TRUE, minSize = 2)) %>% 
        lmFit(., design) %>% 
        eBayes(.) %>% 
        topTable(., coef = 2, number = Inf) %>% 
        mutate(sign = ifelse(adj.P.Val <= 0.001, '***', ifelse(adj.P.Val <= 0.01, '**', ifelse(adj.P.Val <= 0.05, '*', 'ns.')))) %>% 
        tibble::rownames_to_column('Metabolism_pathway') 
    fit1 = gsva(gsvaParam(as.matrix(e12021), gs_path, maxDiff = TRUE)) %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('geo_accession') %>% 
        pivot_longer(!geo_accession,
               names_to = "Metabolism_pathway",
               values_to = 'value') %>% 
        inner_join(gsvaWoetzel, by = 'Metabolism_pathway') %>% 
        mutate(names = paste0(.$Metabolism_pathway, '(', .$sign, ')')) %>% 
        dplyr::rename('sample' = 'geo_accession') %>% 
        inner_join(Woetzelc, by = 'sample') %>% 
        arrange(status) %>% 
        na.omit() 
    fit1$sample = factor(fit1$sample, levels = unique(fit1$sample))
    fit1$Metabolism_pathway = factor(fit1$Metabolism_pathway, levels = c("Aerobic respiration and respiratory electron transport", "Biological oxidations", "Integration of energy metabolism", "Metabolism of amino acids and derivatives", "Metabolism of carbohydrates and carbohydrate derivatives", "Metabolism of lipids", "Metabolism of nucleotides", "Metabolism of vitamins and cofactors"))
    plot1 = fit1 %>% 
        ggplot(aes(x = sample, y = rev(names), fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_gradient(low = "#FBF8B422", high = "#FBB463") +
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        scale_y_discrete(position = "left") +
        theme_minimal() +
        labs(x = "", y = "") +
        guides(fill = guide_legend("Value")) + 
        theme(legend.position = "right",
              legend.text = element_text(size = 14),
              panel.grid = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 14),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title = element_text(face = 'bold', size = 16))
    plot0 = fit1 %>% 
        mutate(y = c("group")) %>% 
        mutate(state = c(rep('Healthy control', 160), rep('Rheumatoid arthritis', 264))) %>% 
        ggplot(aes(x = sample, y = y, fill = state)) +
        geom_tile(aes(fill = state)) +
        scale_fill_manual(values = color) +
        scale_y_discrete(position = "left") +
        theme_minimal() +
        xlab(NULL) + ylab(NULL) +
        labs(fill = 'state', title = 'Woetzel\'s study') + 
        guides(fill = guide_legend("Disease")) + 
        theme(panel.border = element_blank(),
              panel.grid = element_blank(),
              legend.text = element_text(size = 14),
              axis.ticks = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 14),
              axis.text.x = element_blank(),
              plot.title = element_text(face = 'bold', size = 20),
              legend.position = "right")
    plot0 %>% insert_bottom(plot1, height = 8)
#### Cibersortx analysis
    c89408 = data.table::fread("~/Studio/01MetaRA/CIBERSORTx_GSE89048_Results.txt", header = T)[c(1:28, 67:218), ] %>% 
        select(c(1:23)) %>% 
        dplyr::rename('sample' = 'Mixture') %>% 
        pivot_longer(!sample, names_to = 'celltype', values_to = 'value') %>% 
        mutate(group = c(rep('Healthy control', 22*28), rep('Rheumatoid arthritis', 22*152))) %>% 
        ggplot(aes(x = celltype, y = value, fill = group)) +
            geom_boxplot(width = 0.5, position = position_dodge(0.75)) +
            theme_classic2() +
            stat_compare_means(label = "p.signif", method = "wilcox.test") +
            # facet_wrap(~ symbol, nrow = 1) +
            # ggsignif::geom_signif(aes(x = symbol, y = value, group = group), map_signif_level = F, parse = T, y_position = 9.5, comparisons = list(c('Rheumatoid arthritis', 'Healthy control')), size = 0.5, textsize = 2.5, vjust = 0.1, tip_length = 0.01) +
            labs(title = "GSE89408", y = "Estimated proportion", x = "") +
            scale_fill_manual(values = color) + 
            theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 8),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 18, face = 'bold'))
    c93272 = data.table::fread("~/Studio/01MetaRA/CIBERSORTx_GSE93272_Results.txt", header = T) %>% 
        select(c(1:23)) %>% 
        tibble::column_to_rownames('Mixture') %>% 
        t() %>% 
        as.data.frame() %>% 
        select(drugnaive$geo_accession) %>% 
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('geo_accession') %>% 
        pivot_longer(!geo_accession, names_to = 'celltype', values_to = 'value') %>% 
        inner_join(drugnaive[, c('geo_accession', 'DISEASE')], by = 'geo_accession') %>% 
        mutate(group = ifelse(DISEASE == 'RA', 'Rheumatoid arthritis', 'Healthy control')) %>% 
        ggplot(aes(x = celltype, y = value, fill = group)) +
            geom_boxplot(width = 0.5, position = position_dodge(0.75)) +
            theme_classic2() +
            stat_compare_means(label = "p.signif", method = "wilcox.test") +
            # facet_wrap(~ symbol, nrow = 1) +
            # ggsignif::geom_signif(aes(x = symbol, y = value, group = group), map_signif_level = F, parse = T, y_position = 9.5, comparisons = list(c('Rheumatoid arthritis', 'Healthy control')), size = 0.5, textsize = 2.5, vjust = 0.1, tip_length = 0.01) +
            labs(title = "GSE93272", y = "Estimated proportion", x = "") +
            scale_fill_manual(values = color) + 
            theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 8),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 18, face = 'bold'))
    c89408 / c93272 + plot_layout(guides = 'collect')
    ggsave('~/Studio/01MetaRA/F001a.jpg', width = 11.5, height = 10)
    cWoetzel = data.table::fread("~/Studio/01MetaRA/CIBERSORTx_Woetzel_Results.txt", header = T) %>% 
        select(c(1:23)) %>% 
        dplyr::rename('sample' = 'Mixture') %>% 
        pivot_longer(!sample, names_to = 'celltype', values_to = 'value') %>% 
        inner_join(read.csv('~/Studio/01MetaRA/Woetzel/Totalclinic.csv') %>% select(c('sample', 'status')), by = 'sample') %>% 
        mutate(group = ifelse(status == '1', 'Rheumatoid arthritis', 'Normal control')) %>% 
        ggplot(aes(x = celltype, y = value, fill = group)) +
            geom_boxplot(width = 0.5, position = position_dodge(0.75)) +
            theme_classic2() +
            stat_compare_means(label = "p.signif", method = "wilcox.test") +
            # facet_wrap(~ symbol, nrow = 1) +
            # ggsignif::geom_signif(aes(x = symbol, y = value, group = group), map_signif_level = F, parse = T, y_position = 9.5, comparisons = list(c('Rheumatoid arthritis', 'Healthy control')), size = 0.5, textsize = 2.5, vjust = 0.1, tip_length = 0.01) +
            labs(title = "Woetzel\'s study", y = "Estimated proportion", x = "") +
            scale_fill_manual(values = color) + 
            theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 8),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 18, face = 'bold'))
    require(patchwork)
    df = data.table::fread("~/Studio/01MetaRA/CIBERSORTx_GSE89048_Results.txt", header = T)[c(1:28, 67:218), ] %>% 
        select(c(1:23)) %>% 
        dplyr::rename('sample' = 'Mixture') %>% 
        pivot_longer(!sample, names_to = 'celltype', values_to = 'value') %>% 
        mutate(group = c(rep('Healthy control', 22*28), rep('Rheumatoid arthritis', 22*152)))
    x0 = data.frame()
    for(i in unique(df$celltype)){
        x1 = df %>% filter(celltype == i)
        x2 = wilcox.test(value ~ group, data = x1)
        x3 = data.frame(type = i, pvalue = x2$'p.value', statistic = x2$statistic)
        x0 = rbind(x3, x0)
    }
    df89408 <- x0 %>% 
        mutate(group = 'GSE89408')
    df = data.table::fread("~/Studio/01MetaRA/CIBERSORTx_GSE93272_Results.txt", header = T) %>% 
        select(c(1:23)) %>% 
        dplyr::rename('geo_accession' = 'Mixture') %>% 
        pivot_longer(!geo_accession, names_to = 'celltype', values_to = 'value') %>% 
        inner_join(drugnaive[, c('geo_accession', 'DISEASE')], by = 'geo_accession') %>% 
        mutate(group = ifelse(DISEASE == 'RA', 'Rheumatoid arthritis', 'Healthy control'))
    x0 = data.frame()
    for(i in unique(df$celltype)){
        x1 = df %>% filter(celltype == i)
        x2 = wilcox.test(value ~ group, data = x1)
        x3 = data.frame(type = i, pvalue = x2$'p.value', statistic = x2$statistic)
        x0 = rbind(x3, x0)
    }
    df93272 <- x0 %>% 
        mutate(group = 'GSE93272')
    require(ggvenn)
    a = list('GSE89408' = (df89408 %>% filter(pvalue <= 0.05))$type, 
             'GSE93272' = (df93272 %>% filter(pvalue <= 0.05))$type)
    ggvenn(a, stroke_size = 0, set_name_color = "black", set_name_size = 11, show_outside = c("none"),
          fill_color = color, text_size = 10, auto_scale = FALSE)            # draw three-set venn
    ggsave('~/Studio/01MetaRA/F001b.jpg', width = 10, height = 12)
    celltypeselected = (intersect((df89408 %>% filter(pvalue <= 0.05))$type, (df93272 %>% filter(pvalue <= 0.05))$type))[1:5]
#### corrplot 
    require(GSVA)
    g89408c = getGEO(filename = ("~/Studio/01MetaRA/GSE89408_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        pData(.)
    G89408 = data.table::fread('~/Studio/01MetaRA/GSE89408_GEO_count_matrix_rename.txt', header = T)[, c(1, 2:29, 68:219)] %>% 
    tibble::column_to_rownames('V1')
    G89408 = log2 (edgeR::cpm(G89408) +1)
    a89408_1 = gsva(gsvaParam(as.matrix(G89408), gs_path, maxDiff = TRUE)) %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('Mixture') %>% 
        select(-c("Vitamin cofactor", "Lipid", "Nucleotide", "Carbohydrate")) 
    a89408 = data.table::fread("~/Studio/01MetaRA/CIBERSORTx_GSE89048_Results.txt", header = T)[c(1:28, 67:218), ] %>% 
        select(c(1:23)) %>% 
        tibble::column_to_rownames('Mixture') %>% 
        select(celltypeselected) %>% 
        tibble::rownames_to_column('Mixture') %>% 
        merge(a89408_1, by = 'Mixture')
    cor_matrix = cor(a89408[, -1])
    p.mat = corrplot::cor.mtest	(a89408[, -1], conf.level = 0.95)$p
    fit.p = as.data.frame(p.mat[c(1:5), c(6:8)]) %>% tibble::rownames_to_column('symbol') %>% 
        pivot_longer(!symbol,
               names_to = "group",
               values_to = 'pvalue')
    fit1 = as.data.frame(cor_matrix[c(1:5), c(6:8)]) %>% 
        tibble::rownames_to_column('symbol') %>% 
        pivot_longer(!symbol,
               names_to = "group",
               values_to = 'value')  %>% 
        cbind(fit.p[, -c(1:2)]) %>% 
        as.data.frame() %>% 
        mutate(sig = ifelse(pvalue < 0.001, '***', ifelse(pvalue < 0.01, '**', ifelse(pvalue < 0.05, '*', 'ns.')))) %>% 
        mutate(text = paste0(round(value, 2), '(', sig, ')')) %>% 
        mutate(x = c('GSE89408')) %>% 
        mutate(y = paste0(x, group))
    e93272 = read.table('~/Studio/01MetaRA/GSE93272_series.txt', header = T) %>% 
        tibble::column_to_rownames('symbol')
    a93272_1 = gsva(gsvaParam(as.matrix(e93272), gs_path, maxDiff = TRUE)) %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('Mixture') %>% 
        select(-c("Vitamin cofactor", "Lipid", "Nucleotide", "Carbohydrate")) 
    a93272 = data.table::fread("~/Studio/01MetaRA/CIBERSORTx_GSE93272_Results.txt", header = T) %>% 
        select(c(1:23)) %>% 
        tibble::column_to_rownames('Mixture') %>% 
        select(celltypeselected) %>% 
        tibble::rownames_to_column('Mixture') %>% 
        merge(a93272_1, by = 'Mixture')
    cor_matrix = cor(a93272[, -1])
    p.mat = corrplot::cor.mtest(a93272[, -1], conf.level = 0.95)$p
    fit.p = as.data.frame(p.mat[c(1:5), c(6:8)]) %>% tibble::rownames_to_column('symbol') %>% 
        pivot_longer(!symbol,
               names_to = "group",
               values_to = 'pvalue')
    fit2 = as.data.frame(cor_matrix[c(1:5), c(6:8)]) %>% 
        tibble::rownames_to_column('symbol') %>% 
        pivot_longer(!symbol,
               names_to = "group",
               values_to = 'value')  %>% 
        cbind(fit.p[, -c(1:2)]) %>% 
        as.data.frame() %>% 
        mutate(sig = ifelse(pvalue < 0.001, '***', ifelse(pvalue < 0.01, '**', ifelse(pvalue < 0.05, '*', 'ns.')))) %>% 
        mutate(text = paste0(round(value, 2), '(', sig, ')')) %>% 
        mutate(x = c('GSE93272')) %>% 
        mutate(y = paste0(x, group))
    fit0 = rbind(fit1, fit2)
    fit0$symbol = factor(fit0$symbol, levels = rev(unique(fit0$symbol)))
    fit0$y = factor(fit0$y, levels = rev(unique(fit0$y)))
    plot1 = fit0 %>% 
        ggplot(aes(x = symbol, y = reorder(y, rev(x)), fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_gradient2(low = "#D4E5F4", mid = "#FFFEEE55", high = "#FBB9BA") +
        geom_text(aes(label = text), size = 3) + 
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        scale_y_discrete(labels = rev(c('Amino acid', 'Energy', 'TCA cycle', 'Amino acid', 'Energy', 'TCA cycle'))) +
        theme_minimal() +
        labs(title = "Correction heatmap") +
        guides(fill = guide_legend(title = 'Correlation\ncoefficient')) + 
        theme(legend.position = "right",
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 10),
              panel.grid = element_blank(),
              axis.text.y = element_text(size = 12),
              axis.ticks = element_blank(),
              plot.title = element_text(face = 'bold', size = 20),
              axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
              axis.title = element_text(face = 'bold', size = 0))
    ggsave('~/Studio/01MetaRA/F001d.jpg', width = 11, height = 4)
    plot0 = fit0 %>% 
        mutate(z = 'group') %>% 
        ggplot(aes(x = z, y = x, fill = x)) +
        geom_tile(aes(fill = x)) +
        scale_fill_manual(values = c("#B7B7EB", "#FFEBAD")) +
        scale_y_discrete(position = "left") +
        theme_minimal() +
        labs(fill = '', title = '', x = '', y = '') + 
        guides(fill = 'none') + 
        theme(panel.border = element_blank(),
              panel.grid = element_blank(),
              legend.text = element_text(size = 14),
              axis.ticks = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 0),
              axis.text.x = element_blank(),
              plot.title = element_text(face = 'bold', size = 20),
              legend.position = "right")
    ggsave('~/Studio/01MetaRA/F001e.jpg', width = 1, height = 4)
## Figure2 Key genes enrolled
#### limma analysis
    eWoetzel = read.table('~/Studio/01MetaRA/Woetzel/Totalexpression.txt', header =  T)
    cWoetzel = read.table('~/Studio/01MetaRA/Woetzel/Totalclinic.txt', header =  T) %>% 
        mutate(state = ifelse(status == '0', 'NC', 'RA'))
    state = factor(cWoetzel$state, levels = c('NC', 'RA'))
    design = model.matrix(~ 0 + state)
    colnames(design) = c("NC", "RA")
    contrast.matrix = makeContrasts(RA - NC, levels = colnames(design))
    fit = eWoetzel %>% 
        lmFit(., design) %>% 
        contrasts.fit(., contrast.matrix) %>% 
        eBayes(.)
    sig = topTable(fit, n = Inf, adjust = "fdr") %>% 
        na.omit() %>% 
        dplyr::rename('log2FC'="logFC") %>% 
        mutate(group = ifelse(adj.P.Val <= 0.05 & log2FC > 1, 'Up', ifelse(adj.P.Val <= 0.05 & log2FC < - 1, 'Down', 'Not'))) %>% 
        tibble::rownames_to_column('symbol') %>% 
        mutate(log10P = -log10(adj.P.Val)) %>% 
        arrange(group)
    x2 = sig %>% filter(group == 'Up' | group == 'Down')
    require(ggpubr)
    require(ggsci)
    plot2 = ggscatter(sig, "log2FC", "log10P",
        combine = F, merge = T, color = "group", shape = 20, size = 1,
        point = TRUE, font.label = 20, palette = c("#D4E5F4", '#9E9E9E', "#FBB9BA")) +
        geom_vline(xintercept = c(-1, 1), colour = "black", linetype = "dashed") +
        geom_hline(yintercept = -log10(0.05), colour = "black", linetype = "dashed") +
        labs(title = "Woetzel\'s cohort", x = 'log2(FoldChange)', y = 'log10(adj.p.value)') + 
        theme_classic2() +
        guides(color = guide_legend(title = "")) +
        theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 1),
                  axis.text.y = element_text(size = 12),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 12),
                  legend.position = c(0.9, 0.9),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 18, face = 'bold'))
    require(patchwork)
    plot1 + plot2 + plot_layout(guides = 'collect')
    ggsave('~/Studio/01MetaRA/F004.jpg', width = 7, height = 4)
#### Metabolism_related genes from KEGG
    require(KEGGREST)
    require(stringr)
    org = keggList('organism')    
    head(org)
    org[str_detect(org[,3],"human"),]
    hsa_path = keggLink("pathway", "hsa")
    length(unique (names (hsa_path)))
    length(unique(hsa_path))
    meta2 = unique(hsa_path)[grepl('hsa00',unique(hsa_path))]
    hsa_info = lapply(meta2, keggGet)  
    gene_symbol = unlist(lapply(hsa_info, function (x) {
            g = x [[1]]$GENE
            str_split(g[seq(2, length(g), by=2)], ';', simplify = T)[, 1]}))
    gene_ID = unlist(lapply(hsa_info, function(x) {
            g = x [[1]]$GENE
            paste0("hsa:", g[seq(1, length(g), by=2)])}))
    nm = unlist(lapply(hsa_info, function(x) x[[1]]$NAME))
    genes = unlist(lapply(hsa_info, function(x) {
            g = x[[1]]$GENE
            paste(str_split(g[seq(2,length(g),by=2)],';' ,simplify = T)[,1],collapse =';')}))
    meta_genes = data.frame(hsa = meta2, name = nm, genes = genes) 
    genelist = data.frame(gene_ID, gene_symbol)
    genelist = genelist[!duplicated(genelist$gene_symbol),]
    length(genelist$gene_symbol)
    keggGENES = data.frame(symbol = c(unique(genelist$gene_symbol)))
    write.table(keggGENES, '~/Studio/01MetaRA/keggGENES.txt', row.names = T)
#### WGCNA analysis for GSE89408
    require(WGCNA)
    reactome = read.csv('~/Studio/01MetaRA/MetabolismReactome/Reactome.csv')
    metagene = data.frame(symbol = unique(reactome$symbol)) ## 2176 genes
    w89408 = data.table::fread('~/Studio/01MetaRA/GSE89408_GEO_count_matrix_rename.txt', header = T)[, c(1, 2:29, 68:219)] %>% tibble::column_to_rownames('V1') 
    w89408 = log2 (edgeR:: cpm (w89408) +1) %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('symbol') %>% 
        inner_join(metagene, by = 'symbol') %>% 
        tibble::column_to_rownames('symbol') %>% 
        t() %>% 
        as.matrix()
    gsg = goodSamplesGenes(w89408, verbose = 3) 
    gsg$allOK
    enableWGCNAThreads(nThreads = 40) # 设置线程数
    sft = pickSoftThreshold(w89408, powerVector = 1:30, # 尝试 1 到 20
                         networkType = "unsigned")
    sft$powerEstimate
    fig_power1 = ggplot(data = sft$fitIndices, aes(x = Power,y = SFT.R.sq)) +
        geom_point(color = '#FBB9BA') +
        ggrepel::geom_text_repel(aes(label = Power)) +
        geom_hline(aes(yintercept = 0.85), color = '#E7A9C5') +
        labs(title = 'Scale independence',
            x = 'Soft Threshold (power)',
            y = 'Scale Free Topology Model Fit,signed R^2') +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_fill_manual(values = color) + 
        theme(legend.position = "right",
                  axis.text.x = element_text(size = 10, angle = 0, hjust = 1),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 10),
                  axis.title = element_text(size = 12, face = 'bold'),
                  plot.title = element_text(size = 14, face = 'bold'))
    fig_power2 = ggplot(data = sft$fitIndices, aes(x = Power, y = mean.k.)) +
        geom_point(color = '#FBB9BA') +
        ggrepel::geom_text_repel(aes(label = Power)) +
        geom_hline(aes(yintercept = 100), color = '#E7A9C5') +
        labs(title = 'Mean connectivity',
            x = 'Soft Threshold (power)',
            y = 'Mean Connectivity') +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_fill_manual(values = color) + 
        theme(legend.position = "right",
                  axis.text.x = element_text(size = 10, angle = 0, hjust = 1),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 10),
                  axis.title = element_text(size = 12, face = 'bold'),
                  plot.title = element_text(size = 14, face = 'bold'))
    fig_power1 + fig_power2
    ggsave('~/Studio/01MetaRA/F002a.jpg', width = 11.75, height = 4)
#### WCGNA Net
    net = blockwiseModules(w89408,  # 输入数据
        corType = "pearson", # 相关系数算法，pearson|bicor
        power = sft$powerEstimate, # 前面得到的 soft power
        networkType = "unsigned", # unsigned | signed | signed hybrid
        TOMType = "unsigned", # none | unsigned | signed   # 计算 TOM 矩阵
        saveTOMs = TRUE,
        saveTOMFileBase = "blockwiseTOM",
        minModuleSize = 30, # 聚类并划分模块  deepSplit = 0, # 0|1|2|3|4, 值越大得到的模块就越多越小，不设置则使用默认值
        mergeCutHeight = 0.25,  #  对距离小于mergeCutHeight的模块进行合并
        numericLabels = FALSE, # 以数字命名模块
        nThreads = 0, # 0 则使用所有可用线程
        maxBlockSize = 100000 # 需要大于基因的数目
    )
    table(net$colors)    
    wgcna89408 = data.frame(gene_id = names(net$colors), module = net$colors)
    pdf(file = "~/Studio/01MetaRA/module_visualization.pdf", width = 12.2, height = 5.5)
    plotDendroAndColors(dendro = net$dendrograms[[1]], 
        colors = net$colors,
        groupLabels = "Module colors",
        dendroLabels = FALSE, hang = 0.03,
        addGuide = TRUE, guideHang = 0.05)
    dev.off()
#### Quantify module similarity by eigengene correlation
    moduleColors = (as.data.frame(net$colors))$"net$colors"
    MEs = moduleEigengenes(w89408, moduleColors)$eigengenes
    MET = orderMEs(MEs)
    pdf(file = "~/Studio/01MetaRA/eigengenes_trait_relationship.pdf", width = 6, height = 6)
    par(cex = 0.9)
    plotEigengeneNetworks(MET,"", marDendro = c(0, 4, 1, 2), 
                      marHeatmap = c(3, 4, 1, 2), cex.lab = 0.8, xLabelsAngle = 90)
    dev.off()
#### Find the relationships between modules and traits
    datatrail = data.frame(group = c(rep('Healthy control', 28), rep('Rheumatoid arthritis', 152)))
    datatrail$group = factor(datatrail$group, levels = c('Healthy control', 'Rheumatoid arthritis'))
    design = model.matrix(~ 0 + datatrail$group)
    colnames(design) = levels(datatrail$group)
    nGenes = ncol(w89408)
    nSamples = nrow(w89408)
    moduleTraitCor = cor(orderMEs(moduleEigengenes(w89408, moduleColors)$eigengenes), design, use = "p")
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
#### Display the correlation values within a heatmap plot
    pdf(file="~/Studio/01MetaRA/module_trait_relationship.pdf", width = 4, height = 10)
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")", sep = "")
    dim(textMatrix) = dim(moduleTraitCor)
    par(mar = c(6, 8.5, 3, 3))
    labeledHeatmap(Matrix = moduleTraitCor,
                   xLabels = colnames(design),
                   yLabels = names(orderMEs(moduleEigengenes(w89408, moduleColors)$eigengenes)),
                   ySymbols = names(orderMEs(moduleEigengenes(w89408, moduleColors)$eigengenes)),
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.6,
                   zlim = c(-1,1),
                   main = paste("Module-trait relationships"))
    dev.off()
#### Correction with metabolism
    require(GSVA)
    m89408_1 = gsva(gsvaParam(as.matrix(G89408), gs_path, maxDiff = TRUE)) %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>%
        select(c("TCA cycle", "Energy", "Amino acid"))
    moduleTraitCor_1 = cor(orderMEs(moduleEigengenes(w89408, moduleColors)$eigengenes), as.matrix(m89408_1), use = 'p') %>%
        as.data.frame() %>%  
        tibble::rownames_to_column('module') 
    moduleTraitPvalue_1 = corPvalueStudent(cor(orderMEs(moduleEigengenes(w89408, moduleColors)$eigengenes), as.matrix(m89408_1), use = 'p'), nSamples) %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('module')
    moduleTraitCor_2 = moduleTraitCor %>% 
        as.data.frame() %>%  
        tibble::rownames_to_column('module') %>% 
        inner_join(moduleTraitCor_1, by = 'module')
    moduleTraitPvalue_2 = moduleTraitPvalue %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('module') %>% 
        inner_join(moduleTraitPvalue_1, by = 'module')
#### ggplot2 correlation values
    Pvalue = as.data.frame(moduleTraitPvalue_2) %>% 
        pivot_longer(!module,
               names_to = "group",
               values_to = 'Pvalue') %>% 
        arrange(module)
    fit1 = as.data.frame(moduleTraitCor_2) %>% 
        pivot_longer(!module,
               names_to = "group",
               values_to = 'value') %>% 
        arrange(module) %>% 
        cbind(Pvalue[, -c(1:2)]) %>% 
        mutate(sig = ifelse(Pvalue < 0.001, '***', ifelse(Pvalue < 0.01, '**', ifelse(Pvalue < 0.05, '*', 'ns.')))) %>% 
        mutate(text = paste0(round(value, 2), '(', sig, ')'))
    fit1$module = factor(fit1$module, levels = unique(fit1$module))
    fit1$group = factor(fit1$group, levels = c('Healthy control', 'Rheumatoid arthritis', "TCA cycle", "Energy", "Amino acid"))
    plot1 = fit1 %>% 
        ggplot(aes(x = group, y = module, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_gradient2(low = "#D4E5F4", mid = "#FFFEEE55", high = "#FBB9BA") +
        geom_text(aes(label = text), size = 3) + 
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        # scale_x_discrete(position = "bottom") +
        theme_minimal() +
        labs(y = "") +
        guides(fill = guide_legend(title = 'Correlation\ncoefficient')) + 
        theme(legend.position = "right",
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 10),
              panel.grid = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 0),
              axis.ticks = element_blank(),
              axis.text.x = element_text(face = 'bold', size = 12, angle = 45, vjust = 1, hjust = 1),
              axis.title = element_text(face = 'bold', size = 0))
    plot0 = fit1 %>% 
        mutate(y = 'Module') %>% 
        mutate(names = paste0(.$module, '(', Pvalue, ')')) %>% 
        ggplot(aes(x = y, y = module, fill = module)) +
        geom_tile(aes(fill = module)) +
        scale_fill_manual(values = rev(c('#FFEBAD','#C2E3D0', '#E7A9C5', '#9e9e9e', '#85BF67', '#AA8984', '#A8CBDF'))) +
        scale_y_discrete(position = "left") +
        theme_minimal() +
        xlab(NULL) + ylab(NULL) +
        labs(fill = 'module', title = '') + 
        guides(fill = 'none') + 
        theme(panel.border = element_blank(),
              panel.grid = element_blank(),
              legend.text = element_text(size = 14),
              axis.ticks = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 14),
              axis.text.x = element_blank(),
              plot.title = element_text(face = 'bold', size = 20),
              legend.position = "right")
    require(aplot)
    plot0 %>% insert_right(plot1, width = 9)  
    ggsave('~/Studio/01MetaRA/F004.jpg', width = 5.4, height = 6.4)
#### turquoise selected
    turquoise = wgcna89408 %>% filter(module == 'turquoise')
    write.table(turquoise, "~/Studio/01MetaRA/turquoise.txt")
#### limma analysis
    require(limma)
    require(GEOquery)
    g89408c = getGEO(filename = ("~/Studio/01MetaRA/GSE89408_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        pData(.)
    G89408 = data.table::fread('~/Studio/01MetaRA/GSE89408_GEO_count_matrix_rename.txt', header = T)[, c(1, 2:29, 68:219)] %>% tibble::column_to_rownames('V1') 
    G89408 = log2 (edgeR::cpm(G89408) +1)
    status = factor(c(rep('HC', 28), rep('RA', 152)), levels = c('HC', 'RA'))
    design = model.matrix(~ 0 + status)
    colnames(design) = c("HC", "RA")
    contrast.matrix = makeContrasts(RA - HC, levels = colnames(design))
    fit = G89408 %>% 
        lmFit(., design) %>% 
        contrasts.fit(., contrast.matrix) %>% 
        eBayes(.)
    sig = topTable(fit, n = Inf, adjust = "fdr") %>% 
        na.omit() %>% 
        dplyr::rename('log2FC'="logFC") %>% 
        mutate(group = ifelse(adj.P.Val <= 0.05 & log2FC > 1, 'Up', ifelse(adj.P.Val <= 0.05 & log2FC < - 1, 'Down', 'Not'))) %>% 
        tibble::rownames_to_column('symbol') %>% 
        mutate(log10P = -log10(adj.P.Val)) %>% 
        arrange(group)
    x1 = sig %>% filter(group == 'Up' | group == 'Down')
    write.csv(x1, '~/Studio/01MetaRA/x1.csv')
    require(ggpubr)
    require(ggsci)
    plot1 = ggscatter(sig, "log2FC", "log10P",
        combine = F, merge = T, color = "group", shape = 20, size = 1,
        point = TRUE, font.label = 20, palette = c("#D4E5F4", '#9E9E9E', "#FBB9BA")) +
        geom_vline(xintercept = c(-1, 1), colour = "black", linetype = "dashed") +
        geom_hline(yintercept = -log10(0.05), colour = "black", linetype = "dashed") +
        labs(title = "GSE89408", x = 'log2(FoldChange)', y = '-log10(adj.p.value)') + 
        theme_classic2() +
        guides(color = guide_legend(title = "")) +
        theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 1),
                  axis.text.y = element_text(size = 12),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 10),
                  legend.position = 'right',
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 18, face = 'bold'))
    ggsave('~/Studio/01MetaRA/F004a.jpg', width = 3.6, height = 3.2)
#### Venn plot
    require(ggvenn)
    TCA = gs_th %>% filter(Pathway == 'TCA cycle')
    turquoise = read.table("~/Studio/01MetaRA/turquoise.txt", header = T)
    a = list(DEGs = x1$symbol, TCA_cycle = unique(TCA$symbol), turquoise_module = turquoise$gene_id)
    ggvenn::ggvenn(a, stroke_size = 0, set_name_color = "black", set_name_size = 6, show_outside = c("none"),
          fill_color = color, text_size = 3.6, auto_scale = FALSE)
    ggsave('~/Studio/01MetaRA/F004b.jpg', width = 6, height = 4) 
    gene60 = tibble(symbol = intersect(intersect(x1$symbol, TCA$symbol), turquoise$gene_id))
    write.table(gene60, '~/Studio/01MetaRA/gene60.txt', row.names = FALSE)
    write.csv(gene60, '~/Studio/01MetaRA/gene60.csv')
#### Lasso analysis
    set.seed(1234)
    gene60 = read.table('~/Studio/01MetaRA/gene60.txt', header =T)
    require(glmnet)
    require(GEOquery)
    turquoise = read.table("~/Studio/01MetaRA/turquoise.txt", header = TRUE)
    g89408c = getGEO(filename = ("~/Studio/01MetaRA/GSE89408_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        pData(.)
    G89408 = data.table::fread('~/Studio/01MetaRA/GSE89408_GEO_count_matrix_rename.txt', header = T)[, c(1, 2:29, 68:219)] %>% 
    tibble::column_to_rownames('V1')
    G89408 = log2 (edgeR::cpm(G89408) +1)
    l89408 = G89408 %>% 
        t()  %>% 
        as.data.frame() %>% 
        select(gene60$symbol) %>% 
        mutate(status = c(rep('0', 28), rep('1', 152))) %>% 
        select(status, everything())
    X = as.matrix(l89408[,-c(1)])         #this would be another way of defining X
    Y = l89408[,"status"] == "1"              #makes the outcome binary
    cv.model = cv.glmnet(x = X, y = Y, family = "binomial", alpha = 1)  #alpha=1 is lasso
    l.min = cv.model$lambda.min
    plot(cv.model)
    lasso.model = glmnet(x = X, y = Y, family = "binomial", alpha=1, lambda = l.min)
    lasso.model$beta
    head(cv.model)
    plot(cv.model$glmnet.fit, "lambda", label = FALSE)
    cv.model$glmnet.fit
    lasso_result = data.frame(symbol = unlist(lasso.model$beta@Dimnames[1]), beta = paste(lasso.model$beta)) %>% filter(beta != 0)
    line1 = as.data.frame(summary(cv.model$glmnet.fit$beta))
    line2 = data.frame(df = cv.model$glmnet.fit$df, loglambda = log(cv.model$glmnet.fit$lambda)) %>% 
        tibble::rowid_to_column('j')
    line = merge(line1, line2, by = 'j') %>% 
        mutate(group = as.factor(i)) %>% 
        ggplot(aes(x = loglambda, y = x, group = group, color = group)) +
            geom_point(size = 0) +
            geom_line() + 
            theme_classic2() +
            scale_colour_manual(values = color) + 
            geom_hline(yintercept = 0, colour = "9e9e9e", linetype = "dashed", alpha = .8) +
            labs(x = 'Log(λ)', y = 'Coefficients', title = '') +
            theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 8),
                  legend.position = 'none',
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 18, face = 'bold'))
    fit = data.frame(loglambda = log(cv.model$lambda), 
                       cv = cv.model$cvm, 
                       up = cv.model$cvup, 
                       down = cv.model$cvlo, 
                       nm = cv.model$nzero,
                       num = c(1:length(cv.model$lambda)))
    plot1 = fit %>% 
        ggplot(aes(x = loglambda, y = cv)) +
            geom_errorbar(aes(ymax = up, ymin = down), color = "#D6E2E2", width = 0.05, linewidth = 0.25) +
            geom_point(color = '#FBB9BA', size = 0.5, shape = 20) +
            theme_classic2() +
            scale_colour_manual(values = color) + 
            geom_vline(xintercept = c(fit[which.min(fit$cv),]$loglambda, fit[which(fit$num == (fit[which.min(fit$cv), ]$num - 18)), ]$loglambda), colour = "9e9e9e", linetype = "dashed", alpha = .8) +
            labs(x = 'Log(λ)', y = 'Binominal Deviance', title = '') +
            theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 8),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 18, face = 'bold'))
    require(patchwork)
    line + plot1
    ggsave('~/Studio/01MetaRA/F004c.jpg', width = 6.4, height = 3) 
#### WGCNA analysis for Woetzel
    require(WGCNA)
    reactome = read.csv('~/Studio/01MetaRA/MetabolismReactome/Reactome.csv')
    metagene = data.frame(symbol = unique(reactome$symbol))
    eWoetzel = read.table('~/Studio/01MetaRA/Woetzel/Totalexpression.txt', header =  T)
    cWoetzel = read.table('~/Studio/01MetaRA/Woetzel/Totalclinic.txt', header =  T) %>% 
        mutate(state = ifelse(status == '0', 'NC', 'RA'))
    wWoetzel = eWoetzel %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('symbol') %>% 
        inner_join(metagene, by = 'symbol') %>% 
        tibble::column_to_rownames('symbol') %>% 
        t() %>% 
        as.matrix()
    gsg = goodSamplesGenes(wWoetzel, verbose = 3) 
    gsg$allOK
    enableWGCNAThreads(nThreads = 40) # 设置线程数
    sft = pickSoftThreshold(wWoetzel, powerVector = 1:30, # 尝试 1 到 20
                         networkType = "unsigned")
    sft$powerEstimate
    fig_power1 = ggplot(data = sft$fitIndices, aes(x = Power,y = SFT.R.sq)) +
        geom_point(color = '#FBB9BA') +
        ggrepel::geom_text_repel(aes(label = Power)) +
        geom_hline(aes(yintercept = 0.85), color = '#E7A9C5') +
        labs(title = 'Scale independence',
            x = 'Soft Threshold (power)',
            y = 'Scale Free Topology Model Fit,signed R^2') +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_fill_manual(values = color) + 
        theme(legend.position = "right",
                  axis.text.x = element_text(size = 10, angle = 0, hjust = 1),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 10),
                  axis.title = element_text(size = 12, face = 'bold'),
                  plot.title = element_text(size = 14, face = 'bold'))
    fig_power2 = ggplot(data = sft$fitIndices, aes(x = Power, y = mean.k.)) +
        geom_point(color = '#FBB9BA') +
        ggrepel::geom_text_repel(aes(label = Power)) +
        geom_hline(aes(yintercept = 100), color = '#E7A9C5') +
        labs(title = 'Mean connectivity',
            x = 'Soft Threshold (power)',
            y = 'Mean Connectivity') +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_fill_manual(values = color) + 
        theme(legend.position = "right",
                  axis.text.x = element_text(size = 10, angle = 0, hjust = 1),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 10),
                  axis.title = element_text(size = 12, face = 'bold'),
                  plot.title = element_text(size = 14, face = 'bold'))
    fig_power1 + fig_power2
    ggsave('~/Studio/01MetaRA/F002b.jpg', width = 11.75, height = 4)
#### WCGNA Net
    net = blockwiseModules(wWoetzel,  # 输入数据
        corType = "pearson", # 相关系数算法，pearson|bicor
        power = sft$powerEstimate, # 前面得到的 soft power
        networkType = "unsigned", # unsigned | signed | signed hybrid
        TOMType = "unsigned", # none | unsigned | signed   # 计算 TOM 矩阵
        saveTOMs = TRUE,
        saveTOMFileBase = "blockwiseTOM",
        minModuleSize = 30, # 聚类并划分模块  deepSplit = 0, # 0|1|2|3|4, 值越大得到的模块就越多越小，不设置则使用默认值
        mergeCutHeight = 0.25,  #  对距离小于mergeCutHeight的模块进行合并
        numericLabels = FALSE, # 以数字命名模块
        nThreads = 0, # 0 则使用所有可用线程
        maxBlockSize = 100000 # 需要大于基因的数目
    )
    table(net$colors)
    wgcna_result = data.frame(gene_id = names(net$colors), module = net$colors)
    pdf(file = "~/Studio/01MetaRA/module_visualization1.pdf", width = 12.2, height = 5.5)
    plotDendroAndColors(dendro = net$dendrograms[[1]], 
        colors = net$colors,
        groupLabels = "Module colors",
        dendroLabels = FALSE, hang = 0.03,
        addGuide = TRUE, guideHang = 0.05)
    dev.off()
#### Quantify module similarity by eigengene correlation
    moduleColors = (as.data.frame(net$colors))$"net$colors"
    MEs = moduleEigengenes(wWoetzel, moduleColors)$eigengenes
    MET = orderMEs(MEs)
    pdf(file = "~/Studio/01MetaRA/eigengenes_trait_relationship1.pdf", width = 6, height = 6)
    par(cex = 0.9)
    plotEigengeneNetworks(MET,"", marDendro = c(0, 4, 1, 2), 
                      marHeatmap = c(3, 4, 1, 2), cex.lab = 0.8, xLabelsAngle = 90)
    dev.off()
#### Find the relationships between modules and traits
    datatrail = data.frame(group = cWoetzel$state)
    datatrail$group = factor(datatrail$group, levels = c('NC', 'RA'))
    design = model.matrix(~ 0 + datatrail$group)
    colnames(design) = levels(datatrail$group)
    nGenes = ncol(wWoetzel)
    nSamples = nrow(wWoetzel)
    moduleTraitCor = cor(orderMEs(moduleEigengenes(wWoetzel, moduleColors)$eigengenes), design, use = "p")
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
#### Display the correlation values within a heatmap plot
    pdf(file="~/Studio/01MetaRA/module_trait_relationship1.pdf", width = 4, height = 10)
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")", sep = "")
    dim(textMatrix) = dim(moduleTraitCor)
    par(mar = c(6, 8.5, 3, 3))
    labeledHeatmap(Matrix = moduleTraitCor,
                   xLabels = colnames(design),
                   yLabels = names(orderMEs(moduleEigengenes(wWoetzel, moduleColors)$eigengenes)),
                   ySymbols = names(orderMEs(moduleEigengenes(wWoetzel, moduleColors)$eigengenes)),
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.6,
                   zlim = c(-1,1),
                   main = paste("Module-trait relationships"))
    dev.off()
#### Correction with metabolism
    require(GSVA)
    mWoetzel_1 = gsva(gsvaParam(as.matrix(eWoetzel), gs_path, maxDiff = TRUE)) %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>%
        select(c("Energy", "Amino acid"))
    moduleTraitCor_1 = cor(orderMEs(moduleEigengenes(wWoetzel, moduleColors)$eigengenes), as.matrix(mWoetzel_1), use = 'p') %>%
        as.data.frame() %>%  
        tibble::rownames_to_column('module') 
    moduleTraitPvalue_1 = corPvalueStudent(cor(orderMEs(moduleEigengenes(wWoetzel, moduleColors)$eigengenes), as.matrix(mWoetzel_1), use = 'p'), nSamples) %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('module')
    moduleTraitCor_2 = moduleTraitCor %>% 
        as.data.frame() %>%  
        tibble::rownames_to_column('module') %>% 
        inner_join(moduleTraitCor_1, by = 'module') %>% 
        dplyr::rename('Healthy control' = 'NC', 'Rheumatoid arthritis' = 'RA')
    moduleTraitPvalue_2 = moduleTraitPvalue %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('module') %>% 
        inner_join(moduleTraitPvalue_1, by = 'module')
#### ggplot2 correlation values
    Pvalue = as.data.frame(moduleTraitPvalue_2) %>% 
        pivot_longer(!module,
               names_to = "group",
               values_to = 'Pvalue') %>% 
        arrange(module)
    fit1 = as.data.frame(moduleTraitCor_2) %>% 
        pivot_longer(!module,
               names_to = "group",
               values_to = 'value') %>% 
        arrange(module) %>% 
        cbind(Pvalue[, -c(1:2)]) %>% 
        mutate(sig = ifelse(Pvalue < 0.001, '***', ifelse(Pvalue < 0.01, '**', ifelse(Pvalue < 0.05, '*', 'ns.')))) %>% 
        mutate(text = paste0(round(value, 2), '(', sig, ')'))
    fit1$module = factor(fit1$module, levels = unique(fit1$module))
    fit1$group = factor(fit1$group, levels = c('Healthy control', 'Rheumatoid arthritis', "Energy", "Amino acid"))
    plot1 = fit1 %>% 
        ggplot(aes(x = group, y = module, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_gradient2(low = "#D4E5F4", mid = "#FFFEEE55", high = "#FBB9BA") +
        geom_text(aes(label = text), size = 2.5) + 
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        # scale_x_discrete(position = "bottom") +
        theme_minimal() +
        labs(y = "") +
        guides(fill = guide_legend(title = 'Correlation\ncoefficient')) + 
        theme(legend.position = "right",
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 10),
              panel.grid = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 0),
              axis.ticks = element_blank(),
              axis.text.x = element_text(face = 'bold', size = 12, angle = 45, vjust = 1, hjust = 1),
              axis.title = element_text(face = 'bold', size = 0))
    plot0 = fit1 %>% 
        mutate(y = 'Module') %>% 
        mutate(names = paste0(.$module, '(', Pvalue, ')')) %>% 
        ggplot(aes(x = y, y = module, fill = module)) +
        geom_tile(aes(fill = module)) +
        scale_fill_manual(values = rev(c('#FFEBAD','#C2E3D0', '#E7A9C5', '#E3BBED', '#FBB9BA', '#F0988C', '#9e9e9e', '#85BF67', '#AA8984', '#A8CBDF', '#2d2d2d'))) +
        scale_y_discrete(position = "left") +
        theme_minimal() +
        xlab(NULL) + ylab(NULL) +
        labs(fill = 'module', title = '') + 
        guides(fill = 'none') + 
        theme(panel.border = element_blank(),
              panel.grid = element_blank(),
              legend.text = element_text(size = 14),
              axis.ticks = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 14),
              axis.text.x = element_blank(),
              plot.title = element_text(face = 'bold', size = 20),
              legend.position = "right")
    require(aplot)
    plot0 %>% insert_right(plot1, width = 6)  
    ggsave('~/Studio/01MetaRA/F005a.jpg', width = 4.5, height = 8)
#### ROC for wgcna
    z89408 = as.data.frame(orderMEs(moduleEigengenes(w89408, moduleColors)$eigengenes)) %>% 
        tibble::rownames_to_column('sample') %>% 
        mutate(status = c(rep('0', 28), rep('1', 95))) 
    require(pROC)
    require(plotROC)
    x0 = data.frame()
    for(i in 2:8){
        z89408$status = factor(z89408$status, levels = c('1', '0'))
        x1 = roc(response = z89408$status, predictor = z89408[, i])
        x2 = round(auc(response = z89408$status, predictor = z89408[, i]), 4)
        x3 = data.frame(symbol = names(z89408)[i], value = x2)
        x0 = rbind(x3, x0)}
    wa89408 = x0
    y1 = z89408[, - 9] %>% 
        pivot_longer(!sample,
            names_to = "Module",
            values_to = 'value')  %>% 
        mutate(status = c(rep('0', 28*7), rep('1', 95*7))) %>%  
        group_split(Module)
    names(y1) = sort(colnames(z89408)[2:8])
    map_function <- function(y1){
        pROC::roc(
        data = y1,
        response = status,
        predictor = value)}
    roc.list = map(y1, map_function)
    plot1 = pROC::ggroc(roc.list, alpha = 1, size = 0.75, legacy = T) +
        geom_abline(linetype = 3, intercept = 0) +
        labs(title = "GSE89408", color = 'Module') +
        theme_classic() +
        scale_color_manual(values = rev(c('#FFEBAD','#C2E3D0', '#E7A9C5', '#9e9e9e', '#85BF67', '#AA8984', '#A8CBDF'))) +
        theme(axis.text = element_text(size = 10, vjust = 0.5),
              # legend.position = 'left',
              legend.title = element_text(size = 14, face = 'bold'),
              legend.text = element_text(size = 12),
              axis.title = element_text(size = 14, face = 'bold'),
              plot.title = element_text(size = 14, face = 'bold'))    
    zWoetzel = as.data.frame(orderMEs(moduleEigengenes(wWoetzel, moduleColors)$eigengenes)) %>% 
        tibble::rownames_to_column('sample') %>% 
        inner_join(as.data.frame(cWoetzel[, 1:2]), by = 'sample') %>% 
        mutate(group = ifelse(status == 1, '1', '0')) %>% 
        select(-c('status'))
    require(pROC)
    require(plotROC)
    x0 = data.frame()
    for(i in 2:12){
        zWoetzel$group = factor(zWoetzel$group, levels = c('1', '0'))
        x1 = roc(response = zWoetzel$group, predictor = zWoetzel[, i])
        x2 = round(auc(response = zWoetzel$group, predictor = zWoetzel[, i]), 4)
        x3 = data.frame(symbol = names(zWoetzel)[i], value = x2)
        x0 = rbind(x3, x0)}
    waWoetzel = x0
    y1 = zWoetzel[, - 13] %>% 
        pivot_longer(!sample,
            names_to = "Module",
            values_to = 'value')  %>% 
        inner_join(cWoetzel %>% select(c('sample', 'status')), by = 'sample') %>% 
        mutate(group = ifelse(status == 1, '1', '0')) %>% 
        select(-c('status')) %>%  
        group_split(Module)
    names(y1) = sort(colnames(zWoetzel)[2:12])
    map_function <- function(y1){
        pROC::roc(
        data = y1,
        response = group,
        predictor = value)}
    roc.list = map(y1, map_function)
    plot2 = pROC::ggroc(roc.list, alpha = 1, size = 0.75, legacy = T) +
        geom_abline(linetype = 3, intercept = 0) +
        labs(title = "Woetzel\'s cohort", color = 'Module') +
        theme_classic() +
        scale_color_manual(values = rev(c('#FFEBAD','#C2E3D0', '#E7A9C5', '#E3BBED', '#FBB9BA', '#F0988C', '#9e9e9e', '#85BF67', '#AA8984', '#A8CBDF', '#2d2d2d'))) +
        theme(axis.text = element_text(size = 10, vjust = 0.5),
              # legend.position = 'left',
              legend.title = element_text(size = 14, face = 'bold'),
              legend.text = element_text(size = 12),
              axis.title = element_text(size = 14, face = 'bold'),
              plot.title = element_text(size = 14, face = 'bold'))    
    require(patchwork)
    plot1 + plot2
    ggsave('~/Studio/01MetaRA/F006.jpg', width = 7.2, height = 4)
    write.table(wa89408, '~/Studio/01MetaRA/wgcna89408.txt')
    write.table(waWoetzel, '~/Studio/01MetaRA/wgcnawoetzel.txt')
#### GSE89408 Differential expression analysis
    require(GEOquery)
    key = data.frame(symbol = factor(c("COX7C", "COX7B", "PDK1", "PDK3", "COX11", "NDUFAF4", "NDUFB8", "NDUFS1", "SUCLG1", "PDHX", "SUCLA2"), levels = c("COX7C", "COX7B", "PDK1", "PDK3", "COX11", "NDUFAF4", "NDUFB8", "NDUFS1", "SUCLG1", "PDHX", "SUCLA2")))
    g89408c = getGEO(filename = ("~/Studio/01MetaRA/GSE89408_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        pData(.)
    G89408 = data.table::fread('~/Studio/01MetaRA/GSE89408_GEO_count_matrix_rename.txt', header = T)[, c(1, 2:29, 68:219)] %>% 
    tibble::column_to_rownames('V1')
    G89408 = log2 (edgeR::cpm(G89408) +1)
    df89408 = G89408 %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('symbol') %>% 
        inner_join(key, by = 'symbol') %>% 
        tibble::column_to_rownames('symbol') %>% 
        # select(- c('beta')) %>% 
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('gse')  %>% 
        pivot_longer(!gse,
            names_to = "symbol",
            values_to = 'value') %>% 
        mutate(group = factor(c(rep('Healthy control', 28*11), rep('Rheumatoid arthritis', 152*11)), levels = c('Healthy control', 'Rheumatoid arthritis')))
    plot3 = df89408 %>% 
        ggplot(aes(x = symbol, y = value, fill = group)) + 
            geom_boxplot(width = 0.5, position = position_dodge(0.75), size = 0.1, outlier.size = 0.25, outlier.alpha = 0.5) +
            theme_classic2() +
            stat_compare_means(label = "p.signif", method = "wilcox.test", bracket.size = 0.3) +
            # facet_wrap(~ symbol, nrow = 1) +
            # ggsignif::geom_signif(aes(x = symbol, y = value, group = group), map_signif_level = F, parse = T, y_position = 9.5, comparisons = list(c('Rheumatoid arthritis', 'Healthy control')), size = 0.5, textsize = 2.5, vjust = 0.1, tip_length = 0.01) +
            labs(title = "GSE89408", y = "Gene expression", x = "") +
            scale_fill_manual(values = color) + 
            theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 8),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 18, face = 'bold'))
#### ROC
    require(plotROC) 
    require(pROC)
    r89408 = G89408 %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('symbol') %>% 
        inner_join(key, by = 'symbol') %>% 
        tibble::column_to_rownames('symbol') %>% 
        # select(- c('beta')) %>% 
        t() %>% 
        as.data.frame() %>% 
        select(c(COX11, COX7B, COX7C, NDUFAF4, NDUFB8, NDUFS1, PDHX, PDK1, PDK3, SUCLA2, SUCLG1)) %>% 
        mutate(group = c(rep('Healthy control', 28), rep('Rheumatoid arthritis', 152)))
    x0 = data.frame()
    for(i in 1:11){
        r89408$gse = factor(r89408$group, levels = c('Rheumatoid arthritis', 'Healthy control'))
        x1 = roc(response = r89408$group, predictor = r89408[, i])
        x2 = round(auc(response = r89408$group, predictor = r89408[, i]), 4)
        x3 = data.frame(symbol = names(r89408)[i], value = x2)
        x0 = rbind(x3, x0)}
    aaa89408 = x0
    y89408 = r89408[, - 12] %>% 
        pivot_longer(!gse,
            names_to = "Meta",
            values_to = 'value') %>% 
        mutate(group = c(rep('Healthy control', 28*11), rep('Rheumatoid arthritis', 152*11))) %>% 
        group_split(Meta)
    str(y89408)
    map_function <- function(y89408){
        pROC::roc(
        data = y89408,
        response = gse,
        predictor = value)}
    roc.list = map(y89408, map_function)
    names(roc.list) = colnames(r89408)[1:11]
    plot4 = pROC::ggroc(roc.list, alpha = 1, size = 0.5, legacy = T) +
        geom_abline(linetype = 3, intercept = 0) +
        labs(title = "GSE89408", color = 'symbol') +
        theme_classic() +
        scale_color_manual(values = color) +
        theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 10),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 18, face = 'bold'))
    require(patchwork)    
    plot3 + plot4 + plot_layout(design = 'AAAAABBB')
    ggsave('~/Studio/01MetaRA/F004x.jpg', width = 11.5, height = 4)
#### Difference analysis for GSE89408 2e
    fit2 = as.data.frame(G89408) %>% 
        tibble::rownames_to_column('symbol') %>% 
        inner_join(key, by = 'symbol') %>% 
        # select(- c('beta')) %>% 
        tibble::column_to_rownames('symbol') %>% 
        t() %>% 
        as.data.frame() %>% 
        mutate_all(as.numeric) %>% 
        mutate(group = c(rep('Healthy control', 28), rep('Early rheumatoid arthritis', 57), rep('Established rheumatoid arthritis', 95))) 
    fit2$group = factor(fit2$group, levels = c('Healthy control', 'Early rheumatoid arthritis', 'Established rheumatoid arthritis'))
    plot2 = fit2 %>% 
        pivot_longer(!group, 
            names_to = 'symbol', 
            values_to = 'value') %>% 
            ggplot(aes(x = symbol, y = value, fill = group)) + 
            geom_boxplot(width = 0.5, position = position_dodge(0.75), size = 0.1, outlier.size = 0.25, outlier.alpha = 0.5) +
            theme_classic2() +
            stat_compare_means(label = "p.signif", method = "anova", bracket.size = 0.2, size = 2.5) +
            # facet_wrap(~ symbol, nrow = 1) +
            # ggsignif::geom_signif(aes(x = symbol, y = value, group = group), map_signif_level = F, parse = T, y_position = 9.5, comparisons = list(c('Rheumatoid arthritis', 'Healthy control')), size = 0.5, textsize = 2.5, vjust = 0.1, tip_length = 0.01) +
            labs(title = "GSE894048", y = "Gene expression", x = "") +
            scale_fill_manual(values = color, labels = c('Healthy control', 'Early RA', 'Established RA')) + 
            guides(fill = guide_legend(title = '')) +
            theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 9),
                  legend.position = 'right',
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 18, face = 'bold'))
    ggsave('~/Studio/01MetaRA/F004y.jpg', width = 5.5, height = 3.6)
#### ROC 2e
    require(plotROC) 
    require(pROC)
    r89408 = G89408 %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('symbol') %>% 
        inner_join(key, by = 'symbol') %>% 
        tibble::column_to_rownames('symbol') %>% 
        # select(- c('beta')) %>% 
        t() %>% 
        as.data.frame() %>% 
        select(c(COX11, COX7B, COX7C, NDUFAF4, NDUFB8, NDUFS1, PDHX, PDK1, PDK3, SUCLA2, SUCLG1)) %>% 
        mutate(group = c(rep('Healthy control', 28), rep('Early rheumatoid arthritis', 57), rep('Established rheumatoid arthritis', 95))) 
    x0 = data.frame()
    for(i in 1:11){
        r89408ea = r89408 %>% filter(group == 'Healthy control'| group == 'Early rheumatoid arthritis')
        r89408ea$group = factor(r89408ea$group, levels = c('Early rheumatoid arthritis', 'Healthy control'))
        x1 = roc(response = r89408ea$group, predictor = r89408ea[, i])
        x2 = round(auc(response = r89408ea$group, predictor = r89408ea[, i]), 4)
        x3 = data.frame(symbol = names(r89408ea)[i], value = x2)
        x0 = rbind(x3, x0)}
    aab89408 = x0
    y89408 = (r89408 %>% filter(group == 'Healthy control'| group == 'Early rheumatoid arthritis')) %>% 
        pivot_longer(!group,
            names_to = "Meta",
            values_to = 'value') %>% 
        # mutate(group = c(rep('Healthy control', 28*11), rep('Early rheumatoid arthritis', 57*11))) %>% 
        group_split(Meta)
    str(y89408)
    map_function <- function(y89408){
        pROC::roc(
        data = y89408,
        response = group,
        predictor = value)}
    roc.list = map(y89408, map_function)
    names(roc.list) = colnames(r89408)[1:11]
    plot4a = pROC::ggroc(roc.list, alpha = 1, size = 0.75, legacy = T) +
        geom_abline(linetype = 3, intercept = 0) +
        labs(title = "Early RA", color = 'symbol') +
        theme_classic() +
        scale_color_manual(values = color) +
        guides(color = guide_legend(title = '')) +
        theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 9),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 18, face = 'bold'))
    x0 = data.frame()
    for(i in 1:11){
        r89408ea = r89408 %>% filter(group == 'Healthy control'| group == 'Established rheumatoid arthritis')
        r89408ea$group = factor(r89408ea$group, levels = c('Established rheumatoid arthritis', 'Healthy control'))
        x1 = roc(response = r89408ea$group, predictor = r89408ea[, i])
        x2 = round(auc(response = r89408ea$group, predictor = r89408ea[, i]), 4)
        x3 = data.frame(symbol = names(r89408ea)[i], value = x2)
        x0 = rbind(x3, x0)}
    aac89408 = x0
    y89408 = (r89408 %>% filter(group == 'Healthy control'| group == 'Established rheumatoid arthritis')) %>% 
        pivot_longer(!group,
            names_to = "Meta",
            values_to = 'value') %>% 
        #mutate(group = c(rep('Healthy control', 28*11), rep('Established rheumatoid arthritis', 95*11))) %>% 
        group_split(Meta)
    str(y89408)
    map_function <- function(y89408){
        pROC::roc(
        data = y89408,
        response = group,
        predictor = value)}
    roc.list = map(y89408, map_function)
    names(roc.list) = colnames(r89408)[1:11]
    plot4b = pROC::ggroc(roc.list, alpha = 1, size = 0.75, legacy = T) +
        geom_abline(linetype = 3, intercept = 0) +
        labs(title = "Established RA", color = 'symbol') +
        theme_classic() +
        scale_color_manual(values = color) +
        guides(color = guide_legend(title = '')) +
        theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 9),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 18, face = 'bold'))
    plot4a + plot4b + plot_layout(guides = 'collect')
    ggsave('~/Studio/01MetaRA/F004s.jpg', width = 6.3, height = 3.3)
    merge(merge(aaa89408 %>% dplyr::rename('RA' = 'value'), aab89408 %>% dplyr::rename('Early_RA' = 'value'), by = 'symbol'), aac89408 %>% dplyr::rename('Established_RA' = 'value'), by = 'symbol') %>% write.csv(., "~/Studio/01MetaRA/89408AUC.csv")
#### corplot
    require(GSVA)
    key = data.frame(symbol = c("COX7C", "COX7B", "PDK1", "PDK3", "COX11", "NDUFAF4", "NDUFB8", "NDUFS1", "SUCLG1", "PDHX", "SUCLA2"))
    g89408c = getGEO(filename = ("~/Studio/01MetaRA/GSE89408_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        pData(.)
    G89408 = data.table::fread('~/Studio/01MetaRA/GSE89408_GEO_count_matrix_rename.txt', header = T)[, c(1, 2:29, 68:219)] %>% 
    tibble::column_to_rownames('V1')
    G89408 = log2 (edgeR::cpm(G89408) +1)
    z89408_1 = as.data.frame(G89408) %>% 
        tibble::rownames_to_column('symbol') %>% 
        inner_join(key, by = 'symbol') %>% 
        # select(- c('beta')) %>% 
        tibble::column_to_rownames('symbol') %>% 
        t() %>% 
        as.data.frame() %>% 
        mutate_all(as.numeric)
    z89408_2 = gsva(gsvaParam(as.matrix(G89408), gs_path, maxDiff = TRUE)) %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        select('TCA cycle')
    z89408_3 = data.table::fread("~/Studio/01MetaRA/CIBERSORTx_GSE89048_Results.txt", header = T)[c(1:28, 67:218), ] %>% 
        select(c(1:23)) %>% 
        dplyr::rename('sample' = 'Mixture') %>% 
        select(celltypeselected) %>% 
        mutate(group = c(rep('Healthy control', 28), rep('Rheumatoid arthritis', 152)))
    z89408_0 = cbind(z89408_1, z89408_2, z89408_3) %>% 
        select(- c('group')) 
    cor_matrix = cor(z89408_0)
    p.mat = corrplot::cor.mtest(z89408_0, conf.level = 0.95)$p
    fit89408_1 = cor_matrix[1:11,12:17] %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('symbol') %>% 
        pivot_longer(!symbol, names_to = 'module', values_to = 'value') %>% 
        mutate(module = str_remove(module, "^value\\."))
    fit89408_2 = p.mat[1:11,12:17] %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('symbol') %>% 
        pivot_longer(!symbol, names_to = 'module', values_to = 'p') %>% 
        mutate(module = str_remove(module, "^value\\."))
    fit89408_0 = cbind(fit89408_1, pvalue = fit89408_2$p)
    fit89408_0$symbol = factor(fit89408_0$symbol, levels = rev(c("COX11", "COX7B", "COX7C", "NDUFAF4", "NDUFB8", "NDUFS1", "PDHX", "PDK1", "PDK3", "SUCLA2", "SUCLG1")))
    fit89408_0$module = factor(fit89408_0$module, levels = c("TCA cycle", "Neutrophils", "Monocytes", "NK cells resting", "T cells regulatory (Tregs)", "T cells CD4 memory activated"))
    plot1 = fit89408_0 %>% 
        mutate(sig = ifelse(pvalue < 0.001, '***', ifelse(pvalue < 0.01, '**', ifelse(pvalue < 0.05, '*', 'ns.')))) %>% 
        mutate(text = paste0(round(value, 2), '(', sig, ')'))  %>% 
        ggplot(aes(x = module, y = symbol, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_gradient2(low = "#D4E5F4", mid = "#FFFEEE", high = "#FBB9BA") +
        geom_text(aes(label = text), size = 3) + 
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        # scale_x_discrete(position = "bottom") +
        theme_minimal() +
        labs(title = "Correction heatmap") +
        guides(fill = guide_legend(title = 'Correlation\ncoefficient')) + 
        theme(legend.position = "right",
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 10),
              panel.grid = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 12),
              axis.ticks = element_blank(),
              plot.title = element_text(face = 'bold', size = 20),
              axis.text.x = element_text(face = 'bold', size = 12, angle = 45, vjust = 1, hjust = 1),
              axis.title = element_text(face = 'bold', size = 0))
    ggsave('~/Studio/01MetaRA/F004z.jpg', width = 11.5, height = 6)
## Figure3 Functional enrichment of Key genes
#### Functional enrichment
    GO = read.table('~/Studio/01MetaRA/11GO.txt', sep = '\t', header = T) %>% 
        mutate(names = substr(.$Term, start = 12, nchar(.$Term))) %>% 
        mutate(num = c(1:dim(.)[1])) %>% 
        mutate(x = str_sub(.$Category, end = - 8, start = 8)) %>% 
        ggplot(aes(x = reorder(names, rev(num)), y = Count, fill = x)) +
            geom_bar(stat='identity') +
            labs(title = "GO analysis", color = 'x', x = '', y = 'Numbers of enriched') +
            theme_classic2() +
            scale_y_continuous(expand = c(0, 0), breaks = seq(0, 12, 4)) +
            scale_fill_manual(values = color) +
            guides(fill = guide_legend(title = '')) +
            coord_flip() + 
            theme(axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5),
                  axis.text.y = element_text(size = 8, angle = 0),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 8),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 18, face = 'bold'))
    read.table('~/Studio/01MetaRA/11GO.txt', sep = '\t', header = T) %>% write.csv('~/Studio/01MetaRA/11GOx.csv')
    KEGG = read.table('~/Studio/01MetaRA/11pathway.txt', sep = '\t', header = T) %>% 
        mutate(num = c(1:dim(.)[1])) %>% 
        mutate(x = str_sub(.$Category, end = - 9, start = 0))
    x1 = KEGG$Term[1:15] %>% 
        substr(., start = 10, nchar(.))
    x2 = KEGG$Term[16:29] %>% 
        substr(., start = 15, nchar(.))
    pKEGG = KEGG %>% 
        mutate(names = c(x1, x2)) %>% 
        dplyr::rename('enrichmentRatio' = 'X.') %>% 
        ggplot(aes(x = reorder(names, rev(num)), y = enrichmentRatio, fill = x)) +
            geom_bar(stat='identity') +
            theme_classic2() +
            scale_y_continuous(expand = c(0, 0)) +
            scale_fill_manual(values = color) +
            # scale_fill_gradient2(low = "#D4E5F4", mid = "#FFFEEE", high = "#FBB9BA") +
            guides(fill = guide_legend(title = '')) +
            labs(title = "Pathway analysis", x = '', y = 'enrichmentRatio') +
            coord_flip() + 
            theme(axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5),
                  axis.text.y = element_text(size = 8, angle = 0),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 8),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 18, face = 'bold'))
    read.table('~/Studio/01MetaRA/11pathway.txt', sep = '\t', header = T) %>% write.csv('~/Studio/01MetaRA/11pathwayx.csv')
    GO + pKEGG
    ggsave('~/Studio/01MetaRA/F005.jpg', width = 11.8, height = 5)
#### GSEA analysis 
    head(x1) #### result from limma analysis
    x1 %>% 
        inner_join(key, by = 'symbol') %>% 
        select(c('symbol', 'log2FC')) %>% 
        write.table(., '~/Studio/01MetaRA/gesa.txt', sep = '\t', row.names = F)
    gesa = data.table::fread('~/Studio/01MetaRA/11gsea.txt', header = T) %>% 
        mutate(x = c(1:dim(.)[1])) %>% 
        ggplot(aes(x = reorder(Description, rev(x)), y = NES, fill = Type)) +
            geom_bar(stat='identity') + 
            theme_classic2() +
            # scale_y_continuous(expand = c(0, 0)) +
            scale_fill_manual(values = color) +
            guides(fill = guide_legend(title = '')) +
            labs(title = "GSEA", x = '', y = 'NES') +
            scale_x_discrete(labels = c("Pyruvate metabolism\n and Citric Acid (TCA) cycle", "Respiratory electron transport","Respiratory electron transport, \nATP synthesis by chemiosmotic coupling, \nand heat production by uncoupling proteins", "Metabolic pathways", "Oxidative phosphorylation", "Thermogenesis", "Diabetic cardiomyopathy")) +
            coord_flip() + 
            theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
                  axis.text.y = element_text(size = 12, angle = 0),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 8),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 18, face = 'bold'))
    ggsave('~/Studio/01MetaRA/F005c.jpg', width = 6, height = 5)
    read.table("~/Studio/01MetaRA/11GSEA.txt", sep = "\t", header = T) %>% write.csv(., "~/Studio/01MetaRA/11GSEA.csv")
#### Enrichr
    enrichr = data.table::fread('~/Studio/01MetaRA/enrichr.txt', header = T, sep = '\t') %>% 
        filter(Adjusted_p_value < 0.05) %>% 
        tidyr::separate_rows(gene, sep = ', ')
    unique(enrichr$gene)
    plot01 = enrichr %>% 
        mutate(symbol = factor(gene, levels = c("SUCLG1", "SUCLA2", "PDK3", "PDK1", "PDHX", "NDUFS1", "NDUFB8", "NDUFAF4", "COX7C", "COX7B", "COX11"))) %>% 
        dplyr::rename('score' = 'Combined score') %>% 
        ggplot(aes(x = reorder(Name, Index), y = symbol)) + 
        geom_point(aes(size = score, color = Adjusted_p_value)) +
        # geom_tile(aes(fill = Adjusted_p_value)) +
        scale_color_gradient(high = "#D4E5F4", low = "#FBB9BA") +
        scale_size_continuous(range = c(2, 8)) + 
        # geom_text(aes(label = score), size = 3) + 
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        # scale_x_discrete(position = "bottom") +
        theme_minimal() +
        labs(y = "", title = 'DSigDB') +
        guides(fill = guide_legend(title = 'Adjusted\np-value')) + 
        theme(legend.position = "right",
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 10),
              panel.grid = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 12),
              axis.ticks = element_blank(),
              axis.text.x = element_text(face = 'bold', size = 12, angle = 45, vjust = 1, hjust = 1),
              axis.title = element_text(face = 'bold', size = 0),
              plot.title = element_text(size = 18, face = 'bold'))
    data.table::fread('~/Studio/01MetaRA/enrichr.txt', header = T, sep = '\t') %>% 
        filter(Adjusted_p_value < 0.05) %>% 
        write.csv(., "~/Studio/01MetaRA/enrichr2.csv")
    require(patchwork)
    plot01 + plot_layout(design = '#AAAAAAAAAAAAAAAAAAAAAAAAAA')
    ggsave('~/Studio/01MetaRA/F010a.jpg', width = 12, height = 5.5)
#### newAUC
    read.table('~/Studio/01MetaRA/newAUC.txt') %>% write.csv('~/Studio/01MetaRA/newAUC.csv')
    metagene = data.frame(symbol = unique(reactome$symbol)) %>% write.csv('~/Studio/01MetaRA/metagene.csv')
#### Turquoise 
    read.table('~/Studio/01MetaRA/turquoise.txt') %>% write.csv('~/Studio/01MetaRA/turquoise.csv')
## Figure4 Machine learning and correction
#### Multi-Machine learning 
    set.seed(1234)
    key = c("COX7C", "COX7B", "PDK1", "PDK3", "COX11", "NDUFAF4", "NDUFB8", "NDUFS1", "SUCLG1", "PDHX", "SUCLA2")
    library(tidymodels)  # 核心框架
    library(discrim)     # 判别分析
    library(bonsai)      # 增强树模型
    library(baguette)    # 装袋法
    library(rules)       # 规则模型
    m89408 = G89408 %>% 
        t() %>% 
        as.data.frame() %>% 
        dplyr::select(key) %>% 
        mutate(status = factor(c(rep('0', 28), rep('1', 152)), levels = c('0', '1'))) %>% 
        dplyr::select(status, everything())
    rec =  recipe(status ~ ., data = m89408) %>%
        step_normalize(all_numeric_predictors()) %>%   # 标准化数值
        step_zv(all_predictors()) %>%                 # 移除零方差特征
        step_corr(all_numeric_predictors(), threshold = 0.9)
    xgb_mod = boost_tree() %>% set_engine("xgboost") %>% set_mode("classification")   
    dt_mod = decision_tree() %>% set_engine("rpart") %>% set_mode("classification")
    logistic_mod =  logistic_reg() %>% set_engine('glm')
    neural_net = mlp() %>% set_engine('nnet') %>% set_mode('classification')
    naivebayes_mod = naive_Bayes() %>% set_engine('naivebayes')
    kknn_mod = nearest_neighbor() %>% set_engine('kknn') %>% set_mode('classification')
    rf_mod = rand_forest() %>% set_engine('ranger') %>% set_mode('classification')
    svm_mod = svm_rbf() %>% set_engine('kernlab') %>% set_mode('classification')
    c50_spec = C5_rules() %>% set_engine("C5.0") %>% set_mode("classification")
#### lasso model
    set.seed(1234)
    require(glmnet)
    l89408 = G89408 %>% 
        t()  %>% 
        as.data.frame() %>% 
        dplyr::select(key$symbol) %>% 
        mutate(status = c(rep('0', 28), rep('1', 152))) %>% 
        dplyr::select(status, everything())
    X = as.matrix(l89408[,-c(1)])         #this would be another way of defining X
    Y = l89408[,"status"] == "1"              #makes the outcome binary
    cv.model = cv.glmnet(x = X, y = Y, family = "binomial", alpha = 1)  #alpha=1 is lasso
    preds = predict(cv.model, newx = X, type = 'response')
    lasso.model = glmnet(x = X, y = Y, family = "binomial", alpha = 1, lambda = cv.model$lambda.min)
    lasso.model$beta
    write.csv(data.frame(symbol = unlist(lasso.model$beta@Dimnames[1]), beta = paste(lasso.model$beta)), "~/Studio/01MetaRA/lasso.csv")
    assess.glmnet(lasso.model, newx = X, newy = Y )
    lasso = data.frame(.pred_1 = preds, status = l89408$status, Machine_learning = 'Lasso logistic regression: 0.9995', model = 'lasso') %>% dplyr::rename('.pred_1' = 'lambda.1se')
    lasso$status = factor(lasso$status, levels = c('0', '1'))
#### workflow
    wf = workflow_set(preproc=list(rec),          
                       models=list(xgb = xgb_mod, 
                       dt = dt_mod, 
                       log = logistic_mod, 
                       nb = naivebayes_mod, 
                       nnet = neural_net, 
                       knn = kknn_mod, 
                       rf = rf_mod, 
                       svm = svm_mod, 
                       C50 = c50_spec))
    folds = bootstraps(m89408, 1000)
    ctr = control_resamples(save_pred = TRUE)
    wf_res = wf %>%  workflow_map("fit_resamples", #重采样拟合
               resamples = folds, # 重采样 
               control = ctr)
    auc = rank_results(wf_res,rank_metric = "roc_auc") %>% 
        filter(.metric=="roc_auc") %>% 
        dplyr::select(model, mean) 
    write.csv(auc, '~/Studio/01MetaRA/auc11.csv')
    fit = collect_predictions(wf_res) %>%           
        dplyr::group_by(model) %>%   
        dplyr::select(c('.pred_1', 'status', 'model')) %>% 
        mutate('Machine_learning' = case_when(model == 'boost_tree' ~ 'Boost tree: 0.9723',
                                              model == 'decision_tree' ~ 'Decision tree: 0.8565',
                                              model == 'C5_rules' ~ 'C5.0 rule: 0.8502',
                                              model == 'logistic_reg' ~ 'Logistic regression: 0.9655',
                                              model == 'naive_Bayes' ~ 'Naive Bayes: 0.9568',
                                              model == 'mlp' ~ 'Multilayer perceptron: 0.9615',
                                              model == 'nearest_neighbor' ~ 'Nearest neighbor: 0.9407',
                                              model == 'rand_forest' ~ 'Random forest: 0.9766',
                                              model == 'svm_rbf' ~ 'SVM-RFE: 0.9835'))
    fit = rbind(fit, lasso)
    plot0 = fit %>% 
        ggplot(aes(m = .pred_1, d = as.integer(status), color = Machine_learning)) + 
            geom_roc(labels = F, pointsize = 0.3, size = 0.5) +    
            theme_classic2() + 
            scale_color_manual(values = color) +
            labs(x= "1-Specificity", y= "Sensitivity", title = 'Machine learning') + 
            geom_abline(linetype = 3, intercept = 0) +
            guides(fill = guide_legend(title = 'Machine learning')) +
            theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 12),
                  legend.position = 'right',
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 16, face = 'bold'))
    ggsave('~/Studio/01MetaRA/F006.jpg', width = 7, height = 4)
#### Lasso model plot
    plot5 = data.frame(symbol = unlist(lasso.model$beta@Dimnames[1]), beta = paste(lasso.model$beta)) %>% 
        mutate(x = as.numeric(beta)) %>% 
        arrange(x) %>% 
        ggplot(aes(x = reorder(symbol, - x), y = x, fill = x)) +
        geom_bar(stat='identity') +
        labs(title = "Lasso logistic regression", x = "Gene symbol", y = "Beta value", fill = 'Beta value') + 
        scale_fill_gradient(low = "#D4E5F4",  high = "#FBB9BA") + 
        # scale_y_continuous(expand = c(0,0)) + 
        theme_classic() + 
        coord_flip() +
        guides(fill=guide_legend(reverse = T)) + 
        theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 10),
                  legend.position = 'right',
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 16, face = 'bold'))
    ggsave('~/Studio/01MetaRA/F007.jpg', width = 4.5, height = 4)
#### Difference analysis
#### GSE89408 analyssi
    require(GEOquery)
    g89408c = getGEO(filename = ("~/Studio/01MetaRA/GSE89408_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        pData(.)
    G89408 = data.table::fread('~/Studio/01MetaRA/GSE89408_GEO_count_matrix_rename.txt', header = T)[, c(1, 2:29, 68:219)] %>% 
    tibble::column_to_rownames('V1')
    G89408 = log2 (edgeR::cpm(G89408) +1)
    key = data.frame(symbol = c("COX7C", "COX7B", "PDK1", "PDK3", "COX11", "NDUFAF4", "NDUFB8", "NDUFS1", "SUCLG1", "PDHX", "SUCLA2")) 
    fit1 = as.data.frame(G89408) %>% 
        tibble::rownames_to_column('symbol') %>% 
        inner_join(key, by = 'symbol') %>% 
        # select(- c('beta')) %>% 
        tibble::column_to_rownames('symbol') %>% 
        t() %>% 
        as.data.frame() %>% 
        mutate_all(as.numeric) %>% 
        mutate(lasso = COX7C * 2.3558 + COX7B * 1.7983 + PDK1 * 2.7770 + PDK3 * 1.8816 + COX11 * -1.1158 + NDUFAF4 * -0.0719 + NDUFB8 * -0.7949 + NDUFS1 * -0.8306 + SUCLG1 * -0.2759 + PDHX * -1.3711 + SUCLA2 * -0.0634) %>% 
        mutate(group = c(rep('Healthy control', 28), rep('Early RA', 57), rep('Established RA', 95))) 
    fit1$group = factor(fit1$group, levels = c('Healthy control', 'Early RA', 'Established RA'))
    plot1 = fit1 %>% 
        select(c('lasso', 'group')) %>% 
        Rmisc::summarySE(., measurevar = c('lasso'), groupvars = c('group')) %>% 
        ggplot(aes(x = group, y = lasso, fill = group)) + 
            geom_bar(stat='identity') +
            labs(title = "GSE89408", x = "", y = "Lasso value", fill = "") + 
            geom_errorbar(aes(ymin = lasso - se, ymax = lasso + se), position = "dodge", width = 0.25) +
            scale_fill_manual(values = color) + 
            scale_y_continuous(expand = c(0,0), limits = c(0, 40)) +
            geom_signif(aes(x = group, y = lasso), data = fit1, map_signif_level = F, parse = T, comparisons = list(c('Early RA', 'Healthy control'), c('Established RA', 'Healthy control')), size = 0.5, textsize = 2.5, vjust = 0.1, tip_length = 0.01, y_position = c(32, 36, 32), test = "wilcox.test") +
            theme_classic() +
            theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 10),
                  legend.position = 'none',
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 16, face = 'bold'))
    fit1 = as.data.frame(G89408) %>% 
        tibble::rownames_to_column('symbol') %>% 
        inner_join(key, by = 'symbol') %>% 
        # select(- c('beta')) %>% 
        tibble::column_to_rownames('symbol') %>% 
        t() %>% 
        as.data.frame() %>% 
        mutate_all(as.numeric) %>% 
        mutate(lasso = COX7C * 2.3558 + COX7B * 1.7983 + PDK1 * 2.7770 + PDK3 * 1.8816 + COX11 * -1.1158 + NDUFAF4 * -0.0719 + NDUFB8 * -0.7949 + NDUFS1 * -0.8306 + SUCLG1 * -0.2759 + PDHX * -1.3711 + SUCLA2 * -0.0634) %>% 
        mutate(group = c(rep('Healthy control', 28), rep('Rheumatoid arthritis', 152))) 
    fit1$group = factor(fit1$group, levels = c('Healthy control', 'Rheumatoid arthritis'))
    plot1a = fit1 %>% 
        select(c('lasso', 'group')) %>% 
        Rmisc::summarySE(., measurevar = c('lasso'), groupvars = c('group')) %>% 
        ggplot(aes(x = group, y = lasso, fill = group)) + 
            geom_bar(stat='identity') +
            labs(title = "GSE89408", x = "", y = "Lasso value", fill = "") + 
            geom_errorbar(aes(ymin = lasso - se, ymax = lasso + se), position = "dodge", width = 0.25) +
            scale_fill_manual(values = color) + 
            scale_y_continuous(expand = c(0,0), limits = c(0, 38)) +
            geom_signif(annotations = c('< 2.2e-16'), map_signif_level = F, parse = F, comparisons = list(c('Rheumatoid arthritis', 'Healthy control')), size = 0.5, textsize = 2.5, vjust = 0.1, tip_length = 0.02, test = "wilcox.test", y_position = c(34)) +
            theme_classic() +
            theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 10),
                  legend.position = 'none',
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 16, face = 'bold'))
#### GSE77298 analysis
    require(GEOquery)
    G77298c = getGEO(filename = ("~/Studio/01MetaRA/GSE77298_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        pData(.) 
    G77298e = getGEO(filename = ("~/Studio/01MetaRA/GSE77298_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        exprs(.)
    if (length(getGEO(filename = ("~/Studio/01MetaRA/GSE77298_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F)) > 1) {
        idx = grep("GPL570", attr(getGEO(filename = ("~/Studio/01MetaRA/GSE77298_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F), "names"))} else {
        idx = 1}
    qx = as.numeric(quantile(G77298e, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    LogC = (qx[5] > 100) ||
        (qx[6] - qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) {
        G77298e[which(G77298e <= 0)] = NaN
        G77298e = log2(G77298e + 1)
        print("log2 transform finished")} else {
        print("log2 transform not needed")}
    fit2 = G77298e %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('probe_id') %>% 
        inner_join(AnnoProbe::idmap('GPL570'), by = 'probe_id') %>% 
        select('symbol', everything()) %>% 
        mutate(rowMean = rowMeans(.[grep('GSM', names(.))])) %>%
        dplyr::arrange(desc(rowMean)) %>%
        dplyr::rename('symbol' = 'symbol') %>% 
        dplyr::distinct(symbol, .keep_all = T) %>%
        dplyr::select(-c('rowMean', 'probe_id')) %>% 
        tibble::column_to_rownames('symbol') %>% 
        t() %>% 
        as.data.frame() %>% 
        select(key$symbol)  %>% 
        mutate(lasso = COX7C * 2.3558 + COX7B * 1.7983 + PDK1 * 2.7770 + PDK3 * 1.8816 + COX11 * -1.1158 + NDUFAF4 * -0.0719 + NDUFB8 * -0.7949 + NDUFS1 * -0.8306 + SUCLG1 * -0.2759 + PDHX * -1.3711 + SUCLA2 * -0.0634) %>% 
        tibble::rownames_to_column('geo_accession') %>% 
        inner_join(G77298c[, c('geo_accession', 'disease state:ch1')], by = 'geo_accession') %>% 
        dplyr::rename('group' = 'disease state:ch1')
    plot2 = fit2 %>% 
        select(c('lasso', 'group')) %>% 
        Rmisc::summarySE(., measurevar = c('lasso'), groupvars = c('group')) %>% 
        ggplot(aes(x = group, y = lasso, fill = group)) + 
            geom_bar(stat='identity') +
            labs(title = "GSE77298", x = "", y = "Lasso value", fill = "") + 
            geom_errorbar(aes(ymin = lasso - se, ymax = lasso + se), position = "dodge", width = 0.25) +
            scale_fill_manual(values = color) + 
            scale_y_continuous(expand = c(0,0), limits = c(0, 60)) +
            geom_signif(aes(x = group, y = lasso), data = fit2, map_signif_level = F, parse = T, comparisons = list(c('Healthy control;synovium', 'Rheumatoid arthritis;synovium')), size = 0.5, textsize = 2.5, vjust = 0.1, tip_length = 0.03, test = "wilcox.test") +
            scale_x_discrete(labels = c("Healthy control", "Rheumatoid arthritis")) +
            theme_classic() +
            theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 10),
                  legend.position = 'none',
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 16, face = 'bold'))
#### GSE93272 analysis
    e93272 = read.table('~/Studio/01MetaRA/GSE93272_series.txt', header = T) %>% 
        tibble::column_to_rownames('symbol') %>% 
        select(drugnaive$geo_accession) %>% 
        t() %>% 
        as.data.frame() %>% 
        select(key$symbol) %>% 
        mutate(lasso = COX7C * 2.3558 + COX7B * 1.7983 + PDK1 * 2.7770 + PDK3 * 1.8816 + COX11 * -1.1158 + NDUFAF4 * -0.0719 + NDUFB8 * -0.7949 + NDUFS1 * -0.8306 + SUCLG1 * -0.2759 + PDHX * -1.3711 + SUCLA2 * -0.0634) %>% 
        tibble::rownames_to_column('geo_accession') %>% 
        inner_join(drugnaive, by = 'geo_accession') %>% 
        mutate(group = ifelse(DISEASE == 'RA', 'Rheumatoid arthritis', 'Healthy control'))
    plot3 = e93272 %>% 
        Rmisc::summarySE(., measurevar = c('lasso'), groupvars = c('group')) %>% 
        ggplot(aes(x = group, y = lasso, fill = group)) + 
        geom_bar(stat='identity') +
        labs(title = "GSE93272", x = "", y = "Lasso value", fill = "") + 
        geom_errorbar(aes(ymin = lasso - se, ymax = lasso + se), position = "dodge", width = 0.25) +
        scale_fill_manual(values = color) + 
        scale_y_continuous(expand = c(0,0), limits = c(0, 50)) +
        geom_signif(aes(x = group, y = lasso), data = e93272, map_signif_level = F, parse = T, comparisons = list(c('Rheumatoid arthritis', 'Healthy control')), size = 0.5, textsize = 2.5, vjust = 0.1, tip_length = 0.025, y_position = c(45), test = "wilcox.test") +
        theme_classic() + 
        #coord_flip() +
        guides(fill = guide_legend(reverse = T)) + 
        theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 10),
                  legend.position = 'none',
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 16, face = 'bold'))
#### GSE45291 analysis
    require(GEOquery)
    G45291c = getGEO(filename = ("~/Studio/01MetaRA/GSE45291_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        pData(.) 
    G45291e = getGEO(filename = ("~/Studio/01MetaRA/GSE45291_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        exprs(.)
    if (length(getGEO(filename = ("~/Studio/01MetaRA/GSE45291_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F)) > 1) {
        idx = grep("GPL13158", attr(getGEO(filename = ("~/Studio/01MetaRA/GSE45291_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F), "names"))} else {
        idx = 1}
    qx = as.numeric(quantile(G45291e, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    LogC = (qx[5] > 100) ||
        (qx[6] - qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) {
        G45291e[which(G45291e <= 0)] = NaN
        G45291e = log2(G45291e + 1)
        print("log2 transform finished")} else {
        print("log2 transform not needed")}
    GPL13158 = data.table::fread('~/Studio/01MetaRA/GPL13158-5065.txt', header = T) %>% 
        select(c('Gene Symbol', 'ID'))
    fit3 = G45291e %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('ID') %>% 
        inner_join(GPL13158, by = 'ID') %>% 
        select('Gene Symbol', everything()) %>% 
        mutate(rowMean = rowMeans(.[grep('GSM', names(.))])) %>%
        dplyr::arrange(desc(rowMean)) %>%
        dplyr::rename('symbol' = 'Gene Symbol') %>% 
        dplyr::distinct(symbol, .keep_all = T) %>%
        dplyr::select(-c('rowMean', 'ID')) %>% 
        tibble::column_to_rownames('symbol')
    RA45291 = G45291c %>% 
        dplyr::rename("dis" = "disease:ch1") %>% 
        filter(dis != "SLE (Systemic LUPUS Erythomatosus)") %>% 
        select(c('geo_accession', 'dis'))
    fit3a = fit3 %>% 
        select(RA45291$geo_accession) %>% 
        t() %>% 
        as.data.frame() %>% 
        select(key$symbol) %>% 
        mutate(lasso = COX7C * 2.3558 + COX7B * 1.7983 + PDK1 * 2.7770 + PDK3 * 1.8816 + COX11 * -1.1158 + NDUFAF4 * -0.0719 + NDUFB8 * -0.7949 + NDUFS1 * -0.8306 + SUCLG1 * -0.2759 + PDHX * -1.3711 + SUCLA2 * -0.0634) %>% 
        tibble::rownames_to_column('geo_accession') %>% 
        inner_join(RA45291, by = 'geo_accession') %>% 
        mutate(group = ifelse(dis == 'Control', 'Healthy control', 'Rheumatoid arthritis'))
    plot4 = fit3a %>% 
        Rmisc::summarySE(., measurevar = c('lasso'), groupvars = c('group')) %>% 
        ggplot(aes(x = group, y = lasso, fill = group)) + 
        geom_bar(stat='identity') +
        labs(title = "GSE45291", x = "", y = "Lasso value", fill = "") + 
        geom_errorbar(aes(ymin = lasso - se, ymax = lasso + se), position = "dodge", width = 0.25) +
        scale_fill_manual(values = color) + 
        scale_y_continuous(expand = c(0,0), limits = c(0, 70)) +
        geom_signif(comparisons = list(c('Rheumatoid arthritis', 'Healthy control')), annotations = c(3.358e-05), y_position = c(65), size = 0.5, textsize = 2.5, vjust = 0.1, tip_length = 0.4) +
        theme_classic() + 
        #coord_flip() +
        guides(fill = guide_legend(reverse = T)) + 
        theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 10),
                  legend.position = 'none',
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 16, face = 'bold'))
#### combined plot
    require(patchwork)
    plot1a + plot2 + plot4 + plot_layout(ncol = 3)
    ggsave('~/Studio/01MetaRA/F006a.jpg', width = 7, height = 4.5)
    plot1 + plot3 
    ggsave('~/Studio/01MetaRA/F006b.jpg', width = 4, height = 4.5)
#### ROC
    require(plotROC) 
    require(pROC)
    fit89408 = fit1 %>% 
        select(c('lasso', 'group')) %>% 
        mutate(gse = 'GSE89408')
    fit77298 = fit2 %>% 
        select(c('lasso', 'group')) %>% 
        mutate(group = ifelse(group == 'Healthy control;synovium', 'Healthy control', 'Rheumatoid arthritis')) %>% 
        mutate(gse = 'GSE77298')
    fit93272 = e93272 %>% 
        select(c('lasso', 'group')) %>% 
        mutate(gse = 'GSE93272')
    fit45291 = fit3a %>% 
        select(c('lasso', 'group')) %>% 
        mutate(gse = 'GSE45291')
    fit0 = rbind(fit89408, fit77298, fit93272, fit45291) %>% 
        mutate(group = factor(group, levels = c('Rheumatoid arthritis', 'Healthy control')))
    x0 = data.frame()
    for(i in unique(fit0$gse)){
        p1 = fit0 %>% filter(gse == i)
        x1 = roc(response = p1$group, predictor = p1$lasso, data = p1)
        x2 = round(auc(response = p1$group, predictor = p1$lasso), 4)
        x3 = data.frame(symbol = i, value = x2)
        x0 = rbind(x3, x0)}
    write.csv(x0, "~/Studio/01MetaRA/acc.csv")
    y = fit0 %>% 
        group_split(gse)
    str(y)
    map_function <- function(y){
        pROC::roc(
        data = y,
        response = group,
        predictor = lasso)}
    roc.list = map(y, map_function)
    names(roc.list) = c("GSE45291", "GSE77298", "GSE89408", "GSE93272")
    plot6 = pROC::ggroc(roc.list, alpha = 1, size = 0.75, legacy = T) +
        geom_abline(linetype = 3, intercept = 0) +
        labs(title = "GEO datasets", color = 'gse') +
        theme_classic() +
        scale_color_manual(values = color) +
        guides(color = guide_legend(title = '')) +
        theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 10),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 16, face = 'bold'))
    ggsave('~/Studio/01MetaRA/F006d.jpg', width = 4, height = 3.5)
    require(plotROC) 
    require(pROC)
    r89408 = G89408 %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('symbol') %>% 
        inner_join(key, by = 'symbol') %>% 
        tibble::column_to_rownames('symbol') %>% 
        # select(- c('beta')) %>% 
        t() %>% 
        as.data.frame() %>% 
        select(c(COX11, COX7B, COX7C, NDUFAF4, NDUFB8, NDUFS1, PDHX, PDK1, PDK3, SUCLA2, SUCLG1)) %>% 
        mutate(lasso = COX7C * 2.3558 + COX7B * 1.7983 + PDK1 * 2.7770 + PDK3 * 1.8816 + COX11 * -1.1158 + NDUFAF4 * -0.0719 + NDUFB8 * -0.7949 + NDUFS1 * -0.8306 + SUCLG1 * -0.2759 + PDHX * -1.3711 + SUCLA2 * -0.0634) %>% 
        mutate(group = c(rep('Healthy control', 28), rep('Early rheumatoid arthritis', 57), rep('Established rheumatoid arthritis', 95)))
    r89408_1 = r89408 %>% 
        select(c("lasso", "group")) %>% 
        filter(group == 'Healthy control'| group == 'Early rheumatoid arthritis') %>% 
        mutate(case = c('Early RA'))
    r89408_2 = r89408 %>% 
        select(c("lasso", "group")) %>% 
        filter(group == 'Healthy control'| group == 'Established rheumatoid arthritis') %>% 
        mutate(case = c('Established RA'))
    y89408 = rbind(r89408_1, r89408_2) %>% 
        group_split(case)
    map_function <- function(y89408){
        pROC::roc(
        data = y89408,
        response = group,
        predictor = lasso)}
    roc.list = map(y89408, map_function)
    names(roc.list) = c("Early RA", "Established RA")
    plot7 = pROC::ggroc(roc.list, alpha = 1, size = 0.75, legacy = T) +
        geom_abline(linetype = 3, intercept = 0) +
        labs(title = "GEO datasets", color = 'gse') +
        theme_classic() +
        scale_color_manual(values = color) +
        guides(color = guide_legend(title = '')) +
        theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 10),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 16, face = 'bold'))
    ggsave('~/Studio/01MetaRA/F006e.jpg', width = 4.5, height = 3.5)
#### corrplot
    require(GSVA)
    key = data.frame(symbol = sort(c("COX7C", "COX7B", "PDK1", "PDK3", "COX11", "NDUFAF4", "NDUFB8", "NDUFS1", "SUCLG1", "PDHX", "SUCLA2")))
    g89408c = getGEO(filename = ("~/Studio/01MetaRA/GSE89408_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        pData(.)
    G89408 = data.table::fread('~/Studio/01MetaRA/GSE89408_GEO_count_matrix_rename.txt', header = T)[, c(1, 2:29, 68:219)] %>% 
    tibble::column_to_rownames('V1')
    G89408 = log2 (edgeR::cpm(G89408) +1)
    z89408_1 = as.data.frame(G89408) %>% 
        tibble::rownames_to_column('symbol') %>% 
        inner_join(key, by = 'symbol') %>% 
        # select(- c('beta')) %>% 
        tibble::column_to_rownames('symbol') %>% 
        t() %>% 
        as.data.frame() %>% 
        mutate_all(as.numeric) %>% 
        mutate(lasso = COX7C * 2.3558 + COX7B * 1.7983 + PDK1 * 2.7770 + PDK3 * 1.8816 + COX11 * -1.1158 + NDUFAF4 * -0.0719 + NDUFB8 * -0.7949 + NDUFS1 * -0.8306 + SUCLG1 * -0.2759 + PDHX * -1.3711 + SUCLA2 * -0.0634)
    z89408_2 = gsva(gsvaParam(as.matrix(G89408), gs_path, maxDiff = TRUE)) %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        select('TCA cycle')
    z89408_3 = data.table::fread("~/Studio/01MetaRA/CIBERSORTx_GSE89048_Results.txt", header = T)[c(1:28, 67:218), ] %>% 
        select(c(1:23)) %>% 
        dplyr::rename('sample' = 'Mixture') %>% 
        select(celltypeselected) %>% 
        mutate(group = c(rep('Healthy control', 28), rep('Rheumatoid arthritis', 152)))
    z89408_0 = cbind(z89408_1, z89408_2, z89408_3) %>% 
        select(- c('group')) %>% 
        select('lasso', everything())
    cor_matrix = cor(z89408_0)
    p.mat = corrplot::cor.mtest(z89408_0, conf.level = 0.95)$p
    fit89408 = cbind(value = as.data.frame(cor_matrix[,1]), pvalue = as.data.frame(p.mat[,1]))
    colnames(fit89408) = c('value', 'pvalue')
    fit89408 = as.data.frame(fit89408[-1,]) %>% mutate(group = c('GSE89408')) %>% tibble::rownames_to_column('symbol')
    e93272_1 = read.table('~/Studio/01MetaRA/GSE93272_series.txt', header = T) %>% 
        tibble::column_to_rownames('symbol') %>% 
        select(drugnaive$geo_accession) %>% 
        t() %>% 
        as.data.frame() %>% 
        select(key$symbol) %>% 
        mutate(lasso = COX7C * 2.3558 + COX7B * 1.7983 + PDK1 * 2.7770 + PDK3 * 1.8816 + COX11 * -1.1158 + NDUFAF4 * -0.0719 + NDUFB8 * -0.7949 + NDUFS1 * -0.8306 + SUCLG1 * -0.2759 + PDHX * -1.3711 + SUCLA2 * -0.0634) %>% 
        tibble::rownames_to_column('geo_accession') %>% 
        inner_join(drugnaive, by = 'geo_accession') %>% 
        mutate(group = ifelse(DISEASE == 'RA', 'Rheumatoid arthritis', 'Healthy control')) %>% 
        select(c(1:13)) %>% 
        tibble::column_to_rownames('geo_accession') 
    e93272_2 = gsva(gsvaParam(as.matrix(read.table('~/Studio/01MetaRA/GSE93272_series.txt', header = T) %>% tibble::column_to_rownames('symbol')), gs_path, maxDiff = TRUE)) %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        select('TCA cycle')
    e93272_3 = data.table::fread("~/Studio/01MetaRA/CIBERSORTx_GSE93272_Results.txt", header = T) %>% 
        select(c(1:23)) %>% 
        tibble::column_to_rownames('Mixture') %>% 
        select(celltypeselected) %>% 
        tibble::rownames_to_column('geo_accession') %>% 
        inner_join(drugnaive[, c('geo_accession', 'DISEASE')], by = 'geo_accession')
    e93272_0 = cbind(e93272_1, e93272_2, e93272_3) %>% 
        select(c('DISEASE', 'lasso'), everything()) %>% 
        select(-c('geo_accession'))
    cor_matrix = cor(e93272_0[, -1])
    p.mat = corrplot::cor.mtest(e93272_0[, -1], conf.level = 0.95)$p
    fit93272 = cbind(value = as.data.frame(cor_matrix[,1]), pvalue = as.data.frame(p.mat[,1]))
    colnames(fit93272) = c('value', 'pvalue')
    fit93272 = as.data.frame(fit93272[-1,]) %>% mutate(group = c('GSE93272')) %>% tibble::rownames_to_column('symbol')
    fit1 = rbind(fit89408, fit93272)
    fit1$symbol = factor(fit1$symbol, levels = rev(c("COX11", "COX7B", "COX7C", "NDUFAF4", "NDUFB8", "NDUFS1",  "PDHX", "PDK1", "PDK3", "SUCLA2", "SUCLG1", "TCA cycle", "Neutrophils", "Monocytes", "NK cells resting", "T cells regulatory (Tregs)", "T cells CD4 memory activated")))
    plot1 = fit1 %>% 
        mutate(sig = ifelse(pvalue < 0.001, '***', ifelse(pvalue < 0.01, '**', ifelse(pvalue < 0.05, '*', 'ns.')))) %>% 
        mutate(text = paste0(round(value, 2), '(', sig, ')'))  %>% 
        ggplot(aes(x = group, y = symbol, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_gradient2(low = "#D4E5F4", mid = "#FFFEEE", high = "#FBB9BA") +
        geom_text(aes(label = text), size = 3) + 
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        # scale_x_discrete(position = "bottom") +
        theme_minimal() +
        labs(title = "") +
        guides(fill = guide_legend(title = 'Correlation\ncoefficient')) + 
        theme(legend.position = "right",
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 10),
              panel.grid = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 0),
              axis.ticks = element_blank(),
              plot.title = element_text(face = 'bold', size = 20),
              axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
              axis.title = element_text(face = 'bold', size = 0))
    plot0 = fit1 %>% 
        mutate(y = 'group') %>% 
        ggplot(aes(x = y, y = symbol, fill = symbol)) +
        geom_tile(aes(fill = symbol)) +
        scale_fill_manual(values = rev(c("#C2E3D0", "#C2E3D0", "#C2E3D0", "#C2E3D0", "#C2E3D0", "#C2E3D0", "#C2E3D0", "#C2E3D0", "#C2E3D0", "#C2E3D0", "#C2E3D0", "#F6E4D0", "#A8CBDF", "#A8CBDF", "#A8CBDF", "#A8CBDF", "#A8CBDF", "#A8CBDF"))) +
        scale_y_discrete(position = "left") +
        theme_minimal() +
        xlab(NULL) + ylab(NULL) +
        labs(fill = 'module', title = '') + 
        guides(fill = 'none') + 
        theme(panel.border = element_blank(),
              panel.grid = element_blank(),
              legend.text = element_text(size = 14),
              axis.ticks = element_blank(),
              axis.text.y = element_text(size = 12),
              axis.text.x = element_blank(),
              plot.title = element_text(face = 'bold', size = 20),
              legend.position = "right")
    require(aplot)
    plot0 %>% insert_right(plot1, width = 6)  
    ggsave('~/Studio/01MetaRA/F007a.jpg', width = 5, height = 8.5)
#### Analysis plot
    key = data.frame(symbol = c("COX7C", "COX7B", "PDK1", "PDK3", "COX11", "NDUFAF4", "NDUFS1", "NDUFB8", "SUCLG1", "PDHX", "SUCLA2"))
    de93272 = read.table('~/Studio/01MetaRA/GSE93272_series.txt', header = T) %>% 
        tibble::column_to_rownames('symbol') %>% 
        select(drugnaive$geo_accession) %>% 
        t() %>% 
        as.data.frame() %>% 
        select(key$symbol) %>% 
        mutate(LASSO_model = COX7C * 2.3558 + COX7B * 1.7983 + PDK1 * 2.7770 + PDK3 * 1.8816 + COX11 * -1.1158 + NDUFAF4 * -0.0719 + NDUFB8 * -0.7949 + NDUFS1 * -0.8306 + SUCLG1 * -0.2759 + PDHX * -1.3711 + SUCLA2 * -0.0634) %>% 
        # select(-c(1:11)) %>% 
        tibble::rownames_to_column('geo_accession') %>% 
        inner_join(drugnaive, by = 'geo_accession') %>% 
        mutate(group = ifelse(DISEASE == 'RA', 'Rheumatoid arthritis', 'Healthy control')) %>% 
        tibble::column_to_rownames('geo_accession')  %>% 
        filter(DISEASE != 'HC') %>% 
        # select(-c(2:6)) %>% 
        select(-c('title', 'group')) %>% 
        as.data.frame() %>% 
        select(-c(13:17, 27, 30, 36, 37))
    cor_matrix = cor(de93272)
    p.mat = corrplot::cor.mtest(de93272, conf.level = 0.95)$p
    cort93272 = data.frame(cor_matrix[1:12, 13:30]) %>% 
        tibble::rownames_to_column('symbol') %>% 
        pivot_longer(!symbol, names_to = 'module', values_to = 'value') %>% 
        mutate(module = str_remove(module, "^value\\."))
    pv93272 = data.frame(p.mat[1:12, 13:30]) %>% 
        tibble::rownames_to_column('symbol') %>% 
        pivot_longer(!symbol, names_to = 'module', values_to = 'value') %>% 
        mutate(module = str_remove(module, "^value\\."))
    cort93272$symbol = factor(cort93272$symbol, levels = rev(c("COX11", "COX7B", "COX7C", "NDUFAF4", "NDUFB8", "NDUFS1",  "PDHX", "PDK1", "PDK3", "SUCLA2", "SUCLG1", "LASSO_model")))
    plot0 = cort93272 %>% 
        mutate(pvalue = pv93272$value) %>% 
        mutate(sig = ifelse(pvalue < 0.001, '***', ifelse(pvalue < 0.01, '**', ifelse(pvalue < 0.05, '*', 'ns.')))) %>% 
        mutate(text = paste0(round(value, 2), '(', sig, ')')) %>% 
        # mutate(y = 'Lasso value') %>% 
        # slice(-1) %>% 
        ggplot(aes(x = module, y = symbol, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_gradient2(low = "#D4E5F4", mid = "#FFFEEE", high = "#FBB9BA") +
        geom_text(aes(label = text), size = 2.5) + 
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        # scale_x_discrete(position = "bottom") +
        theme_minimal() +
        labs(title = "GSE93272") +
        guides(fill = guide_legend(title = 'Correlation\ncoefficient')) + 
        theme(legend.position = "right",
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 10),
              panel.grid = element_blank(),
              axis.text.y = element_text(size = 12),
              axis.ticks = element_blank(),
              plot.title = element_text(face = 'bold', size = 16),
              axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
              axis.title = element_text(face = 'bold', size = 0))
    ggsave('~/Studio/01MetaRA/F007b.jpg', width = 11.5, height = 6)
#### Drug response
    require(GEOquery)
    g93272c = getGEO(filename = ("~/Studio/01MetaRA/GSE93272_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        pData(.) %>% 
        mutate(patients = substr(.$characteristics_ch1, start = 16, nchar(.$characteristics_ch1))) 
    g93272c$title = gsub('Whole blood from rheumatoid arthritis', '', g93272c$title)
    g93272c$title = gsub('Whole blood from healthy control', '', g93272c$title)
    f93272 = getGEO(filename = ("~/Studio/01MetaRA/GSE93272_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        exprs(.) %>% 
        as.data.frame() %>% 
        select(- drugnaive$geo_accession) %>% 
        tibble::rownames_to_column('probe_id') %>% 
        inner_join(AnnoProbe::idmap('GPL570'), by = 'probe_id') %>% 
        select('symbol', everything()) %>% 
        mutate(rowMean = rowMeans(.[grep('GSM', names(.))])) %>%
        dplyr::arrange(desc(rowMean)) %>%
        dplyr::rename('symbol' = 'symbol') %>% 
        dplyr::distinct(symbol, .keep_all = T) %>%
        dplyr::select(-c('rowMean', 'probe_id')) %>% 
        tibble::column_to_rownames('symbol') %>% 
        t() %>% 
        as.data.frame() %>% 
        select(key$symbol) %>% 
        tibble::rownames_to_column('geo_accession') %>% 
        inner_join(g93272c[, c('title', 'geo_accession')], by = 'geo_accession') %>% 
        mutate(sample_id = substr(.$title, 2, 6)) %>% 
        select(-c('title')) %>% 
        mutate(lasso = COX7C * 2.7741 + COX7B * 1.1632 + PDK1 * 2.8873 + PDK3 * 1.9084 + COX11 * -1.5769 + NDUFAF4 * -0.0832 + SUCLG1 * -1.2878 + PDHX * -1.2707 + SUCLA2 * -0.1232)
    drugresponse = read.table("~/Studio/01MetaRA/drug.response_metadata.txt", header = T) %>% 
        dplyr::rename('patients' = 'INDIVIDUAL.ID') %>% 
        mutate(title = paste0('(', SAMPLE.ID, ')')) %>% 
        mutate(sample_id = substr(.$SAMPLE.ID, 1, 5))
    dy93272 = f93272 %>% 
        inner_join(drugresponse, by = 'sample_id') %>% 
        filter(DISEASE != 'HC')
    (dy93272[, c(2:12, 14)])
    cor_matrix = cor((dy93272[, c(2:12, 14)]))
    p.mat = corrplot::cor.mtest((dy93272[, c(2:12, 14)]), conf.level = 0.95)$p
    p93272 = cbind(value = as.data.frame(cor_matrix[,12]), pvalue = as.data.frame(p.mat[,12]))
    colnames(p93272) = c('value', 'pvalue')
    plot01 = as.data.frame(p93272[-12,]) %>% tibble::rownames_to_column('symbol') %>% 
        mutate(group = c('Lasso value')) %>% 
        mutate(sig = ifelse(pvalue < 0.001, '***', ifelse(pvalue < 0.01, '**', ifelse(pvalue < 0.05, '*', 'ns.')))) %>% 
        mutate(text = paste0(round(value, 2), '(', sig, ')'))  %>% 
        ggplot(aes(x = symbol, y = group, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_gradient2(low = "#D4E5F4", mid = "#FFFEEE", high = "#FBB9BA") +
        geom_text(aes(label = text), size = 3) + 
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        # scale_x_discrete(position = "bottom") +
        theme_minimal() +
        labs(title = "GSE93272") +
        guides(fill = guide_legend(title = 'Correlation\ncoefficient')) + 
        theme(legend.position = "right",
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 10),
              panel.grid = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 10),
              axis.ticks = element_blank(),
              plot.title = element_text(face = 'bold', size = 20),
              axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
              axis.title = element_text(face = 'bold', size = 0))
    ggsave('~/Studio/01MetaRA/F007c.jpg', width = 11.5, height = 2)
    dz93272 = (dy93272[, c(14, 22:43)])
    da93272 = (dz93272[, -c(20:22)] %>% mutate_all(as.numeric))
    x0 = data.frame()
    for(i in 2:20){
        y = cor.test(da93272[, 1], da93272[, i])
        y0 = data.frame(pvalue = y$p.value, value = y$estimate, text = names(da93272)[i])
        x0 = rbind(y0, x0)
    }
    cor.test((dy93272[, c('lasso', 'RF')] %>% na.omit() %>% mutate_all(as.numeric))$lasso, (dy93272[, c('lasso', 'RF')] %>% na.omit() %>% mutate_all(as.numeric))$RF, data = (dy93272[, c('lasso', 'RF')] %>% na.omit() %>% mutate_all(as.numeric)))
    cor.test((dy93272[, c('lasso', 'ANA')] %>% na.omit() %>% mutate_all(as.numeric))$lasso, (dy93272[, c('lasso', 'ANA')] %>% na.omit() %>% mutate_all(as.numeric))$ANA, data = (dy93272[, c('lasso', 'ANA')] %>% na.omit() %>% mutate_all(as.numeric)))
    cor.test((dy93272[, c('lasso', 'CCP')] %>% na.omit() %>% mutate_all(as.numeric))$lasso, (dy93272[, c('lasso', 'CCP')] %>% na.omit() %>% mutate_all(as.numeric))$CCP, data = (dy93272[, c('lasso', 'CCP')] %>% na.omit() %>% mutate_all(as.numeric)))
    di93272 = x0 %>% 
        rbind(., data.frame(text = 'RF', value = 0.2217174, pvalue = 0.2306)) %>% 
        rbind(., data.frame(text = 'ANA', value = -0.0598058, pvalue = 0.7493)) %>% 
        rbind(., data.frame(text = 'CCP', value = -0.1892562, pvalue = 0.3079)) %>% 
        dplyr::rename('symbol' = 'text') %>% 
        mutate(group = c('Lasso value')) %>% 
        mutate(sig = ifelse(pvalue < 0.001, '***', ifelse(pvalue < 0.01, '**', ifelse(pvalue < 0.05, '*', 'ns.')))) %>% 
        mutate(text = paste0(round(value, 2), '(', sig, ')'))  %>% 
        ggplot(aes(x = symbol, y = group, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_gradient2(low = "#D4E5F4", mid = "#FFFEEE", high = "#FBB9BA") +
        geom_text(aes(label = text), size = 2.5) + 
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        # scale_x_discrete(position = "bottom") +
        theme_minimal() +
        labs(title = "GSE93272") +
        guides(fill = guide_legend(title = 'Correlation\ncoefficient')) + 
        theme(legend.position = "right",
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 10),
              panel.grid = element_blank(),
              axis.text.y = element_text(size = 12),
              axis.ticks = element_blank(),
              plot.title = element_text(face = 'bold', size = 20),
              axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
              axis.title = element_text(face = 'bold', size = 0))
    ggsave('~/Studio/01MetaRA/F007d.jpg', width = 11.5, height = 3)
    naive_lasso = e93272 %>% filter(DISEASE != 'HC') %>% select('lasso')
    response_lasso = dy93272 %>% filter(DISEASE != 'HC') %>% select('lasso')
    t.test(naive_lasso$lasso, response_lasso$lasso)
    head(drugresponse)
    cor.test(dy93272$DRUG.TREATMENT.DAYS, dy93272$lasso)
    plot02 = dy93272 %>% 
        ggplot(aes(x = DRUG.TREATMENT.DAYS, y = lasso)) + 
        geom_point(color = '#9E9E9E', size = 0.5) + 
        geom_smooth(method="lm", se=T, color="#FBB9BA") + 
        theme_classic2() + 
        ggplot2::annotate("text", label = paste0("R=", round(-0.1272263, 2), ",P-value=", round(0.08106, 3)), size = 3, x = 200, y = 50 ) + 
        labs(title = "Correction analysis", , x = 'Drug treatment days', y = 'Lasso value') +
        theme(axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5),
              axis.text.y = element_text(size = 8, angle = 0),
              legend.title = element_text(size = 12, face = 'bold'),
              legend.text = element_text(size = 8),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 18, face = 'bold'))
    ggsave('~/Studio/01MetaRA/F007e.jpg', width = 3.5, height = 3.6)
    t.test(lasso ~ group, data = e77298)
#### Corrplot
    require(GSVA)
    key = data.frame(symbol = c("COX7C", "COX7B", "PDK1", "PDK3", "COX11", "NDUFAF4", "NDUFS1", "NDUFB8", "SUCLG1", "PDHX", "SUCLA2"))
    g89408c = getGEO(filename = ("~/Studio/01MetaRA/GSE89408_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        pData(.)
    G89408 = data.table::fread('~/Studio/01MetaRA/GSE89408_GEO_count_matrix_rename.txt', header = T)[, c(1, 2:29, 68:219)] %>% 
    tibble::column_to_rownames('V1')
    G89408 = log2 (edgeR::cpm(G89408) +1)
    z89408_1 = as.data.frame(G89408) %>% 
        tibble::rownames_to_column('symbol') %>% 
        inner_join(key, by = 'symbol') %>% 
        # select(- c('beta')) %>% 
        tibble::column_to_rownames('symbol') %>% 
        t() %>% 
        as.data.frame() %>% 
        mutate_all(as.numeric) %>% 
        mutate(lasso = COX7C * 2.3558 + COX7B * 1.7983 + PDK1 * 2.7770 + PDK3 * 1.8816 + COX11 * -1.1158 + NDUFAF4 * -0.0719 + NDUFB8 * -0.7949 + NDUFS1 * -0.8306 + SUCLG1 * -0.2759 + PDHX * -1.3711 + SUCLA2 * -0.0634)
    z89408_2 = gsva(gsvaParam(as.matrix(G89408), gs_path, maxDiff = TRUE)) %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        select('TCA cycle')
    z89408_3 = data.table::fread("~/Studio/01MetaRA/CIBERSORTx_GSE89048_Results.txt", header = T)[c(1:28, 67:218), ] %>% 
        select(c(1:23)) %>% 
        dplyr::rename('sample' = 'Mixture') %>% 
        select(celltypeselected) %>% 
        mutate(group = c(rep('Healthy control', 28), rep('Rheumatoid arthritis', 152)))
    z89408_0 = cbind(z89408_1, z89408_2, z89408_3) %>% 
        select(- c('group')) %>% 
        select('lasso', everything())
    cor_matrix = cor(z89408_0)
    p.mat = corrplot::cor.mtest(z89408_0, conf.level = 0.95)$p
    fit89408 = cbind(value = as.data.frame(cor_matrix[,1]), pvalue = as.data.frame(p.mat[,1]))
    colnames(fit89408) = c('value', 'pvalue')
    fit89408 = as.data.frame(fit89408[-1,]) %>% mutate(group = c('GSE89408')) %>% tibble::rownames_to_column('symbol')
    e93272_1 = read.table('~/Studio/01MetaRA/GSE93272_series.txt', header = T) %>% 
        tibble::column_to_rownames('symbol') %>% 
        select(drugnaive$geo_accession) %>% 
        t() %>% 
        as.data.frame() %>% 
        select(key$symbol)  %>% 
        tibble::rownames_to_column('geo_accession') %>% 
        inner_join(drugnaive, by = 'geo_accession') %>% 
        mutate(group = ifelse(DISEASE == 'RA', 'Rheumatoid arthritis', 'Healthy control')) %>% 
        select(c(1:10)) %>% 
        tibble::column_to_rownames('geo_accession') 
    e93272_2 = gsva(gsvaParam(as.matrix(read.table('~/Studio/01MetaRA/GSE93272_series.txt', header = T) %>% tibble::column_to_rownames('symbol')), gs_path, maxDiff = TRUE)) %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        select('TCA cycle')
    e93272_3 = data.table::fread("~/Studio/01MetaRA/CIBERSORTx_GSE93272_Results.txt", header = T) %>% 
        select(c(1:23)) %>% 
        tibble::column_to_rownames('Mixture') %>% 
        select(celltypeselected) %>% 
        tibble::rownames_to_column('geo_accession') %>% 
        inner_join(drugnaive[, c('geo_accession', 'DISEASE')], by = 'geo_accession')
    e93272_0 = cbind(e93272_1, e93272_2, e93272_3) %>% 
        select(c('DISEASE'), everything()) %>% 
        select(-c('geo_accession'))
    cor_matrix = cor(e93272_0[, -1])
    p.mat = corrplot::cor.mtest(e93272_0[, -1], conf.level = 0.95)$p
    fit93272_1 = cor_matrix[1:9,10:16] %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('symbol') %>% 
        pivot_longer(!symbol, names_to = 'module', values_to = 'value') %>% 
        mutate(module = str_remove(module, "^value\\."))
    fit93272_2 = p.mat[1:9,10:16] %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('symbol') %>% 
        pivot_longer(!symbol, names_to = 'module', values_to = 'p') %>% 
        mutate(module = str_remove(module, "^value\\."))
    fit93272_0 = cbind(fit93272_1, pvalue = fit93272_2$p)
    fit93272_0$symbol = factor(fit93272_0$symbol, levels = rev(c("COX11", "COX7B", "COX7C", "NDUFAF4", "PDHX", "PDK1", "PDK3", "SUCLA2", "SUCLG1")))
    fit93272_0$module = factor(fit93272_0$module, levels = c("TCA cycle", "Neutrophils", "Monocytes", "NK cells resting", "T cells regulatory (Tregs)", "T cells CD4 memory activated", "T cells CD4 memory resting"))
    fit1 = rbind(fit89408_0, fit93272_0) %>% 
        mutate(gse = c(rep('GSE89408', 63), rep('GSE93272', 63))) %>% 
        mutate(name = paste0(gse, module))
    plot1 = fit1 %>% 
        mutate(sig = ifelse(pvalue < 0.001, '***', ifelse(pvalue < 0.01, '**', ifelse(pvalue < 0.05, '*', 'ns.')))) %>% 
        mutate(text = paste0(round(value, 2), '(', sig, ')'))  %>% 
        ggplot(aes(x = name, y = symbol, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_gradient2(low = "#D4E5F4", mid = "#FFFEEE", high = "#FBB9BA") +
        geom_text(aes(label = sig), size = 3) + 
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        # scale_x_discrete(position = "bottom") +
        theme_minimal() +
        labs(title = "") +
        guides(fill = guide_legend(title = 'Correlation\ncoefficient')) + 
        theme(legend.position = "right",
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 10),
              panel.grid = element_blank(),
              axis.ticks = element_blank(),
              plot.title = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = 12),
              axis.title = element_blank())
    plot0 = fit1 %>% 
        mutate(y = 'group') %>% 
        ggplot(aes(x = name, y = y, fill = gse)) +
        geom_tile(aes(fill = gse)) +
        scale_fill_manual(values = color) +
        scale_y_discrete(position = "left") +
        theme_minimal() +
        xlab(NULL) + ylab(NULL) +
        labs(fill = 'module', title = '') + 
        guides(fill = 'none') + 
        theme(panel.border = element_blank(),
              panel.grid = element_blank(),
              legend.text = element_text(size = 14),
              axis.ticks = element_blank(),
              axis.text.y = element_text(size = 12),
              axis.text.x = element_blank(),
              plot.title = element_text(face = 'bold', size = 20),
              legend.position = "right")
    plot2 = fit1 %>% 
        mutate(y = 'module') %>% 
        ggplot(aes(x = reorder(name, module), y = y, fill = module)) +
        geom_tile(aes(fill = module)) +
        scale_fill_manual(values = c("#AA8984", "#BFB1D0", "#85BF67", 
              "#CCCC99", "#DCD7C1", "#A8CBDF", "#F0988C")) +
        scale_y_discrete(position = "left") +
        scale_x_discrete(labels = rep(c("Monocytes", "Neutrophils", "NK cells resting", "T cells CD4 memory activated", "T cells CD4 memory resting", "T cells regulatory (Tregs)", "TCA cycle"), 2)) +
        theme_minimal() +
        xlab(NULL) + ylab(NULL) +
        labs(fill = 'module', title = '') + 
        guides(fill = 'none') + 
        theme(panel.border = element_blank(),
              panel.grid = element_blank(),
              legend.text = element_text(size = 14),
              axis.ticks = element_blank(),
              axis.text.y = element_text(size = 12),
              axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
              plot.title = element_blank(),
              legend.position = "right")    
    require(aplot)
    plot0 %>% insert_bottom(plot1, height = 15)  %>% insert_bottom(plot2, height = 1)  
    ggsave('~/Studio/01MetaRA/F007x.jpg', width = 4.5, height = 11)
#### nomogram
    set.seed(1234)
    require(rms)
    require(Hmisc)
    key = data.frame(symbol = c("COX7C", "COX7B", "PDK1", "PDK3", "COX11", "NDUFAF4", "NDUFB8", "NDUFS1", "SUCLG1", "PDHX", "SUCLA2")) 
    data = as.data.frame(G89408) %>% 
        tibble::rownames_to_column('symbol') %>% 
        inner_join(key, by = 'symbol') %>% 
        # select(- c('beta')) %>% 
        tibble::column_to_rownames('symbol') %>% 
        t() %>% 
        as.data.frame() %>% 
        mutate_all(as.numeric) %>% 
        mutate(lasso = COX7C * 2.3558 + COX7B * 1.7983 + PDK1 * 2.7770 + PDK3 * 1.8816 + COX11 * -1.1158 + NDUFAF4 * -0.0719 + NDUFB8 * -0.7949 + NDUFS1 * -0.8306 + SUCLG1 * -0.2759 + PDHX * -1.3711 + SUCLA2 * -0.0634) %>% 
        mutate(group = factor(c(rep('HC', 28), rep('RA', 152)), levels = c('HC', 'RA'))) 
    write.csv(d, "~/Studio/01MetaRA/1.csv")
    dd <- datadist(data) 
    options(datadist = "dd")    
    model_penalized <- lrm(group ~ ., data = data, penalty = 0.1, x = TRUE, y = TRUE)
    nom <- nomogram(model_penalized, 
                fun = plogis,        # 逻辑回归使用plogis函数
                fun.at = seq(0.1, 0.9, by = 0.1),  # 预测概率的范围
                funlabel = "Probability of RA",
                lp = FALSE,         # 不显示线性预测值
                abbrev = FALSE,     # 不缩写变量名
                minlength = 1)   
    plot(nom)
    p1 = rbind(data.frame(points = nom$COX11$points, symbol = c("COX11"), line = nom$COX11$COX11),
          data.frame(points = nom$COX7B$points, symbol = c("COX7B"), line = nom$COX7B$COX7B),
          data.frame(points = nom$COX7C$points, symbol = c("COX7C"), line = nom$COX7C$COX7C),
          data.frame(points = nom$NDUFAF4$points, symbol = c("NDUFAF4"), line = nom$NDUFAF4$NDUFAF4),
          data.frame(points = nom$NDUFB8$points, symbol = c("NDUFB8"), line = nom$NDUFB8$NDUFB8),
          data.frame(points = nom$NDUFS1$points, symbol = c("NDUFS1"), line = nom$NDUFS1$NDUFS1),
          data.frame(points = nom$PDHX$points, symbol = c("PDHX"), line = nom$PDHX$PDHX),
          data.frame(points = nom$PDK1$points, symbol = c("PDK1"), line = nom$PDK1$PDK1),
          data.frame(points = nom$PDK3$points, symbol = c("PDK3"), line = nom$PDK3$PDK3),
          data.frame(points = nom$SUCLA2$points, symbol = c("SUCLA2"), line = nom$SUCLA2$SUCLA2),
          data.frame(points = nom$SUCLG1$points, symbol = c("SUCLG1"), line = nom$SUCLG1$SUCLG1),
          data.frame(points = nom$lasso$points, symbol = c("Lasso_model"), line = nom$lasso$lasso),
          data.frame(points = c(seq(0, 100, by = 10)), symbol = c("Total_point"), line = c(seq(0, 100, by = 10)))) %>% 
        mutate(symbol = factor(symbol, levels = rev(c("COX11", "COX7B", "COX7C", "NDUFAF4", "NDUFB8", "NDUFS1", "PDHX", "PDK1", "PDK3", "SUCLA2", "SUCLG1", "Lasso_model", "Total_point"))))
    write.csv(p1, "~/Studio/01MetaRA/2.csv")
    key_points <- p1[1:140, ] %>% group_by(symbol) %>%
        filter(line == min(line) | line == max(line) | line == median(line)) %>% 
        rbind(data.frame(points = c(seq(0, 100, by = 10)), symbol = c("Total_point"), line = c(seq(0, 100, by = 10)))) %>% 
        as.data.frame() %>% 
        filter(symbol != 'Lasso_model') %>% 
        rbind(data.frame(points = c(seq(0, 100, by = 25)), symbol = c("Lasso_model"), line = c(seq(5, 45, by = 10)))) %>% 
        rbind(data.frame(points = c(10.474420), symbol = c("NDUFS1"), line = c(4.5))) %>% 
        rbind(data.frame(points = c(17.827083), symbol = c("NDUFB8"), line = c(3.0))) %>% 
        rbind(data.frame(points = c(31.704436), symbol = c("PDK1"), line = c(5.0))) 
    plot001 = p1 %>% 
        ggplot(aes(x = points, y = symbol, color = symbol)) +  # 固定y=1
        geom_point(size = 2, shape = '|') +
        geom_line() +
        geom_text(data = key_points, aes(label = round(line, 1)),  vjust = -1, size = 2.5) +  # 显示y值标签
        # facet_wrap(~ symbol, ncol = 1) +
        # scale_y_continuous(breaks = NULL) +
        scale_color_manual(values = c('#4E4E4E', '#91D1C2', '#4DBBD5', "#FFEBAD", "#A5AAF9", "#BFB1D0", "#85BF67", 
              "#EAB883", "#DCD7C1", "#A8CBDF", "#F0988C", "#849184", "#A5CCC7")) +
        labs(x = "Total Points", y = NULL) +
        theme_minimal() +
        # scale_x_continuous(breaks = seq(0, 100, by = 10), limits = c(0, max(p1$points) * 1.1), expand = expansion(mult = c(0.02, 0.05))) +
        theme(strip.text = element_blank(),
              #panel.spacing = unit(1, "lines"),
              legend.position = "none", 
              axis.ticks = element_blank(),  # 隐藏所有刻度线
              axis.line.x = element_blank(),
              axis.title.x = element_blank(),  # 移除y轴标签
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = 10),
              panel.grid.major.y = element_blank(),  # 移除水平网格线
              panel.grid.major.x = element_blank(), 
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank())
    ggsave('~/Studio/01MetaRA/F011.jpg', width = 7, height = 3.6)
## Figure5 scRNA seq analysis
#### scRNA preparation
    ## remove.packages(c("Seurat", "Seuratobject"))
    ## install.packages('Seurat', repos = c( 'https://satijalab.r-universe.dev'), lib = "/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SeuratV4")
    ## .libPaths(c("/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SeuratV4", .libPaths()))
    ## packageVersion('Seurat') # 4.4.0
    ## packageVersion("Seuratobject")
    .libPaths(c("/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/Library", .libPaths()))
    library(Seurat)
    scra1 = CreateSeuratObject(counts = Read10X('~/Studio/01MetaRA/SDY2213/scRNA//RA1')[["Gene Expression"]], project = 'RA1', min.features = 200, min.cells = 5)
    scra2 = CreateSeuratObject(counts = Read10X('~/Studio/01MetaRA/SDY2213/scRNA//RA2')[["Gene Expression"]], project = 'RA2', min.features = 200, min.cells = 5)
    scra3 = CreateSeuratObject(counts = Read10X('~/Studio/01MetaRA/SDY2213/scRNA//RA3')[["Gene Expression"]], project = 'RA3', min.features = 200, min.cells = 5)
    scra3a = CreateSeuratObject(counts = Read10X('~/Studio/01MetaRA/SDY2213/scRNA//RA3A'), project = 'RA3A', min.features = 200, min.cells = 5)
    scra3b = CreateSeuratObject(counts = Read10X('~/Studio/01MetaRA/SDY2213/scRNA//RA3B'), project = 'RA3B', min.features = 200, min.cells = 5)
    scra4 = CreateSeuratObject(counts = Read10X('~/Studio/01MetaRA/SDY2213/scRNA//RA4')[["Gene Expression"]], project = 'RA4', min.features = 200, min.cells = 5)
    scra4a = CreateSeuratObject(counts = Read10X('~/Studio/01MetaRA/SDY2213/scRNA//RA4A')[["Gene Expression"]], project = 'RA4A', min.features = 200, min.cells = 5)
    scra4b = CreateSeuratObject(counts = Read10X('~/Studio/01MetaRA/SDY2213/scRNA//RA4B')[["Gene Expression"]], project = 'RA4B', min.features = 200, min.cells = 5)
    scra4c = CreateSeuratObject(counts = Read10X('~/Studio/01MetaRA/SDY2213/scRNA//RA4C')[["Gene Expression"]], project = 'RA4C', min.features = 200, min.cells = 5)
    scra5 = CreateSeuratObject(counts = Read10X('~/Studio/01MetaRA/SDY2213/scRNA//RA5')[["Gene Expression"]], project = 'RA5', min.features = 200, min.cells = 5)
    scra5a = CreateSeuratObject(counts = Read10X('~/Studio/01MetaRA/SDY2213/scRNA//RA5A')[["Gene Expression"]], project = 'RA5A', min.features = 200, min.cells = 5)
    scra5b = CreateSeuratObject(counts = Read10X('~/Studio/01MetaRA/SDY2213/scRNA//RA5B')[["Gene Expression"]], project = 'RA5B', min.features = 200, min.cells = 5)
    scra5c = CreateSeuratObject(counts = Read10X('~/Studio/01MetaRA/SDY2213/scRNA//RA5C')[["Gene Expression"]], project = 'RA5C', min.features = 200, min.cells = 5)
    scra6a = CreateSeuratObject(counts = Read10X('~/Studio/01MetaRA/SDY2213/scRNA//RA6A'), project = 'RA6A', min.features = 200, min.cells = 5)
    scra6b = CreateSeuratObject(counts = Read10X('~/Studio/01MetaRA/SDY2213/scRNA//RA6B'), project = 'RA6B', min.features = 200, min.cells = 5)
    sce_all = merge(scra1, y = c(scra2, scra3, scra3a, scra3b, scra4, scra4a, scra4b, scra4c, scra5, scra5a, scra5b, scra5c, scra6a, scra6b), add.cell.ids = c("RA1", "RA2", "RA3", "RA3A", "RA3B", "RA4", "RA4A", "RA4B", "RA4C", "RA5", "RA5A", "RA5B", "RA5C", "RA6A", "RA6B"), project = "RA Synovial") 
    save(sce_all, '~/Studio/01MetaRA/scRNASDY2213.RData')
#### GSE109449 datasets
    .libPaths(c("/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SeuratV4", .libPaths()))
    packageVersion('Seurat') # 4.4.0
    packageVersion("SeuratObject")
    require(Seurat)
    require(tidyverse)
    meta <- data.table::fread("~/Studio/01MetaRA/GSE109449/GSE109449_singlecell_rnaseq_metadata.tsv.gz",
                   sep = "\t", header = T) %>% 
                   filter(disease == 'RA')
    data <- data.table::fread("~/Studio/01MetaRA/GSE109449/GSE109449_singlecell_rnaseq_gene_counts.tsv.gz",
                   sep = "\t", header = T) %>%
            tibble::column_to_rownames(var = "ID_REF") %>% 
            select(meta$sample_name)
    colnames(meta) <- make.names(colnames(meta))
    rownames(meta) <- meta$sample_name
    sc <- CreateSeuratObject(counts = data, 
                        min.cells = 3, 
                        min.features = 0, 
                        meta.data = meta)
    hg = read.table("~/Studio/01MetaRA/gene_to_ENSG.txt", header = T, sep = "\t")[, c('gene_id', 'gene_name')] %>% 
        distinct(gene_id, .keep_all = TRUE)
    fit02 = data.frame(gene_id = sc@assays$RNA@counts@Dimnames[[1]])
    fit01 = fit02 %>% 
        left_join(hg, by = 'gene_id') %>% 
        mutate(symbol = ifelse(is.na(.$gene_name), .$gene_id, .$gene_name)) %>% 
        mutate(gene = make.unique(symbol))
    str(sc)
    sc@assays$RNA@counts@Dimnames[[1]] = fit01$gene
    sc@assays$RNA@data@Dimnames[[1]] = fit01$gene
    G109449 = sc
    rm(sc)
    G109449 = NormalizeData(G109449, normalization.method = "LogNormalize", scale.factor = 10000)
    G109449 = FindVariableFeatures(G109449)
    G109449 = ScaleData(G109449, features = rownames(G109449))
    G109449 <- RunPCA(G109449, features = VariableFeatures(object = G109449))
    VizDimLoadings(G109449, dims = 1:2, reduction = "pca")
    G109449 <- FindNeighbors(G109449, reduction = "pca", dims = 1:30, verbose = F)
    G109449 <- FindClusters(G109449, verbose = F, resolution = 0.5)
    G109449 <- RunUMAP(G109449, dims = 1:30)
    G109449 <- RunTSNE(G109449, dims = 1:30)  
    DimPlot(G109449, reduction = "umap")
    require(SingleR)
    load("~/Studio/01MetaRA/HumanPrimaryCellAtlas_refMonaco_human.RData")
    load("~/Studio/01MetaRA/HumanPrimaryCellAtlas_hpca.se_human.RData")
    G109449_for_SingleR <- GetAssayData(G109449, slot = "data")
    G109449.ref <- SingleR(test = G109449_for_SingleR, ref = hpca.se, labels = hpca.se$label.main)
    identical(rownames(G109449.ref), colnames(G109449))
    unique(G109449.ref$labels)
    G109449@meta.data$celltype <- G109449.ref$labels    
    save(G109449, file = '~/Studio/01MetaRA/G109449.Rdata')
    load(file = '~/Studio/01MetaRA/G109449.Rdata')
    a1 = as.data.frame(G109449@active.ident) %>% 
        dplyr::rename('ident' = 'G109449@active.ident') %>% 
        mutate(ident = ifelse(ident == '0', 'Fibroblast1', ifelse(ident == '1', 'Fibroblast2', 'Fibroblast3')))
    G109449@meta.data$ident = a1$ident
    cell_type_med <- as.data.frame(G109449@reductions$umap@cell.embeddings) %>% 
        mutate(labels = G109449@meta.data$ident) %>% 
        group_by(labels) %>%
        summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
    label = data.frame(symbol = c("SUCLG1", "SUCLA2", "PDK3", "PDK1", "PDHX", "NDUFS1", "NDUFB8", "NDUFAF4", "COX7C", "COX7B", "COX11"))
    plot01 = as.data.frame(G109449@reductions$umap@cell.embeddings) %>% 
        mutate(labels = G109449@meta.data$ident) %>% 
        ggplot(aes(x = UMAP_1, y = UMAP_2, color = labels)) +
        geom_point(size = 1.5, alpha = 1, shape = 16) +
        theme_classic() + 
        scale_color_manual(values = color) +
        labs(title = "GSE109449") + 
        guides(color = guide_legend(title = '', override.aes = list(size = 2))) +
        ggrepel::geom_text_repel(aes(label = labels), fontface = "bold", size = 3, data = cell_type_med, point.padding=unit(0.5, "lines")) + 
        theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 1, vjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 8),
                  legend.position = 'none',
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 16, face = 'bold'))
    ggsave("~/Studio/01MetaRA/F013a.jpg", width = 4, height = 4)
    plot02 = DotPlot(G109449, features = sort(c("SUCLG1", "SUCLA2", "PDK3", "PDK1", "PDHX", "NDUFS1", "NDUFB8", "NDUFAF4", "COX7C", "COX7B", "COX11" )), assay='RNA', group.by = 'ident')
    exp_mat = plot02$data %>% 
        select(-pct.exp, -avg.exp) %>%  
        pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
        as.data.frame()
    percent_mat = plot02$data %>% 
        select(-avg.exp, -avg.exp.scaled) %>%  
        pivot_wider(names_from = id, values_from = pct.exp) %>% 
        as.data.frame() %>% 
        pivot_longer(!features.plot, names_to = 'cell', values_to = 'value')
    plot013 = exp_mat %>% 
        pivot_longer(!features.plot, names_to = 'cell', values_to = 'module') %>% 
        cbind(percent_mat[, 3]) %>% 
        mutate(symbol = factor(.$features.plot, levels = rev(unique(.$features.plot)))) %>% 
        ggplot(aes(x = cell, y = symbol)) +
        # geom_tile(aes(fill = module)) + 
        geom_point(aes(size = value, color = value)) +
        scale_color_gradient(low = "#D4E5F4", high = "#F39B7F") +
        # geom_text(aes(label = text), size = 2.5) + 
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        # scale_x_discrete(position = "bottom") +
        theme_minimal() +
        labs(title = "Cell type") +
        scale_size(guide = 'none') +
        # guides(color = guide_legend(title = 'Percent')) + 
        theme(legend.position = "right",
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 8),
              panel.grid = element_blank(),
              axis.text.y = element_text(size = 12),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title = element_text(face = 'bold', size = 0),
              plot.title = element_text(size = 16, face = 'bold'))
    plot014 = exp_mat %>% 
        pivot_longer(!features.plot, names_to = 'cell', values_to = 'module') %>% 
        cbind(percent_mat[, 3]) %>% 
        mutate(symbol = factor(.$features.plot, levels = rev(unique(.$features.plot)))) %>% 
        mutate(y = 'x') %>% 
        ggplot(aes(x = cell, y = y, fill = cell)) +
        geom_tile(aes(fill = cell)) +
        scale_fill_manual(values = color) +
        scale_x_discrete(position = "bottom") +
        theme_minimal() +
        xlab(NULL) + ylab(NULL) +
        labs(title = '') + 
        guides(fill = 'none') + 
        theme(panel.border = element_blank(),
              panel.grid = element_blank(),
              legend.text = element_blank(),
              axis.ticks = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
              plot.title = element_blank(),
              legend.position = "none")
    require(aplot)
    plot014 %>% insert_top(plot013, height = 20)  
    ggsave("~/Studio/01MetaRA/F013b.jpg", width = 3.3, height = 4.8)
    gs_th = readxl::read_excel('~/Studio/01MetaRA/MetabolicPathways.xlsx', sheet = 1)
    gs_path = split(gs_th$symbol, list(gs_th$Pathway))
    require(GSVA)
    gsva109449 = gsva(gsvaParam(as.matrix(G109449@assays$RNA@data), gs_path, maxDiff = TRUE))
    gsvax = as.data.frame(gsva109449)
    key = data.frame(symbol = factor(c("COX7C", "COX7B", "PDK1", "PDK3", "COX11", "NDUFAF4", "NDUFB8", "NDUFS1", "SUCLG1", "PDHX", "SUCLA2"), levels = c("COX7C", "COX7B", "PDK1", "PDK3", "COX11", "NDUFAF4", "NDUFB8", "NDUFS1", "SUCLG1", "PDHX", "SUCLA2")))
    str(G109449@assays$RNA@data)
    z1 <- as.data.frame(G109449@assays$RNA@data)[rownames(G109449@assays$RNA@data) %in% key$symbol, ]
    cbz = rbind(z1, gsvax) %>% t() %>% as.data.frame()
    cor_matrix = cor(cbz)
    p.mat = corrplot::cor.mtest(cbz, conf.level = 0.95)$p
    fitc = as.data.frame(cor_matrix[1:11, 12:18]) %>% 
        tibble::rownames_to_column('symbol') %>% 
        pivot_longer(!symbol, names_to = 'module', values_to = 'value')
    fitc1 = as.data.frame(p.mat[1:11, 12:18]) %>% 
        tibble::rownames_to_column('symbol') %>% 
        pivot_longer(!symbol, names_to = 'module', values_to = 'pvalue')
    fitc$symbol = factor(fitc$symbol, levels = rev(c("COX11", "COX7B", "COX7C", "NDUFAF4", "NDUFB8", "NDUFS1",  "PDHX", "PDK1", "PDK3", "SUCLA2", "SUCLG1")))
    fitc$module = factor(fitc$module, levels = rev(c("Amino acid", "Carbohydrate", "Energy",  "Lipid", "Nucleotide", "TCA cycle", "Vitamin cofactor")))
    plot015 = fitc %>% 
        mutate(pvalue = fitc1$pvalue) %>% 
        mutate(sig = ifelse(pvalue < 0.001, '***', ifelse(pvalue < 0.01, '**', ifelse(pvalue < 0.05, '*', 'ns.')))) %>% 
        mutate(text = paste0(round(value, 2), '(', sig, ')'))  %>% 
        ggplot(aes(x = module, y = symbol)) +
        # geom_point(aes(size = value, color = pvalue)) +
        geom_tile(aes(fill = value)) + 
        scale_fill_gradient2(low = "#D4E5F4", mid = "#FFEBAD44", high = "#FBB9BA") +
        geom_text(aes(label = sig), size = 2.5) + 
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        # scale_x_discrete(position = "bottom") +
        theme_minimal() +
        labs(title = "Correction heatmap") +
        guides(fill = guide_legend(title = 'Correlation\ncoefficient')) + 
        theme(legend.position = "right",
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 10),
              panel.grid = element_blank(),
              axis.text.y = element_text(size = 12),
              axis.ticks = element_blank(),
              plot.title = element_text(face = 'bold', size = 16),
              axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
              axis.title = element_text(face = 'bold', size = 0))
    ggsave("~/Studio/01MetaRA/F013c.jpg", width = 4, height = 4.8)
    identical(rownames(G109449@reductions$umap@cell.embeddings), colnames(gsvax))
    all_values <- unlist(as.data.frame(t(gsvax)) %>% select('Amino acid', 'Carbohydrate', 'Energy', 'Lipid', 'Nucleotide', 'TCA cycle', 'Vitamin cofactor'))
    common_limits <- range(all_values, na.rm = TRUE)
    create_plot <- function(data, title, colname) {
        data %>% 
        select(all_of(colname)) %>% 
        dplyr::rename('AC' = all_of(colname)) %>% 
        tibble::rownames_to_column('name') %>% 
        inner_join(as.data.frame(G109449@reductions$umap@cell.embeddings) %>% 
               tibble::rownames_to_column('name'), by = 'name') %>% 
        ggplot(aes(x = UMAP_1, y = UMAP_2)) +
        geom_point(aes(color = AC), size = 1) +
        theme_classic() + 
        scale_color_gradient(low = "#D4E5F4",  high = "#FBB9BA", limits = common_limits) +
        labs(title = title) + 
        guides(color = guide_legend(title = '', override.aes = list(size = 2))) +
        theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            legend.title = element_text(size = 12, face = 'bold'),
            legend.text = element_text(size = 10),
            legend.position = 'right',
            axis.title = element_text(size = 10),
            plot.title = element_text(size = 16, face = 'bold'))}
    fit011 <- create_plot(as.data.frame(t(gsvax)), "Amino acid", "Amino acid")
    fit012 <- create_plot(as.data.frame(t(gsvax)), "Carbohydrate", "Carbohydrate")
    fit013 <- create_plot(as.data.frame(t(gsvax)), "Energy", "Energy")
    fit014 <- create_plot(as.data.frame(t(gsvax)), "Lipid", "Lipid")
    fit015 <- create_plot(as.data.frame(t(gsvax)), "Nucleotide", "Nucleotide")
    fit016 <- create_plot(as.data.frame(t(gsvax)), "TCA cycle", "TCA cycle")
    fit017 <- create_plot(as.data.frame(t(gsvax)), "Vitamin cofactor", "Vitamin cofactor")
    require(patchwork)
    combined_plot <- fit011 + fit012 + fit013 + fit014 + fit015 + fit016 + fit017 + plot_layout(ncol = 4, guides = 'collect')
    ggsave('~/Studio/01MetaRA/F013d.jpg', width = 11.5, height = 5) 
#### Binvignat analysis
    .libPaths(c("/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SeuratV4", .libPaths()))
    packageVersion('Seurat') # 4.4.0
    packageVersion("SeuratObject")
    require(Seurat)
    ## load("/Users/linqi/Studio/01MetaRA/BinvignatscRNA.Rdata") ## V4a
    require(tidyverse)
    ## library(scKidney)
    counts <- GetAssayData(sce, slot = "counts") ## ve slot; v5 layer
    str(counts)
    Binvignat <- CreateSeuratObject(counts)
    hg = read.table("~/Studio/01MetaRA/gene_to_ENSG.txt", header = T, sep = "\t")[, c('gene_id', 'gene_name')] %>% 
        distinct(gene_id, .keep_all = TRUE)
    fit02 = data.frame(gene_id = Binvignat@assays$RNA@counts@Dimnames[[1]])
    fit01 = fit02 %>% 
        left_join(hg, by = 'gene_id') %>% 
        mutate(symbol = ifelse(is.na(.$gene_name), .$gene_id, .$gene_name)) %>% 
        mutate(gene = make.unique(symbol))
    str(Binvignat)
    Binvignat@assays$RNA@counts@Dimnames[[1]] = fit01$gene
    Binvignat@assays$RNA@data@Dimnames[[1]] = fit01$gene
    Binvignat[["percent.mt"]] <- PercentageFeatureSet(Binvignat, pattern = "^MT-")
    VlnPlot(Binvignat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, raster = F)
    plot1 <- FeatureScatter(Binvignat, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(Binvignat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    require(patchwork)
    plot1 + plot2 + plot_layout(guides = 'collect')
    ggsave('~/Studio/01MetaRA/F009.jpg', width = 11.5, height = 4.5)
    metadata <- sce@meta.data
    identical(rownames(metadata), colnames(Binvignat))
    Binvignat = NormalizeData(Binvignat, normalization.method = "LogNormalize", scale.factor = 10000)
    Binvignat = FindVariableFeatures(Binvignat)
    Binvignat = ScaleData(Binvignat, features = rownames(Binvignat))
    save(Binvignat, file = "~/Studio/01MetaRA/BinvignatX.Rdata")
    metadata <- sce@meta.data %>% write.table(., file = "~/Studio/01MetaRA/BinvignatXmeta.txt", sep = "\t")
    metadata = read.table("~/Studio/01MetaRA/BinvignatXmeta.txt", header = T, sep = "\t")
    metadata1 = Binvignat@meta.data %>% 
        cbind(metadata[, -c(1:3)])
    identical(rownames(metadata1), colnames(Binvignat))
    Binvignat@meta.data <- metadata1
    require(Seurat)
    Binvignat <- RunPCA(Binvignat, features = VariableFeatures(object = Binvignat))
    rm(list = ls())
    VizDimLoadings(Binvignat, dims = 1:2, reduction = "pca")
    Binvignat <- FindNeighbors(Binvignat, reduction = "pca", dims = 1:30, verbose = F)
    Binvignat <- FindClusters(Binvignat, verbose = F, resolution = 0.5)
    Binvignat <- RunUMAP(Binvignat, dims = 1:30)
    Binvignat <- RunTSNE(Binvignat, dims = 1:30)  
    DimPlot(Binvignat, reduction = "umap")
    all_markers <- FindAllMarkers(Binvignat, only.pos = TRUE, min.pct= 0.1, logfc.threshold= 0.25)
    head(all_markers)
    write.table(all_markers, "~/Studio/01MetaRA/2.txt", sep = '\t')
    require(SingleR)
    load("~/Studio/01MetaRA/HumanPrimaryCellAtlas_refMonaco_human.RData")
    Binvignat_for_SingleR <- GetAssayData(Binvignat, slot = "data")
    Binvignat.ref <- SingleR(test = Binvignat_for_SingleR, ref = refMonaco, labels =refMonaco$label.main)
    identical(rownames(Binvignat.ref), colnames(Binvignat))
    Binvignat@meta.data$labels <- Binvignat.ref$labels
    save(Binvignat, file = "~/Studio/01MetaRA/BinvignatX.Rdata")
    load(file = "~/Studio/01MetaRA/BinvignatX.Rdata")
    FeaturePlot(Binvignat, features = sort(c("SUCLG1", "SUCLA2", "PDK3", "PDK1", "PDHX", "NDUFS1", "NDUFB8", "NDUFAF4", "COX7C", "COX7B", "COX11" )), raster = F)
    ggsave("~/Studio/01MetaRA/F009a.jpg", width = 11.5, height = 7.5)
    unique(Binvignat@meta.data$cell_type)
    filtered_celltypes <- gsub(" alpha-beta|, terminally differentiated", "", 
                      Binvignat@meta.data$cell_type)
    unique(filtered_celltypes)
    cell_type_med <- as.data.frame(Binvignat@reductions$umap@cell.embeddings) %>% 
        mutate(labels = filtered_celltypes) %>% 
        group_by(labels) %>%
        summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
    label = data.frame(symbol = c("SUCLG1", "SUCLA2", "PDK3", "PDK1", "PDHX", "NDUFS1", "NDUFB8", "NDUFAF4", "COX7C", "COX7B", "COX11"))
    plot02a = VlnPlot(Binvignat, features = sort(c("SUCLG1", "SUCLA2", "PDK3", "PDK1", "PDHX", "NDUFS1", "NDUFB8", "NDUFAF4", "COX7C", "COX7B", "COX11" )), raster = FALSE)
    plot01 = as.data.frame(Binvignat@reductions$umap@cell.embeddings) %>% 
        mutate(labels = filtered_celltypes) %>% 
        ggplot(aes(x = UMAP_1, y = UMAP_2, color = labels)) +
        geom_point(size = 0.005, alpha = 0.5) +
        theme_classic() + 
        scale_color_manual(values = color) +
        labs(title = "Binvignat et al. study") + 
        guides(color = guide_legend(title = '', override.aes = list(size = 2))) +
        ggrepel::geom_text_repel(aes(label = labels), fontface = "bold", size = 2, data = cell_type_med, point.padding=unit(0.5, "lines")) + 
        theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 1, vjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 8),
                  legend.position = 'none',
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 16, face = 'bold'))
    ggsave("~/Studio/01MetaRA/F009b.jpg", width = 4, height = 4)
    plot02a = VlnPlot(Binvignat, features = sort(c("SUCLG1", "SUCLA2", "PDK3", "PDK1", "PDHX", "NDUFS1", "NDUFB8", "NDUFAF4", "COX7C", "COX7B", "COX11" )), raster = FALSE)
    ggsave("~/Studio/01MetaRA/F009.jpg", width = 11.5, height = 8)
    Binvignat@meta.data$filtered_celltypes = filtered_celltypes
    plot02 = DotPlot(Binvignat, features = sort(c("SUCLG1", "SUCLA2", "PDK3", "PDK1", "PDHX", "NDUFS1", "NDUFB8", "NDUFAF4", "COX7C", "COX7B", "COX11" )), assay='RNA', group.by = 'filtered_celltypes')
    exp_mat = plot02$data %>% 
        select(-pct.exp, -avg.exp) %>%  
        pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
        as.data.frame()
    percent_mat = plot02$data %>% 
        select(-avg.exp, -avg.exp.scaled) %>%  
        pivot_wider(names_from = id, values_from = pct.exp) %>% 
        as.data.frame() %>% 
        pivot_longer(!features.plot, names_to = 'cell', values_to = 'value')
    plot03 = exp_mat %>% 
        pivot_longer(!features.plot, names_to = 'cell', values_to = 'module') %>% 
        cbind(percent_mat[, 3]) %>% 
        mutate(symbol = factor(.$features.plot, levels = rev(unique(.$features.plot)))) %>% 
        ggplot(aes(x = cell, y = symbol)) +
        # geom_tile(aes(fill = module)) + 
        geom_point(aes(size = value, color = value)) +
        scale_color_gradient(low = "#D4E5F4", high = "#F39B7F") +
        # geom_text(aes(label = text), size = 2.5) + 
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        # scale_x_discrete(position = "bottom") +
        theme_minimal() +
        labs(title = "Cell type") +
        scale_size(guide = 'none') +
        # guides(color = guide_legend(title = 'Percent')) + 
        theme(legend.position = "right",
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 8),
              panel.grid = element_blank(),
              axis.text.y = element_text(size = 12),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title = element_text(face = 'bold', size = 0),
              plot.title = element_text(size = 16, face = 'bold'))
    plot04 = exp_mat %>% 
        pivot_longer(!features.plot, names_to = 'cell', values_to = 'module') %>% 
        cbind(percent_mat[, 3]) %>% 
        mutate(symbol = factor(.$features.plot, levels = rev(unique(.$features.plot)))) %>% 
        mutate(y = 'x') %>% 
        ggplot(aes(x = cell, y = y, fill = cell)) +
        geom_tile(aes(fill = cell)) +
        scale_fill_manual(values = color) +
        scale_x_discrete(position = "bottom") +
        theme_minimal() +
        xlab(NULL) + ylab(NULL) +
        labs(title = '') + 
        guides(fill = 'none') + 
        theme(panel.border = element_blank(),
              panel.grid = element_blank(),
              legend.text = element_blank(),
              axis.ticks = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
              plot.title = element_blank(),
              legend.position = "none")
    require(aplot)
    plot04 %>% insert_top(plot03, height = 20)  
    ggsave("~/Studio/01MetaRA/F011a.jpg", width = 7.5, height = 6)
    x1 = as.factor(Binvignat@meta.data$label)
    names(x1) = colnames(Binvignat)
    str(x1)
    Binvignat@active.ident = x1
    head(Binvignat@active.ident)
#### percent
    cbinvignat = data.frame(disease = Binvignat@meta.data$disease, cell = Binvignat@meta.data$filtered_celltypes) %>% table(.) %>% chisq.test(.)
    plota2 = data.frame(disease = Binvignat@meta.data$disease, cell = Binvignat@meta.data$filtered_celltypes) %>% 
    unnest(cols = c(disease, cell)) %>% 
    dplyr::count(disease, cell) %>%
    group_by(disease) %>%
    mutate(percentage = n / sum(n) * 100) %>%
    ungroup() %>% 
    ggplot(aes(x = disease, y = percentage, fill = cell)) +
        geom_bar(stat = "identity", position = "fill", ) +
        geom_text(aes(label = sprintf("%.1f%%", percentage)), 
                position = position_fill(vjust = 0.5), size = 2.5) +
        scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
        coord_cartesian(ylim = c(0, 1.1)) + 
        scale_fill_manual(values = color) + 
        annotate("text", x = 1.5, y = 1.02, size = 3, label = "p-value < 2.2e-16", vjust = -0.5) + 
        theme_classic() +
        labs(title = "Binvignat et al. study", x = "", y = "Percentage", fill = '') + 
        scale_x_discrete(labels = c('Normal', 'Rheumatoid arthritis')) +
        theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.text.y = element_text(size = 10),
              legend.title = element_text(size = 12, face = 'bold'),
              legend.text = element_text(size = 8),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16, face = 'bold'))
    ggsave("~/Studio/01MetaRA/F009f.jpg", width = 4.8, height = 6.4)
    disease_med <- as.data.frame(Binvignat@reductions$umap@cell.embeddings) %>% 
        mutate(labels = Binvignat@meta.data$disease) %>% 
        group_by(labels) %>%
        summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
    plot03 = as.data.frame(Binvignat@reductions$umap@cell.embeddings) %>% 
        mutate(labels = Binvignat@meta.data$disease) %>% 
        ggplot(aes(x = UMAP_1, y = UMAP_2, color = labels)) +
        geom_point(size = 0.005, alpha = 0.5) +
        theme_classic() + 
        scale_color_manual(values = color) +
        labs(title = "Binvignat et al. study") + 
        guides(color = guide_legend(title = '', override.aes = list(size = 2))) +
        ggrepel::geom_text_repel(aes(label = labels), fontface = "bold", size = 4, data = disease_med, point.padding = unit(0.5, "lines"), segment.colour="black") + 
        theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 1, vjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 10),
                  legend.position = 'none',
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 16, face = 'bold'))
    ggsave("~/Studio/01MetaRA/F009s.jpg", width = 4, height = 3.9)
#### GSVA
    .libPaths(c("/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SeuratV4", .libPaths()))
    packageVersion('Seurat') # 4.4.0
    packageVersion("SeuratObject")
    require(Seurat)
    require(tidyverse)
    require(UCell)
    require(irGSEA)
    gs_th = readxl::read_excel('~/Studio/01MetaRA/MetabolicPathways.xlsx', sheet = 1)
    gs_path = split(gs_th$symbol, list(gs_th$Pathway))
    load('~/Studio/01MetaRA/BinvignatX.Rdata')
    require(GSVA)
    x1 = as.factor(Binvignat@meta.data$disease)
    names(x1) = colnames(Binvignat)
    str(x1)
    Binvignat@active.ident = x1
    gsvaB = gsva(gsvaParam(as.matrix(Binvignat@assays$RNA@data), gs_path, maxDiff = TRUE))
    save(gsvaB, file = 'x7.Rdata')
    load('~/Studio/01MetaRA/x7.Rdata')
    str(gsvaB)
    str(Binvignat)
    identical(colnames(gsvaB), rownames(Binvignat@meta.data))
    dis = data.frame(disease = Binvignat@meta.data$disease, labels = Binvignat@meta.data$filtered_celltypes) 
    gsva = as.data.frame(gsvaB) %>% 
        t() %>% 
        as.data.frame() %>% 
        cbind(dis) %>% 
        mutate(disease = ifelse(disease == 'normal', 'NC', 'RA'))
    require(limma)
    g = factor(gsva$disease, levels = c('NC', 'RA'))
    design = model.matrix( ~ 1 + g)
    gsvaBinvignat = t(gsva[, 1:7]) %>% 
        as.data.frame() %>% 
        lmFit(., design) %>% 
        eBayes(.) %>% 
        topTable(., coef = 2, number = Inf) %>% 
        mutate(sign = ifelse(adj.P.Val <= 0.001, '***', ifelse(adj.P.Val <= 0.01, '**', ifelse(adj.P.Val <= 0.05, '*', 'ns.')))) %>% 
        tibble::rownames_to_column('Meta')
    require(ggpubr)
    require(ggsci)
    plot005 = gsvaBinvignat %>% 
        mutate(log10P = - log10(adj.P.Val)) %>% 
        mutate(group = ifelse(adj.P.Val <= 0.05 & logFC > 0, 'Up', ifelse(adj.P.Val <= 0.05 & logFC < 0, 'Down', 'Not'))) %>% 
        ggplot(aes(x = logFC, y = log10P, color = group)) +
        geom_point(aes()) + 
        scale_color_manual(values = c("#B7B7EB", '#9E9E9E', "#FBB9BA")) +
        ggrepel::geom_text_repel(aes(label = Meta), size = 3, hjust = 0, vjust = -1) + 
        geom_vline(xintercept = c(0), colour = "black", linetype = "dashed") +
        geom_hline(yintercept = -log10(0.05), colour = "black", linetype = "dashed") +
        labs(title = "GSVA", x = 'log2(FoldChange)', y = '-log10(adj.p.value)') + 
        theme_classic2() +
        guides(color = guide_legend(title = "")) +
        theme(axis.text.x = element_text(size = 10, face = 'bold'),
              axis.text.y = element_text(size = 10),
              legend.position = 'none',
              legend.title = element_text(size = 12, face = 'bold'),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16, face = 'bold'))
    ggsave('~/Studio/01MetaRA/F012a.jpg', width = 2.8, height = 4.2)
#### UMAP for GSVA
    identical(rownames(Binvignat@reductions$umap@cell.embeddings), rownames(gsva))
    all_values <- unlist(gsva %>% select('Amino acid', 'Carbohydrate', 'Energy', 'Lipid', 'Nucleotide', 'TCA cycle', 'Vitamin cofactor'))
    common_limits <- range(all_values, na.rm = TRUE)
    create_plot <- function(data, title, colname) {
        data %>% 
        select(all_of(colname)) %>% 
        dplyr::rename('AC' = all_of(colname)) %>% 
        tibble::rownames_to_column('name') %>% 
        inner_join(as.data.frame(Binvignat@reductions$umap@cell.embeddings) %>% 
               tibble::rownames_to_column('name'), by = 'name') %>% 
        ggplot(aes(x = UMAP_1, y = UMAP_2)) +
        geom_point(aes(color = AC), size = 0.005) +
        theme_classic() + 
        scale_color_gradient(low = "#D4E5F4",  high = "#FBB9BA", limits = common_limits) +
        labs(title = title) + 
        guides(color = guide_legend(title = '', override.aes = list(size = 2))) +
        theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            legend.title = element_text(size = 12, face = 'bold'),
            legend.text = element_text(size = 10),
            legend.position = 'right',
            axis.title = element_text(size = 10),
            plot.title = element_text(size = 16, face = 'bold'))}
    fit011 <- create_plot(gsva, "Amino acid", "Amino acid")
    fit012 <- create_plot(gsva, "Carbohydrate", "Carbohydrate")
    fit013 <- create_plot(gsva, "Energy", "Energy")
    fit014 <- create_plot(gsva, "Lipid", "Lipid")
    fit015 <- create_plot(gsva, "Nucleotide", "Nucleotide")
    fit016 <- create_plot(gsva, "TCA cycle", "TCA cycle")
    fit017 <- create_plot(gsva, "Vitamin cofactor", "Vitamin cofactor")
    require(patchwork)
    combined_plot <- fit011 + fit012 + fit013 + fit014 + fit015 + fit016 + fit017 + plot_layout(ncol = 4, guides = 'collect')
    ggsave('~/Studio/01MetaRA/F012b.jpg', width = 9, height = 4.5) 
#### corplot
    str(Binvignat)
    key = data.frame(symbol = factor(c("COX7C", "COX7B", "PDK1", "PDK3", "COX11", "NDUFAF4", "NDUFB8", "NDUFS1", "SUCLG1", "PDHX", "SUCLA2"), levels = c("COX7C", "COX7B", "PDK1", "PDK3", "COX11", "NDUFAF4", "NDUFB8", "NDUFS1", "SUCLG1", "PDHX", "SUCLA2")))
    str(Binvignat@assays$RNA@data)
    z1 <- as.data.frame(Binvignat@assays$RNA@data)[rownames(Binvignat@assays$RNA@data) %in% key$symbol, ]
    cbz = t(z1) %>% 
        as.data.frame() %>% 
        cbind(., gsva)
    cor_matrix = cor(cbz[, c(1:18)])
    p.mat = corrplot::cor.mtest(cbz[, c(1:18)], conf.level = 0.95)$p
    fitc = as.data.frame(cor_matrix[1:11, 12:18]) %>% 
        tibble::rownames_to_column('symbol') %>% 
        pivot_longer(!symbol, names_to = 'module', values_to = 'value') %>% 
        mutate(module = str_remove(module, "^value\\."))
    fitc1 = as.data.frame(p.mat[1:11, 12:18]) %>% 
        tibble::rownames_to_column('symbol') %>% 
        pivot_longer(!symbol, names_to = 'module', values_to = 'pvalue') %>% 
        mutate(module = str_remove(module, "^value\\."))    
    fitc$symbol = factor(fitc$symbol, levels = rev(c("COX11", "COX7B", "COX7C", "NDUFAF4", "NDUFB8", "NDUFS1",  "PDHX", "PDK1", "PDK3", "SUCLA2", "SUCLG1")))
    fitc$module = factor(fitc$module, levels = rev(c("Amino acid", "Carbohydrate", "Energy",  "Lipid", "Nucleotide", "TCA cycle", "Vitamin cofactor")))
    plot003 = fitc %>% 
        mutate(pvalue = fitc1$pvalue) %>% 
        mutate(sig = ifelse(pvalue < 0.001, '***', ifelse(pvalue < 0.01, '**', ifelse(pvalue < 0.05, '*', 'ns.')))) %>% 
        mutate(text = paste0(round(value, 2), '(', sig, ')'))  %>% 
        ggplot(aes(x = module, y = symbol)) +
        # geom_point(aes(size = value, color = pvalue)) +
        geom_tile(aes(fill = value)) + 
        scale_fill_gradient2(low = "#D4E5F4", mid = "#FFEBAD44", high = "#FBB9BA") +
        geom_text(aes(label = sig), size = 2.5) + 
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        # scale_x_discrete(position = "bottom") +
        theme_minimal() +
        labs(title = "Correction heatmap") +
        guides(fill = guide_legend(title = 'Correlation\ncoefficient')) + 
        theme(legend.position = "right",
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 10),
              panel.grid = element_blank(),
              axis.text.y = element_text(size = 12),
              axis.ticks = element_blank(),
              plot.title = element_text(face = 'bold', size = 16),
              axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
              axis.title = element_text(face = 'bold', size = 0))
    ggsave("~/Studio/01MetaRA/F010e.jpg", width = 4, height = 4.8)
#### Difference analysis
    plot007 = t(z1) %>% 
        as.data.frame() %>% 
        mutate(disease = Binvignat@meta.data$disease) %>% 
        mutate(disease = ifelse(disease == 'normal', 'Normal', 'Rheumatoid arthritis')) %>% 
        pivot_longer(!disease, names_to = 'symbol', values_to = 'value') %>% 
        ggplot(aes(x = symbol, y = value, fill = disease)) +
        geom_violin() +
        # geom_boxplot(alpha = 0.7, width = 0.5) +
        scale_fill_manual(values = color) +
        ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test", size = 4) +
        theme_classic() + 
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        # scale_x_discrete(position = "bottom") +
        labs(title = "Difference analysis", y = "Gene expression", x = "") +
        guides(fill = guide_legend(title = '')) + 
        theme(legend.position = "right",
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 10),
              panel.grid = element_blank(),
              axis.text.y = element_text(size = 12),
              axis.ticks = element_blank(),
              plot.title = element_text(face = 'bold', size = 16),
              axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
              axis.title = element_text(size = 12))
    ggsave("~/Studio/01MetaRA/F010b.jpg", width = 7.5, height = 4.4)
## Figure 6 COX7C analysis
#### Binvignat Umap for COX7C
    .libPaths(c("/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SeuratV4", .libPaths()))
    packageVersion('Seurat') # 4.4.0
    packageVersion("SeuratObject")
    require(Seurat)
    ## load("/Users/linqi/Studio/01MetaRA/BinvignatscRNA.Rdata") ## V4a
    require(tidyverse)
    load(file = '~/Studio/01MetaRA/BinvignatX.Rdata')
    str(Binvignat)
    x1 = as.factor(Binvignat@meta.data$labels)
    names(x1) = colnames(Binvignat)
    str(x1)
    Binvignat@active.ident = x1
    head(Binvignat@active.ident)
    key = data.frame(symbol = factor(c("COX7C", "COX7B", "PDK1", "PDK3", "COX11", "NDUFAF4", "NDUFB8", "NDUFS1", "SUCLG1", "PDHX", "SUCLA2"), levels = c("COX7C", "COX7B", "PDK1", "PDK3", "COX11", "NDUFAF4", "NDUFB8", "NDUFS1", "SUCLG1", "PDHX", "SUCLA2")))
    z1 <- as.data.frame(Binvignat@assays$RNA@data)[rownames(Binvignat@assays$RNA@data) %in% key$symbol, ]
    z2 = as.data.frame(t(z1)) %>% select('COX7C') %>% cbind(disease = Binvignat@meta.data$disease, labels = Binvignat@meta.data$filtered_celltypes) %>% 
        group_by(labels, disease) %>%
        summarise(mean_COX7C = mean(COX7C, na.rm = TRUE),
                  sd_COX7C = sd(COX7C, na.rm = TRUE),
                  n_cells = n(),
                  se_COX7C = sd_COX7C / sqrt(n_cells),
                  .groups = "drop") %>% 
        as.data.frame()
    plotb2 = as.data.frame(Binvignat@reductions$umap@cell.embeddings) %>% 
        cbind(COX7C = as.data.frame(t(z1)) %>% select('COX7C')) %>% 
        ggplot(aes(x = UMAP_1, y = UMAP_2, color = COX7C)) +
        geom_point(aes(color = COX7C), size = 0.005, alpha = 0.8) +
        theme_classic() + 
        scale_color_gradient2(low = "#D4E5F4", mid = "#FFFEEE55", high = "#FBB9BA") +
        labs(title = "Binvignat et al. study") + 
        guides(color = guide_legend(title = '', override.aes = list(size = 2))) +
        ## ggrepel::geom_text_repel(aes(label = label), fontface = "bold", size = 4, data = cell_type_med, point.padding=unit(0.5, "lines")) + 
        theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 1, vjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 10),
                  # legend.key = element_rect(fill = "white", colour = "black"),
                  legend.position = 'right',
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 16, face = 'bold'))
    ggsave("~/Studio/01MetaRA/F009r.jpg", width = 4, height = 3.5)
    load(file = '~/Studio/01MetaRA/x7.Rdata')
    identical(colnames(gsvaB), rownames(Binvignat@meta.data))
    dis = data.frame(disease = Binvignat@meta.data$disease, labels = Binvignat@meta.data$filtered_celltypes) 
    gsva = as.data.frame(gsvaB) %>% 
        t() %>% 
        as.data.frame() %>% 
        cbind(dis) %>% 
        mutate(disease = ifelse(disease == 'normal', 'NC', 'RA'))
    y1 = z1 %>% 
        filter(rownames(.) %in% 'COX7C') %>% 
        t() %>% 
        as.data.frame() %>% 
        cbind(gsva)
    plot004 = y1 %>% 
        ggplot(aes(x = labels, y = COX7C, fill = disease)) +
        geom_violin() +
        # geom_boxplot(alpha = 0.7, width = 0.5) +
        ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test", size = 3) +
        theme_classic() +
        # coord_flip() + 
        scale_fill_manual(values = color, labels = c('Normal', 'Rheumatoid arthritis')) +
        labs(title = "Difference analysis of COX7C expression", x = "", y = "Gene expression", fill = "") + 
        theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.text.y = element_text(size = 10),
              legend.position = 'right',
              legend.title = element_text(size = 12, face = 'bold'),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16, face = 'bold'))
    require(patchwork)
    plot004 + plot_layout(design = '##AAAAAAAAAAAAAAAAAAAA')
    ggsave("~/Studio/01MetaRA/F011e.jpg", width = 8, height = 5.5)
    x0 = data.frame()
    for(j in unique(y1$labels)) {
        y01 <- y1 %>% filter(labels == j)
        if(nrow(y01) >= 3) {
            for(i in 2:8) {
            x1 <- cor.test(y01[, 1], y01[, i])
            x2 <- data.frame(cor = x1$estimate, pvalue = x1$p.value, meta = colnames(y01)[i], cell = j, stringsAsFactors = FALSE)
        x0 <- rbind(x0, x2)}
    } else {
    warning(paste("Not enough samples for cell type:", j, 
                 "- Skipping (n =", nrow(y01), ")"))}}
    plot008 = x0 %>% 
        mutate(sig = ifelse(pvalue < 0.001, '***', ifelse(pvalue < 0.01, '**', ifelse(pvalue < 0.05, '*', 'ns.')))) %>% 
        mutate(text = paste0(round(cor, 2), '(', sig, ')')) %>% 
        ggplot(aes(x = cell, y = meta)) + 
        geom_tile(aes(fill = cor)) + 
        scale_fill_gradient2(low = "#D4E5F4", mid = "#FFEBAD44", high = "#FBB9BA") +
        geom_text(aes(label = sig), size = 3) + 
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        # scale_x_discrete(position = "bottom") +
        theme_minimal() +
        labs(title = "Binvignat et al. study") +
        guides(fill = guide_legend(title = 'Correlation\ncoefficient')) + 
        theme(legend.position = "right",
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 10),
              panel.grid = element_blank(),
              axis.text.y = element_text(size = 12),
              axis.ticks = element_blank(),
              plot.title = element_text(face = 'bold', size = 16),
              axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
              axis.title = element_text(face = 'bold', size = 0))
    ggsave("~/Studio/01MetaRA/F011b.jpg", width = 8.5, height = 4.5)
#### GSE109449
    load(file = '~/Studio/01MetaRA/G109449.Rdata')
    plot005 = as.data.frame(G109449@assays$RNA@data) %>% 
        filter(rownames(.) %in% 'COX7C') %>% 
        t() %>% 
        as.data.frame() %>% 
        cbind(as.data.frame(G109449@reductions$umap@cell.embeddings)) %>% 
        ggplot(aes(x = UMAP_1, y = UMAP_2, color = COX7C)) +
        geom_point(aes(color = COX7C), alpha = 0.6) +
        theme_classic() + 
        scale_color_gradient2(low = "#D4E5F4", mid = "#FFEBAD", high = "#FBB9BA") +
        labs(title = "GSE109449") + 
        guides(color = guide_legend(title = '', override.aes = list(size = 2))) +
        ## ggrepel::geom_text_repel(aes(label = label), fontface = "bold", size = 4, data = cell_type_med, point.padding=unit(0.5, "lines")) + 
        theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 1, vjust = 1),
                  axis.text.y = element_text(size = 10),
                  legend.title = element_text(size = 12, face = 'bold'),
                  legend.text = element_text(size = 10),
                  # legend.key = element_rect(fill = "white", colour = "black"),
                  legend.position = 'right',
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 16, face = 'bold'))
    require(patchwork)
    plotb2 / plot005
    ggsave("~/Studio/01MetaRA/F011c.jpg", width = 3.5, height = 5.5)
    ggsave("~/Studio/01MetaRA/F013e.jpg", width = 4, height = 3.5)
    plot006 = as.data.frame(G109449@assays$RNA@data) %>% 
        filter(rownames(.) %in% 'COX7C') %>% 
        t() %>% 
        as.data.frame() %>% 
        cbind(as.data.frame(t(gsvax)))
    x0 = data.frame()
    for(i in c(2:8)){
        x1 = cor.test(plot006[, 1], plot006[, i])
        x2 = data.frame(cor = x1$estimate, pvalue = x1$p.value, meta = colnames(plot006)[i])
        x0 = rbind(x2, x0)
    }
    plot007 = x0 %>% 
        mutate(group = 'COX7C') %>% 
        mutate(sig = ifelse(pvalue < 0.001, '***', ifelse(pvalue < 0.01, '**', ifelse(pvalue < 0.05, '*', 'ns.')))) %>% 
        mutate(text = paste0(round(cor, 2), '(', sig, ')')) %>% 
        ggplot(aes(x = group, y = meta)) + 
        geom_tile(aes(fill = cor)) + 
        scale_fill_gradient2(low = "#D4E5F4", mid = "#FFEBAD44", high = "#FBB9BA") +
        geom_text(aes(label = sig), size = 3) + 
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        # scale_x_discrete(position = "bottom") +
        theme_minimal() +
        labs(title = "GSE109449") +
        guides(fill = guide_legend(title = 'Correlation\ncoefficient')) + 
        theme(legend.position = "right",
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 10),
              panel.grid = element_blank(),
              axis.text.y = element_text(size = 12),
              axis.ticks = element_blank(),
              plot.title = element_text(face = 'bold', size = 16),
              axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
              axis.title = element_text(face = 'bold', size = 0))
    ggsave("~/Studio/01MetaRA/F013f.jpg", width = 3, height = 2.75)
#### GSE93776
    require(GEOquery)
    g93776c = getGEO(filename = ("~/Studio/01MetaRA/GSE93776_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        pData(.)
    g93776e = getGEO(filename = ("~/Studio/01MetaRA/GSE93776_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F) %>% 
        exprs(.) %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('probe_id') %>% 
        inner_join(AnnoProbe::idmap('GPL570'), by = 'probe_id') %>% 
        select('symbol', everything()) %>% 
        mutate(rowMean = rowMeans(.[grep('GSM', names(.))])) %>%
        dplyr::arrange(desc(rowMean)) %>%
        dplyr::rename('symbol' = 'symbol') %>% 
        dplyr::distinct(symbol, .keep_all = T) %>%
        dplyr::select(-c('rowMean', 'probe_id')) %>% 
        tibble::column_to_rownames('symbol')
    fit = g93776e %>% 
        t() %>% 
        as.data.frame() %>% 
        select(key$symbol)
    plota7 = fit %>% 
        tibble::rownames_to_column('geo_accession') %>% 
        inner_join(g93776c[, c('geo_accession', 'cell type:ch1', 'disease:ch1')], by = 'geo_accession') %>% 
        dplyr::rename('cell' = 'cell type:ch1', 'disease' = 'disease:ch1') %>% 
        filter(cell != 'B.naive') %>% 
        filter(cell != 'B.memory') %>% 
        mutate(disease = ifelse(disease == 'healthy control', 'Healthy control', 'Rheumatoid arthritis')) %>% 
        # Rmisc::summarySE(., measurevar = 'COX7C', groupvars = c('disease', 'cell')) %>% 
        ggplot(aes(x = cell, y = lasso, fill = disease)) +
        geom_boxplot(alpha = 0.7, width = 0.5) +
        ggpubr::stat_compare_means(label = "p.format", method = "t.test", size = 3) +
        theme_classic() +
        scale_fill_manual(values = color) +
        labs(title = "GSE93776", x = "", y = "Gene expression") + 
        theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.text.y = element_text(size = 10),
              legend.position = 'right',
              legend.title = element_text(size = 12, face = 'bold'),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16, face = 'bold'))
    plota7 = g93776e['COX7C', ] %>% 
        t() %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('geo_accession') %>% 
        inner_join(g93776c[, c('geo_accession', 'cell type:ch1', 'disease:ch1')], by = 'geo_accession') %>% 
        dplyr::rename('cell' = 'cell type:ch1', 'disease' = 'disease:ch1') %>% 
        filter(cell != 'B.naive') %>% 
        filter(cell != 'B.memory') %>% 
        mutate(disease = ifelse(disease == 'healthy control', 'Healthy control', 'Rheumatoid arthritis')) %>% 
        # Rmisc::summarySE(., measurevar = 'COX7C', groupvars = c('disease', 'cell')) %>% 
        ggplot(aes(x = cell, y = COX7C, fill = disease)) +
        geom_violin(alpha = 0.7, linewidth = 0.25) +
        ggpubr::stat_compare_means(label = "p.format", method = "t.test", size = 3) +
        theme_classic() +
        scale_fill_manual(values = color) +
        labs(title = "GSE93776", x = "", y = "Gene expression", fill = "") + 
        theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.text.y = element_text(size = 10),
              legend.position = 'right',
              legend.title = element_text(size = 12, face = 'bold'),
              legend.text = element_text(size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16, face = 'bold'))
    require(patchwork)
    plota7 + plot_layout(design = '#AAAAAAAAAAAAAAAAAA')
    ggsave("~/Studio/01MetaRA/F009h.jpg", width = 12, height = 4)
#### Drug Enrichr
    plot09 = data.table::fread('~/Studio/01MetaRA/COX7CDrug.txt', header = T, sep = '\t') %>% 
        dplyr::rename('score' = 'Combined score') %>% 
        mutate(symbol = c('COX7C')) %>% 
        ggplot(aes(x = reorder(Name, Index), y = symbol)) + 
        # geom_point(aes(size = score, color = Adjusted_p_value)) +
        geom_tile(aes(fill = Adjusted_p_value), color = "white") +
        scale_fill_gradient(high = "#D4E5F4", low = "#FBB9BA") +
        scale_size_continuous(range = c(2, 8)) + 
        geom_text(aes(label = score), size = 3) + 
        # scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        # scale_x_discrete(position = "bottom") +
        theme_minimal() +
        labs(y = "", title = "DSigDB Drug Enrichment (COX7C)") +
        guides(fill = guide_legend(title = 'Adjusted\np-value')) + 
        theme(legend.position = "bottom",
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 10),
              panel.grid = element_blank(),
              axis.text.y = element_text(size = 12),
              axis.ticks = element_blank(),
              axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
              axis.title = element_text(face = 'bold', size = 0),
              plot.title = element_text(size = 18, face = 'bold'))
    data.table::fread('~/Studio/01MetaRA/COX7CDrug.txt', header = T, sep = '\t') %>% 
        write.csv(., "~/Studio/01MetaRA/enrichr3.csv")
    require(patchwork)
    plot09 + plot_layout(design = '#AAAAAAAAAAAAAAAAAAAAAAAAAA')
    ggsave('~/Studio/01MetaRA/F013h.jpg', width = 11.5, height = 3.5)
