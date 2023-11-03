## Packages for this study.
    library(tidyverse)              # Data modulation and plotting.
    library(data.table)             # Data convenience features.
    require(readxl)                 # Reading and writing from excel.
    require(Rmisc)                  # Calculation SE for ploting SE bar.
    require(ggsignif)               # Significance of ploting using ggplot2.
    require(scales)                 # Scale Functions for Visualization.
    color = c("#FA7F6F", "#8ECFC9", "#FFBE7A", "#82B0D2", 
              "#6b798e", "#BEB8DC", "#E7DAD2", "#928187",
              "#C4A5DE", "#CFEAF1", "#F6CAE5", "#9E9E9E")
## Figure 1 S100A8 in DLBCL patients 
#### GSE56315 based on GPL570
    require(GEOquery)
    G56315 = getGEO(filename = ('~/Movies/01S100A8/GSE56315_series_matrix.txt.gz'), GSEMatrix = T, 
            getGPL = F)
    if (length(G56315) > 1) {idx = grep("GPL570", attr(G56315, "names"))} else {idx = 1}
    g56315e = data.frame(exprs(G56315))
    g56315c = pData(G56315)
    qx = as.numeric(quantile(g56315e, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC = (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if(LogC){g56315e[which(g56315e <= 0)] = NaN
            G56315e = log2(g56315e)
            print("log2 transform finished")}else{print("log2 transform not needed")}
    count(g56315c$source_name_ch1)
    GPL570 = getGEO(filename = ('~/Movies/01S100A8/GPL570.annot.gz')) %>% 
            Table(.) %>%                   
            dplyr::select('ID','Gene symbol') %>% 
            filter('Gene symbol' != "---") %>% 
            dplyr::rename('symbol' = 'Gene symbol', 'probe_id' = 'ID')
    s56315 = g56315e %>% 
        rownames_to_column('probe_id') %>%
        merge(GPL570, by = 'probe_id') %>%
        dplyr::select(- 'probe_id')  %>%
        dplyr::select(symbol, everything()) %>% 
        mutate(rowMean = rowMeans(.[grep('GSM', names(.))])) %>% 
        arrange(desc(rowMean)) %>% 
        distinct(symbol, .keep_all = T)  %>%
        dplyr::select(-'rowMean') %>% 
        filter(symbol == 'S100A8') %>% 
        gather() %>% 
        dplyr::rename('geo_accession' = 'key', 'S100A8' = 'value') %>% 
        slice(-1) %>% 
        inner_join(g56315c, by = 'geo_accession') %>% 
        mutate(group = ifelse(source_name_ch1 == 'diagnostic diffuse large B-cell lymphoma sample', 'DLBCL', 'Normal')) %>% 
        mutate_at('S100A8', as.numeric)
    fit = s56315 %>% 
        ggplot(aes(fill = group, y = S100A8, x = group)) +
        geom_boxplot(alpha = 0.8)  +
        geom_jitter(data = s56315, aes(y = S100A8), size = 2, shape = 21,
            stroke = 0, show.legend = FALSE, width = 0.1) +
        geom_signif(map_signif_level = F, parse = T,
            comparisons = list(c("Normal", "DLBCL")), size = 0.5, textsize = 3, vjust = 0.1, test = "wilcox.test") +
        labs(title = "GSE56315", y = "The expression of S100A8", x = "") +
        scale_fill_manual(values = color) + 
        theme_classic() +
        theme(legend.position = "none",
        axis.text = element_text(size = 14, vjust = 0.5),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16))
    ggsave("~/Movies/01S100A8/F_G56315Com.png", width = 3, height = 5)
#### GSE83632 based on GPL5175
    G83632 = getGEO(filename = ("~/Movies/01S100A8/GSE83632_series_matrix.txt.gz"), GSEMatrix = TRUE, getGPL = FALSE)
    if (length(G83632) > 1) {
        idx = grep("GPL5175", attr(G83632, "names"))} else {
        idx = 1}
    g83632e = data.frame(exprs(G83632))
    g83632c = pData(G83632)
    qx = as.numeric(quantile(g83632e, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    LogC = (qx[5] > 100) ||
        (qx[6] - qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) {
        g83632e[which(g83632e <= 0)] = NaN
        g83632e = log2(g83632e)
        print("log2 transform finished")} else {
        print("log2 transform not needed")}
    an5175 = read.table(
        file = "~/Movies/01S100A8/GPL5175_noParents.an.txt",
        header = T, sep = "\t", quote = "", fill = T,
        comment.char = "!") %>%
        select("ProbeName", "GeneSymbols") %>% 
        dplyr::rename("ID" = "ProbeName", "symbol" = "GeneSymbols")
    s83632 = g83632e %>% 
        rownames_to_column('ID') %>%
        merge(an5175, by = 'ID') %>%
        dplyr::select(- 'ID')  %>%
        dplyr::select(symbol, everything()) %>% 
        mutate(rowMean = rowMeans(.[grep('GSM', names(.))])) %>% 
        arrange(desc(rowMean)) %>% 
        distinct(symbol, .keep_all = T)  %>%
        dplyr::select(-'rowMean') %>% 
        filter(symbol == 'S100A8') %>% 
        gather() %>% 
        dplyr::rename('geo_accession' = 'key', 'S100A8' = 'value') %>% 
        slice(-1) %>% 
        inner_join(g83632c, by = 'geo_accession') %>% 
        mutate(group = ifelse(characteristics_ch1.1 == 'condition: DLBCL at diagnosis', 'DLBCL', 'Normal')) %>% 
        mutate_at('S100A8', as.numeric)
    fit = s83632 %>% 
        ggplot(aes(fill = group, y = S100A8, x = group)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(data = s83632, aes(y = S100A8), size = 2, shape = 21,
            stroke = 0, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = F, parse = T,
            comparisons = list(c("Normal", "DLBCL")), size = 0.5, textsize = 3, vjust = 0.1, test = "wilcox.test") +
        labs(title = "GSE83632", y = "The expression of S100A8", x = "") +
        #scale_y_continuous(limits = c(3, 12)) +
        scale_fill_manual(values = color) + theme_classic() +
        theme(legend.position = "none",
        axis.text = element_text(size = 14, vjust = 0.5),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16))
    ggsave("~/Movies/01S100A8/F_G83632Com.png", width = 3, height = 5)
#### ROC curve
    require(plotROC)
    require(pROC)
    r56315 = roc(response = s56315$group, predictor = s56315$S100A8)
    auc = round(auc(s56315$group, s56315$S100A8), 4)
    ci.auc(r56315)
    ggroc(r56315, colour = "#82B0D2", size = 1, legacy = T) +
        geom_abline(linetype = 3, intercept = 0) +
        labs(title = "GSE56315") +
        theme_classic() +
        theme(axis.text = element_text(size = 14, vjust = 0.5),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16)) +
        annotate("text", x = 0.75, y = 0.1, label = "AUC = 1", size = 6)
    ggsave("~/Movies/01S100A8/F_G56315roc.png", width = 4.5, height = 5)
    r83632 = roc(response = s83632$group, predictor = s83632$S100A8)
    auc = round(auc(s83632$group, s83632$S100A8), 4)
    ci.auc(r83632)
    ggroc(r83632, colour = "#82B0D2", size = 1, legacy = T) +
        geom_abline(linetype = 3, intercept = 0) +
        labs(title = "GSE83632") +
        theme_classic() +
        theme(axis.text = element_text(size = 14, vjust = 0.5),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16)) +
        annotate("text", x = 0.75, y = 0.1, label = "AUC = 0.9282", size = 6)
    ggsave("~/Movies/01S100A8/F_G83632roc.png", width = 4.5, height = 5)
#### Volcano plot GSE56315
    g56315l = g56315c %>% 
        mutate(states = ifelse(characteristics_ch1 == 'tissue: diffuse large B-cell lymphoma (DLBCL)', 'DLBCL', 'Normal')) %>%
        select(-characteristics_ch1.1)
    g = g56315l$states
    design = model.matrix(~ -1 + g)
    colnames(design) = c("DLBCL", "Normal")
    contrast.matrix = makeContrasts(DLBCL - Normal, levels = design)
    d56315 = g56315e %>% 
        rownames_to_column('probe_id') %>%
        merge(GPL570, by = 'probe_id') %>%
        dplyr::select(- 'probe_id') %>%
        dplyr::select(symbol, everything()) %>% 
        mutate(rowMean = rowMeans(.[grep('GSM', names(.))])) %>% 
        arrange(desc(rowMean)) %>% 
        distinct(symbol, .keep_all = T)  %>%
        dplyr::select(-'rowMean') %>% 
        column_to_rownames('symbol') %>% 
        lmFit(., design) %>% 
        contrasts.fit(., contrast.matrix) %>% 
        eBayes(.) %>% 
        topTable(., n = Inf, adjust = "fdr") %>% 
        na.omit(.) %>% 
        mutate(group = ifelse(adj.P.Val <= 0.05 & logFC > 1, 'Up', ifelse(adj.P.Val <= 0.05 & logFC < -1, 'Down', 'Not'))) %>% 
        tibble::rownames_to_column("symbol") %>% 
        mutate(log10P = -log10(adj.P.Val)) 
    fit = d56315 %>% 
        ggscatter("logFC", "log10P",
        combine = F, merge = T, color = "group", shape = 20, size = 3,
        point = TRUE, font.label = 20, palette = c("#8ECFC9", "grey", "#FA7F6F")) +
        geom_vline(xintercept = c(-1, 1), colour = "black", linetype = "dashed") +
        geom_hline(yintercept = -log10(0.05), colour = "black", linetype = "dashed") +
        labs(title = "GSE56315", x = 'log2FC') + 
        theme_classic() + 
        #coord_flip() +
        ggrepel::geom_label_repel(aes(label = symbol), data = filter(d56315, symbol == "S100A8")) +
        theme(legend.position=c(.9, .88),
        legend.title = element_blank(),
        axis.text = element_text(size = 14, vjust = 0.5),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16))
    ggsave("~/Movies/01S100A8/F_G56315forrest.png", width = 4, height = 5)
#### Volcano plot GSE83632
    g83632l = g83632c %>% 
        mutate(states = ifelse(source_name_ch1 == "Healthy donor, blood", "Health", "DLBCL")) %>%
        select(-source_name_ch1)
    g = g83632l$states
    design = model.matrix(~ -1 + g)
    colnames(design) = c("DLBCL", "Normal")
    contrast.matrix = makeContrasts(DLBCL - Normal, levels = design)
    d83632 = g83632e %>% 
        rownames_to_column('ID') %>%
        merge(an5175, by = 'ID') %>%
        dplyr::select(- 'ID')  %>%
        dplyr::select(symbol, everything()) %>% 
        mutate(rowMean = rowMeans(.[grep('GSM', names(.))])) %>% 
        arrange(desc(rowMean)) %>% 
        distinct(symbol, .keep_all = T)  %>%
        dplyr::select(-'rowMean') %>% 
        column_to_rownames('symbol') %>% 
        lmFit(., design) %>% 
        contrasts.fit(., contrast.matrix) %>% 
        eBayes(.) %>% 
        topTable(., n = Inf, adjust = "fdr") %>% 
        na.omit(.) %>% 
        mutate(group = ifelse(adj.P.Val <= 0.05 & logFC > 1, 'Up', ifelse(adj.P.Val <= 0.05 & logFC < -1, 'Down', 'Not'))) %>% 
        tibble::rownames_to_column("symbol") %>% 
        mutate(log10P = -log10(adj.P.Val)) 
    fit = d83632 %>% 
        ggscatter("logFC", "log10P",
        combine = F, merge = T, color = "group", shape = 20, size = 3,
        point = TRUE, font.label = 20, palette = c("#8ECFC9", "grey", "#FA7F6F")) +
        geom_vline(xintercept = c(-1, 1), colour = "black", linetype = "dashed") +
        geom_hline(yintercept = -log10(0.05), colour = "black", linetype = "dashed") +
        labs(title = "GSE83632", x = 'log2FC') + 
        theme_classic() + 
        #coord_flip() +
        ggrepel::geom_label_repel(aes(label = symbol), data = filter(d83632, symbol == "S100A8")) +
        theme(legend.position=c(.9, .88),
        legend.title = element_blank(),
        axis.text = element_text(size = 14, vjust = 0.5),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16))
    ggsave("~/Movies/01S100A8/F_G83632forrest.png", width = 4, height = 5)
#### Survival curve of GSE87371
    require(GEOquery)
    require(survival)
    require(survminer)
    G87371 = getGEO(filename = "~/Movies/01S100A8/GSE87371_series_matrix.txt.gz", GSEMatrix = T, 
                getGPL = F)
    if (length(G87371) > 1) {
        idx = grep("GPL570", attr(G87371, "names"))} else {
        idx = 1}
    g87371e = data.frame(exprs(G87371))
    g87371c = pData(G87371)
    qx = as.numeric(quantile(g87371e, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    LogC = (qx[5] > 100) ||
        (qx[6] - qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) {
        g87371e[which(g87371e <= 0)] = NaN
        g87371e = log2(g87371e)
        print("log2 transform finished")} else {
        print("log2 transform not needed")}
    GPL570 = getGEO(filename = ('~/Movies/01S100A8/GPL570.annot.gz')) %>% 
            Table(.) %>%                   
            dplyr::select('ID','Gene symbol') %>% 
            filter('Gene symbol' != "---") %>% 
            dplyr::rename('symbol' = 'Gene symbol', 'probe_id' = 'ID')
    #ids=getIDs('gpl570')
    #head(ids)
    g87371ex = g87371e %>%
        rownames_to_column("probe_id") %>%
        inner_join(GPL570, by = "probe_id") %>%
        dplyr::select(-c("probe_id")) %>%
        dplyr::select(symbol, everything()) %>% 
        mutate(rowMean = rowMeans(.[grep('GSM', names(.))])) %>% 
        arrange(desc(rowMean)) %>% 
        distinct(symbol, .keep_all = T)  %>%
        dplyr::select(-'rowMean')
    # write.csv(g87371ex, "~/Movies/01S100A8/G87371.csv")
    colnames(g87371c) = gsub(':ch1', '', colnames(g87371c))
    s87371 = g87371ex[grep('S100A8', g87371ex$symbol),] %>%
        gather() %>%
        dplyr::rename("symbol" = "key", "S100A8" = "value") %>%
        slice(-1) %>%
        mutate_at("S100A8", as.numeric) %>%
        dplyr::rename("geo_accession" = "symbol") %>%
        inner_join(g87371c, by = "geo_accession") %>% 
        filter(os_time != 'NA') %>% 
        mutate_at(c('os_time', 'cens_os', 'pfs_time', 'cens_pfs'), as.numeric) %>% 
        filter(os_time >= 1) %>% 
        mutate(state = ifelse(S100A8 < median(S100A8), 'low', 'high'))
    fit = survfit(Surv(os_time, cens_os) ~ state, data = s87371) 
    p1 = ggsurvplot(fit,
            pval = TRUE, conf.int = TRUE, 
            risk.table = TRUE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("High S100A8", "Low S100A8"),
            legend.title = "",
            title = "GSE87371", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "hv", 
            palette = color,
            ggtheme = theme_classic())  
    p1$plot + theme(axis.title = element_text(size = 14),
                    plot.title = element_text(size = 16),
                    axis.text = element_text(size = 14, vjust = 0.5),
                    legend.position = c(0.85, 0.95),
                    legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5), color = "transparent")) +
              scale_y_continuous(labels = scales::percent)
    data.survdiff = survdiff(Surv(os_time, cens_os) ~ state, data = s87371)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
    up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[1]+1/data.survdiff$exp[2]))
    low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[1]+1/data.survdiff$exp[2]))
    ggsave("~/Movies/01S100A8/F_G87371OS.png", width = 6, height = 5)
    fit = survfit(Surv(pfs_time, cens_pfs) ~ state, data = s87371) 
    data.survdiff = survdiff(Surv(pfs_time, cens_pfs) ~ state, data = s87371)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
    up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[1]+1/data.survdiff$exp[2]))
    low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[1]+1/data.survdiff$exp[2]))
    p1 = ggsurvplot(fit,
            pval = TRUE, conf.int = TRUE, 
            risk.table = TRUE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("High S100A8", "Low S100A8"),
            legend.title = "",
            title = "GSE87371", 
            xlab = " Time (Month)",
            ylab = 'Progression-Free Survival (%)',
            surv.median.line = "hv", 
            palette = color,
            ggtheme = theme_classic())  
    p1$plot + theme(axis.title = element_text(size = 14),
                    plot.title = element_text(size = 16),
                    axis.text = element_text(size = 14, vjust = 0.5),
                    legend.position = c(0.85, 0.95),
                    legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5), color = "transparent")) +
              scale_y_continuous(labels = scales::percent)
    ggsave("~/Movies/01S100A8/F_G87371PFS.png", width = 6, height = 5)
#### Survival curve of NCICCR_DLBCL
    require(survival)
    require(survminer)
    NCICCR_clinic = fread(file = '~/Movies/01S100A8/NCICCR_DLBCL_clinical.tsv')    
    NCICCR_DLBCL = read.table(file = '~/Movies/01S100A8/NCICCR_DLBCL.txt', header = T)
    NCICCRa = NCICCR_DLBCL %>% 
        as_tibble() %>% 
        dplyr::select(- 'gene_name') %>% 
        column_to_rownames("gene_id")
    NCICCRa = NCICCRa + 1
    NCICCRb = log2(NCICCRa) %>% 
        rownames_to_column('gene_id') %>% 
        inner_join(NCICCR_DLBCL[, 1:2], by = 'gene_id') %>% 
        dplyr::select(-'gene_id') %>% 
        dplyr::select('gene_name', everything()) %>% 
        mutate(rowMean = rowMeans(.[grep('DLBCL', names(.))])) %>% 
        arrange(desc(rowMean)) %>% 
        distinct(gene_name, .keep_all = T)  %>%
        dplyr::select(-'rowMean') %>% 
        na.omit()
    # write.csv(NCICCRb, '~/Movies/01S100A8/NCICCR.csv')
    NCICCR_S100A8 = NCICCRb %>% 
        filter(gene_name == 'S100A8') %>% 
        gather() %>% 
        dplyr::rename('case_submitter_id' = 'key', 'S100A8' = 'value') %>% 
        slice(-1) %>% 
        inner_join(NCICCR_clinic, by = 'case_submitter_id') %>% 
        filter(days_to_last_follow_up != "'--") %>% 
        mutate(os = ifelse(vital_status == 'Dead', 1, 0)) %>% 
        mutate_at(c('days_to_last_follow_up', 'S100A8'), as.numeric) %>% 
        filter(days_to_last_follow_up >= 30) %>% 
        mutate(osmonth = days_to_last_follow_up/30) %>% 
        mutate(state = ifelse(S100A8 < median(S100A8), 'low', 'high'))
    fit = survfit(Surv(osmonth, os) ~ state, data = NCICCR_S100A8) 
    data.survdiff = survdiff(Surv(osmonth, os) ~ state, data = NCICCR_S100A8)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
    up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[1]+1/data.survdiff$exp[2]))
    low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[1]+1/data.survdiff$exp[2]))
    p1 = ggsurvplot(fit,
            pval = TRUE, conf.int = TRUE, 
            risk.table = TRUE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("High S100A8", "Low S100A8"),
            legend.title = "",
            title = "NCICCR_DLBCL", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "hv", 
            palette = color,
            ggtheme = theme_classic())  
    p1$plot + theme(axis.title = element_text(size = 14),
                    plot.title = element_text(size = 16),
                    axis.text = element_text(size = 14, vjust = 0.5),
                    legend.position = c(0.85, 0.95),
                    legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5), color = "transparent")) +
              scale_y_continuous(labels = scales::percent)
    ggsave("~/Movies/01S100A8/F_NCICCROS.png", width = 6, height = 5)
#### Survival curve of GSE31312
    G31312 = getGEO(filename = ("~/Movies/01S100A8/GSE31312_series_matrix.txt.gz"), GSEMatrix = T, 
                getGPL = F)
    if (length(G31312) > 1) {
        idx = grep("GPL570", attr(G31312, "names"))} else {
        idx = 1}
    g31312e = data.frame(exprs(G31312))
    g31312c = pData(G31312)
    qx = as.numeric(quantile(g31312e, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    LogC = (qx[5] > 100) ||
        (qx[6] - qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) {
        g31312e[which(g31312e <= 0)] = NaN
        g31312e = log2(g31312e)
        print("log2 transform finished")} else {
        print("log2 transform not needed")}
    g31312ex = g31312e %>%
        rownames_to_column('probe_id') %>%
        merge(GPL570, by = 'probe_id') %>%
        dplyr::select(- 'probe_id') %>%
        dplyr::select(symbol, everything()) %>% 
        mutate(rowMean = rowMeans(.[grep('GSM', names(.))])) %>% 
        arrange(desc(rowMean)) %>% 
        distinct(symbol, .keep_all = T)  %>%
        dplyr::select(-'rowMean')
    g31312cx = read.csv("~/Movies/01S100A8/GSE31312_clinical_data.csv", header = T) %>% 
        inner_join(g31312c, by = 'title') %>% 
        select('geo_accession', OSmonth = 'OS', 
            OS = 'OScensor', PFSmonth = 'PFS', PFS = 'PFScensor',
            'LDH', IPI = 'IPI.score', Bulky = 'Bulky.cod', 
            AAS = 'AAS.code', 
            'Age', 'Gender', COO = 'gene expression profiling subgroup:ch1',
            extran = 'N.extran', type = 'X3.markers.algoritGEP', 
            ECOG = 'Ecog.code', 'CD10', 'Bcl6', 'FOXP1', 
            'GCET', 'MUM1', 'Response') %>% 
        mutate_at(c('OS', 'PFS', 'Age', 'extran', 'CD10', 'Bcl6', 
        'FOXP1', 'GCET', 'MUM1'), as.numeric) %>% 
        mutate(clinic = case_when(Response == 'SD' ~ '0',
                                  Response == 'PD' ~ '0',
                                  Response == 'CR' ~ '1',
                                  Response == 'PR' ~ '1')) %>%
        mutate(year = ifelse(Age <= 60, 'Young', 'Old')) %>% 
        mutate(treat = 'RCHOP')
    s31312 = g31312ex[grep('S100A8', g31312ex$symbol), ] %>% 
        gather() %>% 
        dplyr::rename('symbol' = 'key', 'S100A8' = 'value') %>% 
        slice(- 1 ) %>% 
        mutate_at('S100A8', as.numeric) %>% 
        dplyr::rename('geo_accession' = 'symbol') %>% 
        inner_join(g31312cx, by = 'geo_accession') %>% 
        mutate(state = ifelse(S100A8 < median(S100A8), 'low', 'high'))
    fit = survfit(Surv(OSmonth, OS) ~ state, data = s31312) 
    data.survdiff = survdiff(Surv(OSmonth, OS) ~ state, data = s31312)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
    up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[1]+1/data.survdiff$exp[2]))
    low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[1]+1/data.survdiff$exp[2]))
    p1 = ggsurvplot(fit,
            pval = TRUE, conf.int = TRUE, 
            risk.table = TRUE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("High S100A8", "Low S100A8"),
            legend.title = "",
            title = "GSE31312", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "hv", 
            palette = color,
            ggtheme = theme_classic())  
    p1$plot + theme(axis.title = element_text(size = 14),
                    plot.title = element_text(size = 16),
                    axis.text = element_text(size = 14, vjust = 0.5),
                    legend.position = c(0.85, 0.95),
                    legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5), color = "transparent")) +
              scale_y_continuous(labels = scales::percent)
    ggsave("~/Movies/01S100A8/F_G31312OS.png", width = 6, height = 5)
    data.survdiff = survdiff(Surv(PFSmonth, PFS) ~ state, data = s31312)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
    up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[1]+1/data.survdiff$exp[2]))
    low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[1]+1/data.survdiff$exp[2]))
    fit = survfit(Surv(PFSmonth, PFS) ~ state, data = s31312) 
    p1 = ggsurvplot(fit,
            pval = TRUE, conf.int = TRUE, 
            risk.table = TRUE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("High S100A8", "Low S100A8"),
            legend.title = "",
            title = "GSE31312", 
            xlab = " Time (Month)",
            ylab = 'Progression-Free Survival (%)',
            surv.median.line = "hv", 
            palette = color,
            ggtheme = theme_classic())  
    p1$plot + theme(axis.title = element_text(size = 14),
                    plot.title = element_text(size = 16),
                    axis.text = element_text(size = 14, vjust = 0.5),
                    legend.position = c(0.85, 0.95),
                    legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5), color = "transparent")) +
              scale_y_continuous(labels = scales::percent)
    ggsave("~/Movies/01S100A8/F_G31312PFS.png", width = 6, height = 5)
## Figure 2 Biological function of DEGs
#### Volcano plot for GSE87371 and NCICCR_DLBCL.
    require(limma)
    str(s87371)
    g = s87371$state
    design = model.matrix(~ -1 + g)
    colnames(design) = c("high", "low")
    contrast.matrix = makeContrasts(high - low, levels = design)
    d87371 = g87371e %>% 
        rownames_to_column('probe_id') %>%
        merge(GPL570, by = 'probe_id') %>%
        dplyr::select(- 'probe_id') %>%
        dplyr::select(symbol, everything()) %>% 
        mutate(rowMean = rowMeans(.[grep('GSM', names(.))])) %>% 
        arrange(desc(rowMean)) %>% 
        distinct(symbol, .keep_all = T)  %>%
        dplyr::select(-'rowMean') %>% 
        column_to_rownames('symbol') %>% 
        select(s87371$geo_accession) %>% 
        lmFit(., design) %>% 
        contrasts.fit(., contrast.matrix) %>% 
        eBayes(.) %>% 
        topTable(., n = Inf, adjust = "fdr") %>% 
        na.omit(.) %>% 
        mutate(group = ifelse(adj.P.Val <= 0.05 & logFC > 1, 'Up', ifelse(adj.P.Val <= 0.05 & logFC < -1, 'Down', 'Not'))) %>% 
        tibble::rownames_to_column("symbol") %>% 
        mutate(log10P = -log10(adj.P.Val)) 
    count(d87371$group)
    fit = d87371 %>% 
        ggscatter("logFC", "log10P",
        combine = F, merge = T, color = "group", shape = 20, size = 2,
        point = TRUE, font.label = 20, palette = c("#8ECFC9", "grey", "#FA7F6F")) +
        geom_vline(xintercept = c(-1, 1), colour = "black", linetype = "dashed") +
        geom_hline(yintercept = -log10(0.05), colour = "black", linetype = "dashed") +
        labs(title = "GSE87371", x = 'log2FC') + 
        theme_classic() + 
        #coord_flip() +
        #ggrepel::geom_label_repel(aes(label = symbol), data = filter(d87371, symbol == "S100A8")) +
        theme(legend.position=c(.9, .88),
        legend.title = element_blank(),
        axis.text = element_text(size = 14, vjust = 0.5),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16))
    ggsave("~/Movies/01S100A8/F_G87371Volano.png", width = 4, height = 4)
    g = NCICCR_S100A8$state
    design = model.matrix(~ -1 + g)
    colnames(design) = c("high", "low")
    contrast.matrix = makeContrasts(high - low, levels = design)
    dNCICCR = NCICCRa %>% 
        as_tibble() %>% 
        select(- 'gene_id') %>% 
        mutate(rowMean = rowMeans(.[grep('DLBCL', names(.))])) %>% 
        arrange(desc(rowMean)) %>% 
        distinct(gene_name, .keep_all = T)  %>%
        dplyr::select(-'rowMean') %>% 
        column_to_rownames('gene_name') %>% 
        select(NCICCR_S100A8$case_submitter_id) %>% 
        lmFit(., design) %>% 
        contrasts.fit(., contrast.matrix) %>% 
        eBayes(.) %>% 
        topTable(., n = Inf, adjust = "fdr") %>% 
        na.omit(.) %>% 
        mutate(group = ifelse(adj.P.Val <= 0.05 & logFC > 1, 'Up', ifelse(adj.P.Val <= 0.05 & logFC < -1, 'Down', 'Not'))) %>% 
        tibble::rownames_to_column("symbol") %>% 
        mutate(log10P = -log10(adj.P.Val)) 
    fit = dNCICCR %>% 
        ggscatter("logFC", "log10P",
        combine = F, merge = T, color = "group", shape = 20, size = 2,
        point = TRUE, font.label = 20, palette = c("#8ECFC9", "grey", "#FA7F6F")) +
        geom_vline(xintercept = c(-1, 1), colour = "black", linetype = "dashed") +
        geom_hline(yintercept = -log10(0.05), colour = "black", linetype = "dashed") +
        labs(title = "NCICCR_DLBCL", x = 'log2FC') + 
        theme_classic() + 
        #coord_flip() +
        #ggrepel::geom_label_repel(aes(label = symbol), data = filter(d87371, symbol == "S100A8")) +
        theme(legend.position=c(.9, .88),
        legend.title = element_blank(),
        axis.text = element_text(size = 14, vjust = 0.5),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16))
    ggsave("~/Movies/01S100A8/F_NCICCRVolano.png", width = 4, height = 4)
#### Intersect Genes
    require(ggvenn)
    d87371a = d87371 %>% 
        filter(group != 'Not') %>% 
        select(symbol)
    dNCICCRa = dNCICCR %>% 
        filter(group != 'Not') %>% 
        select(symbol)
    d87371a = d87371 %>% 
        filter(group == 'Up') %>% 
        select(symbol)
    dNCICCRa = dNCICCR %>% 
        filter(group == 'Up') %>% 
        select(symbol)
    diffgenes1 = intersect(d87371a$symbol, dNCICCRa$symbol)
    d87371a = d87371 %>% 
        filter(group == 'Down') %>% 
        select(symbol)
    dNCICCRa = dNCICCR %>% 
        filter(group == 'Down') %>% 
        select(symbol)
    diffgenes2 = intersect(d87371a$symbol, dNCICCRa$symbol)
    diffgenes = intersect(d87371a$symbol, dNCICCRa$symbol)
    datalist = list("GSE87371" = d87371a$symbol, "NCICCR_DLBCL" = dNCICCRa$symbol)
    ggvenn(datalist, fill_color = color,
        fill_alpha = .7,
        stroke_linetype = "solid",
        set_name_size = 10,
        text_size = 9)  
    ggsave("~/Movies/01S100A8/F_Venn1.png", width = 10, height = 10)
#### GO and KEGG plot
    go1 = read.table(file = "~/Movies/01S100A8/GOBP.txt",header = T,sep = "\t") %>% 
        arrange(., PValue)
    go1 = go1[1:10,]
    go2 = read.table(file = "~/Movies/01S100A8/GOCC.txt", header = T, sep = "\t") %>% 
       arrange(PValue)
    go2 = go2[1:10,]
    go3 = read.table(file = "~/Movies/01S100A8/GOMF.txt", header = T, sep = "\t") %>% 
       arrange(PValue)
    go3 = go3[1:10,]
    go = rbind(go1, go2, go3) %>% 
        as_tibble() %>% 
        mutate(goTerm = c(rep('Biological Process', 10), rep('Cellular Component', 10), rep('Molecular Function', 10))) %>% 
        mutate(go_term_order = factor(as.integer(rownames(.)), labels = Term)) %>% 
        ggplot(aes(x = go_term_order, y = Count, fill = goTerm)) +
        geom_bar(stat='identity') +
        labs(title = "The Most Enriched GO Terms", x = "", y = "The Count of Enrichment", fill = 'GO Terms') + 
        scale_fill_manual(values = color) + 
        theme_classic() + coord_flip() +
        guides(fill=guide_legend(reverse=TRUE)) + 
        theme(legend.position = c(0.9, 0.9),
              legend.title = element_blank(),
              axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 1),
              axis.title = element_text(size = 16),
              plot.title = element_text(size = 16))
    ggsave('~/Movies/01S100A8/F_GO.png', width = 10, height = 6)
    kegg = read.table(file = "~/Movies/01S100A8/KEGG.txt",header = T,sep = "\t") %>% 
        arrange(., PValue) %>%  
        as_tibble() %>% 
        slice(c(1:30)) %>% 
        ggplot(aes(x = Count, y = reorder(Term, - PValue), fill=-log(PValue))) +
        geom_bar(stat='identity') +
        labs(title = "The Most Enriched KEGG Terms", x = "The Count of Enrichment", y = "", fill = "-log(P-value)") + 
        scale_fill_gradient(expression(-log["10"](Pvalue)),low = "#8ECFC9", high = "#FA7F6F") + 
        theme_classic() + 
        #guides(fill=guide_legend(reverse=TRUE)) + 
        theme(legend.position = c(0.85, 0.2),
              axis.text.x = element_text(size = 10, angle = 0, hjust = 1, vjust = 1),
              axis.title = element_text(size = 16),
              plot.title = element_text(size = 16))
    ggsave('~/Movies/01S100A8/F_KEGG.png', width = 8, height = 6)
#### GSEA analysis through WEB-based GEne SeT AnaLysis Toolkit
    d87371g = data.frame(symbol = diffgenes) %>% 
        inner_join(d87371, by = 'symbol') %>% 
        arrange(desc(logFC)) %>% 
        select(c('symbol', 'logFC'))
    write_csv(d87371g, file = '~/Movies/01S100A8/gsea87371.csv')
    fit = read.table(file = '~/Movies/01S100A8/GSEA87371.txt', sep = '\t', header = T) %>% 
        arrange(desc(abs(NES))) %>% 
        group_by(sign(NES))
    fit = fit[1:20,] %>% 
        ggplot(aes(NES, fct_reorder(Description, NES), fill=sign(NES)), showCategory=20) + 
            geom_bar(stat='identity') + 
            scale_fill_gradient(low = "#8ECFC9", high = "#FA7F6F") + 
            theme_classic() + 
            labs(title = "GSE87371", x = "Normalized Enrichment Score", y = "") +
            theme(legend.position= 'none',
                legend.title = element_blank(),
                axis.text = element_text(size = 10, vjust = 0.5),
                axis.title = element_text(size = 14),
                plot.title = element_text(size = 16))
    ggsave('~/Movies/01S100A8/GSEA87371.png', width = 8, height = 6)
    dNCICCRg = data.frame(symbol = diffgenes) %>% 
        inner_join(dNCICCR, by = 'symbol') %>% 
        arrange(desc(logFC)) %>% 
        select(c('symbol', 'logFC'))
    write_csv(dNCICCRg, file = '~/Movies/01S100A8/gseaNCICCR.csv')
    fit = read.table(file = '~/Movies/01S100A8/GSEANCICCR.txt', sep = '\t', header = T) %>% 
        arrange(desc(abs(NES))) %>% 
        group_by(sign(NES))
    fit = fit[1:20, ] %>% 
        ggplot(aes(NES, fct_reorder(Description, NES), fill=sign(NES)), showCategory=20) + 
            geom_bar(stat='identity') + 
            scale_fill_gradient(low = "#8ECFC9", high = "#FA7F6F") + 
            theme_classic() + 
            labs(title = "NCICCR_DLBCL", x = "Normalized Enrichment Score", y = "") +
            theme(legend.position= 'none',
                legend.title = element_blank(),
                axis.text = element_text(size = 10, vjust = 0.5),
                axis.title = element_text(size = 14),
                plot.title = element_text(size = 16))
    ggsave('~/Movies/01S100A8/GSEANCICCR.png', width = 8, height = 6)
## Figure S3 Uninvarite and Multivarite Cox analysis
    genes = data.frame(symbol = c('S100A8', 'S100A9', 'IL6', 'IL1B', 'IFNG', 'MMP1', 'CXCL10', 'CEBPB', 'CXCL8','CXCL1', 'CXCL2'))
    x87371 = g87371e %>%
        rownames_to_column("probe_id") %>%
        inner_join(GPL570, by = "probe_id") %>%
        dplyr::select(-c("probe_id")) %>%
        dplyr::select(symbol, everything()) %>% 
        mutate(rowMean = rowMeans(.[grep('GSM', names(.))])) %>% 
        arrange(desc(rowMean)) %>% 
        distinct(symbol, .keep_all = T)  %>%
        dplyr::select(-'rowMean') %>% 
        inner_join(genes, by = 'symbol') %>% 
        t()
    colnames(x87371) = x87371[1,]
    x87371 = x87371[-1, ] %>% 
        as.data.frame() %>% 
        mutate_all(as.numeric) %>% 
        rownames_to_column('geo_accession') %>% 
        inner_join(g87371c, by = "geo_accession") %>% 
        filter(os_time != 'NA') %>% 
        mutate_at(c('os_time', 'cens_os', 'pfs_time', 'cens_pfs'), as.numeric) %>% 
        filter(os_time >= 1)
    abc = names(x87371)[2:12]
    x = data.frame()
    for(i in abc) {
        expr = x87371[, i]
        cox = coxph(Surv(os_time, cens_os) ~ expr, x87371)
        coxsum = summary(cox)
        x = rbind(x, cbind(coxi = i, 
            HR = round(coxsum$coefficients[, 'exp(coef)'], 3), 
            lower = round(coxsum$conf.int[, 3], 3), 
            upper = round(coxsum$conf.int[, 4], 3),
            'HR(95%CI)' = paste(round(coxsum$coefficients[, 'exp(coef)'], 3), '(', round(coxsum$conf.int[, 3], 3), '-', round(coxsum$conf.int[, 4], 3), ')', sep = ''),
            pvalue = round(coxsum$coefficients[, 'Pr(>|z|)'], 3),
            z = round(coxsum$coefficients[, "z"], 3)))}
    fit = x %>% 
        as_tibble() %>% 
        mutate_at(c('HR','lower','upper'), as.numeric) %>% 
        mutate(log_HR = log(HR)) %>% 
        mutate(log_lower = log(lower)) %>% 
        mutate(log_upper = log(upper)) %>% 
        mutate(log_x = paste(round(log_HR, 2),'(',round(log_lower, 2),',', round(log_upper, 2),')', sep = '')) %>% 
        mutate(iterm = factor(c('S100A8', 'S100A9', 'IL6', 'IL1B', 'IFNG', 'MMP1', 'CXCL10', 'CEBPB', 'CXCL8','CXCL1', 'CXCL2'), levels = c('S100A8', 'S100A9', 'IL6', 'IL1B', 'IFNG', 'MMP1', 'CXCL10', 'CEBPB', 'CXCL8','CXCL1', 'CXCL2')))  %>% 
        ggplot(aes(y = fct_rev(iterm))) + 
            geom_point(aes(x=log_HR), shape=16, size=5, color = '#82B0D2') +
            geom_errorbar(aes(xmin=log_lower, xmax=log_upper), width = 0.2) +
            geom_vline(xintercept = 0, linetype = "dashed") +
            labs(x="Log Hazard Ratio", y="", title = 'GSE87371') +
            #coord_cartesian(ylim=c(1,11), xlim=c(-2, 1.5)) +
            theme_classic() + 
            theme(axis.text.y = element_text(size = 12),
                axis.title.y= element_blank(),
                axis.ticks.y= element_blank(),
                axis.line.y = element_blank())
    ggsave('~/Movies/01S100A8/forrest87371.png', width = 5, height = 4)
    xNCICCR = NCICCRa %>%
        select(-1) %>% 
        dplyr::rename('symbol' = 'gene_name') %>%
        inner_join(genes, by = 'symbol') %>% 
        select('symbol', everything()) %>% 
        t()
    colnames(xNCICCR) = xNCICCR[1,]
    xNCICCR = xNCICCR[-1, ] %>% 
        as.data.frame() %>% 
        mutate_all(as.numeric) %>% 
        rownames_to_column('case_submitter_id') %>% 
        inner_join(NCICCR_clinic, by = 'case_submitter_id') %>% 
        filter(days_to_last_follow_up != "'--") %>% 
        mutate(os = ifelse(vital_status == 'Dead', 1, 0)) %>% 
        mutate_at(c('days_to_last_follow_up'), as.numeric) %>% 
        filter(days_to_last_follow_up >= 30) %>% 
        mutate(osmonth = days_to_last_follow_up/30)
    abc = names(xNCICCR)[2:12]
    x = data.frame()
    for(i in abc) {
        expr = xNCICCR[, i]
        cox = coxph(Surv(osmonth, os) ~ expr, xNCICCR)
        coxsum = summary(cox)
        x = rbind(x, cbind(coxi = i, 
            HR = round(coxsum$coefficients[, 'exp(coef)'], 3), 
            lower = round(coxsum$conf.int[, 3], 3), 
            upper = round(coxsum$conf.int[, 4], 3),
            'HR(95%CI)' = paste(round(coxsum$coefficients[, 'exp(coef)'], 3), '(', round(coxsum$conf.int[, 3], 3), '-', round(coxsum$conf.int[, 4], 3), ')', sep = ''),
            pvalue = round(coxsum$coefficients[, 'Pr(>|z|)'], 3),
            z = round(coxsum$coefficients[, "z"], 3)))}
    fit = x %>% 
        as_tibble() %>% 
        mutate_at(c('HR','lower','upper'), as.numeric) %>% 
        mutate(log_HR = log(HR)) %>% 
        mutate(log_lower = log(lower)) %>% 
        mutate(log_upper = log(upper)) %>% 
        mutate(log_x = paste(round(log_HR, 2),'(',round(log_lower, 2),',', round(log_upper, 2),')', sep = '')) %>% 
        mutate(iterm = factor(c('S100A8', 'S100A9', 'IL6', 'IL1B', 'IFNG', 'MMP1', 'CXCL10', 'CEBPB', 'CXCL8','CXCL1', 'CXCL2'), levels = c('S100A8', 'S100A9', 'IL6', 'IL1B', 'IFNG', 'MMP1', 'CXCL10', 'CEBPB', 'CXCL8','CXCL1', 'CXCL2')))  %>% 
        ggplot(aes(y = fct_rev(iterm))) + 
            geom_point(aes(x=log_HR), shape=16, size=5, color = '#82B0D2') +
            geom_errorbar(aes(xmin=log_lower, xmax=log_upper), width = 0.2) +
            geom_vline(xintercept = 0, linetype = "dashed") +
            labs(x="Log Hazard Ratio", y="", title = 'NCICCR_DLBCL') +
            #coord_cartesian(ylim=c(1,11), xlim=c(-2, 1.5)) +
            theme_classic() + 
            theme(axis.text.y = element_text(size = 12),
                axis.title.y= element_blank(),
                axis.ticks.y= element_blank(),
                axis.line.y = element_blank())
    ggsave('~/Movies/01S100A8/forrestNCICCR.png', width = 5, height = 4)
## Figure 4 S100A8 in single cell analysis
    rm(list=ls())
    options(stringsAsFactors = F)
    options(Seurat.object.assay.version = 'v5')
    require(Seurat)
    require(clustree)
    require(cowplot)
    require(patchwork)
    G182434r = fread("~/Movies/01S100A8/GSE182434_raw_count_matrix.txt.gz", data.table = F) %>% 
        column_to_rownames('Gene')
    G182434c = fread("~/Movies/01S100A8/GSE182434_cell_annotation.txt.gz", data.table = F)
    dim(table(G182434c$Patient, G182434c$CellType))
    sample = c("DLBCL002", "DLBCL007", "DLBCL008", "DLBCL111")
    G182434c = G182434c %>% 
        filter(Patient %in% sample) %>% 
        filter(TumorNormal == 'Tumor')
    dim(G182434c)
    G182434e = G182434r %>% 
        dplyr::select(G182434c$ID)
    seu = CreateSeuratObject(counts = Matrix::Matrix(as.matrix(G182434e)), project = "seurat", min.cells = 3, min.features = 200)
    seu[["percent.mt"]] = PercentageFeatureSet(seu, pattern = "^MT-")  ## QC mitochondria
    p1 = seu@meta.data %>% 
        ggplot(aes(fill = orig.ident, y = nFeature_RNA, x = orig.ident)) +
        geom_violin(alpha = 0.8, color = "#FA7F6F")  +
        geom_jitter(aes(y = nFeature_RNA), size = 0.5, shape = 20,
            stroke = 0, show.legend = FALSE, width = 0.1, color = '#6b798e') +
        labs(title = "nFeature_RNA", y = "", x = "") +
        theme_classic() +
        theme(legend.position = "none",
            axis.text.x = element_text(size = 14),
            axis.title = element_text(size = 14),
            plot.title = element_text(size = 16))
    p2 = seu@meta.data %>% 
        ggplot(aes(fill = orig.ident, y = nCount_RNA, x = orig.ident)) +
        geom_violin(alpha = 0.8, color = "#FA7F6F")  +
        geom_jitter(aes(y = nCount_RNA), size = 0.5, shape = 20,
            stroke = 0, show.legend = FALSE, width = 0.1, color = '#6b798e') +
        labs(title = "nCount_RNA", y = "", x = "") +
        theme_classic() +
        theme(legend.position = "none",
            axis.text.x = element_text(size = 14),
            axis.title = element_text(size = 14),
            plot.title = element_text(size = 16))
    p3 = seu@meta.data %>% 
        ggplot(aes(fill = orig.ident, y = percent.mt, x = orig.ident)) +
        geom_violin(alpha = 0.8, color = "#FA7F6F")  +
        geom_jitter(aes(y = percent.mt), size = 0.5, shape = 20,
            stroke = 0, show.legend = FALSE, width = 0.1, color = '#6b798e') +
        labs(title = "percent.mt", y = "", x = "") +
        theme_classic() +
        theme(legend.position = "none",
            axis.text.x = element_text(size = 14),
            axis.title = element_text(size = 14),
            plot.title = element_text(size = 16))
    p = p1 + p2 + p3
    ggsave('~/Movies/01S100A8/Fig_scRNAqc.jpg', width = 8, height = 8)
    seu = subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
    seu = NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
    seu = FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
    scegenes = VariableFeatures(seu)
    grep('S100A8', scegenes)
    p1 = VariableFeaturePlot(seu, cols = c("#6b798e", "#FA7F6F"), pt.size = 1,
            log = NULL, selection.method = NULL, assay = NULL, raster = NULL, 
            raster.dpi = c(512, 512)) %>% 
            LabelPoints(points = 'S100A8', repel = TRUE) + 
            theme(legend.position = c(0.05, 0.95),
                  axis.title.x = element_blank())
    seu = ScaleData(seu, features = rownames(seu))
    seu = RunPCA(seu, npcs = 50, verbose = FALSE)
    DimPlot(seu, reduction = "pca")
    ElbowPlot(seu)
    seu = FindNeighbors(seu, dims = 1:15) %>% 
        FindClusters(resolution = 0.8) %>% 
        RunUMAP(dims = 1:15) %>% 
        RunTSNE(dims = 1:15)
    DimPlot(seu, reduction = "tsne", label = TRUE)
    top10 = FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
        group_by(cluster) %>%
        top_n(n = 10, wt = avg_log2FC) %>% 
        as.data.frame()
    grep('S100A8', top10$gene)
    require(harmony)     # harmony remove batch
    head(G182434c)
    seu@meta.data = seu@meta.data %>% 
        rownames_to_column('ID') %>% 
        inner_join(G182434c, by = 'ID') %>% 
        column_to_rownames('ID')
    seu = RunHarmony(seu, "Patient")
    seu = seu %>% 
        RunUMAP(reduction = "harmony", dims = 1:15) %>% 
        FindNeighbors(reduction = "harmony", dims = 1:15) %>% 
        FindClusters(resolution = 0.5) %>% 
        identity()
    seu = seu %>% 
        RunTSNE(reduction = "harmony", dims = 1:15) %>% 
        RunUMAP(reduction = "harmony", dims = 1:15)
    p2 = DimPlot(seu, reduction = "tsne", group.by = "ident", pt.size = 0.25, cols = color)
    require(SingleR)
    hpca.se = celldex::HumanPrimaryCellAtlasData()
    seu_for_SingleR = GetAssayData(seu, slot="data")
    seu.hesc = SingleR(test = seu_for_SingleR, ref = hpca.se, labels = hpca.se$label.main)
    table(seu.hesc$labels, seu@meta.data$seurat_clusters)
    seu@meta.data$labels = seu.hesc$labels
    p3 = DimPlot(seu, group.by = c("CellType"), reduction = "tsne", pt.size = 0.25, cols = color)
    p4 = VlnPlot(seu, features = c("S100A8"), group.by = c("CellType"), cols = color) +
                theme(axis.title.x = element_blank())
    p6 = DimPlot(seu, group.by = c("Patient"), reduction = "tsne", pt.size = 0.25, cols = color)
    p5 = FeaturePlot(seu, features = c("S100A8"), reduction = "tsne", cols = 
        c("#E7DAD217", "#FF0000"))
    py = (p3 + p5) / p4
    p = p1 + py + plot_layout(widths = c(1, 3)) + plot_annotation(tag_levels = 'A') & 
        theme(plot.tag = element_text(size = 24))
    ggsave('~/Movies/01S100A8/Fig_seu1s100a8.jpg', width = 15, height = 8)
## Figure 5 Immune infiltration analysis
    require(tidyestimate)
    require(immunedeconv)
    require(patchwork)
    require(viridis)
    require(aplot)
    require(readxl)
    require(stringr)
    T = read_excel('~/Movies/01S100A8/gene_table.xlsx', sheet = 'T') %>% 
        filter(cancer == 'DLBC (n=48)')
    M = read_excel('~/Movies/01S100A8/gene_table.xlsx', sheet = 'M') %>% 
        filter(cancer == 'DLBC (n=48)')
    B = read_excel('~/Movies/01S100A8/gene_table.xlsx', sheet = 'B') %>% 
        filter(cancer == 'DLBC (n=48)')
    immun = rbind(T, M, B) %>% 
        as_tibble()
    group = str_split_fixed(immun$infiltrates, "_", 2)
    colnames(group) = c('cell', 'data')
    immun = cbind(immun, group) %>% 
        as_tibble() %>% 
        arrange(data)
    immun$data = factor(immun$data, levels = c(unique(immun$data)))
    immun$infiltrates = factor(immun$infiltrates, levels = c(unique(immun$infiltrates)))
    immun$guide = c(rep('CIBERSORT', 8), rep('CIBERSORT-ABS', 8), rep('EPIC', 3), rep('MCPcounter', 4), rep('quanTIseq', 5), rep('TIDE', 1), rep('TIMER', 3), rep('xCell', 13))
    immun$guide = factor(immun$guide, levels = c(unique(immun$guide)))
    heat0 = immun %>% 
        mutate(group = c("S100A8")) %>% 
        mutate_at('rho', as.numeric) %>% 
        ggplot(aes(x = infiltrates, y = group, fill = guide)) +
        geom_tile(aes(fill = guide)) +
        scale_fill_manual(values = color) +
        scale_y_discrete(position = "left") +
        theme_bw() +
        labs(x = "", y = "Algorithm", fill = 'Algorithm') +
        guides(fill = guide_legend(nrow = 1, byrow=TRUE)) +
        theme(legend.position = "top",
              panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.ticks = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_blank(),
              axis.title = element_blank())
    heat1 = immun %>% 
        mutate(group = c("S100A8")) %>% 
        mutate_at('rho', as.numeric) %>% 
        ggplot(aes(x = infiltrates, y = group, fill = rho)) +
        geom_tile(aes(fill = rho)) +
        scale_fill_gradient2(high = '#FA7F6F', low = '#82B0D2', mid = '#CFEAF1')+
        theme_minimal() +
        guides(fill = guide_legend("rho")) + 
        theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
              panel.grid = element_blank(),
              axis.ticks = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 14),
              axis.text.x = element_text(angle = 90, size = 12, hjust = 1),
              axis.title = element_blank(),
              legend.position = "right")
    require(patchwork)
    heat0 / heat1
    ggsave('~/Movies/01S100A8/F_immunenca.jpg', width = 20, height = 5) 
    fit = estimate_score(g87371ex, is_affymetrix = F) %>% 
        dplyr::rename('geo_accession' = 'sample') %>% 
        inner_join(s87371, by = 'geo_accession') %>% 
        arrange(state) %>% 
        mutate(status = ifelse(estimate < median(estimate), 'low', 'high')) %>% 
        select(c('stromal', 'immune','estimate', 'state', 'geo_accession', 'status', 'os_time', 'pfs_time', 'cens_os', 'cens_pfs'))
    fit1 = survfit(Surv(os_time, cens_os) ~ status, data = fit) 
    p1 = ggsurvplot(fit1, data = fit, 
            pval = TRUE, conf.int = TRUE, 
            risk.table = TRUE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("High S100A8", "Low S100A8"),
            legend.title = "",
            title = "GSE87371", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "hv", 
            palette = color,
            ggtheme = theme_classic())  
    p1$plot + theme(axis.title = element_text(size = 20),
                    plot.title = element_text(size = 20),
                    axis.text = element_text(size = 20, vjust = 0.5),
                    legend.position = c(0.85, 0.95),
                    legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5), color = "transparent")) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Movies/01S100A8/F_immunesur1.jpg', width = 10, height = 5)
    fit$geo_accession = factor(fit$geo_accession, levels = c(unique(fit$geo_accession)))
    group = c('stromal', 'immune', 'estimate')
    z = data.frame()
    for(i in group){
        y = fit %>% select(c(i, 'state'))
        y1 = wilcox.test(y[,1] ~ y[,2])
        y2 = data.frame(pvalue = y1$p.value, group = i)
        z = rbind(y2, z)}
    ze = z %>% 
        mutate(sign = case_when(pvalue >= 0.05 ~ 'ns',
                                pvalue < 0.05 & pvalue >= 0.01 ~ '*',
                                pvalue < 0.01 & pvalue >= 0.001 ~ '**',
                                pvalue < 0.001 ~ '***',)) %>% 
        select(group, everything()) %>% 
        dplyr::rename('name' = 'group') %>% 
        mutate(names = paste0(.$name, '(', .$sign, ')'))
    p0 = fit %>% as_tibble() %>%
        pivot_longer(cols = colnames(.)[1:3],
               values_to = 'value') %>%
               inner_join(ze, by = 'name') %>% 
        filter(name == 'estimate') %>% 
        ggplot(aes(x = geo_accession, y = names, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_viridis(alpha = 0.5, option = 'F', begin = '1', end = '0') +
        scale_y_discrete(position = "left") +
        theme_bw() +
        labs(x = "", y = "ESTIMATE") +
        guides(fill=guide_legend("ESTIMATE")) + 
        theme(legend.position = "right",
              panel.grid = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.ticks = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 14),
              axis.text.x = element_blank(),
              axis.title = element_blank())
    p1 = fit %>% 
        ggplot(aes(fill = state, y = stromal, x = state)) +
        geom_violin(alpha = 0.8) +
        geom_jitter(aes(y = stromal), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = F, parse = T,
            comparisons = list(c('low', 'high')), size = 0.5, 
            textsize = 3, vjust = 0, step_increase = 0.15, position = "identity") +
        labs(title = "GSE87371", y = "Stromal Score", x = "") +
        scale_x_discrete("", labels = c("low" = "low S100A8", 'high' = 'high S100A8')) +
        #scale_y_continuous(limits = c(3, 10)) +
        scale_fill_manual(values = color) + 
        theme_classic() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 12))
    p2 = fit %>% 
        ggplot(aes(fill = state, y = immune, x = state)) +
        geom_violin(alpha = 0.8) +
        geom_jitter(aes(y = immune), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = F, parse = T,
            comparisons = list(c('low', 'high')), size = 0.5, 
            textsize = 3, vjust = 0, step_increase = 0.15, position = "identity") +
        labs(title = "GSE87371", y = "Immune Score", x = "") +
        scale_x_discrete("", labels = c("low" = "low S100A8", 'high' = 'high S100A8')) +
        #scale_y_continuous(limits = c(3, 10)) +
        scale_fill_manual(values = color) + 
        theme_classic() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 12))
    p3 = fit %>% 
        ggplot(aes(fill = state, y = estimate, x = state)) +
        geom_violin(alpha = 0.8) +
        geom_jitter(aes(y = estimate), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = F, parse = T,
            comparisons = list(c('low', 'high')), size = 0.5, 
            textsize = 3, vjust = 0, step_increase = 0.15, position = "identity") +
        labs(title = "GSE87371", y = "Estimate Score", x = "") +
        scale_x_discrete("", labels = c("low" = "low S100A8", 'high' = 'high S100A8')) +
        #scale_y_continuous(limits = c(3, 10)) +
        scale_fill_manual(values = color) + 
        theme_classic() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 12))
    a87371 = read.table("~/Movies/01S100A8/CIBERSORTx_GSE87371.txt", header = T, sep = "\t") %>% 
        dplyr::rename('geo_accession' = 'Mixture') %>% 
        inner_join(s87371, by = 'geo_accession')
    res1 = data.frame(a87371[, c(1:23, 92)]) %>% 
        pivot_longer(cols = colnames(.)[2:23],
               names_to = "cell.type",
               values_to = 'value') %>% 
        arrange(state) 
    res1$geo_accession = factor(res1$geo_accession, levels = c(unique(res1$geo_accession)))
    count(res1$cell.type)
    x1 = res1 %>% 
        filter(cell.type == 'B.cells.memory' | cell.type == 'B.cells.naive' | cell.type == 'Macrophages.M0' | cell.type == 'Macrophages.M1' | cell.type == 'Macrophages.M2' | cell.type == 'Monocytes' | cell.type ==  'T.cells.CD8') %>% 
        dplyr::rename('cell_type' = 'cell.type')
    group = c('B.cells.memory', 'B.cells.naive', 'Macrophages.M0', 'Macrophages.M1', 'Macrophages.M2', 'Monocytes', 'T.cells.CD8')
    z = data.frame()
    for(i in group){
        y = subset(x1, cell_type == i)
        y1 = wilcox.test(y$value ~ y$state)
        y2 = data.frame(pvalue = y1$p.value, group = i)
        z = rbind(y2, z)}
    z1 = z %>% 
        mutate(sign = case_when(pvalue >= 0.05 ~ 'ns',
                                pvalue < 0.05 & pvalue >= 0.01 ~ '*',
                                pvalue < 0.01 & pvalue >= 0.001 ~ '**',
                                pvalue < 0.001 ~ '***',)) %>% 
        select(group, everything()) %>% 
        dplyr::rename('cell_type' = 'group')
    head(z1)
    fit0 = x1 %>% as_tibble() %>% 
        mutate(y = c("group")) %>% 
        ggplot(aes(x = geo_accession, y = y, fill = state)) +
        geom_tile(aes(fill = state)) +
        scale_fill_manual(values = color)+
        scale_y_discrete(position = "left") +
        theme_minimal() +
        xlab(NULL) + ylab(NULL) +
        labs(fill = 'group') + 
        guides(fill = guide_legend("Levels of S100A8")) + 
        theme(panel.border = element_blank(),
              panel.grid = element_blank(),
              axis.ticks = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 14),
              axis.text.x = element_blank(),
              axis.title = element_blank(),
              legend.position = "right")
    fit1 = x1 %>% as_tibble() %>% 
        inner_join(z1, by = 'cell_type') %>% 
        mutate(names = paste0(.$cell_type, '(', .$sign, ')')) %>% 
        ggplot(aes(x = geo_accession, y = names, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_viridis(alpha = 0.5, option = 'D', begin = '1', end = '0') +
        scale_y_discrete(position = "left") +
        theme_bw() +
        labs(x = "", y = "CIBERSORTx") +
        guides(fill=guide_legend("CIBERSORTx")) + 
        theme(legend.position = "right",
              panel.grid = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank())
    m87371 = g87371ex %>% 
        column_to_rownames('symbol') %>% 
        as.matrix()
    x2 = immunedeconv::deconvolute(m87371, "xcell") %>% 
        as.data.frame() %>% 
        pivot_longer(cols = colnames(.)[2:ncol(.)],
               names_to = "geo_accession",
               values_to = 'value') %>% 
        filter(cell_type == 'B cell' | cell_type == 'B cell memory' | cell_type == 'Class-switched memory B cell' | cell_type == 'B cell naive' | cell_type == 'B cell plasma' | cell_type == 'Macrophage' | cell_type == 'Macrophage M1' | cell_type == 'Macrophage M2' | cell_type == 'Monocyte' | cell_type == 'T cell CD8+' | cell_type == 'T cell CD8+ central memory' | cell_type == 'T cell CD8+ effector memory' | cell_type == 'T cell CD8+ naive') %>% 
        inner_join(s87371[,c('geo_accession', 'state')], by = 'geo_accession') %>% 
        arrange(state)
    x2$geo_accession = factor(x2$geo_accession, levels = c(unique(x2$geo_accession)))
    group = c('B cell', 'B cell memory', 'Class-switched memory B cell', 'B cell naive', 'B cell plasma', 'Macrophage', 'Macrophage M1', 'Macrophage M2', 'Monocyte', 'T cell CD8+', 'T cell CD8+ central memory', 'T cell CD8+ effector memory', 'T cell CD8+ naive')
    z = data.frame()
    for(i in group){
        y = subset(x2, cell_type == i)
        y1 = wilcox.test(y$value ~ y$state)
        y2 = data.frame(pvalue = y1$p.value, group = i)
        z = rbind(y2, z)}
    z2 = z %>% 
        mutate(sign = case_when(pvalue >= 0.05 ~ 'ns',
                                pvalue < 0.05 & pvalue >= 0.01 ~ '*',
                                pvalue < 0.01 & pvalue >= 0.001 ~ '**',
                                pvalue < 0.001 ~ '***',)) %>% 
        select(group, everything()) %>% 
        dplyr::rename('cell_type' = 'group')
    fit2 = x2 %>% as_tibble() %>% 
        inner_join(z2, by = 'cell_type') %>% 
        mutate(names = paste0(.$cell_type, '(', .$sign, ')')) %>% 
        ggplot(aes(x = geo_accession, y = names, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        scale_y_discrete(position = "left") +
        theme_bw() +
        labs(x = "", y = "xCell") +
        guides(fill=guide_legend("xCell")) + 
        theme(legend.position = "right",
              panel.grid = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank())
    x3 = immunedeconv::deconvolute(m87371, "quantiseq") %>% 
        as.data.frame() %>% 
        pivot_longer(cols = colnames(.)[2:ncol(.)],
               names_to = "geo_accession",
               values_to = 'value') %>% 
        filter(cell_type == 'B cell' | cell_type == 'Macrophage M1' | cell_type == 'Macrophage M2' | cell_type == 'Monocyte' | cell_type == 'T cell CD8+' ) %>% 
        inner_join(s87371[,c('geo_accession', 'state')], by = 'geo_accession') %>% 
        arrange(state)
    x3$geo_accession = factor(x3$geo_accession, levels = c(unique(x3$geo_accession)))
    group = c('B cell', 'Macrophage M1', 'Macrophage M2', 'Monocyte', 'T cell CD8+')
    z = data.frame()
    for(i in group){
        y = subset(x3, cell_type == i)
        y1 = wilcox.test(y$value ~ y$state)
        y2 = data.frame(pvalue = y1$p.value, group = i)
        z = rbind(y2, z)}
    z3 = z %>% 
        mutate(sign = case_when(pvalue >= 0.05 ~ 'ns',
                                pvalue < 0.05 & pvalue >= 0.01 ~ '*',
                                pvalue < 0.01 & pvalue >= 0.001 ~ '**',
                                pvalue < 0.001 ~ '***',)) %>% 
        select(group, everything()) %>% 
        dplyr::rename('cell_type' = 'group')
    fit3 = x3 %>% as_tibble() %>% 
        inner_join(z3, by = 'cell_type') %>% 
        mutate(names = paste0(.$cell_type, '(', .$sign, ')')) %>% 
        ggplot(aes(x = geo_accession, y = names, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_viridis(alpha = 0.5, option = 'E', begin = '1', end = '0') +
        scale_y_discrete(position = "left") +
        theme_bw() +
        labs(x = "", y = "quanTIseq") +
        guides(fill=guide_legend("quanTIseq")) + 
        theme(legend.position = "right",
              panel.grid = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank())
    require(MCPcounter)
    probesets = read.table('~/Movies/01S100A8/probesets.txt',
                        sep="\t",
                        stringsAsFactors=FALSE,
                        colClasses="character")
    genes = read.table('~/Movies/01S100A8/genes.txt',
                    sep="\t",
                    stringsAsFactors=FALSE,
                    header=TRUE,
                    colClasses="character",
                    check.names=FALSE)
    x4 = MCPcounter::MCPcounter.estimate(m87371, featuresType = "HUGO_symbols", probesets = probesets, genes = genes) %>% 
        as.data.frame() %>% 
        rownames_to_column('cell_type') %>% 
        pivot_longer(cols = colnames(.)[2:ncol(.)],
               names_to = "geo_accession",
               values_to = 'value') %>% 
        filter(cell_type == 'Monocytic lineage' | cell_type == 'B lineage' | cell_type == 'CD8 T cells') %>% 
        inner_join(s87371[,c('geo_accession', 'state')], by = 'geo_accession') %>% 
        arrange(state)
    x4$geo_accession = factor(x4$geo_accession, levels = c(unique(x4$geo_accession)))
    group = c('Monocytic lineage', 'B lineage', 'CD8 T cells')
    z = data.frame()
    for(i in group){
        y = subset(x4, cell_type == i)
        y1 = wilcox.test(y$value ~ y$state)
        y2 = data.frame(pvalue = y1$p.value, group = i)
        z = rbind(y2, z)}
    z4 = z %>% 
        mutate(sign = case_when(pvalue >= 0.05 ~ 'ns',
                                pvalue < 0.05 & pvalue >= 0.01 ~ '*',
                                pvalue < 0.01 & pvalue >= 0.001 ~ '**',
                                pvalue < 0.001 ~ '***',)) %>% 
        select(group, everything()) %>% 
        dplyr::rename('cell_type' = 'group')
    fit4 = x4 %>% as_tibble() %>% 
        inner_join(z4, by = 'cell_type') %>% 
        mutate(names = paste0(.$cell_type, '(', .$sign, ')')) %>% 
        ggplot(aes(x = geo_accession, y = names, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_viridis(alpha = 0.5, option = 'A', begin = '1', end = '0') +
        scale_y_discrete(position = "left") +
        theme_bw() +
        labs(x = "", y = "MCPcounter") +
        guides(fill=guide_legend("MCPcounter")) + 
        theme(legend.position = "right",
              panel.grid = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank())
    fit0 %>% insert_bottom(p0, height = 1) %>% insert_bottom(fit1, height = 7) %>% insert_bottom(fit2, height = 13) %>% insert_bottom(fit3, height = 5) %>% insert_bottom(fit4, height = 3)
    ggsave('~/Movies/01S100A8/F_immune87371.jpg', width = 20, height = 9) 
    nci = estimate_score(NCICCRb, is_affymetrix = F) %>% 
        dplyr::rename('case_submitter_id' = 'sample') %>% 
        inner_join(NCICCR_S100A8, by = 'case_submitter_id') %>% 
        arrange(state) %>% 
        mutate(status = ifelse(estimate < median(estimate), 'low', 'high')) %>% 
        select(c('stromal', 'immune','estimate', 'status', 'state', 'case_submitter_id', 'os', 'osmonth'))
    nci$case_submitter_id = factor(nci$case_submitter_id, levels = c(unique(nci$case_submitter_id)))
    fit1 = survfit(Surv(osmonth, os) ~ status, data = nci) 
    p1 = ggsurvplot(fit1, data = nci, 
            pval = TRUE, conf.int = TRUE, 
            risk.table = TRUE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("High S100A8", "Low S100A8"),
            legend.title = "",
            title = "NCICCR_DLBCL", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "hv", 
            palette = color,
            ggtheme = theme_classic())  
    p1$plot + theme(axis.title = element_text(size = 20),
                    plot.title = element_text(size = 20),
                    axis.text = element_text(size = 20, vjust = 0.5),
                    legend.position = c(0.85, 0.95),
                    legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5), color = "transparent")) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Movies/01S100A8/F_immunesur2.jpg', width = 10, height = 5)
    group = c('stromal', 'immune', 'estimate')
    z = data.frame()
    for(i in group){
        y = nci %>% select(c(i, 'state'))
        y1 = wilcox.test(y[,1] ~ y[,2])
        y2 = data.frame(pvalue = y1$p.value, group = i)
        z = rbind(y2, z)}
    ze = z %>% 
        mutate(sign = case_when(pvalue >= 0.05 ~ 'ns',
                                pvalue < 0.05 & pvalue >= 0.01 ~ '*',
                                pvalue < 0.01 & pvalue >= 0.001 ~ '**',
                                pvalue < 0.001 ~ '***',)) %>% 
        select(group, everything()) %>% 
        dplyr::rename('name' = 'group') %>% 
        mutate(names = paste0(.$name, '(', .$sign, ')'))
    n0 = nci %>% as_tibble() %>%
        pivot_longer(cols = colnames(.)[1:3],
               values_to = 'value') %>%
               inner_join(ze, by = 'name') %>% 
        filter(name == 'estimate') %>% 
        ggplot(aes(x = case_submitter_id, y = names, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_viridis(alpha = 0.5, option = 'F', begin = '1', end = '0') +
        scale_y_discrete(position = "left") +
        theme_bw() +
        labs(x = "", y = "") +
        guides(fill=guide_legend("ESTIMATE")) + 
        theme(legend.position = "right",
              panel.grid = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.ticks = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 14),
              axis.text.x = element_blank(),
              axis.title = element_blank())
    fit1 = nci %>% 
        ggplot(aes(fill = state, y = stromal, x = state)) +
        geom_violin(alpha = 0.8) +
        geom_jitter(aes(y = stromal), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('low', 'high')), size = 0.5, 
            textsize = 5, vjust = 0.4, step_increase = 0.15, position = "identity") +
        labs(title = "NCICCR_DLBCL", y = "Stromal Score", x = "") +
        scale_x_discrete("", labels = c("low" = "low S100A8", 'high' = 'high S100A8')) +
        #scale_y_continuous(limits = c(3, 10)) +
        scale_fill_manual(values = color) + 
        theme_classic() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 12))
    fit2 = nci %>% 
        ggplot(aes(fill = state, y = immune, x = state)) +
        geom_violin(alpha = 0.8) +
        geom_jitter(aes(y = immune), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('low', 'high')), size = 0.5, 
            textsize = 5, vjust = 0.4, step_increase = 0.15, position = "identity") +
        labs(title = "NCICCR_DLBCL", y = "Immune Score", x = "") +
        scale_x_discrete("", labels = c("low" = "low S100A8", 'high' = 'high S100A8')) +
        #scale_y_continuous(limits = c(3, 10)) +
        scale_fill_manual(values = color) + 
        theme_classic() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 12))
    fit3 = nci %>% 
        ggplot(aes(fill = state, y = estimate, x = state)) +
        geom_violin(alpha = 0.8) +
        geom_jitter(aes(y = estimate), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('low', 'high')), size = 0.5, 
            textsize = 5, vjust = 0.4, step_increase = 0.15, position = "identity") +
        labs(title = "NCICCR_DLBCL", y = "Estimate Score", x = "") +
        scale_x_discrete("", labels = c("low" = "low S100A8", 'high' = 'high S100A8')) +
        #scale_y_continuous(limits = c(3, 10)) +
        scale_fill_manual(values = color) + 
        theme_classic() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 12))
    aNCICCR = read.table("~/Movies/01S100A8/CIBERSORTx_NCICCR.txt", header = T, sep = "\t") %>% 
        dplyr::rename('case_submitter_id' = 'Mixture') %>% 
        inner_join(NCICCR_S100A8, by = 'case_submitter_id')
    res2 = data.frame(aNCICCR[,c(1:23, 187)]) %>% 
        pivot_longer(cols = colnames(.)[2:23],
               names_to = "cell.type",
               values_to = 'value') %>% 
               arrange(state) 
    res2$case_submitter_id = factor(res2$case_submitter_id, levels = c(unique(res2$case_submitter_id)))
    b1 = res2 %>% 
        filter(cell.type == 'B.cells.memory' | cell.type == 'B.cells.naive' | cell.type == 'Macrophages.M0' | cell.type == 'Macrophages.M1' | cell.type == 'Macrophages.M2' | cell.type == 'Monocytes' | cell.type ==  'T.cells.CD8') %>% 
        dplyr::rename('cell_type' = 'cell.type')
    group = c('B.cells.memory', 'B.cells.naive', 'Macrophages.M0', 'Macrophages.M1', 'Macrophages.M2', 'Monocytes', 'T.cells.CD8')
    z = data.frame()
    for(i in group){
        y = subset(b1, cell_type == i)
        y1 = wilcox.test(y$value ~ y$state)
        y2 = data.frame(pvalue = y1$p.value, group = i)
        z = rbind(y2, z)}
    z1 = z %>% 
        mutate(sign = case_when(pvalue >= 0.05 ~ 'ns',
                                pvalue < 0.05 & pvalue >= 0.01 ~ '*',
                                pvalue < 0.01 & pvalue >= 0.001 ~ '**',
                                pvalue < 0.001 ~ '***',)) %>% 
        select(group, everything()) %>% 
        dplyr::rename('cell_type' = 'group')
    nx = b1 %>% as_tibble() %>% 
        mutate(y = c("group")) %>% 
        ggplot(aes(x = case_submitter_id, y = y, fill = state)) +
        geom_tile(aes(fill = state)) +
        scale_fill_manual(values = color)+
        scale_y_discrete(position = "left") +
        theme_minimal() +
        xlab(NULL) + ylab(NULL) +
        labs(fill = 'group') + 
        guides(fill = guide_legend("Levels of S100A8")) + 
        theme(panel.border = element_blank(),
              panel.grid = element_blank(),
              axis.ticks = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 14),
              axis.text.x = element_blank(),
              axis.title = element_blank(),
              legend.position = "right")
    n1 = b1 %>% as_tibble() %>% 
        inner_join(z1, by = 'cell_type') %>% 
        mutate(names = paste0(.$cell_type, '(', .$sign, ')')) %>% 
        ggplot(aes(x = case_submitter_id, y = names, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_viridis(alpha = 0.5, option = 'D', begin = '1', end = '0') +
        scale_y_discrete(position = "left") +
        theme_bw() +
        labs(x = "", y = "CIBERSORTx") +
        guides(fill=guide_legend("CIBERSORTx")) + 
        theme(legend.position = "right",
              panel.grid = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank())
    mnci = NCICCRb %>% 
        column_to_rownames('gene_name') %>% 
        as.matrix()
    n2 = immunedeconv::deconvolute(mnci, "xcell") %>% 
        as.data.frame() %>% 
        pivot_longer(cols = colnames(.)[2:ncol(.)],
               names_to = "case_submitter_id",
               values_to = 'value') %>% 
        filter(cell_type == 'B cell' | cell_type == 'B cell memory' | cell_type == 'Class-switched memory B cell' | cell_type == 'B cell naive' | cell_type == 'B cell plasma' | cell_type == 'Macrophage' | cell_type == 'Macrophage M1' | cell_type == 'Macrophage M2' | cell_type == 'Monocyte' | cell_type == 'T cell CD8+' | cell_type == 'T cell CD8+ central memory' | cell_type == 'T cell CD8+ effector memory' | cell_type == 'T cell CD8+ naive') %>% 
        inner_join(NCICCR_S100A8[,c('case_submitter_id', 'state')], by = 'case_submitter_id') %>% 
        arrange(state)
    n2$case_submitter_id = factor(n2$case_submitter_id, levels = c(unique(n2$case_submitter_id)))
    group = c('B cell', 'B cell memory', 'Class-switched memory B cell', 'B cell naive', 'B cell plasma', 'Macrophage', 'Macrophage M1', 'Macrophage M2', 'Monocyte', 'T cell CD8+', 'T cell CD8+ central memory', 'T cell CD8+ effector memory', 'T cell CD8+ naive')
    z = data.frame()
    for(i in group){
        y = subset(n2, cell_type == i)
        y1 = wilcox.test(y$value ~ y$state)
        y2 = data.frame(pvalue = y1$p.value, group = i)
        z = rbind(y2, z)}
    z2 = z %>% 
        mutate(sign = case_when(pvalue >= 0.05 ~ 'ns',
                                pvalue < 0.05 & pvalue >= 0.01 ~ '*',
                                pvalue < 0.01 & pvalue >= 0.001 ~ '**',
                                pvalue < 0.001 ~ '***',)) %>% 
        select(group, everything()) %>% 
        dplyr::rename('cell_type' = 'group')
    n02 = n2 %>% as_tibble() %>% 
        inner_join(z2, by = 'cell_type') %>% 
        mutate(names = paste0(.$cell_type, '(', .$sign, ')')) %>% 
        ggplot(aes(x = case_submitter_id, y = names, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_viridis(alpha = 0.5, option = 'C', begin = '1', end = '0') +
        scale_y_discrete(position = "left") +
        theme_bw() +
        labs(x = "", y = "xCell") +
        guides(fill=guide_legend("xCell")) + 
        theme(legend.position = "right",
              panel.grid = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank())
    res3 = immunedeconv::deconvolute(mnci, "quantiseq") %>% 
        as.data.frame() %>% 
        pivot_longer(cols = colnames(.)[2:ncol(.)],
               names_to = "case_submitter_id",
               values_to = 'value') 
    n3 = res3 %>% 
        filter(cell_type == 'B cell' | cell_type == 'Macrophage M1' | cell_type == 'Macrophage M2' | cell_type == 'Monocyte' | cell_type == 'T cell CD8+' ) %>% 
        inner_join(NCICCR_S100A8[,c('case_submitter_id', 'state')], by = 'case_submitter_id') %>% 
        arrange(state)
    n3$case_submitter_id = factor(n3$case_submitter_id, levels = c(unique(n3$case_submitter_id)))
    group = c('B cell', 'Macrophage M1', 'Macrophage M2', 'Monocyte', 'T cell CD8+')
    z = data.frame()
    for(i in group){
        y = subset(n3, cell_type == i)
        y1 = wilcox.test(y$value ~ y$state)
        y2 = data.frame(pvalue = y1$p.value, group = i)
        z = rbind(y2, z)}
    z3 = z %>% 
        mutate(sign = case_when(pvalue >= 0.05 ~ 'ns',
                                pvalue < 0.05 & pvalue >= 0.01 ~ '*',
                                pvalue < 0.01 & pvalue >= 0.001 ~ '**',
                                pvalue < 0.001 ~ '***',)) %>% 
        select(group, everything()) %>% 
        dplyr::rename('cell_type' = 'group')
    n03 = n3 %>% as_tibble() %>% 
        inner_join(z3, by = 'cell_type') %>% 
        mutate(names = paste0(.$cell_type, '(', .$sign, ')')) %>% 
        ggplot(aes(x = case_submitter_id, y = names, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_viridis(alpha = 0.5, option = 'E', begin = '1', end = '0') +
        scale_y_discrete(position = "left") +
        theme_bw() +
        labs(x = "", y = "quanTIseq") +
        guides(fill=guide_legend("quanTIseq")) + 
        theme(legend.position = "right",
              panel.grid = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank())
    require(MCPcounter)
    probesets = read.table('~/Movies/01S100A8/probesets.txt',
                        sep="\t",
                        stringsAsFactors=FALSE,
                        colClasses="character")
    genes = read.table('~/Movies/01S100A8/genes.txt',
                    sep="\t",
                    stringsAsFactors=FALSE,
                    header=TRUE,
                    colClasses="character",
                    check.names=FALSE)
    n4 = MCPcounter::MCPcounter.estimate(mnci, featuresType = "HUGO_symbols", probesets = probesets, genes = genes) %>% 
        as.data.frame() %>% 
        rownames_to_column('cell_type') %>% 
        pivot_longer(cols = colnames(.)[2:ncol(.)],
               names_to = "case_submitter_id",
               values_to = 'value') %>% 
        filter(cell_type == 'Monocytic lineage' | cell_type == 'B lineage' | cell_type == 'CD8 T cells') %>% 
        inner_join(NCICCR_S100A8[,c('case_submitter_id', 'state')], by = 'case_submitter_id') %>% 
        arrange(state)
    n4$case_submitter_id = factor(n4$case_submitter_id, levels = c(unique(n4$case_submitter_id)))
    group = c('Monocytic lineage', 'B lineage', 'CD8 T cells')
    z = data.frame()
    for(i in group){
        y = subset(n4, cell_type == i)
        y1 = wilcox.test(y$value ~ y$state)
        y2 = data.frame(pvalue = y1$p.value, group = i)
        z = rbind(y2, z)}
    z4 = z %>% 
        mutate(sign = case_when(pvalue >= 0.05 ~ 'ns',
                                pvalue < 0.05 & pvalue >= 0.01 ~ '*',
                                pvalue < 0.01 & pvalue >= 0.001 ~ '**',
                                pvalue < 0.001 ~ '***',)) %>% 
        select(group, everything()) %>% 
        dplyr::rename('cell_type' = 'group')
    n04 = n4 %>% as_tibble() %>% 
        inner_join(z4, by = 'cell_type') %>% 
        mutate(names = paste0(.$cell_type, '(', .$sign, ')')) %>% 
        ggplot(aes(x = case_submitter_id, y = names, fill = value)) +
        geom_tile(aes(fill = value)) +
        scale_fill_viridis(alpha = 0.5, option = 'A', begin = '1', end = '0') +
        scale_y_discrete(position = "left") +
        theme_bw() +
        labs(x = "", y = "MCPcounter") +
        guides(fill=guide_legend("MCPcounter")) + 
        theme(legend.position = "right",
              panel.grid = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank())
    nx %>% insert_bottom(n0, height = 1) %>% insert_bottom(n1, height = 7) %>% insert_bottom(n02, height = 13) %>% insert_bottom(n03, height = 5) %>% insert_bottom(n04, height = 3)
    ggsave('~/Movies/01S100A8/F_immunenci.jpg', width = 20, height = 9) 
## Figure 7 Drug sensetivity analysis 
#### Immune check point analysis
    require(readxl)
    ckt = read_excel('~/Movies/01S100A8/CKTTD.xlsx')
    head(ckt, 20)
    str(ckt$'Entrez Gene ID')
    require(clusterProfiler)
    require(org.Hs.eg.db)
    ensembls = mapIds(org.Hs.eg.db, keys = ckt$'Entrez Gene ID', keytype = "ENTREZID", column="SYMBOL")
    gene = as.data.frame(ensembls) %>% 
        filter(ensembls != 'NA') %>% 
        rbind(., 'S100A8') %>% 
        mutate(guide = paste0('gene', c(1:nrow(.))))
    x = data.frame()
    symbol = gene$ensembls
    for(i in symbol){
        y1 = g87371ex %>% 
            filter(symbol == i)
        x = rbind(y1, x)
    }    
    e87371 = x %>% t()
    colnames(e87371) = e87371[1, ]
    e87371 = e87371[-1, ] %>% 
        as.data.frame() %>% 
        mutate_all(as.numeric)
    y = data.frame()
    for(i in 2:ncol(e87371)){
        y1 = cor.test(e87371[,1], e87371[,i])
        y2 = cor(e87371[,1], e87371[,i])
        y3 = data.frame(rho = y2, pvalue = y1[]$p.value, gene = colnames(e87371)[i])
        y = rbind(y3, y)
    }
    icp87371 = y %>% 
        arrange(gene) %>% 
        mutate(color = ifelse(rho >= 0.3, '1', ifelse(rho <= -0.3, '-1', '0')))        
    x = data.frame()    
    for(i in symbol){
        y1 = NCICCRb %>% 
            filter(gene_name == i)
        x = rbind(y1, x)
    } 
    eNCICCR = x %>% t()
    colnames(eNCICCR) = eNCICCR[1, ]
    eNCICCR = eNCICCR[-1, ] %>% 
        as.data.frame() %>% 
        mutate_all(as.numeric)
    y = data.frame()
    for(i in 2:ncol(eNCICCR)){
        y1 = cor.test(eNCICCR[,1], eNCICCR[, i])
        y2 = cor(eNCICCR[,1], eNCICCR[, i])
        y3 = data.frame(rho = y2, pvalue = y1[]$p.value, gene = colnames(eNCICCR)[i])
        y = rbind(y3, y)
    }
    icpNCICCR = y %>% 
        arrange(gene) %>% 
        mutate(color = ifelse(rho >= 0.3, '1', ifelse(rho <= -0.3, '-1', '0')))
    icp = icp87371$gene
    iNCICCR = data.frame()
    for(i in icp){
        target = icpNCICCR %>% 
            filter(gene == i)
        iNCICCR = rbind(iNCICCR, target)
    }
    iNCICCR
    ix = rbind(iNCICCR, icp87371) %>% 
        as.data.frame() %>% 
        mutate(group = c(rep('NCICCR_DLBC', 93), rep('GSE87371', 93)))
    ix$group = factor(ix$group, levels = c('NCICCR_DLBC', 'GSE87371'))
    cor = ix %>% 
        mutate(cor = case_when(color == '0' ~ '-3 < rho < 3',
                               color == '1' ~ '3 <= rho',
                               color == '-1' ~ ' rho <= -3',)) %>% 
        ggplot(aes(x = gene, y = group, fill = cor)) + 
        geom_tile(aes(fill = cor)) +
        scale_fill_manual(values = c("#82B0D2", "#E7DAD277", "#FA7F6F")) +
        scale_y_discrete(position = "left") +
        theme_void() +
        geom_text(aes(label = round(rho, digits = 2)), color = 'black', size = 3, angle = 90) +
        guides(fill = guide_legend(nrow = 1, byrow=TRUE)) +
        guides(fill = guide_legend("rho")) +
        theme(legend.position = "top",
              panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.ticks = element_blank(),
              axis.text.y = element_text(face = 'bold', size = 14, hjust = 1),
              axis.text.x = element_text(angle = 90, size =10, hjust = 1),
              axis.title = element_blank())
    ggsave('~/Movies/01S100A8/Fig_cor.jpg', width = 20, height = 2)
#### OncoPredict
    rm(list = ls())  
    options(stringsAsFactors = F)
    require(oncoPredict)
    require(data.table)
    require(gtools)
    require(reshape2)
    require(ggpubr)
    dir = '~/Movies/01S100A8/DataFiles/Training Data/'
    GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
    GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
    GDSC2_Res = exp(GDSC2_Res) 
    testExpr = g87371ex %>% 
            column_to_rownames('symbol') %>% 
            as.matrix()
    t = calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData')
    DrugPre = read.csv('~/Movies/01S100A8/calcPhenotype_Output/DrugPredictions87371.csv', header = T) %>% 
        dplyr::rename('geo_accession' = 'X') %>% 
        dplyr::inner_join(s87371[, c(1:2, ncol(s87371))], by = 'geo_accession') %>% 
        dplyr::select(c('S100A8', 'state'), everything()) %>% 
        dplyr::select( - 'geo_accession')
    df = data.frame()
    drugcancer = colnames(DrugPre)[-1]
    for(i in 3:ncol(DrugPre)){
        y1 = wilcox.test(DrugPre[,i] ~ DrugPre[, 2])
        high = DrugPre[, c(2, i)] %>% 
            filter(state == 'high') %>% 
            summarise(mean = mean(.[, 2]), sd = sd((.[, 2]))) %>% 
            mutate(drug = c(colnames(DrugPre)[i])) %>% 
            mutate(group = c('high'))
        low = DrugPre[, c(2, i)] %>% 
            filter(state == 'low') %>% 
            summarise(mean = mean(.[, 2]), sd = sd((.[, 2]))) %>% 
            mutate(drug = c(colnames(DrugPre)[i])) %>% 
            mutate(group = c('low'))
        diff = (high$mean - low$mean)
        y3 = data.frame(pvalue = y1[]$p.value, drug = colnames(DrugPre)[i], diff = diff)
        y2 = rbind(high, low) %>% 
            inner_join(y3, by = 'drug')
        df = rbind(y2, df)
    }
    df87371 = df %>% 
        filter(pvalue < 0.01) %>% 
        arrange(diff, pvalue) %>% 
        filter(diff < 0)
    testExpr = NCICCRb %>% 
            column_to_rownames('gene_name') %>% 
            as.matrix()
    t = calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData')
    DrugPre = read.csv('~/Movies/01S100A8/calcPhenotype_Output/DrugPredictionsNCI.csv', header = T) %>% 
        dplyr::rename('case_submitter_id' = 'X') %>% 
        dplyr::inner_join(NCICCR_S100A8[, c(1,2,ncol(NCICCR_S100A8))], by = 'case_submitter_id') %>% 
        dplyr::select(c('S100A8', 'state'), everything()) %>% 
        dplyr::select( - 'case_submitter_id')
    df = data.frame()
    drugcancer = colnames(DrugPre)[-1]
    for(i in 3:ncol(DrugPre)){
        y1 = wilcox.test(DrugPre[,i] ~ DrugPre[, 2])
        high = DrugPre[, c(2, i)] %>% 
            filter(state == 'high') %>% 
            summarise(mean = mean(.[, 2]), sd = sd((.[, 2]))) %>% 
            mutate(drug = c(colnames(DrugPre)[i])) %>% 
            mutate(group = c('high'))
        low = DrugPre[, c(2, i)] %>% 
            filter(state == 'low') %>% 
            summarise(mean = mean(.[, 2]), sd = sd((.[, 2]))) %>% 
            mutate(drug = c(colnames(DrugPre)[i])) %>% 
            mutate(group = c('low'))
        diff = (high$mean - low$mean)
        y3 = data.frame(pvalue = y1[]$p.value, drug = colnames(DrugPre)[i], diff = diff)
        y2 = rbind(high, low) %>% 
            inner_join(y3, by = 'drug')
        df = rbind(y2, df)
    }
    dfnci = df %>% 
        filter(pvalue < 0.01) %>% 
        arrange(diff, pvalue) %>% 
        filter(diff < 0)
    drugx = intersect(dfnci$drug, df87371$drug)
    drug = data.frame(drug = str_sub(drugx, end = - 6))
    dr87371 = str_sub(df87371$drug, end = - 6)
    drnci = str_sub(dfnci$drug, end = - 6)
    p87371 = df87371 %>% 
        dplyr::select(-drug) %>% 
        mutate(value = ifelse(group == "high", mean, -mean)) %>% 
        cbind(., dr87371) %>% 
        dplyr::rename('drug' = 'dr87371') %>% 
        inner_join(drug, by = 'drug') %>%
        ggplot(aes(x = drug, y = value,  fill = group)) +
        geom_bar(stat = 'identity') +      
        coord_flip() +  
        # geom_text(aes(label = round(mean, digits = 2),                                  
        # vjust = ifelse(group == "high", -0.5, 1),       
        # hjust = ifelse(group == "high", -0.4, 1)), size = 4) +
        scale_y_continuous(labels = abs, expand = expansion(mult = c(0.1, 0.1))) +
        scale_fill_manual(values = c("#FA7F6F", "#8ECFC9")) +
        theme_classic() +
        labs(title = 'GSE87371', y = 'Imputed Sensitivity Score (mean)', x ='') +
            geom_hline(aes(yintercept = 0), color = 'black', linetype = 'solid') +
            theme(legend.position = "none",
                  axis.text.x = element_text(size = 16, angle = 0, hjust = 1, vjust = 0.5),
                  axis.text.y = element_text(size = 16, angle = 0, hjust = 1, vjust = 0.5),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 16, hjust = 0.5))
    pnci= dfnci %>% 
        dplyr::select(-drug) %>% 
        mutate(value = ifelse(group == "high", mean, -mean)) %>% 
        cbind(., drnci) %>% 
        dplyr::rename('drug' = 'drnci') %>% 
        inner_join(drug, by = 'drug') %>% 
        ggplot(aes(x = drug, y = value,  fill = group)) +
        geom_bar(stat = 'identity') +      
        coord_flip() +  
        # geom_text(aes(label = round(mean, digits = 2),                                  
        # vjust = ifelse(group == "high", -0.5, 1),       
        # hjust = ifelse(group == "high", -0.4, 1)), size = 4) +
        scale_y_continuous(labels = abs, expand = expansion(mult = c(0.1, 0.1))) +
        scale_fill_manual(values = c("#FA7F6F", "#8ECFC9")) +
        theme_classic() +
        labs(title = 'NCICCR_DLBCL', x = '', y = 'Imputed Sensitivity Score (mean)') +
        guides(fill = guide_legend("Status of S100A8")) + 
        geom_hline(aes(yintercept = 0), color = 'black', linetype = 'solid') +
            theme(legend.position = "right",
                  axis.text.x = element_text(size = 16, angle = 0, hjust = 1, vjust = 0.5),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.line.y = element_blank(),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 16, hjust = 0.5))
    require(patchwork)
    p87371 | pnci    
    ggsave('~/Movies/01S100A8/Fig_drug1.jpg', width = 16, height =6)
    plot87371 = df87371 %>% 
        dplyr::select(-drug) %>% 
        mutate(value = ifelse(group == "high", mean, -mean)) %>% 
        cbind(., dr87371) %>% 
        dplyr::rename('drug' = 'dr87371') %>% 
        ggplot(aes(x = drug, y = value,  fill = group)) +
        geom_bar(stat = 'identity') +      
        coord_flip() +  
        # geom_text(aes(label = round(mean, digits = 2),                                  
        # vjust = ifelse(group == "high", -0.5, 1),       
        # hjust = ifelse(group == "high", -0.4, 1)), size = 4) +
        scale_y_continuous(labels = abs, expand = expansion(mult = c(0.1, 0.1))) +
        scale_fill_manual(values = c("#FA7F6F", "#8ECFC9")) +
        theme_classic() +
        labs(title = 'GSE87371', y = 'Imputed Sensitivity Score (mean)', x ='') +
            geom_hline(aes(yintercept = 0), color = 'black', linetype = 'solid') +
            theme(legend.position = "none",
                  axis.text.x = element_text(size = 16, angle = 0, hjust = 1, vjust = 0.5),
                  axis.text.y = element_text(size = 16, angle = 0, hjust = 1, vjust = 0.5),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 16, hjust = 0.5))
    plotnci= dfnci %>% 
        dplyr::select(-drug) %>% 
        mutate(value = ifelse(group == "high", mean, -mean)) %>% 
        cbind(., drnci) %>% 
        dplyr::rename('drug' = 'drnci') %>% 
        ggplot(aes(x = drug, y = value,  fill = group)) +
        geom_bar(stat = 'identity') +      
        coord_flip() +  
        # geom_text(aes(label = round(mean, digits = 2),                                  
        # vjust = ifelse(group == "high", -0.5, 1),       
        # hjust = ifelse(group == "high", -0.4, 1)), size = 4) +
        scale_y_continuous(labels = abs, expand = expansion(mult = c(0.1, 0.1))) +
        scale_fill_manual(values = c("#FA7F6F", "#8ECFC9")) +
        theme_classic() +
        labs(title = 'NCICCR_DLBCL', y = 'Imputed Sensitivity Score (mean)', x ='') +
            geom_hline(aes(yintercept = 0), color = 'black', linetype = 'solid') +
            theme(legend.position = "none",
                  axis.text.x = element_text(size = 16, angle = 0, hjust = 1, vjust = 0.5),
                  axis.text.y = element_text(size = 16, angle = 0, hjust = 1, vjust = 0.5),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(size = 16, hjust = 0.5))
    plot87371 | plotnci
    ggsave('~/Movies/01S100A8/Fig_drug2.jpg', width = 16, height =8)
    require(ggvenn)
    datalist = list("GSE87371" = unique(dr87371), "NCICCR_DLBCL" = unique(drnci))
    ggvenn(datalist, fill_color = color,
        fill_alpha = .7,
        stroke_linetype = "solid",
        set_name_size = 10,
        text_size = 9)  
    ggsave("~/Movies/01S100A8/F_Venn2.png", width = 10, height = 10)
    drug = colnames(DrugPre)[-1]
    for(i in 2:ncol(DrugPre)){
        y1 = cor.test(DrugPre[,1], DrugPre[,i])
        y2 = cor(DrugPre[,1], DrugPre[,i])
        y3 = data.frame(rho = y2, pvalue = y1[]$p.value, gene = colnames(DrugPre)[i])
        df = rbind(y3, df)
    }
    drug = str_sub(df$gene, end = - 6)
    df1 = cbind(df, drug) %>% as_tibble() %>% 
        filter(abs(rho) >= 0.3) %>% 
        arrange(-rho) %>% 
        distinct(drug, .keep_all = T)
    df1$drug = factor(df1$drug, levels = c(df1$drug))
    ncidrug = df1$drug
    df3 = df1 %>% 
        mutate(group = ifelse(rho < 0, 'neg', 'pos')) %>% 
        ggplot(aes(x = drug, y = rho, fill = group)) +
            geom_bar(stat="identity") +
            labs(y = 'rho score', x = '') +
            scale_fill_manual(values = c("#8ECFC9", "#FA7F6F")) +
            geom_text(aes(label = round(rho, 2)), size = 4, parse = TRUE) + 
            theme_classic() +
            labs(title = 'NCICCR_DLBCL') +
            geom_hline(aes(yintercept = 0), color = 'black', linetype = 'solid') +
            theme(legend.position = "none",
              axis.text.x = element_text(size = 16, angle = 90, hjust = 1, vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    require(ggvenn)
    diffdrug = intersect(ncidrug, GSEdrug)
    df2 / df3
    ggsave('~/Movies/01S100A8/Fig_drug1.jpg', width = 16, height =10)
#### CellMiner analysis
    require(readxl)
    require(impute)
    require(limma)
    require(WGCNA)
    require(tidyr)
    dat1 = read_excel(path = "~/Movies/01S100A8/DTP_NCI60_ZSCORE.xlsx", skip = 7)
    colnames(dat1) = dat1[1,]
    dat1 = dat1[-1,-c(67,68)]
    table(dat1$'FDA status')
    dat1 = dat1[dat1$'FDA status' %in% c("FDA approved", "Clinical trial"),]
    dat1 = dat1[,-c(1, 3:6)]
    drugDat = dat1 %>% as.matrix()
    rownames(drugDat) = drugDat[,1]
    drug = drugDat[, 2:ncol(drugDat)]
    dimnames = list(rownames(drug),colnames(drug))
    data = matrix(as.numeric(as.matrix(drug)), nrow=nrow(drug), dimnames=dimnames)
    data = data[, which(colMeans(!is.na(data)) > 0.8)]
    mat = impute.knn(data)
    drug = mat$data
    drug = avereps(drug) %>% 
        t() %>% 
        as.data.frame()
    colnames(drug)[1:12]
    dat2 = read_excel(path = "~/Movies/01S100A8/RNA__RNA_seq_composite_expression.xls", skip = 9)
    colnames(dat2) = dat2[1,]
    expro = dat2[-1,-c(2:6)] %>% 
        column_to_rownames('Gene name d')
    inputgene = c("S100A8")
    gl = intersect(inputgene, row.names(expro))
    exp0 = expro[gl,] %>% 
        t() %>% 
        as.data.frame()
    identical(rownames(exp0),rownames(drug))
    exp1 = exp0[which(rownames(exp0) %in% rownames(drug)), ] %>% as.numeric()
    corTab = cor(drug, exp1, method="pearson")
    corPval = corPvalueStudent(corTab, nSamples = nrow(drug))
    fitercor = cbind(corTab, corPval) %>% 
        as.data.frame() %>% 
        dplyr::rename('Pearson' = 'V1', 'Pval' = 'V2') %>% 
        filter(abs(Pearson) > 0.3) %>% 
        filter(Pval < 0.01) %>% 
        rownames_to_column('drug') %>% 
        arrange(-Pearson, Pval) %>% 
        column_to_rownames('drug')
    fitercor1 = fitercor %>% 
        rownames_to_column('drug')
    # write_csv2(fitercor1, file = '~/Movies/01CISD2/DrugSensetivity.csv')
    require(ggpubr)
    g = c('S100A8')
    for(dr in rownames(fitercor)[1:9]){
        df = data.frame(exp = exp1, dr1 = drug[, dr])
        tit = paste0("R: ",round(fitercor[dr, 1], 3),", P-value < 0.01")
        p = ggplot(data = df, aes(x = exp, y = dr1)) + 
            geom_point(alpha = 0.6, shape = 19, size=3, color="#E6550DFF") + 
            geom_smooth(method = lm, formula = y ~ x,aes(colour = "lm"), size = 1.2,se = T)+
            scale_color_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) +
            theme_bw() +
            labs(y = paste0("Activity z scores of ", dr), x = (paste0("The expression of ", g))) +
            annotate("text", label = tit, size = 4, x = -Inf, y = Inf, hjust = -.2, vjust = 2) +
            theme(legend.position = "none",q
              axis.text.x = element_text(size = 10,vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave(paste0("opFig/",g,"-",dr,"-cor.png"), width = 5, height = 5)}
    for(dr in rownames(fitercor)[1:9]){
        df = data.frame(exp = exp1, dr1 = drug[, dr])
        med = median(df$exp)
        df$group = ifelse(df$exp > med,"High","Low")
        p = ggplot(df, aes(group, dr1, fill= group)) + 
            geom_violin(aes(fill = group),trim = FALSE) +
            geom_signif(comparisons = list(c("High","Low")),
                    step_increase = 0.1,
                    map_signif_level = T,
                    margin_top=0.2,
                    tip_length =0.02,
                    test = "t.test")+
            geom_boxplot(width = 0.1,fill = "white") +
            scale_fill_manual(values = c("#E6550DFF", "#3182BDFF")) +
            theme_bw() +
            labs(y= paste0("Activity z scores of ", dr)) +
            theme(legend.position = "none",
              axis.text.x = element_text(size = 10,vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
      ggsave(filename = paste0("opFig/",g,"-",dr,"-violin.png"),
             width = 5, height = 5)}
    fit = fitercor1 %>% 
        mutate(iterm = factor(rownames(fitercor), levels = rownames(fitercor))) %>% 
        mutate(group = ifelse(Pearson >= 0, '1', '2')) %>% 
        ggplot(aes(x = iterm, y = Pearson, fill = group)) +
            geom_bar(stat="identity") +
            labs(y = 'Pearson Score', x = '') +
            scale_fill_manual(values = color) +
            geom_text(aes(label = round(Pearson, 2)), size = 8, parse = TRUE) + 
            theme_classic() +
            # coord_flip() +
            geom_hline(aes(yintercept = 0), color = 'black', linetype = 'solid') +
            theme(legend.position = "none",
              axis.text.x = element_text(size = 36, angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 24, angle = 0, hjust = 1, vjust = 0.5),
              axis.title = element_text(size = 32),
              plot.title = element_text(size = 16))
    require(patchwork)
    ggsave('~/Movies/01S100A8/Fig_drug.jpg', width = 20, height = 16)
