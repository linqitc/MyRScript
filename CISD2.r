## Packages for this study and Loading workplace
    library(tidyverse) # Data modulation and plotting.
    library(data.table) # Data convenience features.
    library(readxl) # Reading and writing from excel.
    library(Rmisc) # Calculation SE for ploting SE bar.
    library(ggsignif) # Significance of ploting using ggplot2.
    library(scales) # Scale Functions for Visualization.
    library(ggsci) 
    library(limma)
    library(GEOquery)
    library(TCGAbiolinks)
    library(survival)
    library(survminer)
## Figure 2 CISD2 expression analysis
#### Clean a TCGA pheno data
    gtex = read.table("~/Documents/01CISD2/samplepair.txt", header = T, sep = "\t")
    tcga_ref = gtex[, 1:2]
    gtex$type = paste0(gtex$TCGA, "_normal_GTEx")
    gtex$sample_type = "normal"
    gtex = gtex[, c("TCGA", "GTEx", "type", "sample_type")]
    names(gtex)[1:2] = c("tissue", "X_primary_site")
    gp = read.table("~/Documents/01CISD2/GTEx_phenotype.txt", header = T, sep = "\t")
    gtex2tcga = merge(gtex, gp, by = "X_primary_site")
    gtex_data = gtex2tcga[, c(5, 2:4)]
    names(gtex_data)[1] = "sample"
    tcga = read.table("~/Documents/01CISD2/TCGA_phenotype.txt", header = T, sep = "\t")
    tcga = merge(tcga_ref, tcga, by.y = "X_primary_disease", by.x = "Detail", all.y = T)
    tcga = tcga[tcga$sample_type %in% c("Primary Tumor", "Solid Tissue Normal"), ]
    tcga$type = ifelse(tcga$sample_type == "Solid Tissue Normal",
        paste(tcga$TCGA, "normal_TCGA", sep = "_"), paste(tcga$TCGA, "tumor_TCGA", sep = "_"))
    tcga$sample_type = ifelse(tcga$sample_type == "Solid Tissue Normal", "normal", "tumor")
    tcga = tcga[, c(3, 2, 6, 5)]
    names(tcga)[2] = "tissue"
#### Remove samples without tpm data
    gtex_exp = fread("~/Documents/01CISD2/gtex_RSEM_gene_tpm.gz", data.table = F)
    gtexS = gtex_data[gtex_data$sample %in% colnames(gtex_exp)[-1], ]
    tcga_exp = fread("~/Documents/01CISD2/tcga_RSEM_gene_tpm.gz", data.table = F)
    tcgaS = tcga[tcga$sample %in% colnames(tcga_exp)[-1], ]
    tcga_gtex = rbind(tcgaS, gtexS)
    idmap = read.delim("~/Documents/01CISD2/gencode.v23.annotation.gene.probemap", as.is = T)
#### Chance CISD2
    target = "CISD2"
    id = idmap$id[which(idmap$gene == target)] ##  ENSG00000145354
    tcga_data = t(tcga_exp[tcga_exp$sample == id, colnames(tcga_exp) %in% c("sample", tcga_gtex$sample)])
    tcga_data = data.frame(tcga_data[-1, ])
    tcga_data = rownames_to_column(tcga_data, "sample")
    names(tcga_data)[2] = "tpm"
    gtex_data = t(gtex_exp[gtex_exp$sample == id, colnames(gtex_exp) %in% c("sample", tcga_gtex$sample)])
    gtex_data = data.frame(gtex_data[-1, ])
    gtex_data = rownames_to_column(gtex_data, "sample")
    names(gtex_data)[2] = "tpm"
    tmp = rbind(tcga_data, gtex_data)
    exp = merge(tmp, tcga_gtex, by = "sample", all.x = T)
    exp1 = exp[, c("tissue", "sample_type", "tpm")]
    exp1 = arrange(exp, tissue)
    exp1$tpm = as.numeric(exp1$tpm)
    colnames(exp1) = c('sample', 'tpm', 'tissue', 'type', 'Group')
    require(ggplot2)
    require(ggpubr)
    require(RColorBrewer)
    p1 = ggboxplot(exp1,
        x = "tissue", y = "tpm", fill = "Group",
        ylab = 'CISD2 expression level (log2 TPM)',
        xlab = '',
        title = 'TCGA-DLBC + GTEx',
        color = "Group",
        palette = c("#3182BDFF", "#E6550DFF"),
        ggtheme = theme_bw())
    count_N = exp1 %>%
        group_by(tissue, Group) %>%
        tally()
    count_N$n = paste("n =", count_N$n)
    p1 = p1 + geom_text(
        data = count_N, aes(label = n, y = -9, color = Group), position = position_dodge2(0.9),
        size = 3, angle = 90, hjust = 0) +
        theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1.2),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16)) +
        theme(legend.position="right")
    comp = compare_means(tpm ~ Group,
        group.by = "tissue", data = exp1,
        symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")))
    p1 + stat_pvalue_manual(comp,
        x = "tissue", y.position = 7.5,
        label = "p.signif", position = position_dodge(0.8))
    ggsave("~/Documents/01CISD2/F_tgca.png", width = 15, height = 5)
#### TCGA_GTEx compare plot
    require(ggsci)
    ggplot(subset(exp1, tissue == 'DLBC'), aes(fill = Group, y = tpm, x = Group)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(data = subset(exp1, tissue == 'DLBC'), aes(y = tpm), size = 2, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c("normal", "tumor")), size = 0.5, textsize = 10, vjust = 0.3) +
        labs(title = "TCGA-DLBC + GTEx", y = "CISD2 expression level (log2 TPM)", x = "") + 
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF")) +
        theme_bw() + 
        scale_x_discrete(labels = c("normal" = "Normal","tumor" = "DLBCL")) +
        scale_y_continuous(limits = c(-3, 7)) +
        theme(legend.position="none",
        axis.text.x = element_text(size = 14, vjust = 0.5),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16))
    ggsave("~/Documents/01CISD2/F_tcgaCom.png", width = 4, height = 5)
#### TCGA_GTEx ROC plot
    require(plotROC)
    require(pROC)
    tcga_roc = roc(data = subset(exp1, tissue == 'DLBC'), response = Group, predictor = tpm)
    ggroc(tcga_roc, colour = "#E6550DFF", size = 1, legacy = T) +
        geom_abline(linetype = 3, intercept = 0) +
        labs(title = "TCGA-DLBC + GTEx") +
        theme_bw() + 
        theme(axis.text = element_text(size = 14, vjust = 0.5),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16)) +
        annotate("text", x = 0.7, y = 0.1, label = "AUC = 0.8181", size = 5)
    ggsave("~/Documents/01CISD2/F_tcgaroc.png", width = 4.5, height = 5)
#### GSE83632 based on GPL5175
    G83632 = getGEO(filename = ("~/Documents/01CISD2/GSE83632_series_matrix.txt.gz"), GSEMatrix = TRUE, getGPL = FALSE)
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
        file = "~/Documents/01CISD2/GPL5175_noParents.an.txt",
        header = T, sep = "\t", quote = "", fill = T,
        comment.char = "!") %>%
        select("ProbeName", "GeneSymbols") %>% 
        dplyr::rename("ID" = "ProbeName", "symbol" = "GeneSymbols")
    g83632ee = g83632e %>%
        rownames_to_column("ID") %>%
        merge(an5175, by = "ID") %>%
        dplyr::select(-"ID") %>%
        distinct(symbol, .keep_all = T)
    g83632x = tibble::column_to_rownames(g83632ee, "symbol")
    c83632 = g83632x %>%
        tibble::rownames_to_column("symbol") %>%
        filter(symbol == "CISD2") %>%
        column_to_rownames("symbol") %>%
        gather() %>% 
        dplyr::rename(., geo_accession = key, CISD2 = value) %>%
        inner_join(g83632c, by = "geo_accession") %>%
        dplyr::select("geo_accession", "CISD2", "source_name_ch1")  %>% 
        mutate(states = ifelse(source_name_ch1 == "Healthy donor, blood", "Health controls", "DLBCL")) %>%
        dplyr::select(-source_name_ch1)
#### Compare plot.
    require(ggsci)
    c83632$states = factor(c83632$states, levels = c('Health controls', 'DLBCL'))
    ggplot(c83632, aes(fill = states, y = CISD2, x = states)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(data = c83632, aes(y = CISD2), size = 2, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c("Health controls", "DLBCL")), size = 0.5, textsize = 10, vjust = 0.3) +
        labs(title = "GSE83632", y = "The expression of CISD2", x = "") +
        scale_y_continuous(limits = c(3, 12)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF")) + theme_bw() +
        theme(legend.position = "none",
        axis.text = element_text(size = 14, vjust = 0.5),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16))
    ggsave("~/Documents/01CISD2/F_G83632Com.png", width = 4, height = 5)
#### ROC plot
    library(plotROC)
    library(pROC)
    r83632 = roc(response = c83632$states, predictor = c83632$CISD2)
    auc = round(auc(c83632$states, c83632$CISD2), 4)
    ggroc(r83632, colour = "#E6550DFF", size = 1, legacy = T) +
        geom_abline(linetype = 3, intercept = 0) +
        labs(title = "GSE83632") +
        theme_bw() +
        theme(axis.text = element_text(size = 14, vjust = 0.5),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16)) +
        annotate("text", x = 0.75, y = 0.1, label = "AUC = 0.8274", size = 6)
    ggsave("~/Documents/01CISD2/F_G83632roc.png", width = 4.5, height = 5)
#### Limma analysis and Volcano Plot
    g83632l = g83632c %>% 
        mutate(states = ifelse(source_name_ch1 == "Healthy donor, blood", "Health", "DLBCL")) %>%
        select(-source_name_ch1)
    g = g83632l$states
    design = model.matrix(~ -1 + g)
    colnames(design) = c("Health", "DLBCL")
    contrast.matrix = makeContrasts(Health - DLBCL, levels = design)
    g83632ex = tibble::column_to_rownames(g83632ee, "symbol")
    fit = lmFit(g83632ex, design)
    fit2 = contrasts.fit(fit, contrast.matrix)
    fit2 = eBayes(fit2)
    sig = topTable(fit2, n = Inf, adjust = "fdr") # Chance adjust = 'fdr' & n = Inf.
    nrDEG = na.omit(sig) # Get rid of NA.
    d83632 = nrDEG
    head(d83632)
    foldChange = 1
    padj = 0.05 # Selection criteria.
    colnames(d83632)[1] = "log2FC"
    d83632 = within(d83632, {
        group = NA
        group[(d83632$adj.P.Val > padj | d83632$adj.P.Val == "NA") | (d83632$log2FC < foldChange) & d83632$log2FC > -foldChange] = "Not"
        group[(d83632$adj.P.Val <= padj & d83632$log2FC > foldChange)] = "Up"
        group[(d83632$adj.P.Val <= padj & d83632$log2FC < -foldChange)] = "Down"})
    df83632 = d83632 %>%
        filter(group == "Up" | group == "Down")
    d83632 = d83632[order(d83632$adj.P.Val), ]
    d83632 = d83632 %>% tibble::rownames_to_column("symbol")
    d83632$log10P = -log10(d83632$adj.P.Val)
    target = "CISD2"
    require(ggpubr)
    require(ggsci)
    ggscatter(d83632, "log2FC", "log10P",
        combine = F, merge = T, color = "group", shape = 20, size = 5,
        point = TRUE, font.label = 20, palette = c("#3182BDFF", "grey", "#E6550DFF")) +
        geom_vline(xintercept = c(-1, 1), colour = "black", linetype = "dashed") +
        geom_hline(yintercept = -log10(0.05), colour = "black", linetype = "dashed") +
        labs(title = "GSE83632") + theme_bw() + ggrepel::geom_label_repel(aes(label = symbol), data = d83632[1016, ]) +
        coord_flip() +
        theme(legend.position=c(.9, .88),
        legend.title = element_blank(),
        axis.text = element_text(size = 14, vjust = 0.5),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16))
    ggsave("~/Documents/01CISD2/F_G83632forrest.png", width = 7, height = 5)
## Figure 2 Compare DLBCL patients with Health controls
#### PCR calculation
    CISD2 = c(24.355, 25.746, 26.293, 27.277, 24.082, 25.965, 26.387, 27.535, 24.121, 25.520, 26.543, 27.309)
    Actin = c(21.074, 24.645, 21.730, 25.418, 21.059, 24.652, 21.840, 25.770, 21.137, 24.574, 21.949, 25.059)
    group = rep(c('DB', 'SUDHL2', 'GM12878', 'SUDHL4'), 3)
    PCR = tibble(CISD2,Actin, group)
    PCR$dct1=PCR$CISD2 - PCR$Actin  ##目的基因Ct-内存基因Ct，即∆Ct
    PCR$ddct1=PCR$dct1 - mean(PCR$dct1[which(PCR$group == 'GM12878')])  ##∆Ct-对照组Ct均值，即∆∆Ct
    PCR$mrna1=2^-PCR$ddct1  ##取-∆∆Ct的2次放，即-2^∆∆Ct
    PCR$group = factor(PCR$group, levels = c('GM12878', 'DB', 'SUDHL4', 'SUDHL2'))
    ggplot(PCR, aes(y = mrna1, x = group, fill = group)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(data = PCR, aes(y = mrna1), size = 2, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        labs(title = "", y = "The relative expression of CISD2", x = '') +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        geom_signif(map_signif_level = T, color = "black", test = 't.test',
            comparisons = list(c("GM12878", "DB"), c('SUDHL4', 'GM12878'), c('SUDHL2', 'GM12878')), 
            size = 0.5, textsize = 5, vjust = 0.5, step_increase = 0.1, position = "identity") + 
        theme(legend.position = "none",
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 20))
    ggsave("~/Documents/01CISD2/F_PCR.png", width = 4.5, height = 5)
#### WB expression
    expr_CISD2 = c(22.943, 38.404, 64.626, 78.109, 32.643, 45.335, 81.987, 56.43, 33.496, 48.961, 68.217, 84.596)
    expr_actin = c(139.053,130.207,121.232,138.579,152.407,144.749,159.401,156.33,147.668,139.943,146.737,146.737)
    group = rep(c('GM12878', 'DB', 'SUDHL4', 'SUDHL2'), 3)
    expr = tibble(expr_CISD2, expr_actin, group, relative = expr_CISD2/expr_actin)
    require(ggsci)
    expr$group = factor(expr$group, levels = c('GM12878', 'DB', 'SUDHL4', 'SUDHL2'))
    ggplot(expr, aes(fill = group, y = relative, x = group)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(data = expr, aes(y = relative), size = 2, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c("GM12878", "DB"), c('SUDHL4', 'GM12878'), c('SUDHL2', 'GM12878')), 
            size = 0.5, textsize = 5, vjust = 0.3, test = "t.test", step_increase = 0.1) +
        labs(title = "", y = "The relative expression of CISD2", x = "") +
        # scale_y_continuous(limits = c(0, 2.4)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + theme_bw() +
        theme(legend.position = "none",
        axis.text.x = element_text(size = 14, angle = 45,  hjust = 1),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 20))
    ggsave("~/Documents/01CISD2/F_exprWB.png", width = 4, height = 5)
## Figure 3 CISD2 Survival Plot
#### GSE31312 based on GPL570 (Trainning dataset)
    an570 = getGEO(filename = ("~/Documents/01CISD2/GPL570.annot.gz")) %>% 
        Table(.) %>% 
        dplyr::select("ID", "Gene symbol") %>%
        filter("Gene symbol" != "---") %>%
        dplyr::rename("symbol" = "Gene symbol")
    G31312 = getGEO(filename = ("~/Documents/01CISD2/GSE31312_series_matrix.txt.gz"), GSEMatrix = T, getGPL = F)
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
        rownames_to_column("ID") %>%
        inner_join(an570, by = "ID") %>%
        dplyr::select(-"ID") %>%
        dplyr::select(symbol, everything()) %>% 
        aggregate(. ~ symbol, median)
    g31312cx = read.csv("~/Documents/01CISD2/GSE31312_clinical_data.csv", header = T) %>% 
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
    c31312 = g31312ex[grep('CISD2', g31312ex$symbol), ] %>% 
        gather() %>% 
        dplyr::rename('symbol' = 'key', 'CISD2' = 'value') %>% 
        slice(- 1 ) %>% 
        mutate_at('CISD2', as.numeric) %>% 
        dplyr::rename('geo_accession' = 'symbol') %>% 
        inner_join(g31312cx, by = 'geo_accession') %>% 
        filter(OSmonth >= 1) %>% 
        mutate(status = ifelse(CISD2 < median(CISD2), 'low', 'high'))
    fit = survfit(Surv(OSmonth, OS) ~ status, data = c31312)
    data.survdiff <- survdiff(Surv(OSmonth, OS) ~ status, data = c31312)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
    up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
    low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
    p1 = ggsurvplot(fit,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2 High", "CISD2 Low"),
            legend.title = "",
            title = "GSE31312", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    p1$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G31312surALL.png', width = 6, height = 4)
    data.survdiff <- survdiff(Surv(PFSmonth, PFS) ~ status, data = c31312)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
    up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
    low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
    fit = survfit(Surv(PFSmonth, PFS) ~ status, data = c31312)
    p1 = ggsurvplot(fit,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2 High", "CISD2 Low"),
            legend.title = "",
            title = "GSE31312", 
            xlab = " Time (Month)",
            ylab = 'Progression-Free Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    p1$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G31312pfsALL.png', width = 6, height = 4)
    fit = survfit(Surv(OSmonth, OS) ~ status, data = subset(c31312, COO == 'ABC'))
    p1 = ggsurvplot(fit,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2 High", "CISD2 Low"),
            legend.title = "",
            title = "GSE31312 (ABC subtype)", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    p1$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G31312surABC.png', width = 6, height = 4)
    gcb31312 = c31312 %>% filter(COO == 'GCB')
    res.cut = surv_cutpoint(gcb31312, time = "OSmonth", event = "OS", variables = 'CISD2')
    gcb31312 = gcb31312 %>% 
        mutate(state = ifelse(CISD2 < res.cut$cutpoint[1, 1], 'low', 'high'))
    fit = survfit(Surv(OSmonth, OS) ~ state, data = gcb31312)
    p1 = ggsurvplot(fit,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2 High", "CISD2 Low"),
            legend.title = "",
            title = "GSE31312 (GCB subtype)", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    p1$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G31312surGCB.png', width = 6, height = 4)
    d31312 = g31312ex[grep('CISD2', g31312ex$symbol), ] %>% 
        gather() %>% 
        dplyr::rename('symbol' = 'key', 'CISD2' = 'value') %>% 
        slice(- 1 ) %>% 
        mutate_at('CISD2', as.numeric) %>% 
        dplyr::rename('geo_accession' = 'symbol') %>% 
        inner_join(g31312cx, by = 'geo_accession') %>% 
        filter(PFSmonth > 1) %>% 
        mutate(status = ifelse(CISD2 < median(CISD2), 'low', 'high'))
    fit = survfit(Surv(PFSmonth, PFS) ~ status, data = d31312)
    p1 = ggsurvplot(fit,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2 High", "CISD2 Low"),
            legend.title = "",
            title = "GSE31312", 
            xlab = " Time (Month)",
            ylab = 'Progress Free Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    p1$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G31312pfsALL.png', width = 6, height = 4)
    fit = survfit(Surv(PFSmonth, PFS) ~ status, data = subset(d31312, COO == 'ABC'))
    p1 = ggsurvplot(fit,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2 High", "CISD2 Low"),
            legend.title = "",
            title = "GSE31312 (ABC subtype)", 
            xlab = " Time (Month)",
            ylab = 'Progress Free Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    p1$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G31312pfsABC.png', width = 6, height = 4)
    gcd31312 = d31312 %>% filter(COO == 'GCB')
    res.cut = surv_cutpoint(gcd31312, time = "PFSmonth", event = "PFS", variables = 'CISD2')
    gcd31312 = gcd31312 %>% 
        mutate(state = ifelse(CISD2 < res.cut$cutpoint[1, 1], 'low', 'high'))
    fit = survfit(Surv(PFSmonth, PFS) ~ state, data = gcd31312)
    p1 = ggsurvplot(fit,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2 High", "CISD2 Low"),
            legend.title = "",
            title = "GSE31312 (GCB subtype)", 
            xlab = " Time (Month)",
            ylab = 'Progress Free Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    p1$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G31312pfsGCB.png', width = 6, height = 4)
#### GSE32918 based on GPL8432 (Validation dataset)
    an8432 = data.table::fread("~/Documents/01CISD2/GPL8432-11703.txt") %>% 
        select("ID", symbol = "Symbol")
    G32918 <- getGEO(filename = ('~/Documents/01CISD2/GSE32918_series_matrix.txt.gz'), GSEMatrix = T, getGPL= F)
    if (length(G32918) > 1) {idx = grep("GPL8432", attr(G32918, "names"))} else {idx = 1}
    g32918e = data.frame(exprs(G32918))
    g32918c = pData(G32918)
    qx = as.numeric(quantile(g32918e, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC = (qx[5] > 100) ||
      (qx[6]-qx[1] > 50 && qx[2] > 0) ||
      (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if(LogC){g32918e[which(g32918e <= 0)] = NaN
        g32918e = log2(g32918e)
        print("log2 transform finished")}else{print("log2 transform not needed")}
    g32918ex = g32918e %>% 
        rownames_to_column('ID') %>%
        merge(an8432, by = 'ID') %>%
        dplyr::select(- 'ID') %>%
        dplyr::select(symbol, everything())
    colnames(g32918c) = gsub(':ch1', '', colnames(g32918c))
    colnames(g32918c) = gsub(' ', '', colnames(g32918c))
    colnames(g32918c) = gsub('-', '', colnames(g32918c))
    c32918 = g32918ex %>%
        filter(symbol == "CISD2") %>%
        column_to_rownames('symbol') %>%
        gather() %>%
        dplyr::rename('geo_accession' = 'key', 'CISD2' = 'value') %>% 
        inner_join(g32918c, by = 'geo_accession') %>% 
        mutate(followup = as.numeric(case_when(followupstatus == 'Alive' ~ 0, followupstatus == 'Dead' ~ 1))) %>% 
        mutate_at('followupyears', as.numeric) %>% 
        mutate(followmonth = followupyears * 12) %>% 
        mutate(status = ifelse(CISD2 < median(CISD2), 'Low', 'High'))
    data.survdiff <- survdiff(Surv(followmonth, followup) ~ status, data = c32918)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
    up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
    low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
    fit = survfit(Surv(followmonth, followup) ~ status, data = c32918) %>% 
            ggsurvplot(.,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2 High", "CISD2 Low"),
            legend.title = "",
            title = "GSE32918", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G32918OS.png', width = 6, height = 4)
    fit = survfit(Surv(followmonth, followup) ~ status, data = c32918 %>% filter(predictedclass == 'GCB')) %>% 
            ggsurvplot(., 
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2 High", "CISD2 Low"),
            legend.title = "",
            title = "GSE32918 (GCB subtype)", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G32918OSgcb.png', width = 6, height = 4)
    abc32918 = subset(c32918, predictedclass == 'ABC') %>% 
        mutate(state = ifelse(CISD2 < surv_cutpoint(., time = "followmonth", event = "followup", variables = 'CISD2')$cutpoint[1,1], 'Low', 'High')) 
    fit = survfit(Surv(followmonth, followup) ~ state, data = abc32918)
    fit = ggsurvplot(fit, 
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2 High", "CISD2 Low"),
            legend.title = "",
            title = "GSE32918 (ABC subtype)", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G32918OSabc.png', width = 6, height = 4)
#### GSE93984 based on GPL570 (Valiation 2).
    G93984 = getGEO(filename = "~/Documents/01CISD2/GSE93984_series_matrix.txt.gz", GSEMatrix = T, getGPL = F)
    if (length(G93984) > 1) {
        idx = grep("GPL570", attr(G93984, "names"))} else {
        idx = 1}
    g93984e = data.frame(exprs(G93984))
    g93984c = pData(G93984)
    qx = as.numeric(quantile(g93984e, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    LogC = (qx[5] > 100) ||
        (qx[6] - qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) {
        g93984e[which(g93984e <= 0)] = NaN
        G93984e = log2(g93984e)
        print("log2 transform finished")} else {
        print("log2 transform not needed")}
    an570 = getGEO(filename = ("~/Documents/01CISD2/GPL570.annot.gz")) %>% 
        Table(.) %>% 
        dplyr::select("ID", "Gene symbol") %>%
        filter("Gene symbol" != "---") %>%
        dplyr::rename("symbol" = "Gene symbol")
    g93984ex = g93984e %>%
        rownames_to_column("ID") %>%
        inner_join(an570, by = "ID") %>%
        dplyr::select(-"ID") %>%
        dplyr::select(symbol, everything()) %>% 
        aggregate(. ~ symbol, median)
    c93984 = g93984ex[grep("CISD2", g93984ex$symbol), ] %>%
        gather() %>%
        dplyr::rename("symbol" = "key", "CISD2" = "value") %>%
        slice(-1) %>%
        mutate_at("CISD2", as.numeric) %>%
        dplyr::rename("geo_accession" = "symbol") %>%
        inner_join(g93984c, by = "geo_accession") %>%
        mutate(average = median(CISD2), status = ifelse(CISD2 < average, "Low", "High")) %>%
        select(-average)
    colnames(c93984) = gsub(':ch1', '', colnames(c93984)) 
    c93984$progression_free_survival_censor = as.numeric(c93984$progression_free_survival_censor)
    c93984$progression_free_survival = as.numeric(c93984$progression_free_survival)
    c93984 = c93984 %>% 
        mutate(pfs = progression_free_survival * 12) %>% 
        mutate(censor = ifelse(progression_free_survival_censor==0,1,2))
    c93984$censor = gsub(2, 0 , c93984$censor)
    require(survival)
    require(survminer)
    data.survdiff <- survdiff(Surv(progression_free_survival, progression_free_survival_censor) ~ status, data = c93984)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
    up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
    low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
    fit = survfit(Surv(progression_free_survival, progression_free_survival_censor) ~ status, data = c93984)
    p1 = ggsurvplot(fit,
        pval = TRUE, conf.int = FALSE,
        risk.table = FALSE, risk.table.col = "strata",
        linetype = "strata",
        legend = "right",
        legend.labs = c("CISD2 High", "CISD2 Low"),
        legend.title = "",
        title = "GSE93984",
        xlab = " Time (Month)",
        ylab = "Progression-Free Survival (%)",
        surv.median.line = "hv",
        palette = c("#E6550DFF", "#3182BDFF"),
        ggtheme = theme_bw())
    p1$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave("~/Documents/01CISD2/F_G93984OS.png", width = 6, height = 4)
#### Different expression analysis
    b31312 = g31312ex[grep('CISD2', g31312ex$symbol), ] %>% 
        gather() %>% 
        dplyr::rename('symbol' = 'key', 'CISD2' = 'value') %>% 
        slice(- 1) %>% 
        mutate_at('CISD2', as.numeric) %>% 
        dplyr::rename('geo_accession' = 'symbol') %>% 
        inner_join(g31312cx, by = 'geo_accession')
    b31312$COO = gsub('UC', 'UNC', b31312$COO)
    head(b31312)
    fit1 = b31312 %>% 
        ggplot(aes(fill = COO, y = CISD2, x = COO)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(data = b31312, aes(y = CISD2), size = 2, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('ABC', 'GCB'), c('GCB', 'UNC'), c('ABC', 'UNC')), size = 0.5, textsize = 5, vjust = 0, step_increase = 0.1, position = "identity") +
        labs(title = "", y = "The expression of CISD2", x = "") +
        scale_y_continuous(limits = c(0, 1)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#31A354FF")) + theme_bw() +
        theme(legend.position = "none",
              axis.text = element_text(size = 14, vjust = 0.5),
              axis.title = element_text(size = 16),
              plot.title = element_text(size = 20))
    ggsave("~/Documents/01CISD2/F_G31312Com.png", width = 5, height = 6)
    fit2 = c32918 %>% 
        ggplot(aes(y = CISD2, x = predictedclass, fill = predictedclass)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(data = c32918, aes(y = CISD2), size = 2, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        labs(title = "", y = "The expression of CISD2", x = "") +
        scale_y_continuous(limits = c(5, 15)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#31A354FF")) + 
        geom_signif(map_signif_level = TRUE, data = c32918,
            comparisons = list(c('ABC', 'GCB'), c('GCB', 'TypeIII'), c('ABC', 'TypeIII')), size = 0.5, textsize = 5, vjust = 0, step_increase = 0.12, stat = "signif", test = 't.test') +
        scale_x_discrete("", labels = c("ABC" = "ABC","GCB" = "GCB", "TypeIII" = "UNC")) + 
        theme_bw() +
        theme(legend.position = "none",
        axis.text = element_text(size = 14, vjust = 0.5),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 20))
    ggsave("~/Documents/01CISD2/F_G32918Com.png", width = 5, height = 6)
    d31312 = c31312 %>% 
             mutate(IPIscore = case_when(IPI == '0' ~ 'Low',
                                         IPI == '1' ~ 'Low',
                                         IPI == '2' ~ 'Medium-Low',
                                         IPI == '3' ~ 'Medium-High',
                                         IPI == '4' ~ 'High',
                                         IPI == '5' ~ 'High')) %>% 
            filter(IPI != 'UNK')
    fit = summary(aov(CISD2 ~ IPI, data = d31312))
    TukeyHSD(aov(CISD2 ~ IPI, data = d31312))
#### COX regression analysis.
    str(c31312)
    c31312$COO = gsub('UC', 'UNC', c31312$COO)
    cox0 = coxph(Surv(OSmonth, OS) ~ COO, c31312)
    cox1 = coxph(Surv(OSmonth, OS) ~ Age, c31312)
    cox2 = coxph(Surv(OSmonth, OS) ~ Gender, c31312)
    cox3 = coxph(Surv(OSmonth, OS) ~ LDH, subset(c31312, LDH != "UNK"))
    x31312 = c31312 %>% 
        filter(IPI != "UNK") %>% 
        mutate_at('IPI', as.factor)
    cox4 = coxph(Surv(OSmonth, OS) ~ IPI, x31312)
    cox5 = coxph(Surv(OSmonth, OS) ~ Bulky, subset(c31312, Bulky == 'nop' | Bulky == 'Bulky'))
    cox6 = coxph(Surv(OSmonth, OS) ~ AAS, subset(c31312, AAS != "UNK"))
    cox7 = coxph(Surv(OSmonth, OS) ~ extran, c31312)
    cox8 = coxph(Surv(OSmonth, OS) ~ ECOG, c31312)
    cox9 = coxph(Surv(OSmonth, OS) ~ CISD2, c31312)
    y31312 = x31312 %>% 
        filter(Bulky == 'nop' | Bulky == 'Bulky') %>% 
        filter(AAS != "UNK") %>% 
        filter(LDH != "UNK") %>% 
        filter(OSmonth >= 1)
    coxa = c("Age", "Gender", "COO", "extran", "ECOG")
    xa = data.frame()
    for(i in coxa) {
        expr = c31312[, i]
        cox = coxph(Surv(OSmonth, OS) ~ expr, c31312)
        coxsum = summary(cox)
        xa = rbind(xa, cbind(coxi = i, 
            HR = round(coxsum$coefficients[, 'exp(coef)'],2), 
            lower = round(coxsum$conf.int[, 3],2), 
            upper = round(coxsum$conf.int[, 4],2),
            'HR(95%CI)' = paste(round(coxsum$coefficients[, 'exp(coef)'],2), '(', round(coxsum$conf.int[, 3],2), '-', round(coxsum$conf.int[, 4],2), ')', sep = ''),
            pvalue = round(coxsum$coefficients[, 'Pr(>|z|)'],2),
            z = round(coxsum$coefficients[, "z"],2)))}    
    coxb = c("LDH")
    xb = data.frame()
    LDH31312 = c31312 %>% filter(LDH != 'UNK')
    for(i in coxb) {
        expr = LDH31312[, i]
        cox = coxph(Surv(OSmonth, OS) ~ expr, LDH31312)
        coxsum = summary(cox)
        xb = rbind(xb, cbind(coxi = i, 
            HR = round(coxsum$coefficients[, 'exp(coef)'],2), 
            lower = round(coxsum$conf.int[, 3],2), 
            upper = round(coxsum$conf.int[, 4],2),
            'HR(95%CI)' = paste(round(coxsum$coefficients[, 'exp(coef)'],2), '(', round(coxsum$conf.int[, 3],2), '-', round(coxsum$conf.int[, 4],2), ')', sep = ''),
            pvalue = round(coxsum$coefficients[, 'Pr(>|z|)'],2),
            z = round(coxsum$coefficients[, "z"],2)))}
    coxc = c("IPI")
    xc = data.frame()
    IPI31312 = c31312 %>% filter(IPI != 'UNK')
    for(i in coxc) {
        expr = IPI31312[, i]
        cox = coxph(Surv(OSmonth, OS) ~ expr, IPI31312)
        coxsum = summary(cox)
        xc = rbind(xc, cbind(coxi = i, 
            HR = round(coxsum$coefficients[, 'exp(coef)'],2), 
            lower = round(coxsum$conf.int[, 3],2), 
            upper = round(coxsum$conf.int[, 4],2),
            'HR(95%CI)' = paste(round(coxsum$coefficients[, 'exp(coef)'],2), '(', round(coxsum$conf.int[, 3],2), '-', round(coxsum$conf.int[, 4],2), ')', sep = ''),
            pvalue = round(coxsum$coefficients[, 'Pr(>|z|)'],2),
            z = round(coxsum$coefficients[, "z"],2)))}
    coxd = c("AAS")
    xd = data.frame()
    AAS31312 = c31312 %>% filter(AAS != 'UNK')
    for(i in coxd) {
        expr = AAS31312[, i]
        cox = coxph(Surv(OSmonth, OS) ~ expr, AAS31312)
        coxsum = summary(cox)
        xd = rbind(xd, cbind(coxi = i, 
            HR = round(coxsum$coefficients[, 'exp(coef)'],2), 
            lower = round(coxsum$conf.int[, 3],2), 
            upper = round(coxsum$conf.int[, 4],2),
            'HR(95%CI)' = paste(round(coxsum$coefficients[, 'exp(coef)'],2), '(', round(coxsum$conf.int[, 3],2), '-', round(coxsum$conf.int[, 4],2), ')', sep = ''),
            pvalue = round(coxsum$coefficients[, 'Pr(>|z|)'],2),
            z = round(coxsum$coefficients[, "z"],2)))}
    coxe = c("Bulky")
    xe = data.frame()
    Bulky31312 = c31312 %>% filter(Bulky == 'nop' | Bulky == 'Bulky')
    for(i in coxe) {
        expr = Bulky31312[, i]
        cox = coxph(Surv(OSmonth, OS) ~ expr, Bulky31312)
        coxsum = summary(cox)
        xe = rbind(xe, cbind(coxi = i, 
            HR = round(coxsum$coefficients[, 'exp(coef)'],2), 
            lower = round(coxsum$conf.int[, 3],2), 
            upper = round(coxsum$conf.int[, 4],2),
            'HR(95%CI)' = paste(round(coxsum$coefficients[, 'exp(coef)'],2), '(', round(coxsum$conf.int[, 3],2), '-', round(coxsum$conf.int[, 4],2), ')', sep = ''),
            pvalue = round(coxsum$coefficients[, 'Pr(>|z|)'],2),
            z = round(coxsum$coefficients[, "z"],2)))}
    coxf = c("CISD2", "status")
    xf = data.frame()
    for(i in coxf) {
        expr = c31312[, i]
        cox = coxph(Surv(OSmonth, OS) ~ expr, c31312)
        coxsum = summary(cox)
        xf = rbind(xf, cbind(coxi = i, 
            HR = round(coxsum$coefficients[, 'exp(coef)'],2), 
            lower = round(coxsum$conf.int[, 3],2), 
            upper = round(coxsum$conf.int[, 4],2),
            'HR(95%CI)' = paste(round(coxsum$coefficients[, 'exp(coef)'],2), '(', round(coxsum$conf.int[, 3],2), '-', round(coxsum$conf.int[, 4],2), ')', sep = ''),
            pvalue = round(coxsum$coefficients[, 'Pr(>|z|)'],2),
            z = round(coxsum$coefficients[, "z"],2)))}
    xf$pvalue = c('<0.05', '<0.05')  
    x = rbind(xa, xb, xc, xd, xe, xf) %>% 
        as_tibble() %>% 
        mutate_at(c('HR','lower','upper'), as.numeric) %>% 
        mutate(log_HR = log(HR)) %>% 
        mutate(log_lower = log(lower)) %>% 
        mutate(log_upper = log(upper)) %>% 
        mutate(iterm = factor(c('Age', 'Gender', 'ABC vs GCB', 'ABC vs UNC', 'extran', 'ECOG', 'LDH', 'IPI 0 vs 1', 'IPI 0 vs 2', 'IPI 0 vs 3', 'IPI 0 vs 4', 'IPI 0 vs 5', 'AAS', 'Bulky', 'CISD2', 'CISD2_d'), levels = c('Age', 'Gender', 'ABC vs GCB', 'ABC vs UNC', 'extran', 'ECOG', 'LDH', 'IPI 0 vs 1', 'IPI 0 vs 2', 'IPI 0 vs 3', 'IPI 0 vs 4', 'IPI 0 vs 5', 'AAS', 'Bulky', 'CISD2', 'CISD2_d'))) %>% 
        mutate(log_x = paste(round(log_HR, 2),'(',round(log_lower, 2),',', round(log_upper, 2),')', sep = '')) %>% 
        slice(-16) %>% 
        ggplot(aes(y = fct_rev(iterm))) + 
        geom_point(aes(x=log_HR), shape=15, size=4, color = '#3182BDFF') +
        geom_linerange(aes(xmin=log_lower, xmax=log_upper)) +
        geom_vline(xintercept = 0, linetype="dashed") +
        labs(x="Log Hazard Ratio", y="") +
        coord_cartesian(ylim=c(1,15), xlim=c(-0.5, 3)) +
        theme_classic() + 
        theme(axis.text.y = element_text(size = 16),
              axis.title.y= element_blank(),
              axis.ticks.y= element_blank(),
              axis.line.y = element_blank())
    ggsave('~/Documents/01CISD2/forrest1.png', width = 8, height = 8)
    y31312 = c31312 %>% 
        filter(IPI != "UNK") %>% 
        mutate_at('IPI', as.factor) %>% 
        filter(Bulky == 'nop' | Bulky == 'Bulky') %>% 
        filter(AAS != "UNK") %>% 
        filter(LDH != "UNK")
    sum_y = summary(coxph(Surv(OSmonth, OS) ~ Age+Gender+COO+LDH+IPI+AAS+Bulky+extran+ECOG+clinic+CISD2, y31312))
    ytable = cbind(HR = sum_y$coefficients[, 'exp(coef)'],
          lower = sum_y$conf.int[, 3],
          upper = sum_y$conf.int[, 4],
          pvalue = sum_y$coefficients[, 'Pr(>|z|)'],
          z = sum_y$coefficients[, "z"]) %>% 
          as_tibble() %>% 
          mutate(log_HR = log(HR)) %>% 
          mutate(log_lower = log(lower)) %>% 
          mutate(log_upper = log(upper)) %>% 
          mutate(iterm = factor(c('Age', 'Gender', 'ABC vs GCB', 'ABC vs UNC', 'LDH', 'IPI 0 vs 1', 'IPI 0 vs 2', 'IPI 0 vs 3', 'IPI 0 vs 4', 'IPI 0 vs 5', 'AAS', 'Bulky', 'extran', 'ECOG', 'CISD2'), levels = c('Age', 'Gender', 'ABC vs GCB', 'ABC vs UNC', 'extran', 'ECOG', 'LDH', 'IPI 0 vs 1', 'IPI 0 vs 2', 'IPI 0 vs 3', 'IPI 0 vs 4', 'IPI 0 vs 5', 'AAS', 'Bulky', 'CISD2'))) %>% 
          mutate(log_x = paste(round(log_HR, 2),'(',round(log_lower, 2),',', round(log_upper, 2),')', sep = ''))
    y = ytable  %>% 
        ggplot(aes(y = fct_rev(iterm))) + 
        geom_point(aes(x=log_HR), shape=16, size=4, color = "#E6550DFF") +
        geom_linerange(aes(xmin=log_lower, xmax=log_upper)) +
        geom_vline(xintercept = 0, linetype="dashed") +
        labs(x="Log Hazard Ratio", y="") +
        coord_cartesian(ylim=c(1,15), xlim=c(-0.5, 3)) +
        theme_classic() + 
        theme(axis.text.y = element_text(size = 16),
              axis.title.y= element_blank(),
              axis.ticks.y= element_blank(),
              axis.line.y = element_blank())
    y
    ggsave('~/Documents/01CISD2/forrest2.png', width = 8, height = 8)
## Figure 4 LASSO Cox and Enrichment plot
#### GSE117556 based on GPL14951 (Training dataset)
    G117556 <- getGEO(filename = ('~/Documents/01CISD2/GSE117556_series_matrix.txt.gz'), GSEMatrix = TRUE, getGPL= FALSE)
    if (length(G117556) > 1) {idx = grep("GPL14951", attr(G117556, "names"))} else {idx = 1}
    g117556e = data.frame(exprs(G117556))
    g117556c = pData(G117556)
    qx = as.numeric(quantile(g117556e, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC = (qx[5] > 100) ||
      (qx[6]-qx[1] > 50 && qx[2] > 0) ||
      (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if(LogC){g117556e[which(g117556e <= 0)] = NaN
            g117556e = log2(g117556e)
            print("log2 transform finished")}else{print("log2 transform not needed")}
    an14951 = data.table::fread("~/Documents/01CISD2/GPL14951-11332.txt") %>% 
        select("ID", 'symbol'="Symbol")
    g117556ex = g117556e %>% 
        rownames_to_column('ID') %>%
        merge(an14951, by = 'ID') %>%
        dplyr::select(- 'ID') %>%
        dplyr::select(symbol, everything()) %>% 
        mutate(rowMean = rowMeans(.[grep('GSM', names(.))])) %>%
        arrange(desc(rowMean)) %>%
        distinct(symbol, .keep_all = T) %>%
        dplyr::select(-'rowMean') 
    colnames(g117556c) = gsub(':ch1', '', colnames(g117556c))
    colnames(g117556c) = gsub(' ', '', colnames(g117556c))
    g117556cx = read_xlsx('~/Documents/01CISD2/DS_JCO.18.01314-2.xlsx', sheet = 1) %>% 
        dplyr::rename('title' = colnames(.)[1]) %>% 
        inner_join(g117556c, by = 'title')
    c117556 = g117556ex[grep('CISD2', g117556ex$symbol), ] %>% 
        gather() %>% 
        dplyr::rename('symbol' = 'key', 'CISD2' = 'value') %>% 
        slice(- 1 ) %>% 
        mutate_at('CISD2', as.numeric) %>% 
        dplyr::rename('geo_accession' = 'symbol') %>% 
        inner_join(g117556cx, by = 'geo_accession') %>% 
        filter(RESP_ASSESS != 0) %>% 
        filter(RESP_ASSESS != 'Not Evaluable') %>% 
        filter(FU_TIME >= 1)
    fit = survfit(Surv(FU_TIME, DEATH_IND) ~ states, c117556) %>% 
          ggsurvplot(.,
            pval = T, conf.int = FALSE, 
            risk.table = T, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2 High", "CISD2 Low"),
            legend.title = "",
            title = "GSE117556", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw()) 
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G117556surALL.png', width = 6, height = 4)
    abc117556 = c117556 %>% 
          filter(molecularsubtype == 'ABC') 
    fit = survfit(Surv(FU_TIME, DEATH_IND) ~ states, abc117556) %>% 
          ggsurvplot(.,
            pval = T, conf.int = FALSE, 
            risk.table = T, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2 High", "CISD2 Low"),
            legend.title= "",
            title= "GSE117556 (ABC subtype)", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G117556surABC.png', width = 6, height = 4)
    gcb117556 = c117556 %>% 
          filter(molecularsubtype == 'GCB') 
    fit = survfit(Surv(FU_TIME, DEATH_IND) ~ states, gcb117556) %>% 
          ggsurvplot(.,
            pval = T, conf.int = FALSE, 
            risk.table = T, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2 High", "CISD2 Low"),
            legend.title= "",
            title= "GSE117556 (GCB subtype)", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G117556surGCB.png', width = 6, height = 4)
    mhg117556 = c117556 %>% 
          filter(molecularsubtype == 'MHG') 
    fit = survfit(Surv(FU_TIME, DEATH_IND) ~ states, mhg117556) %>% 
          ggsurvplot(.,
            pval = T, conf.int = FALSE, 
            risk.table = T, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2 High", "CISD2 Low"),
            legend.title= "",
            title= "GSE117556 (MHG subtype)", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G117556surMHG.png', width = 6, height = 4)
    unc117556 = c117556 %>% 
          filter(molecularsubtype == 'UNC') 
    fit = survfit(Surv(FU_TIME, DEATH_IND) ~ states, unc117556) %>% 
          ggsurvplot(.,
            pval = T, conf.int = FALSE, 
            risk.table = T, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2 High", "CISD2 Low"),
            legend.title= "",
            title= "GSE117556 (UNC subtype)", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G117556surUNC.png', width = 6, height = 4)
    str(c117556)
    fit1 = c117556 %>% 
        ggplot(aes(fill = molecularsubtype, y = CISD2, x = molecularsubtype)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(data = c117556, aes(y = CISD2), size = 2, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('ABC', 'GCB'), c('GCB', 'UNC'), c('ABC', 'UNC'), c('ABC', 'MHG')), size = 0.5, textsize = 5, vjust = 0, step_increase = 0.15, position = "identity") +
        labs(title = "", y = "The expression of CISD2", x = "") +
        scale_y_continuous(limits = c(5, 18)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#31A354FF", '#7876B1FF')) + theme_bw() +
        theme(legend.position = "none",
              axis.text = element_text(size = 14, vjust = 0.5),
              axis.title = element_text(size = 16),
              plot.title = element_text(size = 20))
    ggsave("~/Documents/01CISD2/F_G117556Com.png", width = 5, height = 6)
    c117556 = c117556 %>% 
        mutate(RESP = case_when(RESP_ASSESS == 'CR' ~ 'CR', RESP_ASSESS == 'CRu' ~ 'CR', RESP_ASSESS == 'PD' ~ 'PD', RESP_ASSESS == 'SD' ~ 'SD', RESP_ASSESS == 'PR' ~ 'PR')) %>% 
        mutate(RESX = case_when(RESP == 'CR' ~ '1', RESP == 'PR' ~ '1', RESP == 'PD' ~ '0', RESP == 'SD' ~ '0')) 
    fit1 = c117556 %>% 
        ggplot(aes(fill = RESX, y = CISD2, x = RESX)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(data = c117556, aes(y = CISD2), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('1', '0')), size = 0.5, textsize = 5, vjust = 0, step_increase = 0.15, position = "identity") +
        labs(title = "", y = "The expression of CISD2", x = "") +
        scale_x_discrete("", labels = c("0" = "Ineffective","1" = 'Effective')) +
        scale_y_continuous(limits = c(5, 17)) +
        scale_fill_manual(values = c("#E6550DFF", "#3182BDFF")) + theme_bw() +
        theme(legend.position = "none",
              axis.text = element_text(size = 14, vjust = 0.5),
              axis.title = element_text(size = 16),
              plot.title = element_text(size = 20))
    ggsave('~/Documents/01CISD2/F_G117556effective.png', width = 3, height = 4)
    c117556$RESP = factor(c117556$RESP, levels = c('CR', 'PR', 'SD', 'PD'))
    fit1 = c117556 %>% 
        ggplot(aes(fill = RESP, y = CISD2, x = RESP)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(data = c117556, aes(y = CISD2), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('PR', 'PD'), c('PR', 'SD'), c('CR', 'SD'), c('CR', 'PD')), size = 0.5, 
            textsize = 5, vjust = 0, step_increase = 0.2, position = "identity") +
        labs(title = "", y = "The expression of CISD2", x = "") +
        #scale_x_discrete("", labels = c("0" = "Ineffective","1" = 'Effective')) +
        scale_y_continuous(limits = c(5, 20)) +
        scale_fill_manual(values = c("#E6550DFF", "#3182BDFF", "#7876B1FF", "#20854EFF")) + theme_bw() +
        theme(legend.position = "none",
              axis.text = element_text(size = 14, vjust = 0.5),
              axis.title = element_text(size = 16),
              plot.title = element_text(size = 20))
    ggsave('~/Documents/01CISD2/F_G117556effective1.png', width = 4, height = 4)
    c117556a = c117556 %>% 
        mutate_at('IPI_SCORE', as.factor) %>% 
        mutate(IPI = case_when(IPI_SCORE == '0' ~ 'Low',
                               IPI_SCORE == '1' ~ 'Low',
                               IPI_SCORE == '2' ~ 'Moderate-Low',
                               IPI_SCORE == '3' ~ 'Moderate-High',
                               IPI_SCORE == '4' ~ 'High',
                               IPI_SCORE == '5' ~ 'High'
                               ))
    c117556a$IPI = factor(c117556a$IPI, levels = c('Low', 'Moderate-Low', 'Moderate-High', 'High'))
    fit1 = c117556a %>% 
        ggplot(aes(fill = IPI, y = CISD2, x = IPI)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(data= c117556a, aes(y = CISD2), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(annotations = c("NS","**","NS"), size = 0.5, textsize = 5, vjust = 0,
        y_position = c(15, 17, 19), xmin = c(1, 1, 1), xmax = c(2, 3, 4)) +
        labs(title = "", y = "The expression of CISD2", x = "") +
        scale_y_continuous(limits = c(5, 20)) +
        scale_fill_manual(values = c("#E6550DFF", "#3182BDFF", "#7876B1FF", "#20854EFF", "#6F99ADFF", "#FFDC91FF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 12, vjust = 0.5, angle = 45),
              axis.title = element_text(size = 16),
              plot.title = element_text(size = 20))
    ggsave("~/Documents/01CISD2/F_G117556IPI.png", width = 4, height = 4.5)
    c117556b = c117556 %>% 
        select(c('RESP', 'states'))
    a1 = count(c117556b %>% filter(states == 'low'))
    a2 = count(c117556b %>% filter(states == 'high'))
    a1$RESP = factor(a1$RESP, levels = c('CR', 'PR', 'PD', 'SD'))
    a2$RESP = factor(a1$RESP, levels = c('CR', 'PR', 'PD', 'SD'))
    chisq.test(a1$freq, a2$freq)
    a = rbind(a1, a2) %>% as_tibble() %>% 
        ggplot(aes(fill=RESP, y= freq, x=states)) +
        geom_bar(position="fill", stat="identity") +
        scale_fill_manual(values = c("#E6550DFF", "#3182BDFF", "#7876B1FF", "#20854EFF", "#6F99ADFF", "#FFDC91FF")) + 
        theme_bw() +
        labs(title = "", y = "The Ratio of treatment effective", x = "") +
        scale_y_continuous(labels = scales::percent) +
        guides(fill=guide_legend(title="RESPONSE")) +
        scale_x_discrete("", labels = c("high" = "CISD2 High", "low" = 'CISD2 Low')) + 
        theme(axis.text.x = element_text(size = 12, vjust = 0.5, angle = 45),
              axis.text.y = element_text(size = 12),
              axis.title = element_text(size = 16),
              plot.title = element_text(size = 20))
    ggsave("~/Documents/01CISD2/F_G117556CR.png", width = 4, height = 4.5)
    fit = survfit(Surv(FU_TIME, DEATH_IND) ~ states, a) %>% 
        ggsurvplot(.,
            pval = T, conf.int = FALSE, 
            risk.table = T, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2 High", "CISD2 Low"),
            legend.title= "",
            title= "GSE117556 (RESPONSE = CR)", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G117556surCR.png', width = 6, height = 4)
#### LASSO-COX analysis.
    require(glmnet)
    CISD2_related = read.csv('~/Documents/01CISD2/string_interactions.csv', header = T, sep = ",")
    c_gene = unique(CISD2_related$node1)
    head(g117556ex$symbol)
    x0 = unique(intersect(g117556ex$symbol, c_gene))
    setdiff(c_gene, x0)
    y0 = c('KIAA0831', 'FXC1', 'C19orf52', 'UCRC', 'TMEM49')
    z0 = c(x0, y0)
    x1 = data.frame()
    for(i in z0){
        xi = g117556ex %>% 
            filter(symbol == i)
        x1 = rbind(x1, xi)}
    colnames(x1)
    z1 = c(7949)
    x2 = data.frame()
    for(i in z1){
       xi = g117556ex[i, ]
       x2 = rbind(x2, xi)}
    x_lasso = rbind(x1, x2) %>% 
        select(c(symbol, c117556$geo_accession)) %>% 
        t() %>% 
        as.data.frame()
    colnames(x_lasso) = x_lasso[1, ]
    x_lasso = x_lasso %>% 
        slice(-1) %>% 
        mutate_all(as.numeric) %>% 
        as.matrix()
    str(x_lasso)
    y = data.matrix(Surv(c117556$FU_TIME, c117556$DEATH_IND))
    x = x_lasso
    str(x)
    fit = glmnet(x,y, family = "cox", alpha = 1)
    fitcv = cv.glmnet(x,y,family="cox", alpha=1,nfolds=10)
    coef.min = coef(fitcv, s = "lambda.min")
    active.min = which(coef.min@i != 0)
    lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
    png('~/Documents/01CISD2/LASSO5.png', width = 2000, heigh = 1600, res = 300)
    plot(fit,xvar="lambda",label=T)
    dev.off()   
    png('~/Documents/01CISD2/LASSO6.png', width = 2000, heigh = 1600, res = 300)
    fitcv = cv.glmnet(x,y,family="cox", alpha=1,nfolds=10)
    plot(fitcv)
    coef(fitcv, s="lambda.min")
    coef(fitcv, s="lambda.1se")
    predict(fit, newx = nx, s = c(0.1, 0.05))
    fitcv$index
    coef(fit, s = 0.05)
    dev.off()
    cisd2_related_gene = c('TIMM10','NDUFB2','BCL2L1','BCL2L11','CISD2','BID','NDUFA9','NDUFS5','NDUFB9','BCL2','NDUFA7','MCL1','PMAIP1','PIK3C3','CYCS','UQCRB','NDUFS1','HRK','UVRAG','BBC3','PIK3R4','CISD3','NDUFB1','NRBF2','NDUFB4','FXC1','TMEM49')
    riskscore = TIMM10 * -0.05249048 + NDUFB2 * -0.12915551 + BCL2L1 * -0.11683190 + BCL2L11 * -0.05368164 + CISD2  * 0.24852093 + BID * -0.06101778 + NDUFA9 * -0.02235329 + NDUFS5 * 0.01761829 + NDUFB9 * 0.24596655 + BCL2  * 0.10566931 + NDUFA7 * 0.04073703 + MCL1  * 0.26730071 + PMAIP1 * 0.12489498 + PIK3C3 * 0.15104140 + CYCS  * 0.14281014 + UQCRB  * 0.14679661 + NDUFS1 * 0.05981878 + HRK * 0.03899325 + UVRAG * -0.26143078 + BBC3 * -0.04408451 + PIK3R4 * 0.23802877 + CISD3 * 0.12013327 + NDUFB1 * -0.11920934 + NRBF2 * 0.07456076 + NDUFB4 * -0.14108068 + FXC1  * 0.03498904 + TMEM49 * -0.3655337 
    x117556 = x_lasso %>% 
        as.data.frame() %>% 
        mutate(risk = TIMM10 * -0.05249048 + NDUFB2 * -0.12915551 + BCL2L1 * -0.11683190 + BCL2L11 * -0.05368164 + CISD2 * 0.24852093 + BID * -0.06101778 + NDUFA9 * -0.02235329 + NDUFS5 * 0.01761829 + NDUFB9 * 0.24596655 + BCL2  * 0.10566931 + NDUFA7 * 0.04073703 + MCL1  * 0.26730071 + PMAIP1 * 0.12489498 + PIK3C3 * 0.15104140 + CYCS  * 0.14281014 + UQCRB  * 0.14679661 + NDUFS1 * 0.05981878 + HRK * 0.03899325 + UVRAG * -0.26143078 + BBC3 * -0.04408451 + PIK3R4 * 0.23802877 + CISD3 * 0.12013327 + NDUFB1 * -0.11920934 + NRBF2 * 0.07456076 + NDUFB4 * -0.14108068 + FXC1  * 0.03498904 + TMEM49 * -0.3655337) %>% 
        rownames_to_column('geo_accession') %>% 
        inner_join(c117556, by = 'geo_accession') %>% 
        dplyr::rename('CISD2' = 'CISD2.x') %>% 
        mutate(state = ifelse(risk < median(risk), 'low', 'high')) %>% 
        mutate(RESP = case_when(RESP_ASSESS == 'CR' ~ 'CR', RESP_ASSESS == 'CRu' ~ 'CR', RESP_ASSESS == 'PD' ~ 'PD', RESP_ASSESS == 'SD' ~ 'SD', RESP_ASSESS == 'PR' ~ 'PR')) %>% 
        mutate(RESX = case_when(RESP == 'CR' ~ '1', RESP == 'PR' ~ '1', RESP == 'PD' ~ '0', RESP == 'SD' ~ '0'))  %>% 
        mutate(doubleexpression = ifelse(expressor_RNA == 'double-expressor', '1', '0'))
#### GO and KEGG plot
    go1 = read.table(file = "~/Documents/01CISD2/GOBP.txt",header = T,sep = "\t") %>% 
        arrange(., PValue)
    go1 = go1[1:10,]
    go2 = read.table(file = "~/Documents/01CISD2/GOCC.txt",header = T,sep = "\t") %>% 
       arrange(PValue)
    go2 = go2[1:10,]
    go3 = read.table(file = "~/Documents/01CISD2/GOMF.txt",header = T,sep = "\t") %>% 
       arrange(PValue)
    go3 = go3[1:8,]
    go = rbind(go1, go2, go3) %>% 
        as_tibble() %>% 
        mutate(goTerm = c(rep('Biological Process', 10), rep('Cellular Component', 10), rep('Molecular Function', 8))) %>% 
        mutate(go_term_order = factor(as.integer(rownames(.)), labels = Term)) %>% 
        ggplot(aes(x = go_term_order, y = Count, fill = goTerm)) +
        geom_bar(stat='identity') +
        labs(title = "The Most Enriched GO Terms", x = "", y = "The Count of Enrichment", fill = 'GO Terms') + 
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() + coord_flip() +
        guides(fill=guide_legend(reverse=TRUE)) + 
        theme(legend.position = "right",
              axis.text.x = element_text(size = 10, angle = 0, hjust = 1, vjust = 1),
              axis.title = element_text(size = 16),
              plot.title = element_text(size = 20))
    ggsave('~/Documents/01CISD2/F_GO.png', width = 12, height = 10)
    kegg = read.table(file = "~/Documents/01CISD2/KEGG.txt",header = T,sep = "\t") %>% 
        arrange(., PValue) %>%  
        as_tibble() 
    kegg = kegg[1:20, ] %>% 
        ggplot(aes(x = Count, y = reorder(Term, Count), fill=-log(PValue))) +
        geom_bar(stat='identity') +
        labs(title = "The Most Enriched KEGG Terms", x = "The Count of Enrichment", y = "") + 
        scale_fill_gradient(expression(-log["10"](Pvalue)),low="#3182BDFF",high="#E6550DFF") + 
        theme_bw() + 
        #guides(fill=guide_legend(reverse=TRUE)) + 
        theme(legend.position = "right",
              axis.text.x = element_text(size = 10, angle = 0, hjust = 1, vjust = 1),
              axis.title = element_text(size = 16),
              plot.title = element_text(size = 20))
    ggsave('~/Documents/01CISD2/F_KEGG.png', width = 12, height = 10)
## Figure 5 Association Between CISD2Risk and Clinical Features
#### Corretion with CISD2Risk based on GSE117556
    head(x117556)
    fit = survfit(Surv(FU_TIME, DEATH_IND) ~ state, data = x117556) %>% 
            ggsurvplot(.,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2Risk High", "CISD2Rsik Low"),
            legend.title = "",
            title = "GSE117556", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G117556R1.png', width = 6, height = 4)
    fit = survfit(Surv(PROG_TIME, PROG_IND) ~ state, data = x117556) %>% 
            ggsurvplot(.,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2Risk High", "CISD2Rsik Low"),
            legend.title = "",
            title = "GSE117556", 
            xlab = " Time (Month)",
            ylab = 'Progression-free Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G117556R100.png', width = 6, height = 4)
    fit2 = x117556 %>% 
        ggplot(aes(fill = molecularsubtype, y = risk, x = molecularsubtype)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(data = x117556, aes(y = risk), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('ABC', 'GCB'), c('ABC', 'MHG'), c('ABC', 'UNC'), c('GCB', 'MHG'), c('GCB', 'UNC'), c('UNC', 'MHG')), 
            size = 0.5, textsize = 5, vjust = 0.5, step_increase = 0.1, position = "identity") +
        labs(title = "GSE117556", y = "CISD2Risk value", x = "") +
        #scale_y_continuous(limits = c(3, 11)) +
        scale_fill_manual(values = c("#E6550DFF", "#3182BDFF", "#31A354FF", '#7876B1FF')) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text = element_text(size = 14, vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave("~/Documents/01CISD2/F_G117556R2.png", width = 3, height = 4.5)
    fit3 = x117556 %>% 
        ggplot(aes(fill = RESX, y = risk, x = RESX)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(data = x117556, aes(y = risk), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('1', '0')), size = 0.5, textsize = 5, vjust = 0.5, step_increase = 0.15, position = "identity") +
        labs(title = "GSE117556", y = "CISD2Risk value", x = "") +
        scale_x_discrete("", labels = c("0" = "Clinical Ineffectiveness","1" = 'Clinical Effectiveness')) +
        #scale_y_continuous(limits = c(3, 11)) +
        scale_fill_manual(values = c("#E6550DFF", "#3182BDFF")) + theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 14, angle = 60, hjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G117556R3.png', width = 2.5, height = 7.5)
    x117556$RESP = factor(x117556$RESP, levels = c('CR', 'PR', 'SD', 'PD'))
    fit4 = x117556 %>% 
        ggplot(aes(fill = RESP, y = risk, x = RESP)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(data = x117556, aes(y = risk), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('PR', 'PD'), c('CR', 'PD')), size = 0.5, 
            textsize = 5, vjust = 0.2, step_increase = 0.15, position = "identity") +
        labs(title = "GSE117556", y = "CISD2Risk value", x = "") +
        #scale_x_discrete("", labels = c("0" = "Ineffective","1" = 'Effective')) +
        scale_y_continuous(limits = c(3, 11)) +
        scale_fill_manual(values = c("#20854EFF", "#7876B1FF", "#3182BDFF", "#E6550DFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text = element_text(size = 14, vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G117556R4.png', width = 3, height = 4.5)
    fit5 = x117556 %>% 
        ggplot(aes(fill = MYC_RNA, y = risk, x = MYC_RNA)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(data = x117556, aes(y = risk), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('high', 'average')), size = 0.5, 
            textsize = 5, vjust = 0.2, step_increase = 0.15, position = "identity") +
        labs(title = "GSE117556", y = "CISD2Risk value", x = "") +
        scale_x_discrete("", labels = c("high" = "High Myc RNA","average" = 'Average Myc RNA')) +
        scale_y_continuous(limits = c(3, 11)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G117556R5.png', width = 2.5, height = 5)
    fit6 = x117556 %>% 
        ggplot(aes(fill = BCL2_RNA, y = risk, x = BCL2_RNA)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(data = x117556, aes(y = risk), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('high', 'average')), size = 0.5, 
            textsize = 5, vjust = 0, step_increase = 0.15, position = "identity") +
        labs(title = "GSE117556", y = "CISD2Risk value", x = "") +
        scale_x_discrete("", labels = c("high" = "High BCL2 RNA","average" = 'Average BCL2 RNA')) +
        scale_y_continuous(limits = c(3, 11)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G117556R6.png', width = 2.5, height = 5)
    fit7 = x117556 %>% 
        ggplot(aes(fill = doubleexpression, y = risk, x = doubleexpression)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(data = x117556, aes(y = risk), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('0', '1')), size = 0.5, 
            textsize = 5, vjust = 0, step_increase = 0.15, position = "identity") +
        labs(title = "GSE117556", y = "CISD2Risk value", x = "") +
        scale_x_discrete("", labels = c("0" = "non-double-expressor", '1' = 'double-expressor')) +
        scale_y_continuous(limits = c(3, 11)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G117556R7.png', width = 2.5, height = 5)
    fit8 = x117556 %>% 
        filter(STAGE != 0) %>% 
        mutate(xyz = case_when(STAGE == 'Stage I' ~ '0', STAGE == 'Stage II' ~ '0',
                               STAGE == 'Stage III' ~ '1', STAGE == 'Stage IV' ~ '1')) %>% 
        ggplot(aes(fill = xyz, y = risk, x = xyz)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(aes(y = risk), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('0', '1')), size = 0.5, 
            textsize = 5, vjust = 0, step_increase = 0.15, position = "identity") +
        labs(title = "GSE117556", y = "CISD2Risk value", x = "") +
        scale_x_discrete("", labels = c("0" = "Stage I~II", '1' = 'Stage III~IV')) +
        scale_y_continuous(limits = c(3, 11)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text = element_text(size = 10, vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G117556R8.png', width = 2.5, height = 4.5)
    fit9 = x117556 %>% 
        mutate(xyz = ifelse(IPI_SCORE <= 2, '0', '1')) %>% 
        ggplot(aes(fill = xyz, y = risk, x = xyz)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(aes(y = risk), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('0', '1')), size = 0.5, 
            textsize = 5, vjust = 0, step_increase = 0.15, position = "identity") +
        labs(title = "GSE117556", y = "CISD2Risk value", x = "") +
        scale_x_discrete("", labels = c("0" = "Low IPI", '1' = 'High IPI')) +
        scale_y_continuous(limits = c(3, 11)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text = element_text(size = 10, vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G117556R9.png', width = 2.5, height = 4.5)
    fit10 = x117556 %>% 
        mutate(AGEx = ifelse(AGE <60, '0', '1')) %>% 
        ggplot(aes(fill = AGEx, y = risk, x = AGEx)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(aes(y = risk), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('0', '1')), size = 0.5, 
            textsize = 5, vjust = 0, step_increase = 0.15, position = "identity") +
        labs(title = "GSE117556", y = "CISD2Risk value", x = "") +
        scale_x_discrete("", labels = c("0" = "Under 60 years old", '1' = 'Above 60 years old')) +
        scale_y_continuous(limits = c(3, 11)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G117556R10.png', width = 2.5, height = 4.5)
    fit11 = x117556 %>% 
        mutate(ECOG1 = ifelse(ECOG >= 2, '1', '0')) %>% 
        mutate_at('ECOG1', as.factor) %>% 
        ggplot(aes(fill = ECOG1, y = risk, x = ECOG1)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(aes(y = risk), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('0', '1')), size = 0.5, 
            textsize = 5, vjust = 0, step_increase = 0.15, position = "identity") +
        labs(title = "GSE117556", y = "CISD2Risk value", x = "") +
        scale_x_discrete("", labels = c("0" = "ECOG < 2", '1' = 'ECOG >=2')) +
        scale_y_continuous(limits = c(3, 12)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text = element_text(size = 10, vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G117556R11.png', width = 2.5, height = 4.5)
    fit12 = x117556 %>% 
        ggplot(aes(y = risk, x = MAX_TUMOUR)) +
        geom_point(size=0.5, shape=16) + 
        geom_smooth(method="lm", se = T) + 
        labs(title = "GSE117556", y = "CISD2Risk value", x = "The size of Tumor (mm)") +
        # scale_x_discrete("", labels = c("0" = "ECOG 0", '1' = 'ECOG 1', '2' = 'ECOG 2')) +
        # scale_y_continuous(limits = c(3, 12)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        stat_cor(method = 'pearson', label.x = 50) +
        theme(legend.position = "none",
              axis.text = element_text(size = 10, vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G117556R12.png', width = 3, height = 4.5)
    fit13 = x117556 %>% 
        mutate(LDH1 = ifelse(LDH < 245, '0', '1')) %>% 
        ggplot(aes(fill = LDH1, y = risk, x = LDH1)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(aes(y = risk), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('0', '1')), size = 0.5, 
            textsize = 5, vjust = 0, step_increase = 0.15, position = "identity") +
        labs(title = "GSE117556", y = "CISD2Risk value", x = "") +
        scale_x_discrete("", labels = c("0" = "normal LDH", '1' = 'raised LDH')) +
        scale_y_continuous(limits = c(3, 11)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10,vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G117556R13.png', width = 2.5, height = 4.5)
    fit15 = x117556 %>% 
        mutate_at('INV_EXTRANODAL_BAS_IND', as.factor) %>% 
        ggplot(aes(fill = INV_EXTRANODAL_BAS_IND, y = risk, x = INV_EXTRANODAL_BAS_IND)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(aes(y = risk), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('0', '1')), size = 0.5, 
            textsize = 5, vjust = 0, step_increase = 0.15, position = "identity") +
        labs(title = "GSE117556", y = "CISD2Risk value", x = "") +
        scale_x_discrete("", labels = c("0" = "non Extranodal", '1' = 'Extranodal')) +
        scale_y_continuous(limits = c(3, 11)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10,vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G117556R15.png', width = 2.5, height = 4.5)
    fit = survfit(Surv(FU_TIME, DEATH_IND) ~ state, data = subset(x117556, molecularsubtype == 'ABC')) %>% 
            ggsurvplot(.,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2Risk High", "CISD2Rsik Low"),
            legend.title = "",
            title = "GSE117556 (ABC Subtype)", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G117556R101.png', width = 6, height = 4)
    fit = survfit(Surv(PROG_TIME, PROG_IND) ~ state,  data = subset(x117556, molecularsubtype == 'ABC')) %>% 
            ggsurvplot(.,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2Risk High", "CISD2Rsik Low"),
            legend.title = "",
            title = "GSE117556 (ABC Subtype)", 
            xlab = " Time (Month)",
            ylab = 'Progression-free Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G117556R102.png', width = 6, height = 4)
    fit = survfit(Surv(FU_TIME, DEATH_IND) ~ state, data = subset(x117556, molecularsubtype == 'GCB')) %>% 
            ggsurvplot(.,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2Risk High", "CISD2Rsik Low"),
            legend.title = "",
            title = "GSE117556 (GCB Subtype)", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G117556R103.png', width = 6, height = 4)
    fit = survfit(Surv(PROG_TIME, PROG_IND) ~ state,  data = subset(x117556, molecularsubtype == 'GCB')) %>% 
            ggsurvplot(.,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2Risk High", "CISD2Rsik Low"),
            legend.title = "",
            title = "GSE117556 (GCB Subtype)", 
            xlab = " Time (Month)",
            ylab = 'Progression-free Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G117556R104.png', width = 6, height = 4)
    fit = survfit(Surv(FU_TIME, DEATH_IND) ~ state, data = subset(x117556, molecularsubtype == 'MHG')) %>% 
            ggsurvplot(.,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2Risk High", "CISD2Rsik Low"),
            legend.title = "",
            title = "GSE117556 (MHG Subtype)", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G117556R105.png', width = 6, height = 4)
    fit = survfit(Surv(PROG_TIME, PROG_IND) ~ state,  data = subset(x117556, molecularsubtype == 'MHG')) %>% 
            ggsurvplot(.,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2Risk High", "CISD2Rsik Low"),
            legend.title = "",
            title = "GSE117556 (MHG Subtype)", 
            xlab = " Time (Month)",
            ylab = 'Progression-free Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G117556R106.png', width = 6, height = 4)
    fit = x117556 %>% 
        select(c('expressor_RNA', 'MYC_RNA', 'BCL2_RNA', 'state', 'risk', 'geo_accession')) %>% 
        mutate(BCL2_expression = ifelse(BCL2_RNA == 'average', 'average BCL2', 'High BCL2')) %>% 
        mutate(Myc_expression = ifelse(MYC_RNA == 'average', 'average Myc', 'High Myc')) %>% 
        select(c('expressor_RNA', 'BCL2_expression', 'Myc_expression', 'state', 'risk', 'geo_accession')) %>% 
        pivot_longer(cols = c('expressor_RNA', 'BCL2_expression', 'Myc_expression'))
    fit$value = factor(fit$value, levels = c('average BCL2', 'High BCL2', 'average Myc', 'High Myc', "non-double-expressor", "double-expressor"))
    fit0 = fit  %>% 
        ggplot(aes(fill = value, y = risk, x = value)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(aes(y = risk), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        labs(title = "GSE117556", y = "CISD2Risk value", x = '') +
        #scale_x_discrete("", labels = c("0" = "Low IPI", '1' = 'High IPI')) +
        scale_y_continuous(limits = c(3, 11)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#3182BDFF", "#E6550DFF", "#3182BDFF", "#E6550DFF")) + 
        theme_bw() + 
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('average BCL2', 'High BCL2'), c('average Myc', 'High Myc'), 
            c("non-double-expressor", "double-expressor")), size = 0.5, 
            textsize = 5, vjust = 0, step_increase = 0, position = "identity") +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10, angle = 60, hjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
        ggsave('~/Documents/01CISD2/F_G117556Ra.png', width = 4, height = 6.5)
    fit = survfit(Surv(FU_TIME, DEATH_IND) ~ molecularsubtype, data = subset(x117556, state == 'high')) %>% 
            ggsurvplot(.,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("ABC subtype", "GCB subtype", "MHG subtype", "UNC subtype"),
            legend.title = "",
            title = "GSE117556 - High CISD2Risk value", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF", "#31A354FF", '#7876B1FF'),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.85)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G117556Raa4.png', width = 6, height = 4)
#### Correction analysis based on GSE181063 (Validation dataset)
    G181063 = getGEO(filename = ('~/Documents/01CISD2/GSE181063_series_matrix.txt.gz'), GSEMatrix = T, getGPL = F)
    if (length(G181063) > 1) {idx = grep("GPL14951", attr(G181063, "names"))} else {idx = 1}
    g181063e = data.frame(exprs(G181063))
    g181063c = pData(G181063)
    qx = as.numeric(quantile(g181063e, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC = (qx[5] > 100) ||
      (qx[6]-qx[1] > 50 && qx[2] > 0) ||
      (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if(LogC){g181063e[which(g181063e <= 0)] = NaN
            g1891063e = log2(g181063e)
            print("log2 transform finished")}else{print("log2 transform not needed")}
    an14951 = data.table::fread("~/Documents/01CISD2/GPL14951-11332.txt") %>% 
        select("ID", 'symbol'="Symbol")
    g181063ex = g181063e %>% 
        rownames_to_column('ID') %>%
        merge(an14951, by = 'ID') %>%
        dplyr::select(- 'ID')  %>%
        dplyr::select(symbol, everything()) %>% 
        mutate(rowMean = rowMeans(.[grep('GSM', names(.))])) %>%
        arrange(desc(rowMean)) %>%
        distinct(symbol, .keep_all = T) %>%
        dplyr::select(-rowMean)
    c181063 <- g181063ex[grep("CISD2", g181063ex$symbol), ] %>%
        gather() %>%
        dplyr::rename("symbol" = "key", "CISD2" = "value") %>%
        slice(-1) %>%
        mutate_at("CISD2", as.numeric) %>%
        dplyr::rename("geo_accession" = "symbol") %>%
        inner_join(g181063c, by = "geo_accession") 
    count(c181063$characteristics_ch1)
    colnames(c181063) = gsub(':ch1', '', colnames(c181063))
    head(c181063)
    x181063a = data.frame()
    cisd2_related_gene = c('TIMM10','NDUFB2','BCL2L1','BCL2L11','CISD2','BID','NDUFA9','NDUFS5','NDUFB9','BCL2','NDUFA7','MCL1','PMAIP1','PIK3C3','CYCS','UQCRB','NDUFS1','HRK','UVRAG','BBC3','PIK3R4','CISD3','NDUFB1','NRBF2','NDUFB4','FXC1','TMEM49')
    for(i in cisd2_related_gene){
        y = g181063ex %>% filter(symbol == i)
        x181063a = rbind(x181063a, y)}
    head(x181063a)
    x181063a = x181063a %>% t()
    colnames(x181063a) = x181063a[1, ]
    x181063 = x181063a %>% 
        as.data.frame() %>% 
        slice(-1) %>% 
        mutate_all(as.numeric) %>%  
        mutate(risk = TIMM10 * -0.05249048 + NDUFB2 * -0.12915551 + BCL2L1 * -0.11683190 + BCL2L11 * -0.05368164 + CISD2 * 0.24852093 + BID * -0.06101778 + NDUFA9 * -0.02235329 + NDUFS5 * 0.01761829 + NDUFB9 * 0.24596655 + BCL2  * 0.10566931 + NDUFA7 * 0.04073703 + MCL1  * 0.26730071 + PMAIP1 * 0.12489498 + PIK3C3 * 0.15104140 + CYCS  * 0.14281014 + UQCRB  * 0.14679661 + NDUFS1 * 0.05981878 + HRK * 0.03899325 + UVRAG * -0.26143078 + BBC3 * -0.04408451 + PIK3R4 * 0.23802877 + CISD3 * 0.12013327 + NDUFB1 * -0.11920934 + NRBF2 * 0.07456076 + NDUFB4 * -0.14108068 + FXC1  * 0.03498904 + TMEM49 * -0.3655337) %>% 
        rownames_to_column('geo_accession') %>% 
        inner_join(c181063, by = 'geo_accession') %>% 
        dplyr::rename('CISD2' = 'CISD2.x') %>% 
        filter(characteristics_ch1=='diagnostic_group: DLBCL') %>% 
        mutate_at(c('os_followup_y', 'os_status'), as.numeric) %>% 
        mutate(osmonth = os_followup_y * 12) %>% 
        filter(osmonth >=1) %>% 
        mutate(state = ifelse(risk < median(risk), 'low', 'high')) 
    fit = survfit(Surv(osmonth, os_status) ~ state, data = x181063) %>% 
            ggsurvplot(.,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2Risk High", "CISD2Rsik Low"),
            legend.title = "",
            title = "GSE181063", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G181063R0.png', width = 6, height = 4)
    fit = survfit(Surv(osmonth, os_status) ~ state, data = subset(x181063, pred_combine == 'ABC')) %>% 
            ggsurvplot(.,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2Risk High", "CISD2Rsik Low"),
            legend.title = "",
            title = "GSE181063 (ABC Subtype)", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G181063R01.png', width = 6, height = 4)
    fit = survfit(Surv(osmonth, os_status) ~ state, data = subset(x181063, pred_combine == 'GCB')) %>% 
            ggsurvplot(.,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2Risk High", "CISD2Rsik Low"),
            legend.title = "",
            title = "GSE181063 (GCB Subtype)", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G181063R02.png', width = 6, height = 4)
    fit = survfit(Surv(osmonth, os_status) ~ state, data = subset(x181063, pred_combine == 'MHG')) %>% 
            ggsurvplot(.,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("CISD2Risk High", "CISD2Rsik Low"),
            legend.title = "",
            title = "GSE181063 (MHG Subtype)", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF"),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.95)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G181063R03.png', width = 6, height = 4)
    fit = survfit(Surv(osmonth, os_status) ~ pred_combine, data = subset(x181063, state == 'high')) %>% 
            ggsurvplot(.,
            pval = TRUE, conf.int = FALSE, 
            risk.table = FALSE, risk.table.col = "strata", 
            linetype = "strata", 
            legend = 'right',
            legend.labs = c("ABC subtype", "GCB subtype", "MHG subtype", "UNC subtype"),
            legend.title = "",
            title = "GSE181063 - High CISD2Risk value", 
            xlab = " Time (Month)",
            ylab = 'Over Survival (%)',
            surv.median.line = "none", 
            palette = c("#E6550DFF", "#3182BDFF", "#31A354FF", '#7876B1FF'),
            ggtheme = theme_bw())  
    fit$plot + theme(axis.title = element_text(size = 16),
                    plot.title = element_text(size = 16)) +
              theme(legend.position = c(0.85, 0.85)) +
              theme(legend.background=element_rect(fill = alpha("white", 0)),
                    legend.key=element_rect(fill = alpha("white", .5))) +
              scale_y_continuous(labels = scales::percent)
    ggsave('~/Documents/01CISD2/F_G181063R04.png', width = 6, height = 4)
    fit2 = x181063 %>% 
        ggplot(aes(fill = pred_combine, y = risk, x = pred_combine)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(data = x181063, aes(y = risk), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('ABC', 'GCB'), c('ABC', 'MHG'), c('ABC', 'UNC'), c('MHG', 'GCB'), c('UNC', 'GCB'), c('MHG', 'UNC')), 
            size = 0.5, textsize = 5, vjust = 0.5, step_increase = 0.1, position = "identity") +
        labs(title = "GSE181063", y = "CISD2Risk value", x = "") +
        #scale_y_continuous(limits = c(3, 11)) +
        scale_fill_manual(values = c("#E6550DFF", "#3182BDFF", "#31A354FF", '#7876B1FF')) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text = element_text(size = 14, vjust = 0.5),
              axis.title = element_text(size = 16),
              plot.title = element_text(size = 16))
    ggsave("~/Documents/01CISD2/F_G181063R2.png", width = 3, height = 4.5)
    count(fit8$Stage)
    fit3 = x181063 %>% 
        filter(Stage != 'NA') %>% 
        mutate(xyz = case_when(Stage == 'I' ~ '0', Stage == 'II' ~ '0',
                               Stage == 'III' ~ '1', Stage == 'IV' ~ '1')) %>% 
        ggplot(aes(fill = xyz, y = risk, x = xyz)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(aes(y = risk), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('0', '1')), size = 0.5, 
            textsize = 5, vjust = 0, step_increase = 0.15, position = "identity") +
        labs(title = "GSE181063", y = "CISD2Risk value", x = "") +
        scale_x_discrete("", labels = c("0" = "Stage I~II", '1' = 'Stage III~IV')) +
        scale_y_continuous(limits = c(3, 10)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text = element_text(size = 10, vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G181063R3.png', width = 2.5, height = 4.5)
    fit4 = x181063 %>% 
        filter(ipi_score != 'NA') %>% 
        mutate(xyz = ifelse(ipi_score <= 2, '0', '1')) %>% 
        ggplot(aes(fill = xyz, y = risk, x = xyz)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(aes(y = risk), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('0', '1')), size = 0.5, 
            textsize = 5, vjust = 0, step_increase = 0.15, position = "identity") +
        labs(title = "GSE181063", y = "CISD2Risk value", x = "") +
        scale_x_discrete("", labels = c("0" = "Low IPI", '1' = 'High IPI')) +
        scale_y_continuous(limits = c(4, 10), breaks=round(seq(4,10,length.out=4),2)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text = element_text(size = 10, vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G181063R4.png', width = 2.5, height = 4.5)
    count(x181063$ldh)
    fit5 = x181063 %>% 
        filter (ldh == 'normal'|ldh == 'raised'|ldh == 'low') %>% 
        mutate(ldh1 = case_when(ldh == 'low' ~ 'normal', ldh == 'normal' ~ 'normal', ldh == 'raised' ~ 'raised')) %>%
        ggplot(aes(fill = ldh1, y = risk, x = ldh1)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(aes(y = risk), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('normal', 'raised')), size = 0.5, 
            textsize = 5, vjust = 0, step_increase = 0.15, position = "identity") +
        labs(title = "GSE181063", y = "CISD2Risk value", x = "") +
        scale_x_discrete("", labels = c("normal" = "normal LDH", 'raised' = 'raised LDH')) +
        scale_y_continuous(limits = c(3.5, 10)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text = element_text(size = 10, vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G181063R5.png', width = 2.5, height = 4.5)
    count(x181063$performance_status_ecog)
    fit6 = x181063 %>% 
        filter(performance_status_ecog != 'NA') %>% 
        mutate(ECOG = ifelse(performance_status_ecog >= 2, '1', '0')) %>% 
        ggplot(aes(fill = ECOG, y = risk, x = ECOG)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(aes(y = risk), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('0', '1')), size = 0.5, 
            textsize = 5, vjust = 0, step_increase = 0.15, position = "identity") +
        labs(title = "GSE181063", y = "CISD2Risk value", x = "") +
        scale_x_discrete("", labels = c("0" = "ECOG < 2", '1' = 'ECOG >=2')) +
        scale_y_continuous(limits = c(3, 10.5)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text = element_text(size = 10, vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G181063R6.png', width = 2.5, height = 4.5)
    fit7 = x181063 %>% 
        filter(num_extranodal != 'NA') %>% 
        mutate_at('num_extranodal', as.numeric) %>% 
        mutate(Extranodal = ifelse(num_extranodal >= 1, '1', '0')) %>% 
        ggplot(aes(fill = Extranodal, y = risk, x = Extranodal)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(aes(y = risk), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('0', '1')), size = 0.5, 
            textsize = 5, vjust = 0, step_increase = 0.15, position = "identity") +
        labs(title = "GSE181063", y = "CISD2Risk value", x = "") +
        scale_x_discrete("", labels = c("0" = "non Extranodal", '1' = 'Extranodal')) +
        scale_y_continuous(limits = c(3, 10.5)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text = element_text(size = 10, vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G181063R7.png', width = 2.5, height = 4.5)
    fit8 = x181063 %>% 
        mutate(wbc != 'NA') %>% 
        mutate_at(c('wbc'), as.numeric) %>% 
        filter(wbc <= 1000) %>% 
        ggplot(aes(y = risk, x = wbc)) +
        geom_point(size=0.5, shape=16) + 
        geom_smooth(method="lm", se = T) + 
        labs(title = "GSE181063", y = "CISD2Risk value", x = expression(WBC (10^9/L))) +
        # scale_x_discrete("", labels = c("0" = "ECOG 0", '1' = 'ECOG 1', '2' = 'ECOG 2')) +
        # scale_y_continuous(limits = c(3, 12)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        stat_cor(method = 'pearson', label.x = 0) +
        theme(legend.position = "none",
              axis.text = element_text(size = 10, vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G181063R8.png', width = 2.5, height = 4.5)
    fit9 = x181063 %>% 
        mutate(hb != 'NA') %>% 
        mutate_at(c('hb'), as.numeric) %>% 
        filter(hb <= 1000) %>% 
        ggplot(aes(y = risk, x = hb)) +
        geom_point(size=0.5, shape=16) + 
        geom_smooth(method="lm", se = T) + 
        labs(title = "GSE181063", y = "CISD2Risk value", x = expression(Hb (g/L))) +
        # scale_x_discrete("", labels = c("0" = "ECOG 0", '1' = 'ECOG 1', '2' = 'ECOG 2')) +
        # scale_y_continuous(limits = c(3, 12)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        stat_cor(method = 'pearson', label.x = 0) +
        theme(legend.position = "none",
              axis.text = element_text(size = 10, vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G181063R9.png', width = 2.5, height = 4.5)
    fit10 = x181063 %>% 
        mutate(plts != 'NA') %>% 
        mutate_at(c('plts'), as.numeric) %>% 
        filter(plts <= 2500) %>% 
        ggplot(aes(y = risk, x = plts)) +
        geom_point(size=0.5, shape=16) + 
        geom_smooth(method="lm", se = T) + 
        labs(title = "GSE181063", y = "CISD2Risk value", x = expression(PLT (10^9/L)/L)) +
        # scale_x_discrete("", labels = c("0" = "ECOG 0", '1' = 'ECOG 1', '2' = 'ECOG 2')) +
        # scale_y_continuous(limits = c(3, 12)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        stat_cor(method = 'pearson', label.x = 0) +
        theme(legend.position = "none",
              axis.text = element_text(size = 10, vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G181063R10.png', width = 2.5, height = 4.5)
    fit11 = x181063 %>% 
        mutate(lymphs != 'NA') %>% 
        mutate_at(c('lymphs'), as.numeric) %>% 
        filter(lymphs <= 100) %>% 
        ggplot(aes(y = risk, x = lymphs)) +
        geom_point(size=0.5, shape=16) + 
        geom_smooth(method="lm", se = T) + 
        labs(title = "GSE181063", y = "CISD2Risk value", x = expression(lymphs (10^9/L))) +
        # scale_x_discrete("", labels = c("0" = "ECOG 0", '1' = 'ECOG 1', '2' = 'ECOG 2')) +
        # scale_y_continuous(limits = c(3, 12)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        stat_cor(method = 'pearson', label.x = 0) +
        theme(legend.position = "none",
              axis.text = element_text(size = 10, vjust = 0.5),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G181063R11.png', width = 2.5, height = 4.5)
    fit = x181063 %>% 
        filter(sample_pre_active_treatment != 'Treatment data not available') %>% 
        ggplot(aes(fill = sample_pre_active_treatment, y = risk, x = sample_pre_active_treatment)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(aes(y = risk), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('No', 'Yes')), size = 0.5, 
            textsize = 5, vjust = 0, step_increase = 0.15, position = "identity") +
        labs(title = "GSE181063", y = "CISD2Risk value", x = "") +
        scale_x_discrete("", labels = c("No" = "No active treatment", 'Yes' = 'Active treatment')) +
        scale_y_continuous(limits = c(3.5, 10)) +
        scale_fill_manual(values = c("#E6550DFF", "#3182BDFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10, hjust = 1, angle = 45),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G181063Rb.png', width = 2.5, height = 6)
#### Correction plot based on GSE117556 and GSE181063
    require(corrplot)    
    colnames(x117556)
    tiff('~/Documents/01CISD2/F_G117556corr.tiff', width = 2000, height = 2500)
    y117556 = x117556 %>% 
        select(c('TIMM10','NDUFB2','BCL2L1','BCL2L11','CISD2','BID','NDUFA9','NDUFS5','NDUFB9','BCL2','NDUFA7','MCL1','PMAIP1','PIK3C3','CYCS','UQCRB','NDUFS1','HRK','UVRAG','BBC3','PIK3R4','CISD3','NDUFB1','NRBF2','NDUFB4','FXC1','TMEM49'))
    corrplot(corr = cor(y117556), type = "lower", 
            col = colorRampPalette(c("#3182BDFF", "#E6550DFF"))(20),
            tl.col = 'black', tl.cex = 3, cl.cex = 2.5)
    dev.off() 
    tiff('~/Documents/01CISD2/F_G181063corr.tiff', width = 2000, height = 2500)
    y181063 = x181063 %>% 
        select(c('TIMM10','NDUFB2','BCL2L1','BCL2L11','CISD2','BID','NDUFA9','NDUFS5','NDUFB9','BCL2','NDUFA7','MCL1','PMAIP1','PIK3C3','CYCS','UQCRB','NDUFS1','HRK','UVRAG','BBC3','PIK3R4','CISD3','NDUFB1','NRBF2','NDUFB4','FXC1','TMEM49'))
    corrplot(corr = cor(y181063), type = "lower", 
            col = colorRampPalette(c("#3182BDFF", "#E6550DFF"))(20),
            tl.col = 'black', tl.cex = 3, cl.cex = 2.5)
    dev.off()
## Figure 6 Immune analysis
#### estimate analysis
    require(tidyestimate)
    fit = estimate_score(g117556ex, is_affymetrix = F) %>% 
        dplyr::rename('geo_accession' = 'sample') %>% 
        inner_join(x117556, by = 'geo_accession')
    fit1 = fit %>% 
        ggplot(aes(fill = state, y = stromal, x = state)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(aes(y = stromal), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('low', 'high')), size = 0.5, 
            textsize = 5, vjust = 0.4, step_increase = 0.15, position = "identity") +
        labs(title = "GSE117556", y = "Stromal Score", x = "") +
        scale_x_discrete("", labels = c("low" = "low CISD2Risk value", 'high' = 'high CISD2Risk value')) +
        #scale_y_continuous(limits = c(3, 10)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G117556Ri1.png', width = 3, height = 5)
    fit2 = fit %>% 
        ggplot(aes(fill = state, y = immune, x = state)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(aes(y = immune), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('low', 'high')), size = 0.5, 
            textsize = 5, vjust = 0.4, step_increase = 0.15, position = "identity") +
        labs(title = "GSE117556", y = "Immune Score", x = "") +
        scale_x_discrete("", labels = c("low" = "low CISD2Risk value", 'high' = 'high CISD2Risk value')) +
        #scale_y_continuous(limits = c(3, 10)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G117556Ri2.png', width = 3, height = 5)
    fit3 = fit %>% 
        ggplot(aes(fill = state, y = estimate, x = state)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(aes(y = estimate), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('low', 'high')), size = 0.5, 
            textsize = 5, vjust = 0.4, step_increase = 0.15, position = "identity") +
        labs(title = "GSE117556", y = "Estimate Score", x = "") +
        scale_x_discrete("", labels = c("low" = "low CISD2Risk value", 'high' = 'high CISD2Risk value')) +
        #scale_y_continuous(limits = c(3, 10)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G117556Ri3.png', width = 3, height = 5)
    fit = estimate_score(g181063ex, is_affymetrix = F) %>% 
        dplyr::rename('geo_accession' = 'sample') %>% 
        inner_join(x181063, by = 'geo_accession')
    fit1 = fit %>% 
        ggplot(aes(fill = state, y = stromal, x = state)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(aes(y = stromal), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('low', 'high')), size = 0.5, 
            textsize = 5, vjust = 0.4, step_increase = 0.15, position = "identity") +
        labs(title = "GSE181063", y = "Stromal Score", x = "") +
        scale_x_discrete("", labels = c("low" = "low CISD2Risk value", 'high' = 'high CISD2Risk value')) +
        #scale_y_continuous(limits = c(3, 10)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G181063Ri1.png', width = 3, height = 5)
    fit2 = fit %>% 
        ggplot(aes(fill = state, y = immune, x = state)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(aes(y = immune), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('low', 'high')), size = 0.5, 
            textsize = 5, vjust = 0.4, step_increase = 0.15, position = "identity") +
        labs(title = "GSE181063", y = "Immune Score", x = "") +
        scale_x_discrete("", labels = c("low" = "low CISD2Risk value", 'high' = 'high CISD2Risk value')) +
        #scale_y_continuous(limits = c(3, 10)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G181063Ri2.png', width = 3, height = 5)
    fit3 = fit %>% 
        ggplot(aes(fill = state, y = estimate, x = state)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(aes(y = estimate), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('low', 'high')), size = 0.5, 
            textsize = 5, vjust = 0.4, step_increase = 0.15, position = "identity") +
        labs(title = "GSE181063", y = "Estimate Score", x = "") +
        scale_x_discrete("", labels = c("low" = "low CISD2Risk value", 'high' = 'high CISD2Risk value')) +
        #scale_y_continuous(limits = c(3, 10)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    ggsave('~/Documents/01CISD2/F_G181063Ri3.png', width = 3, height = 5)
#### CIBERSORTx analysis
    b117556x = read.table("~/Documents/01CISD2/CIBERSORTx_Job5_Results.txt", header = T, sep = "\t") %>% 
        dplyr::rename('geo_accession' = 'Mixture') %>% 
        inner_join(x117556, by = 'geo_accession')
    fit = b117556x %>% 
        ggplot(aes(fill = state, y = Macrophages.M2, x = state)) +
        geom_boxplot(alpha = 0.8) +
        geom_jitter(aes(y = Macrophages.M2), size = 1, shape = 20,
            stroke = 0.01, show.legend = FALSE, width = 0.15) +
        stat_boxplot(geom = "errorbar", width = 0.1) +
        geom_signif(map_signif_level = TRUE,
            comparisons = list(c('low', 'high')), size = 0.5, 
            textsize = 5, vjust = 0.4, step_increase = 0.15, position = "identity") +
        labs(title = "GSE181063", y = "Estimate Score", x = "") +
        scale_x_discrete("", labels = c("low" = "low CISD2Risk value", 'high' = 'high CISD2Risk value')) +
        #scale_y_continuous(limits = c(3, 10)) +
        scale_fill_manual(values = c("#3182BDFF", "#E6550DFF", "#7876B1FF", "#20854EFF")) + 
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16))
    res1 = data.frame(b117556x[,c(1:23,191)]) %>% 
        pivot_longer(cols = colnames(.)[2:23],
               names_to = "cell.type",
               values_to = 'value') %>% 
        ggplot(aes(x = cell.type, y = value, fill = state)) + 
        geom_boxplot(alpha = 0.8) + 
        theme_bw() + 
        labs(x = "", y = "Estimated Proportion", title = 'GSE117556') +
        theme(legend.position = "right") + 
        guides(fill=guide_legend(title="CISD2Risk value")) + 
        scale_fill_manual(values = c("#E6550DFF", "#3182BDFF")) + 
        theme(axis.text.x = element_text(angle=60, hjust = 1, size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16)) +
        stat_compare_means(aes(group = state),label = "p.signif", size= 4, label.y = 0.5, method = "kruskal.test")
    ggsave('~/Documents/01CISD2/F_G117556RR.png', width = 15, height = 5)
    b181063x = read.table("~/Documents/01CISD2/CIBERSORTx_Job6_Results.txt", header = T, sep = "\t") %>% 
        dplyr::rename('geo_accession' = 'Mixture') %>% 
        inner_join(x181063, by = 'geo_accession')
    res1 = data.frame(b181063x[,c(1:23,168)]) %>% 
        pivot_longer(cols = colnames(.)[2:23],
               names_to = "cell.type",
               values_to = 'value') %>% 
        ggplot(aes(x = cell.type, y = value, fill = state)) + 
        geom_boxplot(alpha = 0.8) + 
        theme_bw() + 
        labs(x = "", y = "Estimated Proportion", title = 'GSE181063') +
        theme(legend.position = "right") + 
        guides(fill=guide_legend(title="CISD2Risk value")) + 
        scale_fill_manual(values = c("#E6550DFF", "#3182BDFF")) + 
        theme(axis.text.x = element_text(angle=60, hjust = 1, size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16)) +
        stat_compare_means(aes(group = state), label = "p.signif", size= 4, label.y = 0.64, method = "kruskal.test")
    ggsave('~/Documents/01CISD2/F_G181063RR.png', width = 15, height = 5)
#### MCPcounter analysis
    require(MCPcounter)
    m117556 = g117556ex %>% 
        column_to_rownames('symbol') %>% 
        as.matrix()
    probesets = read.table('~/Documents/01CISD2/probesets.txt',
                        sep="\t",
                        stringsAsFactors=FALSE,
                        colClasses="character")
    genes = read.table('~/Documents/01CISD2/genes.txt',
                    sep="\t",
                    stringsAsFactors=FALSE,
                    header=TRUE,
                    colClasses="character",
                    check.names=FALSE)
    res = MCPcounter.estimate(m117556, featuresType = "HUGO_symbols", 
            probesets = probesets, genes = genes)
    res1 = res %>% 
        as.data.frame() %>% 
        rownames_to_column('cell.type') %>% 
        pivot_longer(cols = colnames(.)[2:929],
               names_to = "geo_accession",
               values_to = 'value') %>% 
        inner_join(x117556[c('geo_accession', 'state')], by = 'geo_accession') %>% 
        ggplot(aes(x = cell.type, y = value, fill = state)) + 
        geom_boxplot(alpha = 0.8) + 
        theme_bw() + 
        labs(x = "", y = "MCP-counter scores", title = 'GSE117556') +
        theme(legend.position = "right") + 
        guides(fill=guide_legend(title="CISD2Risk value")) + 
        scale_fill_manual(values = c("#E6550DFF", "#3182BDFF")) + 
        scale_y_continuous(limits = c(7, 17)) + 
        theme(axis.text.x = element_text(angle=60, hjust = 1, size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16)) +
        stat_compare_means(aes(group = state), label = "p.signif", size= 4, label.y = 16, method = "kruskal.test")
    ggsave('~/Documents/01CISD2/F_G117556MP.png', width = 7.5, height = 5)
    m181063 = g181063ex %>% 
        column_to_rownames('symbol') %>% 
        as.matrix()
    res = MCPcounter.estimate(m181063, featuresType = "HUGO_symbols", 
            probesets = probesets, genes = genes)
    res1 = res %>% 
        as.data.frame() %>% 
        rownames_to_column('cell.type') %>% 
        pivot_longer(cols = colnames(.)[2:1311],
               names_to = "geo_accession",
               values_to = 'value') %>% 
        inner_join(x181063[c('geo_accession', 'state')], by = 'geo_accession') %>% 
        ggplot(aes(x = cell.type, y = value, fill = state)) + 
        geom_boxplot(alpha = 0.8) + 
        theme_bw() + 
        labs(x = "", y = "MCP-counter scores", title = 'GSE181063') +
        theme(legend.position = "right") + 
        guides(fill=guide_legend(title="CISD2Risk value")) + 
        scale_fill_manual(values = c("#E6550DFF", "#3182BDFF")) + 
        scale_y_continuous(limits = c(6, 17)) + 
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 16)) +
        stat_compare_means(aes(group = state), label = "p.signif", size= 4, label.y = 16, method = "kruskal.test")
    ggsave('~/Documents/01CISD2/F_G181063MP.png', width = 7.5, height = 5)
## Figure 7 Plot for Univariate and Multivariate analysis
#### Univariate and Multivariate analysis (GSE117556)
    str(x117556)
    y117556 = x117556 %>% 
        mutate(year = ifelse(AGE <60, 'Younger', 'Older')) %>% 
        mutate(ipi = ifelse(IPI_SCORE < 2, '0', '1')) %>% 
        mutate(LDH1 = ifelse(LDH < 245, '0', '1')) %>% 
        mutate(ECOG1 = ifelse(ECOG >= 2, '1', '0')) %>% 
        mutate(COO = ifelse(molecularcoosubtype == 'GCB', 'GCB', 'non-GCB')) %>% 
        mutate(Molecular_subtype = ifelse(molecularsubtype == 'MHG', 'MHG', 'non-MHG')) %>% 
        select(c('year', 'GENDER', 'COO', 'Molecular_subtype', 'ipi', 'ECOG1', 'LDH1', 
                'INV_EXTRANODAL_BAS_IND', 'expressor_RNA',
                'state', 'risk', 'PROG_IND', 'PROG_TIME', 'DEATH_IND', 'FU_TIME'))
    coxa = names(y117556)[1:10]
    xa = data.frame()
    for(i in coxa) {
        expr = y117556[, i]
        cox = coxph(Surv(FU_TIME, DEATH_IND) ~ expr, y117556)
        coxsum = summary(cox)
        xa = rbind(xa, cbind(coxi = i, 
            HR = round(coxsum$coefficients[, 'exp(coef)'],2), 
            lower = round(coxsum$conf.int[, 3],2), 
            upper = round(coxsum$conf.int[, 4],2),
            'HR(95%CI)' = paste(round(coxsum$coefficients[, 'exp(coef)'],2), '(', round(coxsum$conf.int[, 3],2), '-', round(coxsum$conf.int[, 4],2), ')', sep = ''),
            pvalue = round(coxsum$coefficients[, 'Pr(>|z|)'],2),
            z = round(coxsum$coefficients[, "z"],2)))}    
    stage117556 = x117556 %>% 
            filter(STAGE != 0) %>% 
            mutate(Stage = case_when(STAGE == 'Stage I' ~ '0', STAGE == 'Stage II' ~ '0',
                               STAGE == 'Stage III' ~ '1', STAGE == 'Stage IV' ~ '1')) %>% 
            select('Stage', 'PROG_IND', 'PROG_TIME', 'DEATH_IND', 'FU_TIME')
    coxb = names(stage117556)[1]
    xb = data.frame()
    for(i in coxb) {
        expr = stage117556[, i]
        cox = coxph(Surv(FU_TIME, DEATH_IND) ~ expr, stage117556)
        coxsum = summary(cox)
        xb = rbind(xb, cbind(coxi = i, 
            HR = round(coxsum$coefficients[, 'exp(coef)'],2), 
            lower = round(coxsum$conf.int[, 3],2), 
            upper = round(coxsum$conf.int[, 4],2),
            'HR(95%CI)' = paste(round(coxsum$coefficients[, 'exp(coef)'],2), '(', round(coxsum$conf.int[, 3],2), '-', round(coxsum$conf.int[, 4],2), ')', sep = ''),
            pvalue = round(coxsum$coefficients[, 'Pr(>|z|)'],2),
            z = round(coxsum$coefficients[, "z"],2)))}
    x = rbind(xa[1:6,], xb, xa[7:10,]) %>% 
        as_tibble() %>% 
        mutate_at(c('HR','lower','upper'), as.numeric) %>% 
        mutate(log_HR = log(HR)) %>% 
        mutate(log_lower = log(lower)) %>% 
        mutate(log_upper = log(upper)) %>% 
        mutate(iterm = factor(c('Age: < 60 vs. >= 60', 'Gender: female vs. male', 'GCB: yes vs. no', 'MHG: yes vs. no', 'IPI: < 2 vs. >= 2 ', 'ECOG PS: < 2 vs. >=2', 'Stage: I-II vs. III-IV', 'LDH: normal vs. raised', 'Extranodal: yes vs. no', 'Double-expressor: yes vs. no', 'CISD2Risk: high vs. low'), levels = c('Age: < 60 vs. >= 60', 'Gender: female vs. male', 'GCB: yes vs. no', 'MHG: yes vs. no', 'IPI: < 2 vs. >= 2 ', 'ECOG PS: < 2 vs. >=2', 'Stage: I-II vs. III-IV', 'LDH: normal vs. raised', 'Extranodal: yes vs. no', 'Double-expressor: yes vs. no', 'CISD2Risk: high vs. low'))) %>% 
        mutate(log_x = paste(round(log_HR, 2),'(',round(log_lower, 2),',', round(log_upper, 2),')', sep = '')) %>% 
        ggplot(aes(y = fct_rev(iterm))) + 
        geom_point(aes(x=log_HR), shape=16, size=5, color = '#3182BDFF') +
        geom_linerange(aes(xmin=log_lower, xmax=log_upper)) +
        geom_vline(xintercept = 0, linetype="dashed") +
        labs(x="Log Hazard Ratio", y="") +
        coord_cartesian(ylim=c(1,11), xlim=c(-2, 1.5)) +
        theme_classic() + 
        theme(axis.text.y = element_text(size = 16),
              axis.title.y= element_blank(),
              axis.ticks.y= element_blank(),
              axis.line.y = element_blank())
    ggsave('~/Documents/01CISD2/forrest117556ax.png', width = 8, height = 6)
    z117556 = x117556 %>% 
        filter(STAGE != 0) %>% 
        mutate(Stage = case_when(STAGE == 'Stage I' ~ '0', STAGE == 'Stage II' ~ '0',
                               STAGE == 'Stage III' ~ '1', STAGE == 'Stage IV' ~ '1')) %>% 
        mutate(year = ifelse(AGE <60, 'Younger', 'Older')) %>% 
        mutate(ECOG1 = ifelse(ECOG >= 2, '1', '0')) %>% 
        mutate(ipi = ifelse(IPI_SCORE < 2, '0', '1')) %>% 
        mutate(LDH1 = ifelse(LDH < 245, '0', '1')) %>% 
        mutate(COO = ifelse(molecularcoosubtype == 'GCB', 'GCB', 'non-GCB')) %>% 
        mutate(Molecular_subtype = ifelse(molecularsubtype == 'MHG', 'MHG', 'non-MHG')) %>% 
        select(c('year', 'GENDER', 'COO', 'Molecular_subtype', 'ipi', 'ECOG1', 'Stage', 'LDH1', 
                'INV_EXTRANODAL_BAS_IND', 'expressor_RNA', 
                'state', 'risk', 'PROG_IND', 'PROG_TIME', 'DEATH_IND', 'FU_TIME')) %>% 
        coxph(Surv(FU_TIME, DEATH_IND) ~ year + GENDER + COO + Molecular_subtype + ipi + ECOG1 + Stage + LDH1 + INV_EXTRANODAL_BAS_IND + expressor_RNA + state, .) %>% 
        summary(.)
    ytable = cbind(HR = z117556$coefficients[, 'exp(coef)'],
          lower = z117556$conf.int[, 3],
          upper = z117556$conf.int[, 4],
          pvalue = z117556$coefficients[, 'Pr(>|z|)'],
          z = z117556$coefficients[, "z"]) %>% 
          as_tibble() %>% 
          mutate(log_HR = log(HR)) %>% 
          mutate(log_lower = log(lower)) %>% 
          mutate(log_upper = log(upper)) %>% 
          mutate(iterm = factor(c('Age: < 60 vs. >= 60', 'Gender: female vs. male', 'GCB: yes vs. no', 'MHG: yes vs. no', 'IPI: < 2 vs. >= 2 ', 'ECOG PS: < 2 vs. >=2', 'Stage: I-II vs. III-IV', 'LDH: normal vs. raised', 'Extranodal: yes vs. no', 'Double-expressor: yes vs. no', 'CISD2Risk: high vs. low'), levels = c('Age: < 60 vs. >= 60', 'Gender: female vs. male', 'GCB: yes vs. no', 'MHG: yes vs. no', 'IPI: < 2 vs. >= 2 ', 'ECOG PS: < 2 vs. >=2', 'Stage: I-II vs. III-IV', 'LDH: normal vs. raised', 'Extranodal: yes vs. no', 'Double-expressor: yes vs. no', 'CISD2Risk: high vs. low'))) %>% 
          mutate(log_x = paste(round(log_HR, 2),'(',round(log_lower, 2),',', round(log_upper, 2),')', sep = ''))
    y = ytable  %>% 
        ggplot(aes(y = fct_rev(iterm))) + 
        geom_point(aes(x=log_HR), shape=15, size=5, color = "#E6550DFF") +
        geom_linerange(aes(xmin=log_lower, xmax=log_upper)) +
        geom_vline(xintercept = 0, linetype="dashed") +
        labs(x="Log Hazard Ratio", y="") +
        coord_cartesian(ylim=c(1,11), xlim=c(-2, 1.5)) +
        theme_classic() + 
        theme(axis.text.y = element_text(size = 16),
              axis.title.y= element_blank(),
              axis.ticks.y= element_blank(),
              axis.line.y = element_blank())
    ggsave('~/Documents/01CISD2/forrest117556bx.png', width = 8, height = 6)
    y117556 = x117556 %>% 
        mutate(year = ifelse(AGE <60, 'Younger', 'Older')) %>% 
        mutate(ipi = ifelse(IPI_SCORE < 2, '0', '1')) %>% 
        mutate(LDH1 = ifelse(LDH < 245, '0', '1')) %>% 
        mutate(ECOG1 = ifelse(ECOG >= 2, '1', '0')) %>% 
        mutate(COO = ifelse(molecularcoosubtype == 'GCB', 'GCB', 'non-GCB')) %>% 
        mutate(Molecular_subtype = ifelse(molecularsubtype == 'MHG', 'MHG', 'non-MHG')) %>% 
        select(c('year', 'GENDER', 'COO', 'Molecular_subtype', 'ipi', 'ECOG1', 'LDH1', 
                'INV_EXTRANODAL_BAS_IND', 'expressor_RNA',
                'state', 'risk', 'PROG_IND', 'PROG_TIME', 'DEATH_IND', 'FU_TIME'))
    coxa = names(y117556)[1:10]
    xa = data.frame()
    for(i in coxa) {
        expr = y117556[, i]
        cox = coxph(Surv(PROG_TIME, PROG_IND) ~ expr, y117556)
        coxsum = summary(cox)
        xa = rbind(xa, cbind(coxi = i, 
            HR = round(coxsum$coefficients[, 'exp(coef)'],2), 
            lower = round(coxsum$conf.int[, 3],2), 
            upper = round(coxsum$conf.int[, 4],2),
            'HR(95%CI)' = paste(round(coxsum$coefficients[, 'exp(coef)'],2), '(', round(coxsum$conf.int[, 3],2), '-', round(coxsum$conf.int[, 4],2), ')', sep = ''),
            pvalue = round(coxsum$coefficients[, 'Pr(>|z|)'],2),
            z = round(coxsum$coefficients[, "z"],2)))}    
    stage117556 = x117556 %>% 
            filter(STAGE != 0) %>% 
            mutate(Stage = case_when(STAGE == 'Stage I' ~ '0', STAGE == 'Stage II' ~ '0',
                               STAGE == 'Stage III' ~ '1', STAGE == 'Stage IV' ~ '1')) %>% 
            select('Stage', 'PROG_IND', 'PROG_TIME', 'DEATH_IND', 'FU_TIME')
    coxb = names(stage117556)[1]
    xb = data.frame()
    for(i in coxb) {
        expr = stage117556[, i]
        cox = coxph(Surv(PROG_TIME, PROG_IND) ~ expr, stage117556)
        coxsum = summary(cox)
        xb = rbind(xb, cbind(coxi = i, 
            HR = round(coxsum$coefficients[, 'exp(coef)'],2), 
            lower = round(coxsum$conf.int[, 3],2), 
            upper = round(coxsum$conf.int[, 4],2),
            'HR(95%CI)' = paste(round(coxsum$coefficients[, 'exp(coef)'],2), '(', round(coxsum$conf.int[, 3],2), '-', round(coxsum$conf.int[, 4],2), ')', sep = ''),
            pvalue = round(coxsum$coefficients[, 'Pr(>|z|)'],2),
            z = round(coxsum$coefficients[, "z"],2)))}
    x = rbind(xa[1:6,], xb, xa[7:10,]) %>% 
        as_tibble() %>% 
        mutate_at(c('HR','lower','upper'), as.numeric) %>% 
        mutate(log_HR = log(HR)) %>% 
        mutate(log_lower = log(lower)) %>% 
        mutate(log_upper = log(upper)) %>% 
        mutate(iterm = factor(c('Age: < 60 vs. >= 60', 'Gender: female vs. male', 'GCB: yes vs. no', 'MHG: yes vs. no', 'IPI: < 2 vs. >= 2 ', 'ECOG PS: < 2 vs. >=2', 'Stage: I-II vs. III-IV', 'LDH: normal vs. raised', 'Extranodal: yes vs. no', 'Double-expressor: yes vs. no', 'CISD2Risk: high vs. low'), levels = c('Age: < 60 vs. >= 60', 'Gender: female vs. male', 'GCB: yes vs. no', 'MHG: yes vs. no', 'IPI: < 2 vs. >= 2 ', 'ECOG PS: < 2 vs. >=2', 'Stage: I-II vs. III-IV', 'LDH: normal vs. raised', 'Extranodal: yes vs. no', 'Double-expressor: yes vs. no', 'CISD2Risk: high vs. low'))) %>% 
        mutate(log_x = paste(round(log_HR, 2),'(',round(log_lower, 2),',', round(log_upper, 2),')', sep = '')) %>% 
        ggplot(aes(y = fct_rev(iterm))) + 
        geom_point(aes(x=log_HR), shape=16, size=5, color = '#3182BDFF') +
        geom_linerange(aes(xmin=log_lower, xmax=log_upper)) +
        geom_vline(xintercept = 0, linetype="dashed") +
        labs(x="Log Hazard Ratio", y="") +
        coord_cartesian(ylim=c(1,11), xlim=c(-2, 1.5)) +
        theme_classic() + 
        theme(axis.text.y = element_text(size = 16),
              axis.title.y= element_blank(),
              axis.ticks.y= element_blank(),
              axis.line.y = element_blank())
    ggsave('~/Documents/01CISD2/forrest117556ay.png', width = 8, height = 6)
    z117556 = x117556 %>% 
        filter(STAGE != 0) %>% 
        mutate(Stage = case_when(STAGE == 'Stage I' ~ '0', STAGE == 'Stage II' ~ '0',
                               STAGE == 'Stage III' ~ '1', STAGE == 'Stage IV' ~ '1')) %>% 
        mutate(year = ifelse(AGE <60, 'Younger', 'Older')) %>% 
        mutate(ECOG1 = ifelse(ECOG >= 2, '1', '0')) %>% 
        mutate(ipi = ifelse(IPI_SCORE < 2, '0', '1')) %>% 
        mutate(LDH1 = ifelse(LDH < 245, '0', '1')) %>% 
        mutate(COO = ifelse(molecularcoosubtype == 'GCB', 'GCB', 'non-GCB')) %>% 
        mutate(Molecular_subtype = ifelse(molecularsubtype == 'MHG', 'MHG', 'non-MHG')) %>% 
        select(c('year', 'GENDER', 'COO', 'Molecular_subtype', 'ipi', 'ECOG1', 'Stage', 'LDH1', 
                'INV_EXTRANODAL_BAS_IND', 'expressor_RNA', 
                'state', 'risk', 'PROG_IND', 'PROG_TIME', 'DEATH_IND', 'FU_TIME')) %>% 
        coxph(Surv(PROG_TIME, PROG_IND) ~ year + GENDER + COO + Molecular_subtype + ipi + ECOG1 + Stage + LDH1 + INV_EXTRANODAL_BAS_IND + expressor_RNA + state, .) %>% 
        summary(.)
    ytable = cbind(HR = z117556$coefficients[, 'exp(coef)'],
          lower = z117556$conf.int[, 3],
          upper = z117556$conf.int[, 4],
          pvalue = z117556$coefficients[, 'Pr(>|z|)'],
          z = z117556$coefficients[, "z"]) %>% 
          as_tibble() %>% 
          mutate(log_HR = log(HR)) %>% 
          mutate(log_lower = log(lower)) %>% 
          mutate(log_upper = log(upper)) %>% 
          mutate(iterm = factor(c('Age: < 60 vs. >= 60', 'Gender: female vs. male', 'GCB: yes vs. no', 'MHG: yes vs. no', 'IPI: < 2 vs. >= 2 ', 'ECOG PS: < 2 vs. >=2', 'Stage: I-II vs. III-IV', 'LDH: normal vs. raised', 'Extranodal: yes vs. no', 'Double-expressor: yes vs. no', 'CISD2Risk: high vs. low'), levels = c('Age: < 60 vs. >= 60', 'Gender: female vs. male', 'GCB: yes vs. no', 'MHG: yes vs. no', 'IPI: < 2 vs. >= 2 ', 'ECOG PS: < 2 vs. >=2', 'Stage: I-II vs. III-IV', 'LDH: normal vs. raised', 'Extranodal: yes vs. no', 'Double-expressor: yes vs. no', 'CISD2Risk: high vs. low'))) %>% 
          mutate(log_x = paste(round(log_HR, 2),'(',round(log_lower, 2),',', round(log_upper, 2),')', sep = ''))
    y = ytable  %>% 
        ggplot(aes(y = fct_rev(iterm))) + 
        geom_point(aes(x=log_HR), shape=15, size=5, color = "#E6550DFF") +
        geom_linerange(aes(xmin=log_lower, xmax=log_upper)) +
        geom_vline(xintercept = 0, linetype="dashed") +
        labs(x="Log Hazard Ratio", y="") +
        coord_cartesian(ylim=c(1,11), xlim=c(-2, 1.5)) +
        theme_classic() + 
        theme(axis.text.y = element_text(size = 16),
              axis.title.y= element_blank(),
              axis.ticks.y= element_blank(),
              axis.line.y = element_blank())
    ggsave('~/Documents/01CISD2/forrest117556by.png', width = 8, height = 6)
    coxa = names(y117556)[1:11]
    xa = data.frame()
    for(i in coxa) {
        expr = y117556[, i]
        cox = coxph(Surv(PROG_TIME, PROG_IND) ~ expr, y117556)
        coxsum = summary(cox)
        xa = rbind(xa, cbind(coxi = i, 
            HR = round(coxsum$coefficients[, 'exp(coef)'],2), 
            lower = round(coxsum$conf.int[, 3],2), 
            upper = round(coxsum$conf.int[, 4],2),
            'HR(95%CI)' = paste(round(coxsum$coefficients[, 'exp(coef)'],2), '(', round(coxsum$conf.int[, 3],2), '-', round(coxsum$conf.int[, 4],2), ')', sep = ''),
            pvalue = round(coxsum$coefficients[, 'Pr(>|z|)'],2),
            z = round(coxsum$coefficients[, "z"],2)))}    
    stage117556 = x117556 %>% 
            filter(STAGE != 0) %>% 
            mutate(Stage = case_when(STAGE == 'Stage I' ~ '0', STAGE == 'Stage II' ~ '0',
                               STAGE == 'Stage III' ~ '1', STAGE == 'Stage IV' ~ '1')) %>% 
            select('Stage', 'PROG_IND', 'PROG_TIME', 'DEATH_IND', 'FU_TIME')
    coxb = names(stage117556)[1]
    xb = data.frame()
    for(i in coxb) {
        expr = stage117556[, i]
        cox = coxph(Surv(PROG_TIME, PROG_IND) ~ expr, stage117556)
        coxsum = summary(cox)
        xb = rbind(xb, cbind(coxi = i, 
            HR = round(coxsum$coefficients[, 'exp(coef)'],2), 
            lower = round(coxsum$conf.int[, 3],2), 
            upper = round(coxsum$conf.int[, 4],2),
            'HR(95%CI)' = paste(round(coxsum$coefficients[, 'exp(coef)'],2), '(', round(coxsum$conf.int[, 3],2), '-', round(coxsum$conf.int[, 4],2), ')', sep = ''),
            pvalue = round(coxsum$coefficients[, 'Pr(>|z|)'],2),
            z = round(coxsum$coefficients[, "z"],2)))}
    x = rbind(xa[1:12,], xb, xa[13,]) %>% 
        as_tibble() %>% 
        mutate_at(c('HR','lower','upper'), as.numeric) %>% 
        mutate(log_HR = log(HR)) %>% 
        mutate(log_lower = log(lower)) %>% 
        mutate(log_upper = log(upper)) %>% 
        mutate(iterm = factor(c('Age', 'Gender', 'ABC vs GCB', 'ABC vs UNC', 'Clinical Response', 'IPI', 'ECOG 0 vs ECOG 1', 'ECOG 0 vs ECOG 2', 'Extran', 'Myc RNA', 'BCL2 RNA', 'Double expressor', 'Clinical Stage', 'CISD2Risk value'), levels = c('Age', 'Gender', 'ABC vs GCB', 'ABC vs UNC', 'Clinical Response', 'IPI', 'ECOG 0 vs ECOG 1', 'ECOG 0 vs ECOG 2', 'Extran', 'Myc RNA', 'BCL2 RNA', 'Double expressor', 'Clinical Stage', 'CISD2Risk value'))) %>% 
        mutate(log_x = paste(round(log_HR, 2),'(',round(log_lower, 2),',', round(log_upper, 2),')', sep = '')) %>% 
        ggplot(aes(y = fct_rev(iterm))) + 
        geom_point(aes(x=log_HR), shape=16, size=5, color = '#3182BDFF') +
        geom_linerange(aes(xmin=log_lower, xmax=log_upper)) +
        geom_vline(xintercept = 0, linetype="dashed") +
        labs(x="Log Hazard Ratio", y="") +
        coord_cartesian(ylim=c(1,14), xlim=c(-2.5, 1.5)) +
        theme_classic() + 
        theme(axis.text.y = element_text(size = 16),
              axis.title.y= element_blank(),
              axis.ticks.y= element_blank(),
              axis.line.y = element_blank())
    ggsave('~/Documents/01CISD2/forrest11755603.png', width = 8, height = 8)
    z117556 = x117556 %>% 
        filter(STAGE != 0) %>% 
        mutate(Stage = case_when(STAGE == 'Stage I' ~ '0', STAGE == 'Stage II' ~ '0',
                               STAGE == 'Stage III' ~ '1', STAGE == 'Stage IV' ~ '1')) %>% 
        mutate(year = ifelse(AGE <60, 'Younger', 'Older')) %>% 
        mutate_at('ECOG', as.factor) %>% 
        mutate(xyz = ifelse(IPI_SCORE < 2, '0', '1')) %>% 
        select(c('year', 'GENDER', 'molecularcoosubtype', 'RESX', 'xyz', 'ECOG', 
                'INV_EXTRANODAL_BAS_IND', 'MYC_RNA', 'BCL2_RNA', 'expressor_RNA', 'STAGE',
                'state', 'risk', 'PROG_IND', 'PROG_TIME', 'DEATH_IND', 'FU_TIME')) %>% 
        coxph(Surv(PROG_TIME, PROG_IND) ~ year + GENDER +
              RESX + xyz + ECOG + INV_EXTRANODAL_BAS_IND + 
              MYC_RNA + BCL2_RNA + expressor_RNA + STAGE + state, .) %>% 
        summary(.)
    ytable = cbind(HR = z117556$coefficients[, 'exp(coef)'],
          lower = z117556$conf.int[, 3],
          upper = z117556$conf.int[, 4],
          pvalue = z117556$coefficients[, 'Pr(>|z|)'],
          z = z117556$coefficients[, "z"]) %>% 
          as_tibble() %>% 
          mutate(log_HR = log(HR)) %>% 
          mutate(log_lower = log(lower)) %>% 
          mutate(log_upper = log(upper)) %>% 
          mutate(iterm = factor(c('Age', 'Gender', 'ABC vs GCB', 'ABC vs UNC', 'Clinical Response', 'IPI', 'ECOG 0 vs ECOG 1', 'ECOG 0 vs ECOG 2', 'Extran', 'Myc RNA', 'BCL2 RNA', 'Double expressor', 'Clinical Stage', 'CISD2Risk value'), levels = c('Age', 'Gender', 'ABC vs GCB', 'ABC vs UNC', 'Clinical Response', 'IPI', 'ECOG 0 vs ECOG 1', 'ECOG 0 vs ECOG 2', 'Extran', 'Myc RNA', 'BCL2 RNA', 'Double expressor', 'Clinical Stage', 'CISD2Risk value'))) %>% 
          mutate(log_x = paste(round(log_HR, 2),'(',round(log_lower, 2),',', round(log_upper, 2),')', sep = ''))
    y = ytable  %>% 
        ggplot(aes(y = fct_rev(iterm))) + 
        geom_point(aes(x=log_HR), shape=15, size=5, color = "#E6550DFF") +
        geom_linerange(aes(xmin=log_lower, xmax=log_upper)) +
        geom_vline(xintercept = 0, linetype="dashed") +
        labs(x="Log Hazard Ratio", y="") +
        coord_cartesian(ylim=c(1,14), xlim=c(-2.5, 1.5)) +
        theme_classic() + 
        theme(axis.text.y = element_text(size = 16),
              axis.title.y= element_blank(),
              axis.ticks.y= element_blank(),
              axis.line.y = element_blank())
    ggsave('~/Documents/01CISD2/forrest11755604.png', width = 8, height = 8)
#### Univariate and Multivariate analysis (GSE181063)
    str(x181063)
    x1810631 = x181063 %>% 
        mutate(year = ifelse(age_at_diagnosis <60, 'Younger', 'Older'))  %>% 
        mutate(Gender = ifelse(Sex == 'F', 'Female', 'Male')) %>% 
        select(c('year', 'Gender', 'pred_combine', 
                'state', 'risk', 'osmonth', 'os_status'))
    x1810631$pred_combine = factor(x1810631$pred_combine, levels = c('MHG', 'ABC', 'GCB', 'UNC'))
    coxa = names(x1810631)[1:4]
    xa = data.frame()
    for(i in coxa) {
        expr = x1810631[, i]
        cox = coxph(Surv(osmonth, os_status) ~ expr, x1810631)
        coxsum = summary(cox)
        xa = rbind(xa, cbind(coxi = i, 
            HR = round(coxsum$coefficients[, 'exp(coef)'],2), 
            lower = round(coxsum$conf.int[, 3],2), 
            upper = round(coxsum$conf.int[, 4],2),
            'HR(95%CI)' = paste(round(coxsum$coefficients[, 'exp(coef)'],2), '(', round(coxsum$conf.int[, 3],2), '-', round(coxsum$conf.int[, 4],2), ')', sep = ''),
            pvalue = round(coxsum$coefficients[, 'Pr(>|z|)'],2),
            z = round(coxsum$coefficients[, "z"],2)))}    
    x1810632 = x181063 %>% 
        filter (ldh == 'normal'|ldh == 'raised') %>% 
        select(c('ldh', 'osmonth', 'os_status'))
    coxb = names(x1810632)[1]
    xb = data.frame()
    for(i in coxb) {
        expr = x1810632[, i]
        cox = coxph(Surv(osmonth, os_status) ~ expr, x1810632)
        coxsum = summary(cox)
        xb = rbind(xb, cbind(coxi = i, 
            HR = round(coxsum$coefficients[, 'exp(coef)'],2), 
            lower = round(coxsum$conf.int[, 3],2), 
            upper = round(coxsum$conf.int[, 4],2),
            'HR(95%CI)' = paste(round(coxsum$coefficients[, 'exp(coef)'],2), '(', round(coxsum$conf.int[, 3],2), '-', round(coxsum$conf.int[, 4],2), ')', sep = ''),
            pvalue = round(coxsum$coefficients[, 'Pr(>|z|)'],2),
            z = round(coxsum$coefficients[, "z"],2)))}
    x1810633 = x181063 %>% 
        filter(Stage != 'NA') %>% 
        mutate(xyz = case_when(Stage == 'I' ~ '0', Stage == 'II' ~ '0',
                               Stage == 'III' ~ '1', Stage == 'IV' ~ '1')) %>% 
        select(c('xyz', 'osmonth', 'os_status'))
    coxc = names(x1810633)[1]
    xc = data.frame()
    for(i in coxc) {
        expr = x1810633[, i]
        cox = coxph(Surv(osmonth, os_status) ~ expr, x1810633)
        coxsum = summary(cox)
        xc = rbind(xc, cbind(coxi = i, 
            HR = round(coxsum$coefficients[, 'exp(coef)'],2), 
            lower = round(coxsum$conf.int[, 3],2), 
            upper = round(coxsum$conf.int[, 4],2),
            'HR(95%CI)' = paste(round(coxsum$coefficients[, 'exp(coef)'],2), '(', round(coxsum$conf.int[, 3],2), '-', round(coxsum$conf.int[, 4],2), ')', sep = ''),
            pvalue = round(coxsum$coefficients[, 'Pr(>|z|)'],2),
            z = round(coxsum$coefficients[, "z"],2)))}
    x1810634 = x181063 %>% 
        filter(ipi_score != 'NA') %>% 
        mutate(IPI = ifelse(ipi_score < 2, '0', '1')) %>% 
        select(c('IPI', 'osmonth', 'os_status'))
    coxd = names(x1810634)[1]
    xd = data.frame()
    for(i in coxd) {
        expr = x1810634[, i]
        cox = coxph(Surv(osmonth, os_status) ~ expr, x1810634)
        coxsum = summary(cox)
        xd = rbind(xd, cbind(coxi = i, 
            HR = round(coxsum$coefficients[, 'exp(coef)'],2), 
            lower = round(coxsum$conf.int[, 3],2), 
            upper = round(coxsum$conf.int[, 4],2),
            'HR(95%CI)' = paste(round(coxsum$coefficients[, 'exp(coef)'],2), '(', round(coxsum$conf.int[, 3],2), '-', round(coxsum$conf.int[, 4],2), ')', sep = ''),
            pvalue = round(coxsum$coefficients[, 'Pr(>|z|)'],2),
            z = round(coxsum$coefficients[, "z"],2)))}
    x1810635 = x181063 %>% 
        filter(performance_status_ecog != 'NA') %>% 
        mutate_at('performance_status_ecog', as.numeric) %>% 
        mutate(ECOG = ifelse(performance_status_ecog >= 2, '1', '0')) %>% 
        select(c('ECOG', 'osmonth', 'os_status'))
    coxe = names(x1810635)[1]
    xe = data.frame()
    for(i in coxe) {
        expr = x1810635[, i]
        cox = coxph(Surv(osmonth, os_status) ~ expr, x1810635)
        coxsum = summary(cox)
        xe = rbind(xe, cbind(coxi = i, 
            HR = round(coxsum$coefficients[, 'exp(coef)'],2), 
            lower = round(coxsum$conf.int[, 3],2), 
            upper = round(coxsum$conf.int[, 4],2),
            'HR(95%CI)' = paste(round(coxsum$coefficients[, 'exp(coef)'],2), '(', round(coxsum$conf.int[, 3],2), '-', round(coxsum$conf.int[, 4],2), ')', sep = ''),
            pvalue = round(coxsum$coefficients[, 'Pr(>|z|)'],2),
            z = round(coxsum$coefficients[, "z"],2)))}
    x1810636 = x181063 %>% 
        filter(num_extranodal != 'NA') %>% 
        mutate_at('num_extranodal', as.numeric) %>% 
        mutate(Extranodal = ifelse(num_extranodal >= 1, '1', '0')) %>% 
        select(c('Extranodal', 'osmonth', 'os_status'))
    coxf = names(x1810636)[1]
    xf = data.frame()
    for(i in coxf) {
        expr = x1810636[, i]
        cox = coxph(Surv(osmonth, os_status) ~ expr, x1810636)
        coxsum = summary(cox)
        xf = rbind(xf, cbind(coxi = i, 
            HR = round(coxsum$coefficients[, 'exp(coef)'],2), 
            lower = round(coxsum$conf.int[, 3],2), 
            upper = round(coxsum$conf.int[, 4],2),
            'HR(95%CI)' = paste(round(coxsum$coefficients[, 'exp(coef)'],2), '(', round(coxsum$conf.int[, 3],2), '-', round(coxsum$conf.int[, 4],2), ')', sep = ''),
            pvalue = round(coxsum$coefficients[, 'Pr(>|z|)'],2),
            z = round(coxsum$coefficients[, "z"],2)))}
    x = rbind(xa[1:5,], xb ,xc, xd, xe, xf, xa[6,]) %>% 
        as_tibble() %>% 
        mutate_at(c('HR','lower','upper'), as.numeric) %>% 
        mutate(log_HR = log(HR)) %>% 
        mutate(log_lower = log(lower)) %>% 
        mutate(log_upper = log(upper)) %>% 
        mutate(iterm = factor(c('Age', 'Gender', 'MHG vs ABC', 'MHG vs GCB', 'MHG vs UNC', 'LDH', 'Clinical Stage', 'IPI', 'ECOG 0-1 vs >=2', 'Extranodal', 'CISD2Risk value'), levels = c('Age', 'Gender', 'MHG vs ABC', 'MHG vs GCB', 'MHG vs UNC', 'LDH', 'Clinical Stage', 'IPI', 'ECOG 0-1 vs >=2', 'Extranodal', 'CISD2Risk value'))) %>% 
        mutate(log_x = paste(round(log_HR, 2),'(',round(log_lower, 2),',', round(log_upper, 2),')', sep = '')) %>% 
        ggplot(aes(y = fct_rev(iterm))) + 
        geom_point(aes(x=log_HR), shape=16, size=5, color = '#3182BDFF') +
        geom_linerange(aes(xmin=log_lower, xmax=log_upper)) +
        geom_vline(xintercept = 0, linetype="dashed") +
        labs(x="Log Hazard Ratio", y="") +
        coord_cartesian(ylim=c(1,11), xlim=c(-1.5, 2)) +
        theme_classic() + 
        theme(axis.text.y = element_text(size = 16),
              axis.title.y= element_blank(),
              axis.ticks.y= element_blank(),
              axis.line.y = element_blank())
    ggsave('~/Documents/01CISD2/forrest181063a1.png', width = 8, height = 6)
    z181063 = x181063 %>% 
        mutate(year = ifelse(age_at_diagnosis <60, 'Younger', 'Older'))  %>% 
        mutate(Gender = ifelse(Sex == 'F', 'Female', 'Male')) %>% 
        mutate_at('num_extranodal', as.numeric) %>% 
        mutate(Extranodal = ifelse(num_extranodal >= 1, '1', '0')) %>%
        mutate_at('performance_status_ecog', as.numeric) %>% 
        mutate(ECOG = ifelse(performance_status_ecog >= 2, '1', '0')) %>% 
        mutate(IPI = ifelse(ipi_score < 2, '0', '1')) %>% 
        mutate(xyz = case_when(Stage == 'I' ~ '0', Stage == 'II' ~ '0',
                               Stage == 'III' ~ '1', Stage == 'IV' ~ '1')) %>% 
        mutate(ldh1 = case_when(ldh == 'low' ~ 'normal', ldh == 'normal' ~ 'normal', ldh == 'raised' ~ 'raised')) %>% 
        filter(Stage != 'NA' & IPI != 'NA') 
    z181063$pred_combine = factor(z181063$pred_combine, levels = c('MHG', 'ABC', 'GCB', 'UNC'))
    z181063 = z181063 %>% 
        select(c('year', 'Gender', 'pred_combine', 'ldh1', 'xyz', 
                'IPI', 'ECOG', 'Extranodal', 'state', 'osmonth', 'os_status')) %>% 
        coxph(Surv(osmonth, os_status) ~ year + Gender + pred_combine + ldh1 + xyz +
              IPI + ECOG + Extranodal + state, .) %>% 
        summary(.)
    ytable = cbind(HR = z181063$coefficients[, 'exp(coef)'],
          lower = z181063$conf.int[, 3],
          upper = z181063$conf.int[, 4],
          pvalue = z181063$coefficients[, 'Pr(>|z|)'],
          z = z181063$coefficients[, "z"]) %>% 
          as_tibble() %>% 
          mutate(log_HR = log(HR)) %>% 
          mutate(log_lower = log(lower)) %>% 
          mutate(log_upper = log(upper)) %>% 
          mutate(iterm = factor(c('Age', 'Gender', 'MHG vs ABC', 'MHG vs GCB', 'MHG vs UNC', 'LDH', 'Clinical Stage', 'IPI', 'ECOG 0-1 vs >=2', 'Extranodal', 'CISD2Risk value'), levels = c('Age', 'Gender', 'MHG vs ABC', 'MHG vs GCB', 'MHG vs UNC', 'LDH', 'Clinical Stage', 'IPI', 'ECOG 0-1 vs >=2', 'Extranodal', 'CISD2Risk value'))) %>% 
          mutate(log_x = paste(round(log_HR, 2),'(',round(log_lower, 2),',', round(log_upper, 2),')', sep = ''))
    y = ytable %>% 
        ggplot(aes(y = fct_rev(iterm))) + 
        geom_point(aes(x=log_HR), shape=15, size=5, color = "#E6550DFF") +
        geom_linerange(aes(xmin=log_lower, xmax=log_upper)) +
        geom_vline(xintercept = 0, linetype="dashed") +
        labs(x="Log Hazard Ratio", y="") +
        coord_cartesian(ylim=c(1,11), xlim=c(-1.5, 2)) +
        theme_classic() + 
        theme(axis.text.y = element_text(size = 16),
              axis.title.y= element_blank(),
              axis.ticks.y= element_blank(),
              axis.line.y = element_blank())
    ggsave('~/Documents/01CISD2/forrest181063b1.png', width = 8, height = 6)
## Figure 8 Construciton of nomogram
#### nomogram plot OS in GSE117556
    require(rms)
    png(file = "~/Documents/01CISD2/F_G117556non2x.png", 
    width = 3000, height = 2000, units = "px", bg = "white", res = 250)
    require(rms)
    x117556n = x117556 %>% 
        filter(STAGE != 0) %>% 
        mutate(Stage = case_when(STAGE == 'Stage I' ~ 'stage I-II', STAGE == 'Stage II' ~ 'stage I-II',
                STAGE == 'Stage III' ~ 'stage III-IV', STAGE == 'Stage IV' ~ 'stage III-IV')) %>% 
        mutate(year = ifelse(AGE <60, 'Younger', 'Older')) %>% 
        mutate(ECOG1 = ifelse(ECOG >= 2, 'ECOG PS >= 2', 'ECOG PS < 2')) %>% 
        mutate(ipi = ifelse(IPI_SCORE < 2, 'IPI < 2', 'IPI >= 2')) %>% 
        mutate(LDH1 = ifelse(LDH < 245, 'normal', 'raised')) %>% 
        mutate(COO = ifelse(molecularcoosubtype == 'GCB', 'GCB', 'non-GCB')) %>% 
        mutate(MHG = ifelse(molecularsubtype == 'MHG', 'yes', 'no')) %>% 
        mutate(INV_EXTRANODAL_BAS = ifelse(INV_EXTRANODAL_BAS_IND == 0, 'no', 'yes')) %>% 
        mutate(double = ifelse(expressor_RNA == 'double-expressor', 'yes', 'no')) %>% 
        select(c('year', 'GENDER', 'COO', 'MHG', 'ipi', 'ECOG1', 'Stage', 
                'LDH1', 'INV_EXTRANODAL_BAS', 'double', 
                'state', 'risk', 'PROG_IND', 'PROG_TIME', 'DEATH_IND', 'FU_TIME', 'IPI_SCORE')) %>%
        dplyr::rename('Age' = 'year', 'Gender' = 'GENDER', 'IPI' = 'ipi', 'ECOG_PS' = 'ECOG1', 'Stage' = 'Stage', 'LDH' = 'LDH1', 'Extranodal' = 'INV_EXTRANODAL_BAS', 'Double_expressor' = 'double', 'CISD2Risk' = 'state')
    dd = datadist(x117556n)
    options(datadist="dd")  
    f1 = psm(Surv(FU_TIME, DEATH_IND) ~ Age + Gender + COO + MHG + IPI + ECOG_PS + Stage + LDH + Extranodal + Double_expressor + CISD2Risk,
            x=T, y=T, data = x117556n) 
    surv <- Survival(f1)
    surv1 <- function(x)surv(1*12, lp=x)
    surv2 <- function(x)surv(1*36, lp=x)
    surv3 <- function(x)surv(1*60, lp=x)
    nom1 <- nomogram(f1, fun=list(surv1,surv2,surv3),
               funlabel=c('1 year survival probability',
                          '3 year survival probability',
                          '5 year survival probability'),
               fun.at=c('0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1'))
    plot(nom1, lmgp = .2, cex.axis = .8)
    dev.off()
    rcorrcens(Surv(FU_TIME, DEATH_IND) ~ predict(f1), data =  x117556n)
#### Time-ROC curve OS in GSE117556
    require(timeROC)
    troc_res = timeROC(T = x117556n$FU_TIME,
                       delta = x117556n$DEATH_IND,
                       marker = x117556n$risk,
                       cause = 1,
                       weighting="marginal",
                       times = c(1 * 12, 3 * 12, 5 * 12),
                       ROC = TRUE,
                       iid = TRUE)
    troc_df <- data.frame(
            TP_1year = troc_res$TP[, 1],
            FP_1year = troc_res$FP[, 1],
            TP_3year = troc_res$TP[, 2],
            FP_3year = troc_res$FP[, 2],
            TP_5year = troc_res$TP[, 3],
            FP_5year = troc_res$FP[, 3])
    ggplot(data = troc_df) +
        geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#E6550DFF") +
        geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#3182BDFF") +
        geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#7876B1FF") +
        geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
        theme_bw() +
        scale_x_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1)) +
        scale_y_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1)) +
        annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year  = ", sprintf("%.3f", troc_res$AUC[[1]])), color = "#E6550DFF") +
        annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 year = ", sprintf("%.3f", troc_res$AUC[[2]])), color = "#3182BDFF") +
        annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 year = ", sprintf("%.3f", troc_res$AUC[[3]])), color = "#7876B1FF") +
        labs(x = "False positive rate (%)", y = "True positive rate (%)", title = 'CISD2Risk') +
        theme(axis.text = element_text(size = 8, color = "black", hjust =1),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"))
    ggsave("~/Documents/01CISD2/F_G117556tROCrisk.png", width = 5, height = 5)
    troc_res = timeROC(T = x117556n$FU_TIME,
                       delta = x117556n$DEATH_IND,
                       marker = -predict(f1),
                       cause = 1,
                       weighting="marginal",
                       times = c(1 * 12, 3 * 12, 5 * 12),
                       ROC = TRUE,
                       iid = TRUE)
    troc_df <- data.frame(
            TP_1year = troc_res$TP[, 1],
            FP_1year = troc_res$FP[, 1],
            TP_3year = troc_res$TP[, 2],
            FP_3year = troc_res$FP[, 2],
            TP_5year = troc_res$TP[, 3],
            FP_5year = troc_res$FP[, 3])
    ggplot(data = troc_df) +
        geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#E6550DFF") +
        geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#3182BDFF") +
        geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#7876B1FF") +
        geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
        theme_bw() +
        scale_x_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1)) +
        scale_y_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1)) +
        annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year  = ", sprintf("%.3f", troc_res$AUC[[1]])), color = "#E6550DFF") +
        annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 year = ", sprintf("%.3f", troc_res$AUC[[2]])), color = "#3182BDFF") +
        annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 year = ", sprintf("%.3f", troc_res$AUC[[3]])), color = "#7876B1FF") +
        labs(x = "False positive rate (%)", y = "True positive rate (%)", title = 'nomogram') +
        theme(axis.text = element_text(size = 8, color = "black", hjust =1),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"))
    ggsave("~/Documents/01CISD2/F_G117556tROCnomo.png", width = 5, height = 5)
    troc_res = timeROC(T = x117556n$FU_TIME,
                       delta = x117556n$DEATH_IND,
                       marker = x117556n$IPI_SCORE,
                       cause = 1,
                       weighting="marginal",
                       times = c(1 * 12, 3 * 12, 5 * 12),
                       ROC = TRUE,
                       iid = TRUE)
    troc_df <- data.frame(
            TP_1year = troc_res$TP[, 1],
            FP_1year = troc_res$FP[, 1],
            TP_3year = troc_res$TP[, 2],
            FP_3year = troc_res$FP[, 2],
            TP_5year = troc_res$TP[, 3],
            FP_5year = troc_res$FP[, 3])
    ggplot(data = troc_df) +
        geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#E6550DFF") +
        geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#3182BDFF") +
        geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#7876B1FF") +
        geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
        theme_bw() +
        scale_x_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1)) +
        scale_y_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1)) +
        annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year  = ", sprintf("%.3f", troc_res$AUC[[1]])), color = "#E6550DFF") +
        annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 year = ", sprintf("%.3f", troc_res$AUC[[2]])), color = "#3182BDFF") +
        annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 year = ", sprintf("%.3f", troc_res$AUC[[3]])), color = "#7876B1FF") +
        labs(x = "False positive rate (%)", y = "True positive rate (%)", title = 'IPI') +
        theme(axis.text = element_text(size = 8, color = "black", hjust =1),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"))
    ggsave("~/Documents/01CISD2/F_G117556tROCIPI.png", width = 5, height = 5)
#### Calibrate curve OS in GSE117556
    f2<-psm(Surv(FU_TIME, DEATH_IND) ~ Age + Gender + COO + MHG + IPI + ECOG_PS + Stage + LDH + Extranodal + Double_expressor + CISD2Risk,
            x=T, y=T, data = x117556n, dist='lognormal') 
    cal2<-calibrate(f2, cmethod="KM", method="boot", u= 12, m = 100, B=1000) 
    cal22 = cal2[,c("mean.predicted","KM","std.err")] %>% 
            as_tibble() %>% 
            mutate(min = KM - std.err) %>% 
            mutate(max = KM + std.err)
    cal20 <- data.frame(prob = attr(cal2,"predicted")) 
    f3<-psm(Surv(FU_TIME, DEATH_IND) ~  Age + Gender + COO + MHG + IPI + ECOG_PS + Stage + LDH + Extranodal + Double_expressor + CISD2Risk,
            x=T, y=T, data = x117556n, dist='lognormal') 
    cal3<-calibrate(f3, cmethod="KM", method="boot", u= 36, m = 100, B=1000) 
    cal33 = cal3[,c("mean.predicted","KM","std.err")] %>% 
            as_tibble() %>% 
            mutate(min = KM - std.err) %>% 
            mutate(max = KM + std.err)
    cal30 <- data.frame(prob = attr(cal3,"predicted")) 
    f4<-psm(Surv(FU_TIME, DEATH_IND) ~  Age + Gender + COO + MHG + IPI + ECOG_PS + Stage + LDH + Extranodal + Double_expressor + CISD2Risk,
            x=T, y=T, data = x117556n, dist='lognormal') 
    cal4<-calibrate(f4, cmethod="KM", method="boot", u= 60, m = 100, B=1000) 
    cal44 = cal4[,c("mean.predicted","KM","std.err")] %>% 
            as_tibble() %>% 
            mutate(min = KM - std.err) %>% 
            mutate(max = KM + std.err)
    cal40 <- data.frame(prob = attr(cal4,"predicted"))
    ggplot() +
        geom_line(data=cal22, aes(x=mean.predicted, y=KM), linewidth = 0.5, color = "#E6550DFF") +
        geom_point(data=cal22, aes(x=mean.predicted, y=KM),col="#E6550DFF",size = 2, shape = 15) +
        geom_errorbar(data=cal22, aes(x=mean.predicted, y=KM, ymin = min, ymax = max), color = "#E6550DFF", width = 0.01) +
        geom_line(data=cal33, aes(x=mean.predicted, y=KM), linewidth = 0.5, color = "#3182BDFF") +
        geom_point(data=cal33, aes(x=mean.predicted, y=KM),col="#3182BDFF",size = 2, shape = 16) +
        geom_errorbar(data=cal33, aes(x=mean.predicted, y=KM, ymin = min, ymax = max), color = "#3182BDFF", width = 0.01) +
        geom_line(data=cal44, aes(x=mean.predicted, y=KM), linewidth = 0.5, color = "#7876B1FF") +
        geom_point(data=cal44, aes(x=mean.predicted, y=KM),col="#7876B1FF",size = 2, shape = 17) +
        geom_errorbar(data=cal44, aes(x=mean.predicted, y=KM, ymin = min, ymax = max), color = "#7876B1FF", width = 0.01) +
        geom_abline(slope = 1,intercept = 0,lty=2) + 
        labs(x = "Predited Probability (%)", y = "Observed Probability (%)", title = '') +
        scale_x_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1)) +
        scale_y_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1), sec.axis = sec_axis(~.*400, 
                           name = '', breaks = c(50,100,150,200,250,300,350))) +
        geom_histogram(data = cal20, aes(x=prob, after_stat(count/200)), binwidth = 1/200, alpha = 0.75, fill="#E6550DFF") +
        geom_histogram(data = cal30, aes(x=prob, after_stat(count/200)), binwidth = 1/200, alpha = 0.75, fill="#3182BDFF") +
        geom_histogram(data = cal40, aes(x=prob, after_stat(count/200)), binwidth = 1/200, alpha = 0.75, fill="#7876B1FF") +
        theme_bw() +
        annotate('text',x=0.2,y=0.95,label = '- 1 Year', size = 4, color = "#E6550DFF") +
        annotate('text',x=0.2,y=0.90,label = '- 3 Year', size = 4, color = "#3182BDFF") +
        annotate('text',x=0.2,y=0.85,label = '- 5 Year', size = 4, color = "#7876B1FF") +
        theme(legend.position="top",
        axis.text = element_text(size = 11, color = "black"),
        axis.title.x = element_text(size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
        axis.title.y = element_text(size = 14, color = "black", margin = margin(c(0, 15, 0, 0))))
    ggsave("~/Documents/01CISD2/F_G117556cabayy.png", width = 5.5, height = 5)
#### nomogram plot PFS in GSE117556
    require(rms)
    png(file = "~/Documents/01CISD2/F_G117556non2xx.png", 
    width = 3000, height = 2000, units = "px", bg = "white", res = 250)
    require(rms)
    x117556n = x117556 %>% 
        filter(STAGE != 0) %>% 
        mutate(Stage = case_when(STAGE == 'Stage I' ~ 'stage I-II', STAGE == 'Stage II' ~ 'stage I-II',
                STAGE == 'Stage III' ~ 'stage III-IV', STAGE == 'Stage IV' ~ 'stage III-IV')) %>% 
        mutate(year = ifelse(AGE <60, 'Younger', 'Older')) %>% 
        mutate(ECOG1 = ifelse(ECOG >= 2, 'ECOG PS >= 2', 'ECOG PS < 2')) %>% 
        mutate(ipi = ifelse(IPI_SCORE < 2, 'IPI < 2', 'IPI >= 2')) %>% 
        mutate(LDH1 = ifelse(LDH < 245, 'normal', 'raised')) %>% 
        mutate(COO = ifelse(molecularcoosubtype == 'GCB', 'GCB', 'non-GCB')) %>% 
        mutate(MHG = ifelse(molecularsubtype == 'MHG', 'yes', 'no')) %>% 
        mutate(INV_EXTRANODAL_BAS = ifelse(INV_EXTRANODAL_BAS_IND == 0, 'no', 'yes')) %>% 
        mutate(double = ifelse(expressor_RNA == 'double-expressor', 'yes', 'no')) %>% 
        select(c('year', 'GENDER', 'COO', 'MHG', 'ipi', 'ECOG1', 'Stage', 
                'LDH1', 'INV_EXTRANODAL_BAS', 'double', 
                'state', 'risk', 'PROG_IND', 'PROG_TIME', 'DEATH_IND', 'FU_TIME', 'IPI_SCORE')) %>%
        dplyr::rename('Age' = 'year', 'Gender' = 'GENDER', 'IPI' = 'ipi', 'ECOG_PS' = 'ECOG1', 'Stage' = 'Stage', 'LDH' = 'LDH1', 'Extranodal' = 'INV_EXTRANODAL_BAS', 'Double_expressor' = 'double', 'CISD2Risk' = 'state')
    dd = datadist(x117556n)
    options(datadist="dd")  
    f1 = psm(Surv(PROG_TIME, PROG_IND) ~ Age + Gender + COO + MHG + IPI + ECOG_PS + Stage + LDH + Extranodal + Double_expressor + CISD2Risk,
            x=T, y=T, data = x117556n) 
    surv <- Survival(f1)
    surv1 <- function(x)surv(1*12, lp=x)
    surv2 <- function(x)surv(1*36, lp=x)
    surv3 <- function(x)surv(1*60, lp=x)
    nom1 <- nomogram(f1, fun=list(surv1,surv2,surv3),
               funlabel=c('1 year survival probability',
                          '3 year survival probability',
                          '5 year survival probability'),
               fun.at=c('0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1'))
    plot(nom1, lmgp = .2, cex.axis = .8)
    dev.off()
#### Time-ROC curve PFS in GSE117556
    require(timeROC)
    x117556n$nomo = predict(f1)
    troc_res = timeROC(T = x117556n$PROG_TIME,
                       delta = x117556n$PROG_IND,
                       marker = x117556n$risk,
                       cause = 1,
                       weighting="marginal",
                       times = c(1 * 12, 3 * 12, 5 * 12),
                       ROC = TRUE,
                       iid = TRUE)
    troc_df <- data.frame(
            TP_1year = troc_res$TP[, 1],
            FP_1year = troc_res$FP[, 1],
            TP_3year = troc_res$TP[, 2],
            FP_3year = troc_res$FP[, 2],
            TP_5year = troc_res$TP[, 3],
            FP_5year = troc_res$FP[, 3])
    ggplot(data = troc_df) +
        geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#E6550DFF") +
        geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#3182BDFF") +
        geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#7876B1FF") +
        geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
        theme_bw() +
        scale_x_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1)) +
        scale_y_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1)) +
        annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year  = ", sprintf("%.3f", troc_res$AUC[[1]])), color = "#E6550DFF") +
        annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 year = ", sprintf("%.3f", troc_res$AUC[[2]])), color = "#3182BDFF") +
        annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 year = ", sprintf("%.3f", troc_res$AUC[[3]])), color = "#7876B1FF") +
        labs(x = "False positive rate (%)", y = "True positive rate (%)", title = 'CISD2Risk') +
        theme(axis.text = element_text(size = 8, color = "black", hjust =1),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"))
    ggsave("~/Documents/01CISD2/F_G117556tROCriskp.png", width = 5, height = 5)
    troc_res = timeROC(T = x117556n$PROG_TIME,
                       delta = x117556n$PROG_IND,
                       marker = -predict(f1),
                       cause = 1,
                       weighting="marginal",
                       times = c(1 * 12, 3 * 12, 5 * 12),
                       ROC = TRUE,
                       iid = TRUE)
    troc_df <- data.frame(
            TP_1year = troc_res$TP[, 1],
            FP_1year = troc_res$FP[, 1],
            TP_3year = troc_res$TP[, 2],
            FP_3year = troc_res$FP[, 2],
            TP_5year = troc_res$TP[, 3],
            FP_5year = troc_res$FP[, 3])
    ggplot(data = troc_df) +
        geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#E6550DFF") +
        geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#3182BDFF") +
        geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#7876B1FF") +
        geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
        theme_bw() +
        scale_x_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1)) +
        scale_y_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1)) +
        annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year  = ", sprintf("%.3f", troc_res$AUC[[1]])), color = "#E6550DFF") +
        annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 year = ", sprintf("%.3f", troc_res$AUC[[2]])), color = "#3182BDFF") +
        annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 year = ", sprintf("%.3f", troc_res$AUC[[3]])), color = "#7876B1FF") +
        labs(x = "False positive rate (%)", y = "True positive rate (%)", title = 'nomogram') +
        theme(axis.text = element_text(size = 8, color = "black", hjust =1),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"))
    ggsave("~/Documents/01CISD2/F_G117556tROCnomop.png", width = 5, height = 5)
    troc_res = timeROC(T = x117556n$PROG_TIME,
                       delta = x117556n$PROG_IND,
                       marker = x117556n$IPI_SCORE,
                       cause = 1,
                       weighting="marginal",
                       times = c(1 * 12, 3 * 12, 5 * 12),
                       ROC = TRUE,
                       iid = TRUE)
    troc_df <- data.frame(
            TP_1year = troc_res$TP[, 1],
            FP_1year = troc_res$FP[, 1],
            TP_3year = troc_res$TP[, 2],
            FP_3year = troc_res$FP[, 2],
            TP_5year = troc_res$TP[, 3],
            FP_5year = troc_res$FP[, 3])
    ggplot(data = troc_df) +
        geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#E6550DFF") +
        geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#3182BDFF") +
        geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#7876B1FF") +
        geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
        theme_bw() +
        scale_x_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1)) +
        scale_y_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1)) +
        annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year  = ", sprintf("%.3f", troc_res$AUC[[1]])), color = "#E6550DFF") +
        annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 year = ", sprintf("%.3f", troc_res$AUC[[2]])), color = "#3182BDFF") +
        annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 year = ", sprintf("%.3f", troc_res$AUC[[3]])), color = "#7876B1FF") +
        labs(x = "False positive rate (%)", y = "True positive rate (%)", title = 'IPI') +
        theme(axis.text = element_text(size = 8, color = "black", hjust =1),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"))
    ggsave("~/Documents/01CISD2/F_G117556tROCIPIp.png", width = 5, height = 5)
#### Calibrate curve PFS in GSE117556
    f2<-psm(Surv(PROG_TIME, PROG_IND) ~ Age + Gender + COO + MHG + IPI + ECOG_PS + Stage + LDH + Extranodal + Double_expressor + CISD2Risk,
            x=T, y=T, data = x117556n, dist='lognormal') 
    cal2<-calibrate(f2, cmethod="KM", method="boot", u= 12, m = 100, B=1000) 
    cal22 = cal2[,c("mean.predicted","KM","std.err")] %>% 
            as_tibble() %>% 
            mutate(min = KM - std.err) %>% 
            mutate(max = KM + std.err)
    cal20 <- data.frame(prob = attr(cal2,"predicted")) 
    f3<-psm(Surv(PROG_TIME, PROG_IND) ~  Age + Gender + COO + MHG + IPI + ECOG_PS + Stage + LDH + Extranodal + Double_expressor + CISD2Risk,
            x=T, y=T, data = x117556n, dist='lognormal') 
    cal3<-calibrate(f3, cmethod="KM", method="boot", u= 36, m = 100, B=1000) 
    cal33 = cal3[,c("mean.predicted","KM","std.err")] %>% 
            as_tibble() %>% 
            mutate(min = KM - std.err) %>% 
            mutate(max = KM + std.err)
    cal30 <- data.frame(prob = attr(cal3,"predicted")) 
    f4<-psm(Surv(PROG_TIME, PROG_IND) ~  Age + Gender + COO + MHG + IPI + ECOG_PS + Stage + LDH + Extranodal + Double_expressor + CISD2Risk,
            x=T, y=T, data = x117556n, dist='lognormal') 
    cal4<-calibrate(f4, cmethod="KM", method="boot", u= 60, m = 100, B=1000) 
    cal44 = cal4[,c("mean.predicted","KM","std.err")] %>% 
            as_tibble() %>% 
            mutate(min = KM - std.err) %>% 
            mutate(max = KM + std.err)
    cal40 <- data.frame(prob = attr(cal4,"predicted"))
    ggplot() +
        geom_line(data=cal22, aes(x=mean.predicted, y=KM), linewidth = 0.5, color = "#E6550DFF") +
        geom_point(data=cal22, aes(x=mean.predicted, y=KM),col="#E6550DFF",size = 2, shape = 15) +
        geom_errorbar(data=cal22, aes(x=mean.predicted, y=KM, ymin = min, ymax = max), color = "#E6550DFF", width = 0.01) +
        geom_line(data=cal33, aes(x=mean.predicted, y=KM), linewidth = 0.5, color = "#3182BDFF") +
        geom_point(data=cal33, aes(x=mean.predicted, y=KM),col="#3182BDFF",size = 2, shape = 16) +
        geom_errorbar(data=cal33, aes(x=mean.predicted, y=KM, ymin = min, ymax = max), color = "#3182BDFF", width = 0.01) +
        geom_line(data=cal44, aes(x=mean.predicted, y=KM), linewidth = 0.5, color = "#7876B1FF") +
        geom_point(data=cal44, aes(x=mean.predicted, y=KM),col="#7876B1FF",size = 2, shape = 17) +
        geom_errorbar(data=cal44, aes(x=mean.predicted, y=KM, ymin = min, ymax = max), color = "#7876B1FF", width = 0.01) +
        geom_abline(slope = 1,intercept = 0,lty=2) + 
        labs(x = "Predited Probability (%)", y = "Observed Probability (%)", title = '') +
        scale_x_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1)) +
        scale_y_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1), sec.axis = sec_axis(~.*400, 
                           name = '', breaks = c(50,100,150,200,250,300,350))) +
        geom_histogram(data = cal20, aes(x=prob, after_stat(count/200)), binwidth = 1/200, alpha = 0.75, fill="#E6550DFF") +
        geom_histogram(data = cal30, aes(x=prob, after_stat(count/200)), binwidth = 1/200, alpha = 0.75, fill="#3182BDFF") +
        geom_histogram(data = cal40, aes(x=prob, after_stat(count/200)), binwidth = 1/200, alpha = 0.75, fill="#7876B1FF") +
        theme_bw() +
        annotate('text',x=0.2,y=0.95,label = '- 1 Year', size = 4, color = "#E6550DFF") +
        annotate('text',x=0.2,y=0.90,label = '- 3 Year', size = 4, color = "#3182BDFF") +
        annotate('text',x=0.2,y=0.85,label = '- 5 Year', size = 4, color = "#7876B1FF") +
        theme(legend.position="top",
        axis.text = element_text(size = 11, color = "black"),
        axis.title.x = element_text(size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
        axis.title.y = element_text(size = 14, color = "black", margin = margin(c(0, 15, 0, 0))))
    ggsave("~/Documents/01CISD2/F_G117556cabayyp.png", width = 5.5, height = 5)
#### Nomogram model OS in GSE181063
    require(rms)
    png(file = "~/Documents/01CISD2/F_G181063non1x.png", 
    width = 3000, height = 2000, units = "px", bg = "white", res = 250)
    x181063n = x181063 %>% 
        mutate(year = ifelse(age_at_diagnosis < 60, 'Younger', 'Older'))  %>% 
        mutate(Gender = ifelse(Sex == 'F', 'Female', 'Male')) %>% 
        #mutate_at('num_extranodal', as.numeric) %>% 
        #mutate(Extranodal = ifelse(num_extranodal >= 1, 'extranodal', 'non extrandal')) %>%
        mutate_at('performance_status_ecog', as.numeric) %>% 
        mutate(ECOG = ifelse(performance_status_ecog >= 2, 'PS >= 2', 'PS < 2')) %>% 
        mutate(IPI = ifelse(ipi_score < 2, 'Low IPI', 'High IPI')) %>% 
        mutate(xyz = case_when(Stage == 'I' ~ 'Stage I-II', Stage == 'II' ~ 'Stage I-II',
                               Stage == 'III' ~ 'Stage III-IV', Stage == 'IV' ~ 'Stage III-IV')) %>% 
        mutate(ldh1 = case_when(ldh == 'low' ~ 'normal LDH', ldh == 'normal' ~ 'normal LDH', ldh == 'raised' ~ 'raised LDH')) %>% 
        filter(Stage != 'NA' & IPI != 'NA') %>% 
        select(c('year', 'Gender', 'pred_combine', 'ldh1', 'xyz', 
                'IPI', 'ECOG', 'state', 'osmonth', 'os_status', 'risk')) %>%
        dplyr::rename('Age'='year', 'Molecular_subtype'='pred_combine', 'LDH' = 'ldh1',
                'Clinical_Stage'='xyz', 'CISD2Risk'='state')
    dd = datadist(x181063n)
    options(datadist="dd")  
    f1 = psm(Surv(osmonth, os_status) ~ Age + Gender + Molecular_subtype + LDH + Clinical_Stage + 
            IPI + ECOG + CISD2Risk, x=T, y=T, data = x181063n) 
    surv <- Survival(f1)
    surv1 <- function(x)surv(1*12, lp=x)
    surv2 <- function(x)surv(1*36, lp=x)
    surv3 <- function(x)surv(1*60, lp=x)
    nom1 <- nomogram(f1, fun=list(surv1,surv2,surv3),
               funlabel=c('1-Year survival probability',
                          '3-Years survival probability',
                          '5-Years survival probability'),
               fun.at=c('0.9','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1'))
    plot(nom1, lmgp = .2, cex.axis = .8)
    dev.off()
    rcorrcens(Surv(osmonth, os_status) ~ predict(f1), data =  x181063n)
    require(timeROC)
    troc_res = timeROC(T = x181063n$osmonth,
                       delta = x181063n$os_status,
                       marker = x181063n$risk,
                       cause = 1,
                       weighting="marginal",
                       times = c(1 * 12, 3 * 12, 5 * 12),
                       ROC = TRUE,
                       iid = TRUE)
    troc_df <- data.frame(
            TP_1year = troc_res$TP[, 1],
            FP_1year = troc_res$FP[, 1],
            TP_3year = troc_res$TP[, 2],
            FP_3year = troc_res$FP[, 2],
            TP_5year = troc_res$TP[, 3],
            FP_5year = troc_res$FP[, 3])
    ggplot(data = troc_df) +
        geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#E6550DFF") +
        geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#3182BDFF") +
        geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#7876B1FF") +
        geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
        theme_bw() +
        scale_x_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1)) +
        scale_y_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1)) +
        annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year  = ", sprintf("%.3f", troc_res$AUC[[1]])), color = "#E6550DFF") +
        annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 year = ", sprintf("%.3f", troc_res$AUC[[2]])), color = "#3182BDFF") +
        annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 year = ", sprintf("%.3f", troc_res$AUC[[3]])), color = "#7876B1FF") +
        labs(x = "False positive rate (%)", y = "True positive rate (%)", title = 'CISD2Risk') +
        theme(axis.text = element_text(size = 8, color = "black", hjust =1),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"))
    ggsave("~/Documents/01CISD2/F_G181063tROCriskp.png", width = 5, height = 5)
    troc_res = timeROC(T = x181063n$osmonth,
                       delta = x181063n$os_status,
                       marker = -predict(f1),
                       cause = 1,
                       weighting="marginal",
                       times = c(1 * 12, 3 * 12, 5 * 12),
                       ROC = TRUE,
                       iid = TRUE)
    troc_df <- data.frame(
            TP_1year = troc_res$TP[, 1],
            FP_1year = troc_res$FP[, 1],
            TP_3year = troc_res$TP[, 2],
            FP_3year = troc_res$FP[, 2],
            TP_5year = troc_res$TP[, 3],
            FP_5year = troc_res$FP[, 3])
    ggplot(data = troc_df) +
        geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#E6550DFF") +
        geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#3182BDFF") +
        geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#7876B1FF") +
        geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
        theme_bw() +
        scale_x_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1)) +
        scale_y_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1)) +
        annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year  = ", sprintf("%.3f", troc_res$AUC[[1]])), color = "#E6550DFF") +
        annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 year = ", sprintf("%.3f", troc_res$AUC[[2]])), color = "#3182BDFF") +
        annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 year = ", sprintf("%.3f", troc_res$AUC[[3]])), color = "#7876B1FF") +
        labs(x = "False positive rate (%)", y = "True positive rate (%)", title = 'nomogram') +
        theme(axis.text = element_text(size = 8, color = "black", hjust =1),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"))
    ggsave("~/Documents/01CISD2/F_G181063tROCnomop.png", width = 5, height = 5)
#### Calibrate curve
    f2<-psm(Surv(FU_TIME, DEATH_IND) ~ Age + Gender + Molecular_subtype + LDH + Clinical_Stage + IPI + ECOG + 
            CISD2Risk,
            x=T, y=T, data = x117556n, dist='lognormal') 
    cal2<-calibrate(f2, cmethod="KM", method="boot", u= 12, m = 100, B=1000) 
    cal22 = cal2[,c("mean.predicted","KM","std.err")] %>% 
            as_tibble() %>% 
            mutate(min = KM - std.err) %>% 
            mutate(max = KM + std.err)
    cal20 <- data.frame(prob = attr(cal2,"predicted")) 
    f3<-psm(Surv(FU_TIME, DEATH_IND) ~  Age + Gender + Molecular_subtype + LDH + Clinical_Stage + IPI + ECOG + 
            CISD2Risk,
            x=T, y=T, data = x117556n, dist='lognormal') 
    cal3<-calibrate(f3, cmethod="KM", method="boot", u= 36, m = 100, B=1000) 
    cal33 = cal3[,c("mean.predicted","KM","std.err")] %>% 
            as_tibble() %>% 
            mutate(min = KM - std.err) %>% 
            mutate(max = KM + std.err)
    cal30 <- data.frame(prob = attr(cal3,"predicted")) 
    f4<-psm(Surv(FU_TIME, DEATH_IND) ~  Age + Gender + Molecular_subtype + LDH + Clinical_Stage + IPI + ECOG + 
            CISD2Risk,
            x=T, y=T, data = x117556n, dist='lognormal') 
    cal4<-calibrate(f4, cmethod="KM", method="boot", u= 60, m = 100, B=1000) 
    cal44 = cal4[,c("mean.predicted","KM","std.err")] %>% 
            as_tibble() %>% 
            mutate(min = KM - std.err) %>% 
            mutate(max = KM + std.err)
    cal40 <- data.frame(prob = attr(cal4,"predicted"))
    ggplot() +
        geom_line(data=cal22, aes(x=mean.predicted, y=KM), linewidth = 0.5, color = "#E6550DFF") +
        geom_point(data=cal22, aes(x=mean.predicted, y=KM),col="#E6550DFF",size = 2, shape = 15) +
        geom_errorbar(data=cal22, aes(x=mean.predicted, y=KM, ymin = min, ymax = max), color = "#E6550DFF", width = 0.01) +
        geom_line(data=cal33, aes(x=mean.predicted, y=KM), linewidth = 0.5, color = "#3182BDFF") +
        geom_point(data=cal33, aes(x=mean.predicted, y=KM),col="#3182BDFF",size = 2, shape = 16) +
        geom_errorbar(data=cal33, aes(x=mean.predicted, y=KM, ymin = min, ymax = max), color = "#3182BDFF", width = 0.01) +
        geom_line(data=cal44, aes(x=mean.predicted, y=KM), linewidth = 0.5, color = "#7876B1FF") +
        geom_point(data=cal44, aes(x=mean.predicted, y=KM),col="#7876B1FF",size = 2, shape = 17) +
        geom_errorbar(data=cal44, aes(x=mean.predicted, y=KM, ymin = min, ymax = max), color = "#7876B1FF", width = 0.01) +
        geom_abline(slope = 1,intercept = 0,lty=2) + 
        labs(x = "Predited Probability (%)", y = "Observed Probability (%)", title = 'GSE117556') +
        scale_x_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1)) +
        scale_y_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1), sec.axis = sec_axis(~.*400, 
                           name = '', breaks = c(50,100,150,200,250,300,350))) +
        geom_histogram(data = cal20, aes(x=prob, after_stat(count/200)), binwidth = 1/200, alpha = 0.75, fill="#E6550DFF") +
        geom_histogram(data = cal30, aes(x=prob, after_stat(count/200)), binwidth = 1/200, alpha = 0.75, fill="#3182BDFF") +
        geom_histogram(data = cal40, aes(x=prob, after_stat(count/200)), binwidth = 1/200, alpha = 0.75, fill="#7876B1FF") +
        theme_bw() +
        annotate('text',x=0.2,y=0.95,label = '- 1 Year', size = 4, color = "#E6550DFF") +
        annotate('text',x=0.2,y=0.90,label = '- 3 Year', size = 4, color = "#3182BDFF") +
        annotate('text',x=0.2,y=0.85,label = '- 5 Year', size = 4, color = "#7876B1FF") +
        theme(legend.position="top",
        axis.text = element_text(size = 11, color = "black"),
        axis.title.x = element_text(size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
        axis.title.y = element_text(size = 14, color = "black", margin = margin(c(0, 15, 0, 0))))
    ggsave("~/Documents/01CISD2/F_G117556cabax.png", width = 5.5, height = 5)
    f2<-psm(Surv(osmonth, os_status) ~ Age + Gender + Molecular_subtype + LDH + Clinical_Stage + IPI + ECOG + 
            CISD2Risk, x=T, y=T, data = x181063n, dist='lognormal') 
    cal2<-calibrate(f2, cmethod="KM", method="boot", u= 12, m = 100, B=1000) 
    cal22 = cal2[,c("mean.predicted","KM","std.err")] %>% 
            as_tibble() %>% 
            mutate(min = KM - std.err) %>% 
            mutate(max = KM + std.err)
    cal20 <- data.frame(prob = attr(cal2,"predicted")) 
    f3<-psm(Surv(osmonth, os_status) ~ Age + Gender + Molecular_subtype + LDH + Clinical_Stage + IPI + ECOG + 
            CISD2Risk, x=T, y=T, data = x181063n, dist='lognormal') 
    cal3<-calibrate(f3, cmethod="KM", method="boot", u= 36, m = 100, B=1000) 
    cal33 = cal3[,c("mean.predicted","KM","std.err")] %>% 
            as_tibble() %>% 
            mutate(min = KM - std.err) %>% 
            mutate(max = KM + std.err)
    cal30 <- data.frame(prob = attr(cal3,"predicted")) 
    f4<-psm(Surv(osmonth, os_status) ~ Age + Gender + Molecular_subtype + LDH + Clinical_Stage + IPI + ECOG + 
            CISD2Risk, x=T, y=T, data = x181063n, dist='lognormal') 
    cal4<-calibrate(f4, cmethod="KM", method="boot", u= 60, m = 100, B=1000) 
    cal44 = cal4[,c("mean.predicted","KM","std.err")] %>% 
            as_tibble() %>% 
            mutate(min = KM - std.err) %>% 
            mutate(max = KM + std.err)
    cal40 <- data.frame(prob = attr(cal4,"predicted"))
    ggplot() +
        geom_line(data=cal22, aes(x=mean.predicted, y=KM), linewidth = 0.5, color = "#E6550DFF") +
        geom_point(data=cal22, aes(x=mean.predicted, y=KM),col="#E6550DFF",size = 2, shape = 15) +
        geom_errorbar(data=cal22, aes(x=mean.predicted, y=KM, ymin = min, ymax = max), color = "#E6550DFF", width = 0.01) +
        geom_line(data=cal33, aes(x=mean.predicted, y=KM), linewidth = 0.5, color = "#3182BDFF") +
        geom_point(data=cal33, aes(x=mean.predicted, y=KM),col="#3182BDFF",size = 2, shape = 16) +
        geom_errorbar(data=cal33, aes(x=mean.predicted, y=KM, ymin = min, ymax = max), color = "#3182BDFF", width = 0.01) +
        geom_line(data=cal44, aes(x=mean.predicted, y=KM), linewidth = 0.5, color = "#7876B1FF") +
        geom_point(data=cal44, aes(x=mean.predicted, y=KM),col="#7876B1FF",size = 2, shape = 17) +
        geom_errorbar(data=cal44, aes(x=mean.predicted, y=KM, ymin = min, ymax = max), color = "#7876B1FF", width = 0.01) +
        geom_abline(slope = 1,intercept = 0,lty=2) + 
        labs(x = "Predited Probability (%)", y = "Observed Probability (%)", title = 'GSE181063') +
        scale_x_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1)) +
        scale_y_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0), breaks=round(seq(0,1,length.out=5),1),
                           sec.axis = sec_axis(~.*400, 
                           name = '', breaks = c(50,100,150,200,250,300,350))) +
        geom_histogram(data = cal20, aes(x=prob, after_stat(count/200)), binwidth = 1/200, alpha = 0.75, fill="#E6550DFF") +
        geom_histogram(data = cal30, aes(x=prob, after_stat(count/200)), binwidth = 1/200, alpha = 0.75, fill="#3182BDFF") +
        geom_histogram(data = cal40, aes(x=prob, after_stat(count/200)), binwidth = 1/200, alpha = 0.75, fill="#7876B1FF") +
        theme_bw() +
        annotate('text',x=0.2,y=0.95,label = '- 1 Year', size = 4, color = "#E6550DFF") +
        annotate('text',x=0.2,y=0.90,label = '- 3 Year', size = 4, color = "#3182BDFF") +
        annotate('text',x=0.2,y=0.85,label = '- 5 Year', size = 4, color = "#7876B1FF") +
        theme(legend.position="top",
        axis.text = element_text(size = 11, color = "black"),
        axis.title.x = element_text(size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
        axis.title.y = element_text(size = 14, color = "black", margin = margin(c(0, 15, 0, 0))))
    ggsave("~/Documents/01CISD2/F_G181063cabax.png", width = 5.5, height = 5)
