#### Packages needed in this study.
    library(tidyverse)
    library(meta)
    library(mada)
#### Bibliometrix analysis
    library('bibliometrix')
    biblioshiny()
#### Meta analysis
    x00 = c('study', 'LncRNA', 'TP', 'FP', 'FN', 'TN')
    x01 = c('Ge2019','ANRIL', 96, 84, 42, 56)
    x02 = c('Ming2019', 'NEAT1', 70, 46, 20, 44)
    x03 = c('Du2020', 'PACER', 37, 1, 18, 44)
    x04 = c('Wang2020a', 'LINC00987', 22, 3, 7, 30)
    x05 = c('Wang2020b', 'PVT1', 62, 21, 18, 59)
    x06 = c('Hao2021', 'OIP5-AS1', 54, 10,  8, 45)
    x07 = c('Liu2021', 'CASC2', 39,  4, 11, 46)
    x08 = c('Zhao2021', 'LUCAT1', 67, 40,  3, 30)
    x09 = c('Mao2022', 'PVT1', 269,  7 , 31 , 58)
    x10 = c('Wei2022', 'MEG3', 33, 10, 10, 16)
    x11 = c('Xu2022', 'LINC00599', 23, 2, 7, 28)
    x12 = c('Xiao2023', 'LINC00612', 24, 2, 8, 32)
    x13 = c('Ge2019', 'ANRIL', 122, 70, 14, 70)
    x14 = c('Ming2019', 'NEAT1', 84, 28, 6, 62)
    x15 = c('Qi2019a', 'NR-026690', 55, 19,  1, 16)
    x16 = c('Qi2019b', 'ENST00000447867', 47, 11,  9, 24)
    x17 = c('Liu2020', 'MALAT1', 119, 20, 1, 100)
    x18 = c('Wang2020b', 'PVT1', 67, 5, 13, 75)
    x19 = c('Ge2019', 'ANRIL', 122, 77, 14, 61)
    x20 = c('Ming2019', 'NEAT1', 79, 28, 11, 56)
    x21 = c('Qi2019a', 'NR-026690', 33, 7, 18, 44)
    x22 = c('Qi2019b', 'ENST00000447867', 19, 5, 32, 46)
    x23 = c('Liu2020', 'MALAT1', 77, 20, 43, 100)
    x24 = c('Wang2020', 'PVT1', 40, 8, 40, 72)
    x25 = c('Wu2020', 'Lnc-IL7R', 17, 8, 2, 19)
    x26 = c('Chen2021', 'NEAT1', 103, 15, 19, 163)
    x27 = c('Bamodu2022', 'Lnc-IL7R', 88, 3, 19, 15)
    y = rbind(x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27) %>% as.data.frame() 
    colnames(y) = y[1, ]
    y = y[-1,]
    y$group = c(rep('COPD vs NCs', 12), rep('AECOPD vs NCs', 6), rep('AECOPD vs COPD', 9))
    y$TP = as.numeric(y$TP)
    y$TN = as.numeric(y$TN)
    y$FP = as.numeric(y$FP)
    y$FN = as.numeric(y$FN)
    y$specimen = c('Plasma', 'Plasma', 'Serum', 'Lung tissue', 'PBMCs', 'Blood', 'Serum', 'Blood', 'Serum', 'Blood', 'Lung tissue', 'Blood', 'Plasma', 'Plasma', 'PBMCs', 'PBMCs', 'Plasma', 'PBMCs', 'Plasma', 'Plasma', 'PBMCs', 'PBMCs', 'Plasma', 'PBMCs', 'Plasma', 'Plasma', 'Plasma')
    y$regulated = c('D','U','U','D','U','U', 'D', 'U', 'U','U', 'U', 'D', 'D', 'U', 'U','U', 'U', 'D' ,'D', 'U', 'U','U', 'U', 'D', 'D','U', 'D')
    y1 = y[1:18,]
    y2 = y[19:27,]
    sensitivity_logit1 = metaprop(y1$TP, y1$TP+y1$FN, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = y1$study, byvar = y1$group)
    specificity_logit1 = metaprop(y1$TN, y1$TN+y1$FP, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = y1$study, byvar = y1$group)
    metareg(sensitivity_logit1, y1$specimen, digits= 2)
    metareg(sensitivity_logit1, y1$regulated, digits= 2)
    PLR_logit = metaprop(y1$TP, y1$TP+y1$FP, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = y1$study, byvar = y1$group)
    NLR_logit = metaprop(y1$TN, y1$TN+y1$FN, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = y1$study, byvar = y1$group)
    pdf("Graph_NGPGA.pdf", width = 12, height = 10)
    Graph_NGPGA = meta::forest(sensitivity_logit1, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Sensitivity', layout = "JAMA")
    dev.off()
    pdf("Graph_NGPGB.pdf", width = 12, height = 10)
    Graph_NGPGB = meta::forest(specificity_logit1, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Specificity', layout = "JAMA")
    dev.off()
    DOR_modely1 = metabin(event.e = y1$TP, n.e = y1$TP+y1$FP, event.c = y1$FN, n.c = y1$FN+y1$TN, sm= 'OR', comb.fixed = F, comb.random = T, method = 'Inverse', byvar = y1$group, prediction = F)
    pdf("Graph_NGPGE.pdf", width = 12, height = 10)
    Graph_NGPGE = meta::forest(DOR_modely1, digits= 2, rightcols= c('effect.ci'), studlab = y1$study, xlab= 'Diagnostic Odds Ratio', layout = 'JAMA')
    dev.off()
    sensitivity_logit2 = metaprop(y2$TP, y2$TP+y2$FN, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = y2$study, byvar = y2$group)
    metareg(sensitivity_logit2, y2$specimen, digits= 2)
    metareg(sensitivity_logit2, y2$regulated, digits= 2)
    PLR_logit = metaprop(y2$TP, y2$TP+y2$FP, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = y2$study, byvar = y2$group)
    NLR_logit = metaprop(y2$TN, y2$TN+y2$FN, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = y2$study, byvar = y2$group)
    specificity_logit2 = metaprop(y2$TN, y2$TN+y2$FP, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = y2$study, byvar = y2$group)
    pdf("Graph_NGPGC.pdf", width = 12, height = 10)
    Graph_NGPGC = meta::forest(sensitivity_logit2, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Sensitivity', layout = "JAMA")
    dev.off()
    pdf("Graph_NGPGD.pdf", width = 12, height = 10)
    Graph_NGPGD = meta::forest(specificity_logit2, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Specificity', layout = "JAMA")
    dev.off()
    DOR_modely2 = metabin(event.e = y2$TP, n.e = y2$TP+y2$FP, event.c = y2$FN, n.c = y2$FN+y2$TN, sm= 'OR', comb.fixed = F, comb.random = T, method = 'Inverse', byvar = y2$group, prediction = F, studlab = y2$study) 
    pdf("Graph_NGPGF.pdf", width = 12, height = 10)
    Graph_NGPGF = meta::forest(DOR_modely2, digits= 2, rightcols= c('effect.ci'), xlab= 'Diagnostic Odds Ratio', layout = "JAMA")
    dev.off()
#### SROC curve
    pdf("Graph_NGPG1.pdf", width = 12, height = 10)
    fit.y1 <- reitsma(y1)
    plot(fit.y1, main = 'SROC curve', cex.main = 2, cex.lab = 1.6)
    ROCellipse(y1, add = TRUE, pch = 1, corr = 0, suppress = FALSE, xlim = c(0, 1), ylim = c(0, 1), ellipsecol = "red")
    dev.off()
    pdf("Graph_NGPG2.pdf", width = 12, height = 10)
    fit.y2 <- reitsma(y2)
    plot(fit.y2, main = 'SROC curve', cex.main = 2, cex.lab = 1.6)
    ROCellipse(y2, add = TRUE, pch = 1, corr = 0, suppress = FALSE, xlim = c(0, 1), ylim = c(0, 1), ellipsecol = "red")
    dev.off()
#### Sensitivity analysis
    DOR_y1 = metainf(DOR_modely1)
    pdf("Graph_NGPGG1.pdf", width = 12, height = 10)
    Graph_NGPGG1 = meta::forest(DOR_y1, digits= 2, xlab= 'Diagnostic Odds Ratio', studlab = paste('Omitting', y1$study), layout = 'JAMA', xlim = c(0.5, 30))
    dev.off()
    sensitivity_logit1 = metaprop(y1$TP, y1$TP+y1$FN, comb.fixed = FALSE, comb.random = TRUE, sm = 'PRAW')
    sen_y1 = metainf(sensitivity_logit1)
    pdf("Graph_NGPGG2.pdf", width = 12, height = 10)
    Graph_NGPGG2 = meta::forest(sen_y1, digits= 2, xlab= 'Sensitivity', studlab = paste('Omitting', y1$study), layout = 'JAMA', xlim = c(0.7, 1))
    dev.off()
    specificity_logit1 = metaprop(y1$TN, y1$TN+y1$FP, comb.fixed = FALSE, comb.random = TRUE, sm = 'PRAW')
    spe_y1 = metainf(specificity_logit1)
    pdf("Graph_NGPGG3.pdf", width = 12, height = 10)
    Graph_NGPGG3 = meta::forest(spe_y1, digits= 2, xlab= 'Specificity', studlab = paste('Omitting', y1$study), layout = 'JAMA', xlim = c(0.5, 1))
    dev.off()
    DOR_y2 = metainf(DOR_modely2)
    pdf("Graph_NGPGH1.pdf", width = 12, height = 10)
    Graph_NGPGH1 = meta::forest(DOR_y2, digits= 2, xlab= 'Diagnostic Odds Ratio', studlab = paste('Omitting', y2$study), layout = 'JAMA', xlim = c(0.5, 30))
    dev.off()
    sensitivity_logit1 = metaprop(y2$TP, y2$TP+y2$FN, comb.fixed = FALSE, comb.random = TRUE, sm = 'PRAW')
    sen_y2 = metainf(sensitivity_logit1)
    pdf("Graph_NGPGH2.pdf", width = 12, height = 10)
    Graph_NGPGH2 = meta::forest(sen_y2, digits= 2, xlab= 'Sensitivity', studlab = paste('Omitting', y2$study), layout = 'JAMA', xlim = c(0.5, 1))
    dev.off()
    specificity_logit1 = metaprop(y2$TN, y2$TN+y2$FP, comb.fixed = FALSE, comb.random = TRUE, sm = 'PRAW')
    spe_y2 = metainf(specificity_logit1)
    pdf("Graph_NGPGH3.pdf", width = 12, height = 10)
    Graph_NGPGH3 = meta::forest(spe_y2, digits= 2, xlab= 'Specificity', studlab = paste('Omitting', y2$study), layout = 'JAMA', xlim = c(0.5, 1))
    dev.off()
#### Funnel plot
    pdf("Graph_NGPGI1.pdf", width = 12, height = 10)
    trf<-trimfill(metabin(event.e =y1$TP, n.e =y1$TP+y1$FP, event.c =y1$FN, n.c =y1$FN+y1$TN, sm= 'OR', comb.fixed = F, comb.random = T, method = 'Inverse', prediction = F, studlab =y1$study))
    summary(trf)
    funnel(trf, xlab = "Risk Ratio",
    contour = c(.90,.95,.99), cex.main = 2, cex.lab = 1.6,
    col.contour = c("grey25","grey50","grey75"), bg = "white")
    dev.off()
    metabias(metabin(event.e =y1$TP, n.e =y1$TP+y1$FP, event.c =y1$FN, n.c =y1$FN+y1$TN, sm= 'OR', comb.fixed = F, comb.random = T, method = 'Inverse', prediction = F, studlab =y1$study), method.bias = "linreg")
    pdf("Graph_NGPGI2.pdf", width = 12, height = 10)
    trf<-trimfill(metabin(event.e =y2$TP, n.e =y2$TP+y2$FP, event.c =y2$FN, n.c =y2$FN+y2$TN, sm= 'OR', comb.fixed = F, comb.random = T, method = 'Inverse', prediction = F, studlab =y2$study))
    summary(trf)
    funnel(trf, xlab = "Risk Ratio",
    contour = c(.90,.95,.99), 
    col.contour = c("grey25","grey50","grey75"), bg = "white")
    dev.off()
    metabias(metabin(event.e =y2$TP, n.e =y2$TP+y2$FP, event.c =y2$FN, n.c =y2$FN+y2$TN, sm= 'OR', comb.fixed = F, comb.random = T, method = 'Inverse', prediction = F, studlab =y2$study), method.bias = "linreg")
#### QUADAS-2 checklist
    library(readxl)
    a = read_excel('~/Documents/Paper_LncRNAMeta/QUADAS2.xlsx')
    b = a %>% 
        select(c('Studlab', 'Patient selection', 'Index test', 'Reference standard', 'Flow and Timing')) %>%
        pivot_longer(-Studlab, names_to = 'group', values_to = 'value') %>%
        count(group, value) %>%
        mutate(percent = n/14) %>%
        as.data.frame()
    b$value = as.factor(b$value)
    b$value = factor(b$value, levels = c('low', 'unclear', 'high'))
    b$group = as.factor(b$group)
    b$group = factor(b$group, levels = c('Flow and Timing', 'Reference standard', 'Index test', 'Patient selection'))
    str(b)
    x = ggplot(b,aes(x = group, y = percent, fill = value))+
        geom_bar(stat="identity")+
        scale_y_continuous(labels = scales::percent_format(scale = 100))+ 
        coord_flip()+
        theme_bw()+
        guides(fill=guide_legend(title=NULL)) +
        labs(title = "Risk of Bias", x ="", y = "")+
        theme(axis.text.x=element_text(vjust=1,hjust=1,size=12), axis.text.y=element_text(vjust=1,hjust=1,size=16), plot.title = element_text(face = 'bold', size = 20))
    x + scale_fill_manual(values = c("#0cef4cee","#e7f313ee", "#e75858ee" ))
    ggsave('1.png', width = 10, height = 4)
    b = a %>% 
        select(c('Studlab', 'Patient selection-a', 'Index test-a', 'Reference standard-a')) %>%
        pivot_longer(-Studlab, names_to = 'group', values_to = 'value') %>%
        count(group, value) %>%
        mutate(percent = n/14) %>%
        as.data.frame()
    b$value = as.factor(b$value)
    b$value = factor(b$value, levels = c('low', 'unclear', 'high'))
    b$group = as.factor(b$group)
    b$group = factor(b$group, levels = c('Reference standard-a', 'Index test-a', 'Patient selection-a'))
    str(b)
    x = ggplot(b,aes(x = group, y = percent, fill = value))+
        geom_bar(stat="identity")+
        scale_y_continuous(labels = scales::percent_format(scale = 100))+ 
        scale_x_discrete(labels=c('Reference standard', 'Index test', 'Patient selection'))+
        coord_flip()+
        theme_bw()+
        guides(fill=guide_legend(title=NULL)) +
        labs(title = "Applicability Concerns", x ="", y = "")+
        theme(axis.text.x=element_text(vjust=1,hjust=1,size=12), axis.text.y=element_text(vjust=1,hjust=1,size=16), plot.title = element_text(face = 'bold', size = 20))
    x + scale_fill_manual(values = c("#0cef4cee","#e7f313ee", "#e75858ee" ))
    ggsave('2.png', width = 10, height = 4)
#### Subgroup analysis: specimen
    y11 = subset(y1, group == 'COPD vs NCs')
    sensitivity_logit1 = metaprop(TP, TP+FN, data = y11, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    specificity_logit1 = metaprop(TN, TN+FP, data = y11,comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    metareg(sensitivity_logit1, specimen, digits= 2)
    metareg(specificity_logit1, specimen, digits= 2)
    metareg(sensitivity_logit1, regulated, digits= 2)
    metareg(specificity_logit1, regulated, digits= 2)
    PLR_logit = metaprop(TP, TP+FP, data = y11, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    NLR_logit = metaprop(TN, TN+FN, data = y11, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    pdf("Graph_NGPGA1.pdf", width = 12, height = 10)
    Graph_NGPGA1 = meta::forest(sensitivity_logit1, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Sensitivity', layout = "JAMA")
    dev.off()
    pdf("Graph_NGPGB1.pdf", width = 12, height = 10)
    Graph_NGPGB1 = meta::forest(specificity_logit1, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Specificity', layout = "JAMA")
    dev.off()
    DOR_modely1 = metabin(event.e = TP, n.e = TP+FP, event.c = FN, n.c = FN+TN, sm= 'OR', data = y11, comb.fixed = F, comb.random = T, method = 'Inverse', byvar = specimen, prediction = F)
    metareg(DOR_modely1, specimen, digits= 2)
    DOR_modely1 = metabin(event.e = TP, n.e = TP+FP, event.c = FN, n.c = FN+TN, sm= 'OR', data = y11, comb.fixed = F, comb.random = T, method = 'Inverse', byvar = regulated, prediction = F)
    metareg(DOR_modely1, regulated, digits= 2)
    pdf("Graph_NGPGE1.pdf", width = 12, height = 10)
    Graph_NGPGE1 = meta::forest(DOR_modely1, digits= 2, rightcols= c('effect.ci'), studlab = study, xlab= 'Diagnostic Odds Ratio', layout = 'JAMA')
    dev.off()
    y12 = subset(y1, group == 'AECOPD vs NCs')
    sensitivity_logit12 = metaprop(TP, TP+FN, data = y12, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    specificity_logit12 = metaprop(TN, TN+FP, data = y12,comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    metareg(sensitivity_logit12, specimen, digits= 2)
    metareg(sensitivity_logit12, regulated, digits= 2)
    metareg(specificity_logit12, specimen, digits= 2)
    metareg(specificity_logit12, regulated, digits= 2)
    PLR_logit = metaprop(TP, TP+FP, data = y12, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    NLR_logit = metaprop(TN, TN+FN, data = y12, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    pdf("Graph_NGPGA11.pdf", width = 12, height = 10)
    Graph_NGPGA11 = meta::forest(sensitivity_logit12, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Sensitivity', layout = "JAMA")
    dev.off()
    pdf("Graph_NGPGB11.pdf", width = 12, height = 10)
    Graph_NGPGB11 = meta::forest(specificity_logit12, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Specificity', layout = "JAMA")
    dev.off()
    DOR_modely12 = metabin(event.e = TP, n.e = TP+FP, event.c = FN, n.c = FN+TN, sm= 'OR', data = y12, comb.fixed = F, comb.random = T, method = 'Inverse', byvar = specimen, prediction = F)
    metareg(DOR_modely12, specimen, digits= 2)
    metareg(DOR_modely12, regulated, digits= 2)
    pdf("Graph_NGPGE11.pdf", width = 12, height = 10)
    Graph_NGPGE11 = meta::forest(DOR_modely12, digits= 2, rightcols= c('effect.ci'), studlab = study, xlab= 'Diagnostic Odds Ratio', layout = 'JAMA')
    dev.off()
    sensitivity_logit2 = metaprop(TP, TP+FN, data = y2, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    specificity_logit2 = metaprop(TN, TN+FP, data = y2, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    metareg(sensitivity_logit2, specimen, digits= 2)
    metareg(sensitivity_logit2, regulated, digits= 2)
    metareg(specificity_logit2, specimen, digits= 2)
    metareg(specificity_logit2, regulated, digits= 2)
    PLR_logit = metaprop(TP, TP+FP, data = y2,  comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    NLR_logit = metaprop(TN, TN+FN, data = y2,  comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    pdf("Graph_NGPGC1.pdf", width = 12, height = 10)
    Graph_NGPGC1 = meta::forest(sensitivity_logit2, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Sensitivity', layout = "JAMA")
    dev.off()
    pdf("Graph_NGPGD1.pdf", width = 12, height = 10)
    Graph_NGPGD1 = meta::forest(specificity_logit2, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Specificity', layout = "JAMA")
    dev.off()
    DOR_modely2 = metabin(event.e = TP, n.e = TP+FP, event.c = FN, n.c = FN+TN, data = y2, sm= 'OR', comb.fixed = F, comb.random = T, method = 'Inverse', byvar = specimen, prediction = F, studlab = study) 
    metareg(DOR_modely2, specimen, digits= 2)
    metareg(DOR_modely2, regulated, digits= 2)
    metareg(sensitivity_logit1, specimen + regulated, digits= 2)
    metareg(specificity_logit1, specimen + regulated, digits= 2)
    metareg(DOR_modely1, y1$specimen + y1$regulated, digits= 2)
    metareg(sensitivity_logit2, specimen + regulated, digits= 2)
    metareg(specificity_logit2, specimen + regulated, digits= 2)
    metareg(DOR_modely2, y2$specimen + y2$regulated, digits= 2)
    pdf("Graph_NGPGF1.pdf", width = 12, height = 10)
    Graph_NGPGF1 = meta::forest(DOR_modely2, digits= 2, rightcols= c('effect.ci'), xlab= 'Diagnostic Odds Ratio', layout = "JAMA")
    dev.off()
#### Regression analysis.
    metareg(sensitivity_logit2, y2$specimen, digits= 2)
    metareg(sensitivity_logit1, y1$specimen, digits= 2)
    metareg(DOR_modely1, y1$specimen, digits= 2)
    metareg(sensitivity_logit1, y1$regulated, digits= 2)

#### Subgroup analysis: regulated
    y11 = subset(y1, group == 'COPD vs NCs')
    sensitivity_logit1 = metaprop(TP, TP+FN, data = y11, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = regulated)
    specificity_logit1 = metaprop(TN, TN+FP, data = y11,comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = regulated)
    metareg(sensitivity_logit1, regulated, digits= 2)
    metareg(sensitivity_logit1, regulated, digits= 2)
    PLR_logit = metaprop(TP, TP+FP, data = y11, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = regulated)
    NLR_logit = metaprop(TN, TN+FN, data = y11, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = regulated)
    pdf("Graph_NGPGA2.pdf", width = 12, height = 10)
    Graph_NGPGA2 = meta::forest(sensitivity_logit1, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Sensitivity', layout = "JAMA")
    dev.off()
    pdf("Graph_NGPGB2.pdf", width = 12, height = 10)
    Graph_NGPGB2 = meta::forest(specificity_logit1, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Specificity', layout = "JAMA")
    dev.off()
    DOR_modely1 = metabin(event.e = TP, n.e = TP+FP, event.c = FN, n.c = FN+TN, sm= 'OR', data = y11, comb.fixed = F, comb.random = T, method = 'Inverse', byvar = regulated, prediction = F)
    pdf("Graph_NGPGE2.pdf", width = 12, height = 10)
    Graph_NGPGE2 = meta::forest(DOR_modely1, digits= 2, rightcols= c('effect.ci'), studlab = study, xlab= 'Diagnostic Odds Ratio', layout = 'JAMA')
    dev.off()
    y12 = subset(y1, group == 'AECOPD vs NCs')
    sensitivity_logit12 = metaprop(TP, TP+FN, data = y12, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = regulated)
    specificity_logit12 = metaprop(TN, TN+FP, data = y12,comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = regulated)
    metareg(sensitivity_logit12, regulated, digits= 2)
    metareg(sensitivity_logit12, regulated, digits= 2)
    PLR_logit = metaprop(TP, TP+FP, data = y12, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = regulated)
    NLR_logit = metaprop(TN, TN+FN, data = y12, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = regulated)
    pdf("Graph_NGPGA21.pdf", width = 12, height = 10)
    Graph_NGPGA21 = meta::forest(sensitivity_logit12, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Sensitivity', layout = "JAMA")
    dev.off()
    pdf("Graph_NGPGB21.pdf", width = 12, height = 10)
    Graph_NGPGB11 = meta::forest(specificity_logit12, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Specificity', layout = "JAMA")
    dev.off()
    DOR_modely12 = metabin(event.e = TP, n.e = TP+FP, event.c = FN, n.c = FN+TN, sm= 'OR', data = y12, comb.fixed = F, comb.random = T, method = 'Inverse', byvar = regulated, prediction = F)
    pdf("Graph_NGPGE21.pdf", width = 12, height = 10)
    Graph_NGPGE21 = meta::forest(DOR_modely12, digits= 2, rightcols= c('effect.ci'), studlab = study, xlab= 'Diagnostic Odds Ratio', layout = 'JAMA')
    dev.off()
    sensitivity_logit2 = metaprop(TP, TP+FN, data = y2, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = regulated)
    metareg(sensitivity_logit2, regulated, digits= 2)
    metareg(sensitivity_logit2, regulated, digits= 2)
    PLR_logit = metaprop(TP, TP+FP, data = y2,  comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = regulated)
    NLR_logit = metaprop(TN, TN+FN, data = y2,  comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = regulated)
    specificity_logit2 = metaprop(TN, TN+FP, data = y2, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = regulated)
    pdf("Graph_NGPGC12.pdf", width = 12, height = 10)
    Graph_NGPGC12 = meta::forest(sensitivity_logit2, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Sensitivity', layout = "JAMA")
    dev.off()
    pdf("Graph_NGPGD12.pdf", width = 12, height = 10)
    Graph_NGPGD12 = meta::forest(specificity_logit2, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Specificity', layout = "JAMA")
    dev.off()
    DOR_modely2 = metabin(event.e = TP, n.e = TP+FP, event.c = FN, n.c = FN+TN, data = y2, sm= 'OR', comb.fixed = F, comb.random = T, method = 'Inverse', byvar = regulated, prediction = F, studlab = study) 
    pdf("Graph_NGPGF12.pdf", width = 12, height = 10)
    Graph_NGPGF12 = meta::forest(DOR_modely2, digits= 2, rightcols= c('effect.ci'), xlab= 'Diagnostic Odds Ratio', layout = "JAMA")
    dev.off()
#### multi-lncRNA analysis
    y11_pvt1 = subset(y1, group == 'COPD vs NCs' &LncRNA == 'PVT1')
    sensitivity_logit1 = metaprop(TP, TP+FN, data = y11_pvt1, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    specificity_logit1 = metaprop(TN, TN+FP, data = y11_pvt1,comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    PLR_logit = metaprop(TP, TP+FP, data = y11_pvt1, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    NLR_logit = metaprop(TN, TN+FN, data = y11_pvt1, comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    pdf("Graph_NGPGApvt.pdf", width = 12, height = 10)
    Graph_NGPGApvt = meta::forest(sensitivity_logit1, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Sensitivity', layout = "JAMA")
    dev.off()
    pdf("Graph_NGPGBpvt.pdf", width = 12, height = 10)
    Graph_NGPGBpvt = meta::forest(specificity_logit1, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Specificity', layout = "JAMA")
    dev.off()
    DOR_modely1 = metabin(event.e = TP, n.e = TP+FP, event.c = FN, n.c = FN+TN, sm= 'OR', data = y11_pvt1, comb.fixed = F, comb.random = T, method = 'Inverse', byvar = specimen, prediction = F)
    pdf("Graph_NGPGEpvt.pdf", width = 12, height = 10)
    Graph_NGPGEpvt = meta::forest(DOR_modely1, digits= 2, rightcols= c('effect.ci'), studlab = study, xlab= 'Diagnostic Odds Ratio', layout = 'JAMA')
    dev.off()
    sensitivity_logit2 = metaprop(TP, TP+FN, data = subset(y2, LncRNA == 'NEAT1'), comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    PLR_logit = metaprop(TP, TP+FP, data = subset(y2, LncRNA == 'NEAT1'),  comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    NLR_logit = metaprop(TN, TN+FN, data = subset(y2, LncRNA == 'NEAT1'),  comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    specificity_logit2 = metaprop(TN, TN+FP, data = subset(y2, LncRNA == 'NEAT1'), comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    pdf("Graph_NGPGC1ne.pdf", width = 12, height = 10)
    Graph_NGPGC1ne = meta::forest(sensitivity_logit2, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Sensitivity', layout = "JAMA")
    dev.off()
    pdf("Graph_NGPGD1ne.pdf", width = 12, height = 10)
    Graph_NGPGD1ne = meta::forest(specificity_logit2, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Specificity', layout = "JAMA")
    dev.off()
    DOR_modely2 = metabin(event.e = TP, n.e = TP+FP, event.c = FN, n.c = FN+TN, data = subset(y2, LncRNA == 'NEAT1'), sm= 'OR', comb.fixed = F, comb.random = T, method = 'Inverse', byvar = specimen, prediction = F, studlab = study) 
    pdf("Graph_NGPGF1ne.pdf", width = 12, height = 10)
    Graph_NGPGF1ne = meta::forest(DOR_modely2, digits= 2, rightcols= c('effect.ci'), xlab= 'Diagnostic Odds Ratio', layout = "JAMA")
    dev.off()
    sensitivity_logit2 = metaprop(TP, TP+FN, data = subset(y2, LncRNA == 'Lnc-IL7R'), comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    PLR_logit = metaprop(TP, TP+FP, data = subset(y2, LncRNA == 'Lnc-IL7R'),  comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    NLR_logit = metaprop(TN, TN+FN, data = subset(y2, LncRNA == 'Lnc-IL7R'),  comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    specificity_logit2 = metaprop(TN, TN+FP, data = subset(y2, LncRNA == 'Lnc-IL7R'), comb.fixed = FALSE, comb.random = TRUE, sm = 'PLOGIT', method.ci = 'CP', studlab = study, byvar = specimen)
    pdf("Graph_NGPGC117.pdf", width = 12, height = 10)
    Graph_NGPGC117 = meta::forest(sensitivity_logit2, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Sensitivity', layout = "JAMA")
    dev.off()
    pdf("Graph_NGPGD117.pdf", width = 12, height = 10)
    Graph_NGPGD117 = meta::forest(specificity_logit2, digits= 2, rightcols= c('effect', 'ci'), xlab= 'Specificity', layout = "JAMA")
    dev.off()
    DOR_modely2 = metabin(event.e = TP, n.e = TP+FP, event.c = FN, n.c = FN+TN, data = subset(y2, LncRNA == 'Lnc-IL7R'), sm= 'OR', comb.fixed = F, comb.random = T, method = 'Inverse', byvar = specimen, prediction = F, studlab = study) 
    pdf("Graph_NGPGF117.pdf", width = 12, height = 10)
    Graph_NGPGF117 = meta::forest(DOR_modely2, digits= 2, rightcols= c('effect.ci'), xlab= 'Diagnostic Odds Ratio', layout = "JAMA")
    dev.off()
#### Bioinformatics analysis
    require(tidyverse)
    Subcellular = read.table('Subcellular_Location-2.txt', header = T, sep = '\t')
    head(Subcellular)
    Subcellular$log10 = log(Subcellular$P.value)
    plot = Subcellular %>% 
        ggplot(aes(x = log10, y = Set, color = Sub_Class, size = Count)) +
        geom_point(alpha=0.7) +
        labs(title = "Subcellular Location", y = "", x = "log10(P.value)") + 
        theme_classic() +
        scale_size(range = c(.1, 12)) +
        theme(
        axis.text = element_text(size = 16, vjust = 0.5),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 20))
    ggsave('subcellular.png', width = 10, height = 6)
    Hallmark = read.table('Cancer_Hallmark-2.txt', header = T, sep = '\t')
    head(Hallmark)
    Hallmark$log10 = log(Hallmark$P.value)
    plot = Hallmark %>% 
        ggplot(aes(x = log10, y = Set, size = Count)) +
        geom_point(alpha=0.7) +
        labs(title = "Hallmark", y = "", x = "log10(P.value)") + 
        theme_classic() +
        scale_size(range = c(.1, 12)) +
        theme(
        axis.text = element_text(size = 16, vjust = 0.5),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 20))
    ggsave('Hallmark.png', width = 10, height = 6)
    Binding = read.table('RNA_Binding_Protein-2.txt', header = T, sep = '\t')
    head(Binding)
    Binding$log10 = log(Binding$P.value)
    duplicated(Binding$Set)
    Binding1 = Binding %>% 
        distinct(Set, .keep_all = T) ## 1052 
    str(Binding1)
    library(ggplot2)
    kegg = read.table("KeGGPathway.txt", header=T, sep="\t")  #导入文件
    plot = kegg %>% 
        ggplot(aes(x=Fold.Enrichment, y=Term)) + 
        geom_point(aes(size=Count,color=-1*log10(PValue))) +  #点的大小根据Count数变化，颜色根据P值变化
        scale_colour_gradient(low="blue", high="red") +  #设置图例
        labs(color=expression(-log[10](P.value)),
            title = "KEGG pathway analysis",
            size="Gene number", y ='',
            x="Fold enrichment") +
        scale_size(range = c(.1, 12)) +
        theme_classic() +  #设置背景
        theme(
        axis.text = element_text(size = 16, vjust = 0.5),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 20))
    ggsave(filename="kegg.png", height=8, width=12, dpi=600)  #保存图片
    GO_BP = read.table("GOTERM_BP_DIRECT.txt", header=T, sep="\t")
    GO_BP = GO_BP[order(-GO_BP$Count), ]
    GO_BP10 = GO_BP[1:10,]
    head(GO_BP10)
    GO_CC = read.table("GOTERM_CC_DIRECT.txt", header=T, sep="\t", quote = "")
    GO_CC = GO_CC[order(-GO_CC$Count), ]
    GO_CC10 = GO_CC[1:10,]
    head(GO_CC10)
    GO_MF = read.table("GOTERM_MF_DIRECT.txt", header=T, sep="\t", quote = "")
    GO_MF = GO_MF[order(-GO_MF$Count), ]
    GO_MF10 = GO_MF[1:10,]
    head(GO_MF10)
    GO = rbind(GO_BP10, GO_CC10, GO_MF10) %>% 
        as.data.frame()
    head(GO,20)
    GO$Term = factor(GO$Term, levels = GO$Term)
    GO %>% 
        mutate(group = c(rep('BP',10), rep('CC',10), rep('MF',10))) %>% 
        ggplot(aes(x = factor(Term) , y = Fold.Enrichment, fill = group)) +
        geom_bar(stat="identity", width=0.8) + coord_flip() +
        labs(x ='', title = "The Top 10 Enriched GO Terms", 
            y ="Fold enrichment") + 
        theme_classic() +
        theme(
        axis.text = element_text(size = 12, vjust = 0.5),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 20))
    ggsave('GO.png', width = 12, height = 8)
