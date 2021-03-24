################# Comparison of S. aureus, 2019-nCoV & H1N1  ###################
###############################################################################
# Before executing this file load the 'data_analysis_functions.R',            #
# 'staphAureus_data_analysis.R', '2019_nCoV_data_analysis.R and               #
# 'inflA_H1N1_data_analysis.R' in this Environment.                           #
###############################################################################


# Plot UniProt overview
comp_UP_overview <- matrix(c(c(nrow(usa300_uniProt),
                               nrow(mu50_uniProt),
                               nrow(cov_uniProt),
                               nrow(pur_h1n1_uniProt),
                               nrow(scb_h1n1_uniProt),
                               nrow(stp_h1n1_uniProt),
                               nrow(tex_h1n1_uniProt)),
                             round(c(mean(usa300_len),
                                     mean(mu50_len),
                                     mean(cov_len),
                                     mean(pur_h1n1_len),
                                     mean(scb_h1n1_len),
                                     mean(stp_h1n1_len),
                                     mean(tex_h1n1_len))),
                             c(length(na.omit(usa300_sig)),
                               length(na.omit(mu50_sig)),
                               length(na.omit(cov_sig)),
                               length(na.omit(pur_h1n1_sig)),
                               length(na.omit(scb_h1n1_sig)),
                               length(na.omit(stp_h1n1_sig)),
                               length(na.omit(tex_h1n1_sig))),
                             round(c(mean(na.omit(usa300_sig)),
                                     mean(na.omit(mu50_sig)),
                                     mean(na.omit(cov_sig)),
                                     mean(na.omit(pur_h1n1_sig)),
                                     mean(na.omit(scb_h1n1_sig)),
                                     mean(na.omit(stp_h1n1_sig)),
                                     mean(na.omit(tex_h1n1_sig))))),
                          nrow = 7, ncol = 4)
colnames(comp_UP_overview) <- c("Anz. Sequenzen", "Länge [AS]",
                                "Anz. Signalpeptide","Signallänge [AS]")
rownames(comp_UP_overview) <- c("S. aureus, USA300", "S. aureus, Mu50", "SARS-CoV-2",
                                "H1N1, Puerto Rico / 1934", "H1N1, S. Canterbury / 2000",
                                "H1N1, St.Petersburg / 2006", "H1N1, Texas / 2007")
comp_up_table <- tableGrob(comp_UP_overview, theme = table_theme)
comp_up_title <- textGrob("UniProt Übersicht")
comp_up_subtt <- textGrob(expression(paste("Stämme von ", italic("S. aureus, SARS-Coronaviren, Influenzaviren A"))))
grid.arrange(comp_up_title, comp_up_subtt, comp_up_table, ncol=1,
             heights = unit(c(10,10,65),rep("mm",3)))

# Boxplot of all lengths
layout(mat = matrix(c(1, 2, 1, 2, 3, 3), 
                    nrow = 3, 
                    ncol = 2,
                    byrow=TRUE))
par(mar=c(5.1, 4.1, 0.1, 1.1), tck=0)

boxplot(usa300_len, mu50_len, col = blues9[2], axis=FALSE, 
        names=expression(italic("S. aureus, USA300"), italic("S. aureus, Mu50")))
boxplot(cov_len, col = blues9[2], show.names=TRUE, names=expression(italic("SARS-CoV-2")))
boxplot(pur_h1n1_len, scb_h1n1_len, stp_h1n1_len, tex_h1n1_len, col = blues9[2], 
        axis=FALSE, names = c(expression(italic("H1N1, Puerto Rico / 1934"),
                                         italic("H1N1, S.Canterbury / 2000"),
                                         italic("H1N1, St.Petersburg / 2006"),
                                         italic("H1N1, Texas / 2007"))))
par(mar=c(5.1, 4.1, 4.1, 1.1), tck=NA)


##### Epitopes #####
# Plot prediction overview
comp_Pred_overview_nmp <- matrix(c(round(c(mean(usa300_nmp_cnt),
                                           mean(mu50_nmp_cnt),
                                           mean(cov_nmp_cnt),
                                           mean(pur_h1n1_nmp_cnt),
                                           mean(scb_h1n1_nmp_cnt),
                                           mean(stp_h1n1_nmp_cnt),
                                           mean(tex_h1n1_nmp_cnt))),
                                   round(c(mean(as.numeric(unlist(usa300_nmp_dns))),
                                           mean(na.omit(as.numeric(unlist(mu50_nmp_dns)))),
                                           mean(as.numeric(unlist(cov_nmp_dns))),
                                           mean(as.numeric(unlist(pur_h1n1_nmp_dns))),
                                           mean(as.numeric(unlist(scb_h1n1_nmp_dns))),
                                           mean(as.numeric(unlist(stp_h1n1_nmp_dns))),
                                           mean(as.numeric(unlist(tex_h1n1_nmp_dns)))),
                                         digits = 3),
                                   round(c(mean(na.omit(as.numeric(unlist(usa300_nmp_sgd))))*(mean(usa300_sig)-9),
                                           mean(na.omit(as.numeric(unlist(mu50_nmp_sgd))))*(mean(mu50_sig)-9),
                                           mean(na.omit(as.numeric(unlist(cov_nmp_sgd))))*(mean(cov_sig)-9),
                                           mean(na.omit(as.numeric(unlist(pur_h1n1_nmp_sgd))))*(mean(pur_h1n1_sig)-9),
                                           mean(na.omit(as.numeric(unlist(scb_h1n1_nmp_sgd))))*(mean(scb_h1n1_sig)-9),
                                           mean(na.omit(as.numeric(unlist(stp_h1n1_nmp_sgd))))*(mean(stp_h1n1_sig)-9),
                                           mean(na.omit(as.numeric(unlist(tex_h1n1_nmp_sgd))))*(mean(tex_h1n1_sig)-9))),
                                   round(c(mean(na.omit(as.numeric(unlist(usa300_nmp_sgd)))),
                                           mean(na.omit(as.numeric(unlist(mu50_nmp_sgd)))),
                                           mean(na.omit(as.numeric(unlist(cov_nmp_sgd)))),
                                           mean(na.omit(as.numeric(unlist(pur_h1n1_nmp_sgd)))),
                                           mean(na.omit(as.numeric(unlist(scb_h1n1_nmp_sgd)))),
                                           mean(na.omit(as.numeric(unlist(stp_h1n1_nmp_sgd)))),
                                           mean(na.omit(as.numeric(unlist(tex_h1n1_nmp_sgd))))),
                                         digits = 3)),
                                nrow = 7, ncol = 4)

comp_Pred_overview_spt <- matrix(c(round(c(mean(usa300_spt_cnt),
                                           mean(mu50_spt_cnt),
                                           mean(cov_spt_cnt),
                                           mean(pur_h1n1_spt_cnt),
                                           mean(scb_h1n1_spt_cnt),
                                           mean(stp_h1n1_spt_cnt),
                                           mean(tex_h1n1_spt_cnt))),
                                   round(c(mean(as.numeric(unlist(usa300_spt_dns))),
                                           mean(as.numeric(unlist(mu50_spt_dns))),
                                           mean(as.numeric(unlist(cov_spt_dns))),
                                           mean(as.numeric(unlist(pur_h1n1_spt_dns))),
                                           mean(as.numeric(unlist(scb_h1n1_spt_dns))),
                                           mean(as.numeric(unlist(stp_h1n1_spt_dns))),
                                           mean(as.numeric(unlist(tex_h1n1_spt_dns)))),
                                         digits = 3),
                                   round(c(mean(na.omit(as.numeric(unlist(usa300_spt_sgd))))*(mean(usa300_sig)-9),
                                           mean(na.omit(as.numeric(unlist(mu50_spt_sgd))))*(mean(mu50_sig)-9),
                                           mean(na.omit(as.numeric(unlist(cov_spt_sgd))))*(mean(cov_sig)-9),
                                           mean(na.omit(as.numeric(unlist(pur_h1n1_spt_sgd))))*(mean(pur_h1n1_sig)-9),
                                           mean(na.omit(as.numeric(unlist(scb_h1n1_spt_sgd))))*(mean(scb_h1n1_sig)-9),
                                           mean(na.omit(as.numeric(unlist(stp_h1n1_spt_sgd))))*(mean(stp_h1n1_sig)-9),
                                           mean(na.omit(as.numeric(unlist(tex_h1n1_spt_sgd))))*(mean(tex_h1n1_sig)-9))),
                                   round(c(mean(na.omit(as.numeric(unlist(usa300_spt_sgd)))),
                                           mean(na.omit(as.numeric(unlist(mu50_spt_sgd)))),
                                           mean(na.omit(as.numeric(unlist(cov_spt_sgd)))),
                                           mean(na.omit(as.numeric(unlist(pur_h1n1_spt_sgd)))),
                                           mean(na.omit(as.numeric(unlist(scb_h1n1_spt_sgd)))),
                                           mean(na.omit(as.numeric(unlist(stp_h1n1_spt_sgd)))),
                                           mean(na.omit(as.numeric(unlist(tex_h1n1_spt_sgd))))),
                                         digits = 3),
                                   c(sum(usa300_len>2000),
                                     sum(mu50_len>2000),
                                     sum(cov_len>2000),
                                     sum(pur_h1n1_len>2000),
                                     sum(scb_h1n1_len>2000),
                                     sum(stp_h1n1_len>2000),
                                     sum(tex_h1n1_len>2000))),
                                 nrow = 7, ncol = 5)
colnames(comp_Pred_overview_nmp) <- c("Anz. Epitope", "Epitopdichte", "Anz. Signalepitope", "Signal-Epitopdichte")
colnames(comp_Pred_overview_spt) <- c("Anz. Epitope", "Epitopdichte", "Anz. Signalepitope", "Signal-Epitopdichte", "Gekappte Seq.")
rownames(comp_Pred_overview_nmp) <- c("S. aureus, USA300", "S. aureus, Mu50", "SARS-CoV-2",
                                      "H1N1, Puerto Rico / 1934", "H1N1, S. Canterbury / 2000",
                                      "H1N1, St.Petersburg / 2006", "H1N1, Texas / 2007")
rownames(comp_Pred_overview_spt) <- c("S. aureus, USA300", "S. aureus, Mu50", "SARS-CoV-2",
                                      "H1N1, Puerto Rico / 1934", "H1N1, S. Canterbury / 2000",
                                      "H1N1, St.Petersburg / 2006", "H1N1, Texas / 2007")
comp_nmp_table <- tableGrob(comp_Pred_overview_nmp, theme = table_theme)
comp_spt_table <- tableGrob(comp_Pred_overview_spt, theme = table_theme)
comp_pred_title <- textGrob("Vorhersagen Übersicht")
comp_nmp_subtt <- textGrob(expression("Oben: NetMHCpan, unten: SYFPEITHI"))
grid.arrange(comp_pred_title, comp_nmp_subtt, comp_nmp_table,
             comp_spt_table, ncol=1, heights = unit(c(10,10,60,60), rep("mm", 4)))

# Boxplot about densities
layout(mat = matrix(c(1, 2, 1, 2, 3, 3), 
                    nrow = 3, 
                    ncol = 2,
                    byrow=TRUE))
par(mar=c(5.1, 4.1, 0.1, 1.1), tck=0)

boxplot(usa300_nmp_dns, mu50_nmp_dns, col = blues9[2], axis=FALSE, 
        names=expression(italic("S. aureus, USA300"), italic("S. aureus, Mu50")))
boxplot(cov_nmp_dns, col = blues9[2], show.names=TRUE, names=expression(italic("SARS-CoV-2")))
boxplot(pur_h1n1_nmp_dns, scb_h1n1_nmp_dns, stp_h1n1_nmp_dns, tex_h1n1_nmp_dns, col = blues9[2], 
        axis=FALSE, names = c(expression(italic("H1N1, Puerto Rico / 1934"),
                                         italic("H1N1, S.Canterbury / 2000"),
                                         italic("H1N1, St.Petersburg / 2006"),
                                         italic("H1N1, Texas / 2007"))))
# Boxplot about signal densities
layout(mat = matrix(c(1, 2, 1, 2, 3, 3), 
                    nrow = 3, 
                    ncol = 2,
                    byrow=TRUE))
par(mar=c(5.1, 4.1, 0.1, 1.1), tck=0)

boxplot(usa300_nmp_sgd, mu50_nmp_sgd, col = blues9[2], axis=FALSE, 
        names=expression(italic("S. aureus, USA300"), italic("S. aureus, Mu50")))
boxplot(cov_nmp_sgd, col = blues9[2], show.names=TRUE, names=expression(italic("SARS-CoV-2")))
boxplot(pur_h1n1_nmp_sgd, scb_h1n1_nmp_sgd, stp_h1n1_nmp_sgd, tex_h1n1_nmp_sgd, col = blues9[2], 
        axis=FALSE, names = c(expression(italic("H1N1, Puerto Rico / 1934"),
                                         italic("H1N1, S.Canterbury / 2000"),
                                         italic("H1N1, St.Petersburg / 2006"),
                                         italic("H1N1, Texas / 2007"))))
# Boxplot about epitope counts
layout(mat = matrix(c(1, 2, 1, 2, 3, 3), 
                    nrow = 3, 
                    ncol = 2,
                    byrow=TRUE))
par(mar=c(5.1, 4.1, 0.1, 1.1), tck=0)

boxplot(usa300_nmp_cnt, mu50_nmp_cnt, col = blues9[2], axis=FALSE, 
        names=expression(italic("S. aureus, USA300"), italic("S. aureus, Mu50")))
boxplot(cov_nmp_cnt, col = blues9[2], show.names=TRUE, names=expression(italic("SARS-CoV-2")))
boxplot(pur_h1n1_nmp_cnt, scb_h1n1_nmp_cnt, stp_h1n1_nmp_cnt, tex_h1n1_nmp_cnt, col = blues9[2], 
        axis=FALSE, names = c(expression(italic("H1N1, Puerto Rico / 1934"),
                                         italic("H1N1, S.Canterbury / 2000"),
                                         italic("H1N1, St.Petersburg / 2006"),
                                         italic("H1N1, Texas / 2007"))))
par(mar=c(5.1, 4.1, 4.1, 1.1), tck=NA)

# Plot Correlations
comp_up_cor_overview <- matrix(c(rep("",4),rep(c("N","","S",""),6),
                                 format(c(usa300_lenmas_cor,
                                                     "<1e-318",
                                                    usa300_lensig_cor,
                                                    usa300_lensig_pvl,
                                                    usa300_lencnt_nmp_cor,
                                                    usa300_lencnt_nmp_pvl,
                                                    usa300_lencnt_spt_cor,
                                                    usa300_lencnt_spt_pvl,
                                                    usa300_lendns_nmp_cor,
                                                    usa300_lendns_nmp_pvl,
                                                    usa300_lendns_spt_cor,
                                                    usa300_lendns_spt_pvl,
                                                    usa300_lensgd_nmp_cor,
                                                    usa300_lensgd_nmp_pvl,
                                                    usa300_lensgd_spt_cor,
                                                    usa300_lensgd_spt_pvl,
                                                    usa300_cntdns_nmp_cor,
                                                    usa300_cntdns_nmp_pvl,
                                                    usa300_cntdns_spt_cor,
                                                    usa300_cntdns_spt_pvl,
                                                    usa300_cntsgd_nmp_cor,
                                                    usa300_cntsgd_nmp_pvl,
                                                    usa300_cntsgd_spt_cor,
                                                    usa300_cntsgd_spt_pvl,
                                                    usa300_dnssgd_nmp_cor,
                                                    usa300_dnssgd_nmp_pvl,
                                                    usa300_dnssgd_spt_cor,
                                                    usa300_dnssgd_spt_pvl),
                                       digits = 2),
                                 format(c(mu50_lenmas_cor,
                                                     "<1e-318",
                                                    mu50_lensig_cor,
                                                    mu50_lensig_pvl,
                                                    mu50_lencnt_nmp_cor,
                                                    mu50_lencnt_nmp_pvl,
                                                    mu50_lencnt_spt_cor,
                                                    mu50_lencnt_spt_pvl,
                                                    mu50_lendns_nmp_cor,
                                                    mu50_lendns_nmp_pvl,
                                                    mu50_lendns_spt_cor,
                                                    mu50_lendns_spt_pvl,
                                                    mu50_lensgd_nmp_cor,
                                                    mu50_lensgd_nmp_pvl,
                                                    mu50_lensgd_spt_cor,
                                                    mu50_lensgd_spt_pvl,
                                                    mu50_cntdns_nmp_cor,
                                                    mu50_cntdns_nmp_pvl,
                                                    mu50_cntdns_spt_cor,
                                                    mu50_cntdns_spt_pvl,
                                                    mu50_cntsgd_nmp_cor,
                                                    mu50_cntsgd_nmp_pvl,
                                                    mu50_cntsgd_spt_cor,
                                                    mu50_cntsgd_spt_pvl,
                                                    mu50_dnssgd_nmp_cor,
                                                    mu50_dnssgd_nmp_pvl,
                                                    mu50_dnssgd_spt_cor,
                                                    mu50_dnssgd_spt_pvl),
                                       digits = 2),
                                 format(c(cov_lenmas_cor,
                                                    cov_lenmas_pvl,
                                                    cov_lensig_cor,
                                                    "<1e-318",
                                                    cov_lencnt_nmp_cor,
                                                    cov_lencnt_nmp_pvl,
                                                    cov_lencnt_spt_cor,
                                                    cov_lencnt_spt_pvl,
                                                    cov_lendns_nmp_cor,
                                                    cov_lendns_nmp_pvl,
                                                    cov_lendns_spt_cor,
                                                    cov_lendns_spt_pvl,
                                                    cov_lensgd_nmp_cor,
                                                    cov_lensgd_nmp_pvl,
                                                    cov_lensgd_spt_cor,
                                                    cov_lensgd_spt_pvl,
                                                    cov_cntdns_nmp_cor,
                                                    cov_cntdns_nmp_pvl,
                                                    cov_cntdns_spt_cor,
                                                    cov_cntdns_spt_pvl,
                                                    cov_cntsgd_nmp_cor,
                                                    cov_cntsgd_nmp_pvl,
                                                    cov_cntsgd_spt_cor,
                                                    cov_cntsgd_spt_pvl,
                                                    cov_dnssgd_nmp_cor,
                                                    cov_dnssgd_nmp_pvl,
                                                    cov_dnssgd_spt_cor,
                                                    cov_dnssgd_spt_pvl),
                                       digits = 2),
                                 
                                 format(c(pur_h1n1_lenmas_cor,
                                                    pur_h1n1_lenmas_pvl), digits = 2),
                                 rep("N/A",2),
                                 format(c(pur_h1n1_lencnt_nmp_cor,
                                                    pur_h1n1_lencnt_nmp_pvl,
                                                    pur_h1n1_lencnt_spt_cor,
                                                    pur_h1n1_lencnt_spt_pvl,
                                                    pur_h1n1_lendns_nmp_cor,
                                                    pur_h1n1_lendns_nmp_pvl,
                                                    pur_h1n1_lendns_spt_cor,
                                                    pur_h1n1_lendns_spt_pvl), digits = 2),
                                 rep("N/A",4),
                                 format(c(pur_h1n1_cntdns_nmp_cor,
                                                    pur_h1n1_cntdns_nmp_pvl,
                                                    pur_h1n1_cntdns_spt_cor,
                                                    pur_h1n1_cntdns_spt_pvl), digits = 2),
                                 rep("N/A",8),
                                 
                                 format(c(scb_h1n1_lenmas_cor,
                                                    scb_h1n1_lenmas_pvl), digits = 2),
                                 rep("N/A",2),
                                 format(c(scb_h1n1_lencnt_nmp_cor,
                                                    scb_h1n1_lencnt_nmp_pvl,
                                                    scb_h1n1_lencnt_spt_cor,
                                                    scb_h1n1_lencnt_spt_pvl,
                                                    scb_h1n1_lendns_nmp_cor,
                                                    scb_h1n1_lendns_nmp_pvl,
                                                    scb_h1n1_lendns_spt_cor,
                                                    scb_h1n1_lendns_spt_pvl), digits = 2),
                                 rep("N/A",4),
                                 format(c(scb_h1n1_cntdns_nmp_cor,
                                                    scb_h1n1_cntdns_nmp_pvl,
                                                    scb_h1n1_cntdns_spt_cor,
                                                    scb_h1n1_cntdns_spt_pvl), digits = 2),
                                 rep("N/A",8),
                                 
                                 format(c(stp_h1n1_lenmas_cor,
                                                    stp_h1n1_lenmas_pvl), digits = 2),
                                 rep("N/A",2),
                                 format(c(stp_h1n1_lencnt_nmp_cor,
                                                    stp_h1n1_lencnt_nmp_pvl,
                                                    stp_h1n1_lencnt_spt_cor,
                                                    stp_h1n1_lencnt_spt_pvl,
                                                    stp_h1n1_lendns_nmp_cor,
                                                    stp_h1n1_lendns_nmp_pvl,
                                                    stp_h1n1_lendns_spt_cor,
                                                    stp_h1n1_lendns_spt_pvl), digits = 2),
                                 rep("N/A",4),
                                 format(c(stp_h1n1_cntdns_nmp_cor,
                                                    stp_h1n1_cntdns_nmp_pvl,
                                                    stp_h1n1_cntdns_spt_cor,
                                                    stp_h1n1_cntdns_spt_pvl), digits = 2),
                                 rep("N/A",8),
                                 
                                 format(c(tex_h1n1_lenmas_cor,
                                                    tex_h1n1_lenmas_pvl), digits = 2),
                                 rep("N/A",2),
                                 format(c(tex_h1n1_lencnt_nmp_cor,
                                                    tex_h1n1_lencnt_nmp_pvl,
                                                    tex_h1n1_lencnt_spt_cor,
                                                    tex_h1n1_lencnt_spt_pvl,
                                                    tex_h1n1_lendns_nmp_cor,
                                                    tex_h1n1_lendns_nmp_pvl,
                                                    tex_h1n1_lendns_spt_cor,
                                                    tex_h1n1_lendns_spt_pvl), digits = 2),
                                 rep("N/A",4),
                                 format(c(tex_h1n1_cntdns_nmp_cor,
                                                    tex_h1n1_cntdns_nmp_pvl,
                                                    tex_h1n1_cntdns_spt_cor,
                                                    tex_h1n1_cntdns_spt_pvl), digits = 2),
                                 rep("N/A",8),
                                 
                                 "Korrelation", "P-Wert",
                                 rep("",26)),
                               nrow = 28, ncol = 9)
rownames(comp_up_cor_overview) <- c("Länge ~ Masse","", "Länge ~ Signallänge","",
                                      "Länge ~ Anz. Epitope","",
                                      "Länge ~ Anz. Epitope","",
                                      "Länge ~ Epitopdichte","",
                                      "Länge ~ Epitopdichte","",
                                      "Länge ~ Signal-Epitopdichte","",
                                      "Länge ~ Signal-Epitopdichte","",
                                      "Anz. Epitope ~ Epitopdichte","",
                                      "Anz. Epitope ~ Epitopdichte","",
                                      "Anz. Epitope ~ Signal-Epitopdichte","",
                                      "Anz. Epitope ~ Signal-Epitopdichte","",
                                      "Epitopdichte ~ Signal-Epitopdichte","",
                                      "Epitopdichte ~ Signal-Epitopdichte","")
colnames(comp_up_cor_overview) <- c("","USA300", "Mu50", "SARS-CoV-2", "Puerto Rico",
                                    "S. Canterbury", "St.Petersburg", "Texas","")
comp_up_cor_title <- textGrob("Sequenzkorrelationen mittels Spearman")
comp_up_cor_table <- tableGrob(comp_up_cor_overview, theme = table_theme4)

l <- comp_up_cor_table$layout
for (row in 2:29) {
  comp_up_cor_table$grobs[ which(l$t==row & l$l==2 & l$name=="core-bg")][[1]][["gp"]] <- gpar(fill="white", col="white")
  comp_up_cor_table$grobs[ which(l$t==row & l$l==10 & l$name=="core-bg")][[1]][["gp"]] <- gpar(fill="white", col="white")
  comp_up_cor_table$grobs[ which(l$t==row & l$l==2 & l$name=="core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
  comp_up_cor_table$grobs[ which(l$t==row & l$l==10 & l$name=="core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
}
grid.arrange(comp_up_cor_title, comp_up_cor_table,
             heights = unit(c(10,230), rep("mm",2)))


comp_uniProt_all <- rbind(usa300_uniProt, mu50_uniProt, cov_uniProt, pur_h1n1_uniProt,
                          scb_h1n1_uniProt, stp_h1n1_uniProt, tex_h1n1_uniProt)
comp_uniProt_all$STRAIN <- c(rep("USA300",nrow(usa300_uniProt)),
                             rep("Mu50",nrow(mu50_uniProt)),
                             rep("SARS-CoV-2",nrow(cov_uniProt)),
                             rep("H1N1 Puerto Rico/1934",nrow(pur_h1n1_uniProt)),
                             rep("H1N1 S.Canterbury/2000",nrow(scb_h1n1_uniProt)),
                             rep("H1N1 St.Petersburg/2006",nrow(stp_h1n1_uniProt)),
                             rep("H1N1 Texas/2007",nrow(tex_h1n1_uniProt)))
comp_nmp_all <- rbind.fill(usa300_nmp, mu50_nmp, cov_nmp, pur_h1n1_nmp,
                          scb_h1n1_nmp, stp_h1n1_nmp, tex_h1n1_nmp)
comp_nmp_all$STRAIN <- c(rep("USA300",nrow(usa300_nmp)),
                             rep("Mu50",nrow(mu50_nmp)),
                             rep("SARS-CoV-2",nrow(cov_nmp)),
                             rep("H1N1 Puerto Rico/1934",nrow(pur_h1n1_nmp)),
                             rep("H1N1 S.Canterbury/2000",nrow(scb_h1n1_nmp)),
                             rep("H1N1 St.Petersburg/2006",nrow(stp_h1n1_nmp)),
                             rep("H1N1 Texas/2007",nrow(tex_h1n1_nmp)))
comp_spt_all <- rbind.fill(usa300_spt, mu50_spt, cov_spt, pur_h1n1_spt,
                          scb_h1n1_spt, stp_h1n1_spt, tex_h1n1_spt)
comp_spt_all$STRAIN <- c(rep("USA300",nrow(usa300_spt)),
                             rep("Mu50",nrow(mu50_spt)),
                             rep("SARS-CoV-2",nrow(cov_spt)),
                             rep("H1N1 Puerto Rico/1934",nrow(pur_h1n1_spt)),
                             rep("H1N1 S.Canterbury/2000",nrow(scb_h1n1_spt)),
                             rep("H1N1 St.Petersburg/2006",nrow(stp_h1n1_spt)),
                             rep("H1N1 Texas/2007",nrow(tex_h1n1_spt)))
comp_nmp_cnt_all <- count_epis(comp_nmp_all)
comp_spt_cnt_all <- count_epis(comp_spt_all)

comp_cnts <- function(){
  plot(comp_nmp_cnt_all, comp_spt_cnt_all, xlab = "Anzahl Epitope, NetMHCpan", ylab = "Anzahl Epitope, SYFPEITHI",
       main = expression(paste("Sequenzen der Vergleichsstämme von ", italic("S. aureus, SARS-CoV-2, Influenza A"))), col.main="white")
  abline(coef = c(0,1), col="red")}

# Strain Variations
comp_lenloc <- lm(as.numeric(unlist(comp_uniProt_all$LENGTH))~comp_uniProt_all$STRAIN)
comp_masloc <- lm(as.numeric(unlist(comp_uniProt_all$MASS))~comp_uniProt_all$STRAIN)
comp_sigloc <- lm(as.numeric(unlist(comp_uniProt_all$SIGNAL))~comp_uniProt_all$STRAIN)
comp_cntloc_nmp <- lm(comp_nmp_cnt_all~comp_nmp_all$STRAIN)
comp_dnsloc_nmp <- lm(as.numeric(unlist(comp_nmp_all$DENSITY))~comp_nmp_all$STRAIN)
comp_sgdloc_nmp <- lm(as.numeric(unlist(comp_nmp_all$`SIG DENSITY`))~comp_nmp_all$STRAIN)
comp_cntloc_spt <- lm(comp_spt_cnt_all~comp_spt_all$STRAIN)
comp_dnsloc_spt <- lm(as.numeric(unlist(comp_spt_all$DENSITY))~comp_spt_all$STRAIN)
comp_sgdloc_spt <- lm(as.numeric(unlist(comp_spt_all$`SIG DENSITY`))~comp_spt_all$STRAIN)

par(mfrow=c(3,3))
hist(rstandard(comp_lenloc))
hist(rstandard(comp_masloc))
hist(rstandard(comp_sigloc))
hist(rstandard(comp_cntloc_nmp))
hist(rstandard(comp_dnsloc_nmp))
hist(rstandard(comp_sgdloc_nmp))
hist(rstandard(comp_cntloc_spt))
hist(rstandard(comp_dnsloc_spt))
hist(rstandard(comp_sgdloc_spt))
par(mfrow=c(1,1))

comp_lenloc_fvl <- unlist(summary(aov(as.numeric(unlist(comp_uniProt_all$LENGTH))~comp_uniProt_all$STRAIN)))["F value1"]
comp_lenloc_pvl <- unlist(summary(aov(as.numeric(unlist(comp_uniProt_all$LENGTH))~comp_uniProt_all$STRAIN)))["Pr(>F)1"]
comp_masloc_fvl <- unlist(summary(aov(as.numeric(unlist(comp_uniProt_all$MASS))~comp_uniProt_all$STRAIN)))["F value1"]
comp_masloc_pvl <- unlist(summary(aov(as.numeric(unlist(comp_uniProt_all$MASS))~comp_uniProt_all$STRAIN)))["Pr(>F)1"]
comp_sigloc_fvl <- unlist(summary(aov(as.numeric(unlist(comp_uniProt_all$SIGNAL))~comp_uniProt_all$STRAIN)))["F value1"]
comp_sigloc_pvl <- unlist(summary(aov(as.numeric(unlist(comp_uniProt_all$SIGNAL))~comp_uniProt_all$STRAIN)))["Pr(>F)1"]

comp_cntloc_nmp_fvl <- unlist(summary(aov(comp_nmp_cnt_all~comp_nmp_all$STRAIN)))["F value1"]
comp_cntloc_nmp_pvl <- unlist(summary(aov(comp_nmp_cnt_all~comp_nmp_all$STRAIN)))["Pr(>F)1"]
comp_dnsloc_nmp_fvl <- unlist(summary(aov(as.numeric(unlist(comp_nmp_all$DENSITY))~comp_nmp_all$STRAIN)))["F value1"]
comp_dnsloc_nmp_pvl <- unlist(summary(aov(as.numeric(unlist(comp_nmp_all$DENSITY))~comp_nmp_all$STRAIN)))["Pr(>F)1"]
comp_sgdloc_nmp_fvl <- unlist(summary(aov(as.numeric(unlist(comp_nmp_all$`SIG DENSITY`))~comp_nmp_all$STRAIN)))["F value1"]
comp_sgdloc_nmp_pvl <- unlist(summary(aov(as.numeric(unlist(comp_nmp_all$`SIG DENSITY`))~comp_nmp_all$STRAIN)))["Pr(>F)1"]

comp_cntloc_spt_fvl <- unlist(summary(aov(comp_spt_cnt_all~comp_spt_all$STRAIN)))["F value1"]
comp_cntloc_spt_pvl <- unlist(summary(aov(comp_spt_cnt_all~comp_spt_all$STRAIN)))["Pr(>F)1"]
comp_dnsloc_spt_fvl <- unlist(summary(aov(as.numeric(unlist(comp_spt_all$DENSITY))~comp_spt_all$STRAIN)))["F value1"]
comp_dnsloc_spt_pvl <- unlist(summary(aov(as.numeric(unlist(comp_spt_all$DENSITY))~comp_spt_all$STRAIN)))["Pr(>F)1"]
comp_sgdloc_spt_fvl <- unlist(summary(aov(as.numeric(unlist(comp_spt_all$`SIG DENSITY`))~comp_spt_all$STRAIN)))["F value1"]
comp_sgdloc_spt_pvl <- unlist(summary(aov(as.numeric(unlist(comp_spt_all$`SIG DENSITY`))~comp_spt_all$STRAIN)))["Pr(>F)1"]

comp_loc_cor_overview <- matrix(c(rep("",3), rep(c("N","S"),3),
                                  round(as.numeric(c(comp_lenloc_fvl,
                                                   comp_masloc_fvl,
                                                   comp_sigloc_fvl,
                                                   comp_cntloc_nmp_fvl,
                                                   comp_cntloc_spt_fvl,
                                                   comp_dnsloc_nmp_fvl,
                                                   comp_dnsloc_spt_fvl,
                                                   comp_sgdloc_nmp_fvl,
                                                   comp_sgdloc_spt_fvl)),
                                      digits = 2),
                                format(as.numeric(c(comp_lenloc_pvl,
                                                   comp_masloc_pvl,
                                                   comp_sigloc_pvl,
                                                   comp_cntloc_nmp_pvl,
                                                   comp_cntloc_spt_pvl,
                                                   comp_dnsloc_nmp_pvl,
                                                   comp_dnsloc_spt_pvl,
                                                   comp_sgdloc_nmp_pvl,
                                                   comp_sgdloc_spt_pvl)),
                                      digits = 2)),
                              nrow = 9, ncol = 3)
rownames(comp_loc_cor_overview) <- c("Länge","Masse","Signallänge","Anz. Epitope",
                                   "Anz. Epitope","Epitopdichte","Epitopdichte",
                                   "Signal-Epitopdichte","Signal-Epitopdichte")
colnames(comp_loc_cor_overview) <- c("","F-Wert", "P-Wert")
comp_loc_cor_title1 <- textGrob("Varianzanalyse der Sequenz- und Epitopdaten unter den Stämmen")
comp_loc_cor_title2 <- textGrob(expression(paste("Stämme von ", italic("S. aureus, SARS-Coronaviren, Influenzaviren A"))))
comp_loc_cor_table <- tableGrob(comp_loc_cor_overview, theme = table_theme1)

l <- comp_loc_cor_table$layout
for (row in 2:10) {
  comp_loc_cor_table$grobs[ which(l$t==row & l$l==2 & l$name=="core-bg")][[1]][["gp"]] <- gpar(fill="white", col="white")
  comp_loc_cor_table$grobs[ which(l$t==row & l$l==2 & l$name=="core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
}
grid.arrange(comp_loc_cor_title1, comp_loc_cor_title2, comp_loc_cor_table, heights = unit(c(5,20,80), rep("mm", 3)))

par(mfrow=c(3,3), cex.lab=1.5, cex.axis=1.5)
tuk_plot(TukeyHSD(aov(as.numeric(unlist(comp_uniProt_all$LENGTH))~comp_uniProt_all$STRAIN)), 
         "Länge", "",rev(c("SC-PR","SP-PR","TX-PR","MU-PR","CV-PR","US-PR","SP-SC",
                                    "TX-SC","MU-SC","CV-SC","US-SC","TX-SP","MU-SP","CV-SP",
                                    "US-SP","MU-TX","CV-TX","US-TX","CV-MU","US-MU","US-CV")), 1,1)
tuk_plot(TukeyHSD(aov(as.numeric(unlist(comp_uniProt_all$MASS))~comp_uniProt_all$STRAIN)),
         "Masse", "",rev(c("SC-PR","SP-PR","TX-PR","MU-PR","CV-PR","US-PR","SP-SC",
                                    "TX-SC","MU-SC","CV-SC","US-SC","TX-SP","MU-SP","CV-SP",
                                    "US-SP","MU-TX","CV-TX","US-TX","CV-MU","US-MU","US-CV")), 1,1)
tuk_plot(TukeyHSD(aov(as.numeric(unlist(comp_uniProt_all$SIGNAL))~comp_uniProt_all$STRAIN)),
         "Signallänge", "",rev(c("SC-PR","SP-PR","TX-PR","MU-PR","CV-PR","US-PR","SP-SC",
                                 "TX-SC","MU-SC","CV-SC","US-SC","TX-SP","MU-SP","CV-SP",
                                 "US-SP","MU-TX","CV-TX","US-TX","CV-MU","US-MU","US-CV")), 1,1)
tuk_plot(TukeyHSD(aov(comp_nmp_cnt_all~comp_nmp_all$STRAIN)),
         "Anz. Epitope, NetMHCpan", "",rev(c("SC-PR","SP-PR","TX-PR","MU-PR","CV-PR","US-PR","SP-SC",
                                                      "TX-SC","MU-SC","CV-SC","US-SC","TX-SP","MU-SP","CV-SP",
                                                      "US-SP","MU-TX","CV-TX","US-TX","CV-MU","US-MU","US-CV")), 1,1)
tuk_plot(TukeyHSD(aov(as.numeric(unlist(comp_nmp_all$DENSITY))~comp_nmp_all$STRAIN)),
         "Epitopdichte, NetMHCpan", "",rev(c("SC-PR","SP-PR","TX-PR","MU-PR","CV-PR","US-PR","SP-SC",
                                                      "TX-SC","MU-SC","CV-SC","US-SC","TX-SP","MU-SP","CV-SP",
                                                      "US-SP","MU-TX","CV-TX","US-TX","CV-MU","US-MU","US-CV")), 1,1)
tuk_plot(TukeyHSD(aov(as.numeric(unlist(comp_nmp_all$`SIG DENSITY`))~comp_nmp_all$STRAIN)),
         "Signal-Epitopdichte, NetMHCpan", "",rev(c("SC-PR","SP-PR","TX-PR","MU-PR","CV-PR","US-PR","SP-SC",
                                                    "TX-SC","MU-SC","CV-SC","US-SC","TX-SP","MU-SP","CV-SP",
                                                    "US-SP","MU-TX","CV-TX","US-TX","CV-MU","US-MU","US-CV")), 1,1)
tuk_plot(TukeyHSD(aov(comp_spt_cnt_all~comp_spt_all$STRAIN)),
         "Anz. Epitope, SYFPEITHI", "",rev(c("SC-PR","SP-PR","TX-PR","MU-PR","CV-PR","US-PR","SP-SC",
                                             "TX-SC","MU-SC","CV-SC","US-SC","TX-SP","MU-SP","CV-SP",
                                             "US-SP","MU-TX","CV-TX","US-TX","CV-MU","US-MU","US-CV")), 1,1)

tuk_plot(TukeyHSD(aov(as.numeric(unlist(comp_spt_all$DENSITY))~comp_spt_all$STRAIN)),
         "Epitopdichte, SYFPEITHI", "",rev(c("SC-PR","SP-PR","TX-PR","MU-PR","CV-PR","US-PR","SP-SC",
                                                      "TX-SC","MU-SC","CV-SC","US-SC","TX-SP","MU-SP","CV-SP",
                                                      "US-SP","MU-TX","CV-TX","US-TX","CV-MU","US-MU","US-CV")), 1,1)
tuk_plot(TukeyHSD(aov(as.numeric(unlist(comp_spt_all$`SIG DENSITY`))~comp_spt_all$STRAIN)),
         "Singal-Epitopdichte, SYFPEITHI", "",rev(c("SC-PR","SP-PR","TX-PR","MU-PR","CV-PR","US-PR","SP-SC",
                                                             "TX-SC","MU-SC","CV-SC","US-SC","TX-SP","MU-SP","CV-SP",
                                                             "US-SP","MU-TX","CV-TX","US-TX","CV-MU","US-MU","US-CV")), 1,1)
par(mfrow=c(1,1))





