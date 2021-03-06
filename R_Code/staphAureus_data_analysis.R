####################### Data collection S. aureus #############################
###############################################################################
# Before executing this file load the 'data_analysis_functions.R in this      #
# Environment.                                                                #
###############################################################################

##### Proteins #####

# Read UniProt Files
sa_uniProt_cw <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/example_SA_uniProt.xlsx", sheet = "cell_wall")
sa_uniProt_cp <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/example_SA_uniProt.xlsx", sheet = "cytoplasm")
sa_uniProt_sc <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/example_SA_uniProt.xlsx", sheet = "secreted")
sa_uniProt_cm <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/example_SA_uniProt.xlsx", sheet = "cell_membrane")
sa_uniProt_alllocs <- rbind.fill(sa_uniProt_cw,sa_uniProt_cm,sa_uniProt_cp,sa_uniProt_sc)
sa_uniProt_alllocs$LOCATION=c(rep("cell_wall",nrow(sa_uniProt_cw)),rep("cell_membrane",nrow(sa_uniProt_cm)), 
                              rep("cytoplasm",nrow(sa_uniProt_cp)), rep("secreted",nrow(sa_uniProt_sc)))
sa_all_seqs <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/sa_distinct_uniProt.xlsx", sheet = "All")

# Lengths
sa_len_cw <- as.numeric( unlist( sa_uniProt_cw["LENGTH"]))
sa_len_cm <- as.numeric( unlist( sa_uniProt_cm["LENGTH"]))
sa_len_cp <- as.numeric( unlist( sa_uniProt_cp["LENGTH"]))
sa_len_sc <- as.numeric( unlist( sa_uniProt_sc["LENGTH"]))
sa_len_alllocs <- c( sa_len_cw, sa_len_cp, sa_len_sc)
# Distinct sequences
sa_len_all_seq <- as.numeric( unlist( sa_all_seqs["LENGTH"]))

# Masses
sa_mas_cw <- as.numeric(gsub(",","",unlist(sa_uniProt_cw["MASS"], use.names = FALSE)))
sa_mas_cm <- as.numeric(gsub(",","",unlist(sa_uniProt_cm["MASS"], use.names = FALSE)))
sa_mas_cp <- as.numeric(gsub(",","",unlist(sa_uniProt_cp["MASS"], use.names = FALSE)))
sa_mas_sc <- as.numeric(gsub(",","",unlist(sa_uniProt_sc["MASS"], use.names = FALSE)))
sa_mas_alllocs <- c( sa_mas_cw, sa_mas_cm, sa_mas_cp, sa_mas_sc)
sa_mas_all_seq <- as.numeric( gsub( ",", "", unlist( sa_all_seqs["MASS"], use.names = FALSE)))

# Signal regions
sa_sig_cw_na <- as.numeric( unlist( sa_uniProt_cw["SIGNAL"], use.names = FALSE))
sa_sig_cm_na <- as.numeric( unlist( sa_uniProt_cm["SIGNAL"], use.names = FALSE))
sa_sig_cp_na <- as.numeric( unlist( sa_uniProt_cp["SIGNAL"], use.names = FALSE))
sa_sig_sc_na <- as.numeric( unlist( sa_uniProt_sc["SIGNAL"], use.names = FALSE))
sa_sig_alllocs_na <- c(sa_sig_cw_na, sa_sig_cm_na, sa_sig_cp_na, sa_sig_sc_na)
sa_sig_alllocs <- na.omit(sa_sig_alllocs_na)
sa_sig_all_seq_na <- as.numeric( unlist( sa_all_seqs["SIGNAL"], use.names = FALSE))
sa_sig_all_seq <- na.omit(sa_sig_all_seq_na)

# Plot UniProt overview
sa_UP_overview <- matrix(c(nrow(sa_all_seqs), 
                           nrow(sa_uniProt_cw),
                           nrow(sa_uniProt_cp),
                           nrow(sa_uniProt_sc),
                           nrow(sa_uniProt_cm),
                           round(c(mean(sa_len_all_seq),
                                   mean(sa_len_cw),
                                   mean(sa_len_cp),
                                   mean(sa_len_sc),
                                   mean(sa_len_cm))),
                           length(sa_sig_all_seq),
                           length(na.omit(sa_sig_cw_na)),
                           length(na.omit(sa_sig_cp_na)),
                           length(na.omit(sa_sig_sc_na)),
                           length(na.omit(sa_sig_cm_na)),
                           round(mean(sa_sig_all_seq)),
                           round(mean(sa_sig_cw_na)),
                           "N/A", 
                           round(mean(na.omit(sa_sig_sc_na))),
                           round(mean(na.omit(sa_sig_cm_na)))),
                         nrow = 5, ncol = 4)
colnames(sa_UP_overview) <- c("Anz. Sequenzen", "L�nge [AS]", "Anz. Signalpeptide","Signall�nge [AS]")
rownames(sa_UP_overview) <- c("Proteom", "Zellwand", "Zytoplasma", "Extrazellularraum", "Zellmembran")
sa_up_table <- tableGrob(sa_UP_overview, theme = table_theme)
sa_up_title <- textGrob("UniProt �bersicht")
sa_up_subtt <- textGrob(expression(italic("S. aureus")))
sa_up_table <- gtable_add_grob(sa_up_table, grobs = segmentsGrob(x0 = unit(0,"npc"),
                                                                 y0 = unit(0,"npc"),
                                                                 x1 = unit(1,"npc"),
                                                                 y1 = unit(0,"npc"),
                                                                 gp = gpar(lwd = 1.5)),
                               t = 2, l = 1, r = 5 )
grid.arrange(sa_up_title, sa_up_subtt, sa_up_table, ncol=1,
             heights = unit(c(10,10,40), rep("mm",3)))

# Boxplot of all lengths
par(mar=c(3.1, 4.1, 0.1, 1.1), tck=0)
boxplot(sa_len_all_seq, sa_len_cw,sa_len_cp,sa_len_sc,sa_len_cm, col = blues9[2], axis=FALSE,
        names=expression(paste(italic("S. aureus"),", Proteom"), "Zellwand", "Zytoplasma",
                         "Extrazellularraum", "Zellmembran"))

par(mar=c(5.1, 4.1, 4.1, 1.1), tck=NA)


##### Epitopes #####

# Read NetMHCpan files
sa_nmp_cw <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/example_SA_netMHCPan.xlsx", sheet = "cell_wall")
sa_nmp_cm <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/example_SA_netMHCPan.xlsx", sheet = "cell_membrane")
sa_nmp_cp <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/example_SA_netMHCPan.xlsx", sheet = "cytoplasm")
sa_nmp_sc <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/example_SA_netMHCPan.xlsx", sheet = "secreted")
sa_nmp_alllocs <- rbind.fill(sa_nmp_cw,sa_nmp_cm,sa_nmp_cp,sa_nmp_sc)
sa_nmp_alllocs$LOCATION=c(rep("cell_wall",nrow(sa_nmp_cw)),rep("cell_membrane",nrow(sa_nmp_cm)),
                          rep("cytoplasm",nrow(sa_nmp_cp)),rep("secreted",nrow(sa_nmp_sc)))
sa_nmp_all_seq <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/sa_distinct_netMHCPan.xlsx", sheet = "All")

# Read SYFPEITHI files
sa_spt_cw <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/example_SA_SYFPEITHI.xlsx", sheet = "cell_wall")
sa_spt_cm <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/example_SA_SYFPEITHI.xlsx", sheet = "cell_membrane")
sa_spt_cp <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/example_SA_SYFPEITHI.xlsx", sheet = "cytoplasm")
sa_spt_sc <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/example_SA_SYFPEITHI.xlsx", sheet = "secreted")
sa_spt_alllocs <- rbind.fill(sa_spt_cw,sa_spt_cm,sa_spt_cp,sa_spt_sc)
sa_spt_alllocs$LOCATION=c(rep("cell_wall",nrow(sa_spt_cw)),rep("cell_membrane",nrow(sa_spt_cm)),
                          rep("cytoplasm",nrow(sa_spt_cp)),rep("secreted",nrow(sa_spt_sc)))
sa_spt_all_seq <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/sa_distinct_SYFPEITHI.xlsx", sheet = "All")

# Densities
sa_nmp_dns_cw <- as.numeric(unlist(sa_nmp_cw["DENSITY"]))
sa_nmp_dns_cm <- as.numeric(unlist(sa_nmp_cm["DENSITY"]))
sa_nmp_dns_cp <- as.numeric(unlist(sa_nmp_cp["DENSITY"]))
sa_nmp_dns_sc <- as.numeric(unlist(sa_nmp_sc["DENSITY"]))
sa_nmp_dns_alllocs <- c(sa_nmp_dns_cw, sa_nmp_dns_cm, sa_nmp_dns_cp, sa_nmp_dns_sc)
sa_nmp_dns_all_seq <- as.numeric(unlist(sa_nmp_all_seq["DENSITY"]))

sa_spt_dns_cw <- as.numeric(unlist(sa_spt_cw["DENSITY"]))
sa_spt_dns_cm <- as.numeric(unlist(sa_spt_cm["DENSITY"]))
sa_spt_dns_cp <- as.numeric(unlist(sa_spt_cp["DENSITY"]))
sa_spt_dns_sc <- as.numeric(unlist(sa_spt_sc["DENSITY"]))
sa_spt_dns_alllocs <- unlist(list.append(sa_spt_dns_cw, sa_spt_dns_cm, sa_spt_dns_cp, sa_spt_dns_sc), use.names = FALSE)
sa_spt_dns_all_seq <- as.numeric(unlist(sa_spt_all_seq["DENSITY"]))

# Signal Densities
sa_nmp_sgd_cw <- as.numeric(unlist(sa_nmp_cw["SIG DENSITY"]))
sa_nmp_sgd_cm <- as.numeric(unlist(sa_nmp_cm["SIG DENSITY"]))
sa_nmp_sgd_cp <- as.numeric(unlist(sa_nmp_cp["SIG DENSITY"]))
sa_nmp_sgd_sc <- as.numeric(unlist(sa_nmp_sc["SIG DENSITY"]))
sa_nmp_sgd_alllocs <- as.numeric(unlist(list.append(sa_nmp_sgd_cw, sa_nmp_sgd_cm, sa_nmp_sgd_cp, sa_nmp_sgd_sc), use.names = FALSE))
sa_nmp_sgd_all_seq <- as.numeric(unlist(sa_nmp_all_seq["SIG DENSITY"]))

sa_spt_sgd_cw <- as.numeric(unlist(sa_spt_cw["SIG DENSITY"]))
sa_spt_sgd_cm <- as.numeric(unlist(sa_spt_cm["SIG DENSITY"]))
sa_spt_sgd_cp <- as.numeric(unlist(sa_spt_cp["SIG DENSITY"]))
sa_spt_sgd_sc <- as.numeric(unlist(sa_spt_sc["SIG DENSITY"]))
sa_spt_sgd_alllocs <- as.numeric(unlist(list.append(sa_spt_sgd_cw, sa_spt_sgd_cw, sa_spt_sgd_cp, sa_spt_sgd_sc), use.names = FALSE))
sa_spt_sgd_all_seq <- as.numeric(unlist(sa_spt_all_seq["SIG DENSITY"]))

# Epitope counts
sa_nmp_cnt_cw <- count_epis(sa_nmp_cw)
sa_nmp_cnt_cm <- count_epis(sa_nmp_cm)
sa_nmp_cnt_cp <- count_epis(sa_nmp_cp)
sa_nmp_cnt_sc <- count_epis(sa_nmp_sc)
sa_nmp_cnt_alllocs <- c(sa_nmp_cnt_cw, sa_nmp_cnt_cm, sa_nmp_cnt_cp, sa_nmp_cnt_sc)
sa_nmp_cnt_all_seq <- count_epis(sa_nmp_all_seq)

sa_spt_cnt_cw <- count_epis(sa_spt_cw)
sa_spt_cnt_cm <- count_epis(sa_spt_cm)
sa_spt_cnt_cp <- count_epis(sa_spt_cp)
sa_spt_cnt_sc <- count_epis(sa_spt_sc)
sa_spt_cnt_alllocs <- c(sa_spt_cnt_cw, sa_spt_cnt_cm, sa_spt_cnt_cp, sa_spt_cnt_sc)
sa_spt_cnt_all_seq <- count_epis(sa_spt_all_seq)

sa_cnts <- function(){
        plot(sa_nmp_cnt_all_seq, sa_spt_cnt_all_seq, xlab = "Anzahl Epitope, NetMHCpan", ylab = "Anzahl Epitope, SYFPEITHI",
             main = expression(paste("Sequenzen des Proteoms von ", italic("S. aureus"))), col.main="white")
        abline(coef = c(0,1), col="red")}

# Plot prediction overview
sa_Pred_overview_nmp <- matrix(c(round(c(mean(sa_nmp_cnt_all_seq),
                                         mean(sa_nmp_cnt_cw),
                                         mean(sa_nmp_cnt_cp),
                                         mean(sa_nmp_cnt_sc),
                                         mean(sa_nmp_cnt_cm))),
                                 round(c(mean(na.omit(sa_nmp_dns_all_seq)),
                                         mean(na.omit(sa_nmp_dns_cw)),
                                         mean(na.omit(sa_nmp_dns_cp)),
                                         mean(na.omit(sa_nmp_dns_sc)),
                                         mean(na.omit(sa_nmp_dns_cm))),
                                       digits = 3),
                                 round(c(mean(na.omit(as.numeric(unlist(sa_nmp_sgd_all_seq))))*(mean(sa_sig_all_seq)-9),
                                         mean(na.omit(as.numeric(unlist(sa_nmp_sgd_cw))))*(mean(sa_sig_cw_na)-9),
                                         0,
                                         mean(na.omit(as.numeric(unlist(sa_nmp_sgd_sc))))*(mean(na.omit(sa_sig_sc_na))-9),
                                         mean(na.omit(as.numeric(unlist(sa_nmp_sgd_cm))))*(mean(na.omit(sa_sig_cm_na))-9))),
                                 c(round(mean(na.omit(sa_nmp_sgd_all_seq)), digits = 3),
                                   round(mean(na.omit(sa_nmp_sgd_cw)), digits = 3),
                                   "N/A",
                                   round(mean(na.omit(sa_nmp_sgd_sc)), digits = 3),
                                   round(mean(na.omit(sa_nmp_sgd_cm)), digits = 3))),
                               nrow = 5, ncol = 4)
sa_Pred_overview_spt <- matrix(c(round(c(mean(sa_spt_cnt_all_seq),
                                         mean(sa_spt_cnt_cw),
                                         mean(sa_spt_cnt_cp),
                                         mean(sa_spt_cnt_sc),
                                         mean(sa_spt_cnt_cm))),
                                 round(c(mean(na.omit(sa_spt_dns_all_seq)),
                                         mean(na.omit(sa_spt_dns_cw)),
                                         mean(na.omit(sa_spt_dns_cp)),
                                         mean(na.omit(sa_spt_dns_sc)),
                                         mean(na.omit(sa_spt_dns_cm))),
                                       digits = 3),
                                 round(c(mean(na.omit(as.numeric(unlist(sa_spt_sgd_all_seq))))*(mean(sa_sig_all_seq)-9),
                                         mean(na.omit(as.numeric(unlist(sa_spt_sgd_cw))))*(mean(sa_sig_cw_na)-9),
                                         0,
                                         mean(na.omit(as.numeric(unlist(sa_spt_sgd_sc))))*(mean(na.omit(sa_sig_sc_na))-9),
                                         mean(na.omit(as.numeric(unlist(sa_spt_sgd_cm))))*(mean(na.omit(sa_sig_cm_na))-9))),
                                 c(round(mean(na.omit(sa_spt_sgd_all_seq)), digits = 3),
                                   round(mean(na.omit(sa_spt_sgd_cw)), digits = 3),
                                   "N/A",
                                   round(mean(na.omit(sa_spt_sgd_sc)), digits = 3),
                                   round(mean(na.omit(sa_spt_sgd_cm)), digits = 3)),
                                 c(sum(sa_len_all_seq>2000),
                                   sum(sa_len_cw>2000),
                                   sum(sa_len_cp>2000),
                                   sum(sa_len_sc>2000),
                                   sum(sa_len_cm>2000))),
                               nrow = 5, ncol = 5)
colnames(sa_Pred_overview_nmp) <- c("Anz. Epitope", "Epitopdichte", "Anz. Signalepitope", "Signal-Epitopdichte")
colnames(sa_Pred_overview_spt) <- c("Anz. Epitope", "Epitopdichte", "Anz. Signalepitope", "Signal-Epitopdichte", "Gekappte Seq.")
rownames(sa_Pred_overview_nmp) <- c("Proteom", "Zellwand", "Zytoplasma", "Extrazellularraum", "Zellmembran")
rownames(sa_Pred_overview_spt) <- c("Proteom", "Zellwand", "Zytoplasma", "Extrazellularraum", "Zellmembran")

sa_Pred_theme <- ttheme_minimal(colhead = list(bg_params = list(fill="white")),
                                core = list(bg_params = list(fill=c(blues9[1], blues9[2]))))
sa_nmp_table <- tableGrob(sa_Pred_overview_nmp, theme = sa_Pred_theme)
sa_spt_table <- tableGrob(sa_Pred_overview_spt, theme = sa_Pred_theme)
sa_pred_title <- textGrob("Vorhersagen")
sa_pred_subtt <- textGrob(expression(paste(italic("S. aureus"), ", oben: NetMHCpan, unten: SYFPEITHI")))
sa_nmp_table <- gtable_add_grob(sa_nmp_table, grobs = segmentsGrob(x0 = unit(0,"npc"),
                                                                   y0 = unit(0,"npc"),
                                                                   x1 = unit(1,"npc"),
                                                                   y1 = unit(0,"npc"),
                                                                   gp = gpar(lwd = 1.5)),
                                t = 2,l=5,r=1 )
sa_spt_table <- gtable_add_grob(sa_spt_table, grobs = segmentsGrob(x0 = unit(0,"npc"),
                                                                   y0 = unit(0,"npc"),
                                                                   x1 = unit(1,"npc"),
                                                                   y1 = unit(0,"npc"),
                                                                   gp = gpar(lwd = 1.5)),
                                t = 2,l=6,r=1 )
A <- textGrob("A")
B <- textGrob("B")
space <- textGrob("")
grid.arrange(sa_pred_title, space, sa_pred_subtt, space, sa_nmp_table, A,
             sa_spt_table, B, ncol=2, heights = unit(c(10,10,50,50), rep("mm", 4)),
             widths = unit(c(150,90), rep("mm",2)))

# Boxplots about densities
boxplot(sa_nmp_dns_all_seq,sa_nmp_dns_cw,sa_nmp_dns_cp,sa_nmp_dns_sc,sa_nmp_dns_cm,
        names = c(expression(paste(italic("S. aureus"),", Proteom"), "Zellwand", "Zytoplasma",
                             "Extrazellularraum", "Zellmembran")), col = blues9[2])
# Boxplots about signal densities
boxplot(sa_nmp_sgd_all_seq,sa_nmp_sgd_cw,sa_nmp_sgd_cp,sa_nmp_sgd_sc,sa_nmp_sgd_cm,
        names = c(expression(paste(italic("S. aureus"),", Proteom"), "Zellwand", "Zytoplasma",
                             "Extrazellularraum", "Zellmembran")), col = blues9[2])
# Boxplots about epitope numbers
boxplot(sa_nmp_cnt_all_seq,sa_nmp_cnt_cw,sa_nmp_cnt_cp,sa_nmp_cnt_sc,sa_nmp_cnt_cm,
        names = c(expression(paste(italic("S. aureus"),", Proteom"), "Zellwand", "Zytoplasma",
                             "Extrazellularraum", "Zellmembran")), col = blues9[2])


sa_dns_diff <- matrix(round(c(mean(na.omit(sa_nmp_sgd_all_seq-sa_nmp_dns_all_seq)),
                              mean(na.omit(sa_spt_sgd_all_seq-sa_nmp_dns_all_seq))), digits = 4),
                      nrow = 2)
colnames(sa_dns_diff) <- c("Differenz Epitopdichten")
rownames(sa_dns_diff) <- c("NetMHCpan", "SYFPEITHI")
sa_dns_diff_table <- tableGrob(sa_dns_diff, theme=table_theme1)
grid.arrange(sa_dns_diff_table)




################################# Analysis ####################################

# Correlation Length - Mass
sa_lenmas_cor <- cor.test(sa_len_all_seq, sa_mas_all_seq, method = "spearman")["estimate"]
sa_lenmas_pvl <- cor.test(sa_len_all_seq, sa_mas_all_seq, method = "spearman", exact = FALSE)["p.value"]

# Correlation Length - Signal
sa_lensig_cor <- cor.test(sa_len_all_seq, sa_sig_all_seq_na, method = "spearman")["estimate"]
sa_lensig_pvl <- cor.test(sa_len_all_seq, sa_sig_all_seq_na, method = "spearman", exact = FALSE)["p.value"]

# Correlation Length - Epi count
sa_lencnt_nmp_cor <- cor.test(sa_len_all_seq, sa_nmp_cnt_all_seq, method = "spearman")["estimate"]
sa_lencnt_nmp_pvl <- cor.test(sa_len_all_seq, sa_nmp_cnt_all_seq, method = "spearman", exact = FALSE)["p.value"]
sa_lencnt_spt_cor <- cor.test(sa_len_all_seq, sa_spt_cnt_all_seq, method = "spearman")["estimate"]
sa_lencnt_spt_pvl <- cor.test(sa_len_all_seq, sa_spt_cnt_all_seq, method = "spearman", exact = FALSE)["p.value"]

# Correlation Length - Density
sa_lendns_nmp_cor <- cor.test(sa_len_all_seq, sa_nmp_dns_all_seq, method = "spearman")["estimate"]
sa_lendns_nmp_pvl <- cor.test(sa_len_all_seq, sa_nmp_dns_all_seq, method = "spearman", exact = FALSE)["p.value"]
sa_lendns_spt_cor <- cor.test(sa_len_all_seq, sa_spt_dns_all_seq, method = "spearman")["estimate"]
sa_lendns_spt_pvl <- cor.test(sa_len_all_seq, sa_spt_dns_all_seq, method = "spearman", exact = FALSE)["p.value"]

# Correlation Length - Signal Density
sa_lensgd_nmp_cor <- cor.test(sa_len_all_seq, sa_nmp_sgd_all_seq, method = "spearman")["estimate"]
sa_lensgd_nmp_pvl <- cor.test(sa_len_all_seq, sa_nmp_sgd_all_seq, method = "spearman", exact = FALSE)["p.value"]
sa_lensgd_spt_cor <- cor.test(sa_len_all_seq, sa_spt_sgd_all_seq, method = "spearman")["estimate"]
sa_lensgd_spt_pvl <- cor.test(sa_len_all_seq, sa_spt_sgd_all_seq, method = "spearman", exact = FALSE)["p.value"]

# Correlation Epi Count - Density
sa_cntdns_nmp_cor <- cor.test(sa_nmp_cnt_all_seq, sa_nmp_dns_all_seq, method = "spearman")["estimate"]
sa_cntdns_nmp_pvl <- cor.test(sa_nmp_cnt_all_seq, sa_nmp_dns_all_seq, method = "spearman", exact = FALSE)["p.value"]
sa_cntdns_spt_cor <- cor.test(sa_spt_cnt_all_seq, sa_spt_dns_all_seq, method = "spearman")["estimate"]
sa_cntdns_spt_pvl <- cor.test(sa_spt_cnt_all_seq, sa_spt_dns_all_seq, method = "spearman", exact = FALSE)["p.value"]

# Correlation Epi Count - Signal Density
sa_cntsgd_nmp_cor <- cor.test(sa_nmp_cnt_all_seq, sa_nmp_sgd_all_seq, method = "spearman")["estimate"]
sa_cntsgd_nmp_pvl <- cor.test(sa_nmp_cnt_all_seq, sa_nmp_sgd_all_seq, method = "spearman", exact = FALSE)["p.value"]
sa_cntsgd_spt_cor <- cor.test(sa_spt_cnt_all_seq, sa_spt_sgd_all_seq, method = "spearman")["estimate"]
sa_cntsgd_spt_pvl <- cor.test(sa_spt_cnt_all_seq, sa_spt_sgd_all_seq, method = "spearman", exact = FALSE)["p.value"]

# Correlation Density - Signal Density
sa_dnssgd_nmp_cor <- cor.test(sa_nmp_dns_all_seq, sa_nmp_sgd_all_seq, method = "spearman")["estimate"]
sa_dnssgd_nmp_pvl <- cor.test(sa_nmp_dns_all_seq, sa_nmp_sgd_all_seq, method = "spearman", exact = FALSE)["p.value"]
sa_dnssgd_spt_cor <- cor.test(sa_spt_dns_all_seq, sa_spt_sgd_all_seq, method = "spearman")["estimate"]
sa_dnssgd_spt_pvl <- cor.test(sa_spt_dns_all_seq, sa_spt_sgd_all_seq, method = "spearman", exact = FALSE)["p.value"]

# Correlation Signal Lengt - Signal Density
sa_sigsgd_nmp_cor <- cor.test(sa_sig_all_seq_na, sa_nmp_sgd_all_seq, method = "spearman")["estimate"]
sa_sigsgd_nmp_pvl <- cor.test(sa_sig_all_seq_na, sa_nmp_sgd_all_seq, method = "spearman", exact = FALSE)["p.value"]
sa_sigsgd_spt_cor <- cor.test(sa_sig_all_seq_na, sa_spt_sgd_all_seq, method = "spearman")["estimate"]
sa_sigsgd_spt_pvl <- cor.test(sa_sig_all_seq_na, sa_spt_sgd_all_seq, method = "spearman", exact = FALSE)["p.value"]

# Plot histograms to see normal distribution
par(mfrow=c(3,3))
hist(sa_len_all_seq, main = "", breaks = 100)
hist(sa_mas_all_seq, main = "", breaks = 100)
hist(sa_sig_all_seq_na, main = "", breaks = 100)
hist(sa_nmp_cnt_all_seq, main = "", breaks = 100)
hist(sa_nmp_dns_all_seq, main = "", breaks = 100)
hist(sa_nmp_sgd_all_seq, main = "", breaks = 100)
hist(sa_spt_dns_all_seq, main = "", breaks = 100)
hist(sa_spt_cnt_all_seq, main = "", breaks = 100)
hist(sa_spt_sgd_all_seq, main = "", breaks = 100)
par(mfrow=c(1,1))


# Correlation Type Plots
par(mfrow=c(2,2))
plot(sa_len_all_seq, sa_mas_all_seq, xlab = "Length", ylab = "Mass", 
     main = expression(paste(italic("S. aureus"), " Sequences")))
plot(sa_len_all_seq, sa_sig_all_seq_na, xlab = "Length", ylab = "Signal", 
     main = expression(paste(italic("S. aureus"), " Sequences")))
plot(sa_len_all_seq, sa_nmp_cnt_all_seq, col="blue", ylim = c(0,60), xlab="", ylab="")
par(new=TRUE)
plot(sa_len_all_seq, sa_spt_cnt_all_seq, col="red", ylim = c(0,60), xlab="Sequence Length",
     ylab="Epitopes", main = expression(paste(italic("S. aureus"), " Sequences")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)
plot(sa_nmp_dns_all_seq, sa_nmp_sgd_all_seq, col="blue", ylim = c(0,0.4), xlim = c(0,0.1),
     xlab="", ylab="",main = expression(paste(italic("S. aureus"), " Sequences")))
par(new=TRUE)
plot(sa_spt_dns_all_seq, sa_spt_sgd_all_seq, xlab = "Density", ylab = "Signal Density",
     ylim = c(0,0.4), xlim = c(0,0.1), col = "red")
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch=1, col = c("blue", "red"))
par(new=FALSE)
par(mfrow=c(1,1))

# Plot Sequence Correlations
sa_seq_cor_overview <- matrix(c(rep("",2),rep(c("N","S"),7),
                                round(as.numeric(c(sa_lenmas_cor, sa_lensig_cor,
                                                 sa_lencnt_nmp_cor, sa_lencnt_spt_cor,
                                                 sa_lendns_nmp_cor, sa_lendns_spt_cor,
                                                 sa_lensgd_nmp_cor, sa_lensgd_spt_cor,
                                                 sa_sigsgd_nmp_cor, sa_sigsgd_spt_cor,
                                                 sa_cntdns_nmp_cor, sa_cntdns_spt_cor,
                                                 sa_cntsgd_nmp_cor, sa_cntsgd_spt_cor,
                                                 sa_dnssgd_nmp_cor, sa_dnssgd_spt_cor)),
                                      digits = 4),
                                format(c("<0.05", "<0.05",
                                         "<0.05", "<0.05",
                                         sa_lendns_nmp_pvl, 
                                         sa_lendns_spt_pvl,
                                         "<0.05", "<0.05",
                                         "<0.05", "<0.05",
                                         "<0.05", "<0.05",
                                         sa_cntsgd_nmp_pvl, "<0.05",
                                         "<0.05", "<0.05"), digits = 2)),
                              nrow = 16, ncol = 3)
rownames(sa_seq_cor_overview) <- c("L�nge ~ Masse", "L�nge ~ Signall�nge",
                                   "L�nge ~ Anz. Epitope",
                                   "L�nge ~ Anz. Epitope",
                                   "L�nge ~ Epitopdichte",
                                   "L�nge ~ Epitopdichte",
                                   "L�nge ~ Signal-Epitopdichte",
                                   "L�nge ~ Signal-Epitopdichte",
                                   "Signall�nge ~ Signal-Epitopdichte",
                                   "Signall�nge ~ Signal-Epitopdichte",
                                   "Anz. Epitope ~ Epitopdichte",
                                   "Anz. Epitope ~ Epitopdichte",
                                   "Anz. Epitope ~ Signal-Epitopdichte",
                                   "Anz. Epitope ~ Signal-Epitopdichte",
                                   "Epitopdichte ~ Signal-Epitopdichte",
                                   "Epitopdichte ~ Signal-Epitopdichte")
colnames(sa_seq_cor_overview) <- c("","Korrelation [-1,1]", "P-Wert")
sa_seq_cor_title1 <- textGrob("Sequenzkorrelationen")
sa_seq_cor_title2 <- textGrob(expression(paste(italic("S. aureus"),", alle pwv. Sequenzen")))
sa_seq_cor_table <- tableGrob(sa_seq_cor_overview, theme = table_theme1)

l <- sa_seq_cor_table$layout
for (row in 2:17) {
        sa_seq_cor_table$grobs[ which(l$t==row & l$l==2 & l$name=="core-bg")][[1]][["gp"]] <- gpar(fill="white", col="white")
        sa_seq_cor_table$grobs[ which(l$t==row & l$l==2 & l$name=="core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
}
grid.arrange(sa_seq_cor_title1, sa_seq_cor_title2, sa_seq_cor_table, heights = unit(c(5,20,120), rep("mm", 3)))

# Location Variations
sa_lenloc <- lm(as.numeric(unlist(sa_uniProt_alllocs$LENGTH))~sa_uniProt_alllocs$LOCATION)
sa_masloc <- lm(as.numeric(unlist(sa_uniProt_alllocs$MASS))~sa_uniProt_alllocs$LOCATION)
sa_sigloc <- lm(as.numeric(unlist(sa_uniProt_alllocs$SIGNAL))~sa_uniProt_alllocs$LOCATION)
sa_cntloc_nmp <- lm(sa_nmp_cnt_alllocs~sa_nmp_alllocs$LOCATION)
sa_dnsloc_nmp <- lm(as.numeric(unlist(sa_nmp_alllocs$DENSITY))~sa_nmp_alllocs$LOCATION)
sa_sgdloc_nmp <- lm(as.numeric(unlist(sa_nmp_alllocs$`SIG DENSITY`))~sa_nmp_alllocs$LOCATION)
sa_cntloc_spt <- lm(sa_spt_cnt_alllocs~sa_spt_alllocs$LOCATION)
sa_dnsloc_spt <- lm(as.numeric(unlist(sa_spt_alllocs$DENSITY))~sa_spt_alllocs$LOCATION)
sa_sgdloc_spt <- lm(as.numeric(unlist(sa_spt_alllocs$`SIG DENSITY`))~sa_spt_alllocs$LOCATION)

par(mfrow=c(3,3))
hist(rstandard(sa_lenloc), breaks = 50)
hist(rstandard(sa_masloc), breaks = 50)
hist(rstandard(sa_sigloc), breaks = 50)
hist(rstandard(sa_cntloc_nmp), breaks = 50)
hist(rstandard(sa_dnsloc_nmp), breaks = 50)
hist(rstandard(sa_sgdloc_nmp), breaks = 50)
hist(rstandard(sa_cntloc_spt), breaks = 50)
hist(rstandard(sa_dnsloc_spt), breaks = 50)
hist(rstandard(sa_sgdloc_spt), breaks = 50)
par(mfrow=c(1,1))

sa_lenloc_fvl <- unlist(summary(aov(sa_lenloc)))["F value1"]
sa_lenloc_pvl <- unlist(summary(aov(sa_lenloc)))["Pr(>F)1"]
sa_masloc_fvl <- unlist(summary(aov(sa_masloc)))["F value1"]
sa_masloc_pvl <- unlist(summary(aov(sa_masloc)))["Pr(>F)1"]
sa_sigloc_fvl <- unlist(summary(aov(sa_sigloc)))["F value1"]
sa_sigloc_pvl <- unlist(summary(aov(sa_sigloc)))["Pr(>F)1"]

sa_cntloc_nmp_fvl <- unlist(summary(aov(sa_cntloc_nmp)))["F value1"]
sa_cntloc_nmp_pvl <- unlist(summary(aov(sa_cntloc_nmp)))["Pr(>F)1"]
sa_dnsloc_nmp_fvl <- unlist(summary(aov(sa_dnsloc_nmp)))["F value1"]
sa_dnsloc_nmp_pvl <- unlist(summary(aov(sa_dnsloc_nmp)))["Pr(>F)1"]
sa_sgdloc_nmp_fvl <- unlist(summary(aov(sa_sgdloc_nmp)))["F value1"]
sa_sgdloc_nmp_pvl <- unlist(summary(aov(sa_sgdloc_nmp)))["Pr(>F)1"]

sa_cntloc_spt_fvl <- unlist(summary(aov(sa_cntloc_spt)))["F value1"]
sa_cntloc_spt_pvl <- unlist(summary(aov(sa_cntloc_spt)))["Pr(>F)1"]
sa_dnsloc_spt_fvl <- unlist(summary(aov(sa_dnsloc_spt)))["F value1"]
sa_dnsloc_spt_pvl <- unlist(summary(aov(sa_dnsloc_spt)))["Pr(>F)1"]
sa_sgdloc_spt_fvl <- unlist(summary(aov(sa_sgdloc_spt)))["F value1"]
sa_sgdloc_spt_pvl <- unlist(summary(aov(sa_sgdloc_spt)))["Pr(>F)1"]

sa_loc_cor_overview <- matrix(c("","","",rep(c("N","S"),3),
                                round(as.numeric(c(sa_lenloc_fvl,
                                                 sa_masloc_fvl,
                                                 sa_sigloc_fvl,
                                                 sa_cntloc_nmp_fvl,
                                                 sa_cntloc_spt_fvl,
                                                 sa_dnsloc_nmp_fvl,
                                                 sa_dnsloc_spt_fvl,
                                                 sa_sgdloc_nmp_fvl,
                                                 sa_sgdloc_spt_fvl)),
                                    digits = 1),
                              rep("<0.05",9)),
                              nrow = 9, ncol = 3)
rownames(sa_loc_cor_overview) <- c("L�nge","Masse","Signall�nge","Anz. Epitope",
                                   "Anz. Epitope","Epitopdichte","Epitopdichte",
                                   "Signal-Epitopdichte","Signal-Epitopdichte")
colnames(sa_loc_cor_overview) <- c("","F-Wert", "P-Wert")
sa_loc_cor_title1 <- textGrob("Varianzanalyse der Sequenz- und Epitopdaten unter den Lokationen")
sa_loc_cor_title2 <- textGrob(expression(italic("S. aureus")))
sa_loc_cor_table <- tableGrob(sa_loc_cor_overview, theme = table_theme1)

l <- sa_loc_cor_table$layout
for (row in 2:10) {
        sa_loc_cor_table$grobs[ which(l$t==row & l$l==2 & l$name=="core-bg")][[1]][["gp"]] <- gpar(fill="white", col="white")
        sa_loc_cor_table$grobs[ which(l$t==row & l$l==2 & l$name=="core-fg")][[1]][["gp"]] <- gpar(fontface="bold")
}
grid.arrange(sa_loc_cor_title1, sa_loc_cor_title2, sa_loc_cor_table, heights = unit(c(5,20,80), rep("mm", 3)))

par(mar = c(5.1, 6.1, 4.1, 1.1), mfrow=c(3,3))
tuk_plot(TukeyHSD(aov(as.numeric(unlist(sa_uniProt_alllocs$LENGTH))~sa_uniProt_alllocs$LOCATION)), 
         "L�nge", "", rev(c("Zlw-Zmb", "Ztp-Zmb","Exz-Zmb","Ztp-Zlw","Exz-Zlw","Exz-Ztp")), 1,1)
tuk_plot(TukeyHSD(aov(as.numeric(unlist(sa_uniProt_alllocs$MASS))~sa_uniProt_alllocs$LOCATION)),
         "Masse", "", rev(c("Zlw-Zmb", "Ztp-Zmb","Exz-Zmb","Ztp-Zlw","Exz-Zlw","Exz-Ztp")), 1,1)
tuk_plot(TukeyHSD(aov(as.numeric(unlist(sa_uniProt_alllocs$SIGNAL))~sa_uniProt_alllocs$LOCATION)),
         "Signall�nge", "", rev(c("Zlw-Zmb","Exz-Zmb","Exz-Zlw")), 1,1)
tuk_plot(TukeyHSD(aov(sa_nmp_cnt_alllocs~sa_nmp_alllocs$LOCATION)),
         "Anz. Epitope, NetMHCpan", "", rev(c("Zlw-Zmb", "Ztp-Zmb","Exz-Zmb","Ztp-Zlw","Exz-Zlw","Exz-Ztp")), 1,1)
tuk_plot(TukeyHSD(aov(as.numeric(unlist(sa_nmp_alllocs$DENSITY))~sa_nmp_alllocs$LOCATION)),
         "Epitopdichte, NetMHCpan", "", rev(c("Zlw-Zmb", "Ztp-Zmb","Exz-Zmb","Ztp-Zlw","Exz-Zlw","Exz-Ztp")), 1,1)
tuk_plot(TukeyHSD(aov(as.numeric(unlist(sa_nmp_alllocs$`SIG DENSITY`))~sa_nmp_alllocs$LOCATION)),
         "Signal-Epitopdichte, NetMHCpan", "", rev(c("Zlw-Zmb","Exz-Zmb","Exz-Zlw")), 1,1)
tuk_plot(TukeyHSD(aov(sa_spt_cnt_alllocs~sa_spt_alllocs$LOCATION)),
         "Anz. Epitope, SYFPEITHI", "", rev(c("Zlw-Zmb", "Ztp-Zmb","Exz-Zmb","Ztp-Zlw","Exz-Zlw","Exz-Ztp")), 1,1)
tuk_plot(TukeyHSD(aov(as.numeric(unlist(sa_spt_alllocs$DENSITY))~sa_spt_alllocs$LOCATION)),
         "Epitopdichte, SYFPEITHI", "", rev(c("Zlw-Zmb", "Ztp-Zmb","Exz-Zmb","Ztp-Zlw","Exz-Zlw","Exz-Ztp")), 1,1)
tuk_plot(TukeyHSD(aov(as.numeric(unlist(sa_spt_alllocs$`SIG DENSITY`))~sa_spt_alllocs$LOCATION)),
         "Singal-Epitopdichte, SYFPEITHI", "", rev(c("Zlw-Zmb","Exz-Zmb","Exz-Zlw")), 1,1)
par(mfrow=c(1,1))
par(mar = c(5.1, 4.1, 4.1, 2.1))