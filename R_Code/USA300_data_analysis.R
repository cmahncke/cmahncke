########################### Data collection USA300  ###########################
###############################################################################
# Before executing this file load the 'data_analysis_functions.R' and         #
# in this Environment.                                                        #
###############################################################################

# Read USA300 Files
usa300_uniProt <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/usa300_uniProt.xlsx", sheet="All")

# Length
usa300_len <- as.numeric(unlist(usa300_uniProt["LENGTH"]))

# Mass
usa300_mas <- as.numeric(gsub(",","",unlist(usa300_uniProt["MASS"], use.names = FALSE)))

# Signal regions
usa300_sig_na <- as.numeric( unlist( usa300_uniProt["SIGNAL"]))
usa300_sig <- na.omit(usa300_sig_na)

# Plot UniProt overview
usa300_UP_overview <- matrix(c(nrow(usa300_uniProt),
                            round(mean(usa300_len)), 
                            length(na.omit(usa300_sig)),
                            round(mean(na.omit(usa300_sig)))),
                          nrow = 1, ncol = 4)
colnames(usa300_UP_overview) <- c("Number of Sequences", "Length [AA]", "Number of Signals","Signal Length [AA]")
rownames(usa300_UP_overview) <- c("USA300")
usa300_up_table <- tableGrob(usa300_UP_overview, theme = table_theme)
usa300_up_title <- textGrob("UniProt Overview")
usa300_up_subtt <- textGrob(expression(paste("Strain ", italic("USA300"))))
grid.arrange(usa300_up_title, usa300_up_subtt, usa300_up_table, ncol=1,
             heights = unit(c(10,10,30),rep("mm",3)))





##### Epitopes #####

usa300_nmp <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/usa300_netMHCpan.xlsx", sheet = "All")
usa300_spt <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/usa300_SYFPEITHI.xlsx", sheet = "All")

# Densities
usa300_nmp_dns <- as.numeric(unlist(usa300_nmp["DENSITY"]))
usa300_spt_dns <- as.numeric(unlist(usa300_spt["DENSITY"]))

# Signal Densities
usa300_nmp_sgd <- as.numeric(unlist(usa300_nmp["SIG DENSITY"]))
usa300_spt_sgd <- as.numeric(unlist(usa300_spt["SIG DENSITY"]))

# Epitope counts
usa300_nmp_cnt <- count_epis(usa300_nmp)
usa300_spt_cnt <- count_epis(usa300_spt)

# Plot prediction overview
usa300_Pred_overview_nmp <- matrix(c(round( mean( usa300_nmp_cnt)),
                                   round( mean( as.numeric( unlist( usa300_nmp_dns))),
                                          digits = 3),
                                   round( mean( na.omit( as.numeric( unlist( usa300_nmp_sgd)))),
                                          digits = 3)),
                                 nrow = 1, ncol = 3)

usa300_Pred_overview_spt <- matrix(c(round( mean( usa300_spt_cnt)),
                                   round( mean( as.numeric( unlist( usa300_spt_dns))),
                                          digits = 3),
                                   round( mean( na.omit( as.numeric( unlist( usa300_spt_sgd)))),
                                          digits = 3)),
                                 nrow = 1, ncol = 3)
colnames(usa300_Pred_overview_nmp) <- c("Number of Epitopes", "Density", "Signal Density")
colnames(usa300_Pred_overview_spt) <- c("Number of Epitopes", "Density", "Signal Density")
rownames(usa300_Pred_overview_nmp) <- "USA300"
rownames(usa300_Pred_overview_spt) <- "USA300"
usa300_nmp_table <- tableGrob(usa300_Pred_overview_nmp, theme = table_theme)
usa300_spt_table <- tableGrob(usa300_Pred_overview_spt, theme = table_theme)
usa300_pred_title <- textGrob("Epitopes Overview")
usa300_nmp_subtt <- textGrob(expression(paste(italic("S.aureus"), "-Strain ",
                                            italic("USA300"), ", NetMHCpan")))
usa300_spt_subtt <- textGrob(expression(paste(italic("S.aureus"), "-Strain ",
                                            italic("USA300"), ", SYFPEITHI")))
grid.arrange(usa300_pred_title, usa300_nmp_subtt, usa300_nmp_table, usa300_spt_subtt, usa300_spt_table,
             ncol=1, heights = unit(c(20,1,30,1,30), rep("mm", 5)))

############################### Analysis ##################################

# Correlation Length - Mass
usa300_lenmas_cor <- cor.test(usa300_len, usa300_mas, method = "spearman")["estimate"]
usa300_lenmas_pvl <- cor.test(usa300_len, usa300_mas, method = "spearman", exact = FALSE)["p.value"]

# Correlation Length - Signal
usa300_lensig_cor <- cor.test(usa300_len, usa300_sig_na, method = "spearman")["estimate"]
usa300_lensig_pvl <- cor.test(usa300_len, usa300_sig_na, method = "spearman", exact = FALSE)["p.value"]

# Correlation Length - Epi count
usa300_lencnt_nmp_cor <- cor.test(usa300_len, usa300_nmp_cnt, method = "spearman")["estimate"]
usa300_lencnt_nmp_pvl <- cor.test(usa300_len, usa300_nmp_cnt, method = "spearman", exact = FALSE)["p.value"]
usa300_lencnt_spt_cor <- cor.test(usa300_len, usa300_spt_cnt, method = "spearman")["estimate"]
usa300_lencnt_spt_pvl <- cor.test(usa300_len, usa300_spt_cnt, method = "spearman", exact = FALSE)["p.value"]

# Correlation Length - Density
usa300_lendns_nmp_cor <- cor.test(usa300_len, usa300_nmp_dns, method = "spearman")["estimate"]
usa300_lendns_nmp_pvl <- cor.test(usa300_len, usa300_nmp_dns, method = "spearman", exact = FALSE)["p.value"]
usa300_lendns_spt_cor <- cor.test(usa300_len, usa300_spt_dns, method = "spearman")["estimate"]
usa300_lendns_spt_pvl <- cor.test(usa300_len, usa300_spt_dns, method = "spearman", exact = FALSE)["p.value"]

# Correlation Length - Signal Density
usa300_lensgd_nmp_cor <- cor.test(usa300_len, usa300_nmp_sgd, method = "spearman")["estimate"]
usa300_lensgd_nmp_pvl <- cor.test(usa300_len, usa300_nmp_sgd, method = "spearman", exact = FALSE)["p.value"]
usa300_lensgd_spt_cor <- cor.test(usa300_len, usa300_spt_sgd, method = "spearman")["estimate"]
usa300_lensgd_spt_pvl <- cor.test(usa300_len, usa300_spt_sgd, method = "spearman", exact = FALSE)["p.value"]

# Correlation Epi Count - Density
usa300_cntdns_nmp_cor <- cor.test(usa300_nmp_cnt, usa300_nmp_dns, method = "spearman")["estimate"]
usa300_cntdns_nmp_pvl <- cor.test(usa300_nmp_cnt, usa300_nmp_dns, method = "spearman", exact = FALSE)["p.value"]
usa300_cntdns_spt_cor <- cor.test(usa300_spt_cnt, usa300_spt_dns, method = "spearman")["estimate"]
usa300_cntdns_spt_pvl <- cor.test(usa300_spt_cnt, usa300_spt_dns, method = "spearman", exact = FALSE)["p.value"]

# Correlation Epi Count - Signal Density
usa300_cntsgd_nmp_cor <- cor.test(usa300_nmp_cnt, usa300_nmp_sgd, method = "spearman")["estimate"]
usa300_cntsgd_nmp_pvl <- cor.test(usa300_nmp_cnt, usa300_nmp_sgd, method = "spearman", exact = FALSE)["p.value"]
usa300_cntsgd_spt_cor <- cor.test(usa300_spt_cnt, usa300_spt_sgd, method = "spearman")["estimate"]
usa300_cntsgd_spt_pvl <- cor.test(usa300_spt_cnt, usa300_spt_sgd, method = "spearman", exact = FALSE)["p.value"]


# Correlation Density - Signal Density
usa300_dnssgd_nmp_cor <- cor.test(usa300_nmp_dns, usa300_nmp_sgd, method = "spearman")["estimate"]
usa300_dnssgd_nmp_pvl <- cor.test(usa300_nmp_dns, usa300_nmp_sgd, method = "spearman", exact = FALSE)["p.value"]
usa300_dnssgd_spt_cor <- cor.test(usa300_spt_dns, usa300_spt_sgd, method = "spearman")["estimate"]
usa300_dnssgd_spt_pvl <- cor.test(usa300_spt_dns, usa300_spt_sgd, method = "spearman", exact = FALSE)["p.value"]

# Correlation Type Plots
# Länge~Masse
par(mfrow=c(4,2))
plot(usa300_len, usa300_mas, xlab = "Länge", ylab = "Masse", 
     main = expression(italic("S.aureus, USA300")))

# Länge~Signallänge
plot(usa300_len, usa300_sig_na, xlab = "Länge", ylab = "Signallänge", 
     main = expression(italic("S.aureus, USA300")))

# Länge~Anz. Epitope
plot(usa300_len, usa300_nmp_cnt, col="blue", ylim = c(0,90), xlab="", ylab="")
par(new=TRUE)
plot(usa300_len, usa300_spt_cnt, col="red", ylim = c(0,90), xlab="Länge",
     ylab="Epitope", main = expression(italic("S.aureus, USA300")))
legend("right", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

# Länge~Epitopdichte
plot(usa300_len, usa300_nmp_dns, col="blue", ylim = c(0,0.25), xlab="", ylab="")
par(new=TRUE)
plot(usa300_len, usa300_spt_dns, col="red", ylim = c(0,0.25), xlab="Länge",
     ylab="Epitopdichte", main = expression(italic("S.aureus, USA300")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

# Länge~Signal-Epitopdichte
plot(usa300_len, usa300_nmp_sgd, col="blue", ylim = c(0,0.7), xlab="", ylab="")
par(new=TRUE)
plot(usa300_len, usa300_spt_sgd, col="red", ylim = c(0,0.7), xlab="Länge",
     ylab="Signal-Epitopdichte", main = expression(italic("S.aureus, USA300")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

# Anz.Epitope~Epitopdichte
plot(usa300_nmp_cnt, usa300_nmp_dns, col="blue", ylim = c(0,0.25), xlab="", ylab="")
par(new=TRUE)
plot(usa300_spt_cnt, usa300_spt_dns, col="red", ylim = c(0,0.25), xlab="Länge",
     ylab="Epitopdichte", main = expression(italic("S.aureus, USA300")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

# Anz.Epitope~Signal-Epitopdichte
plot(usa300_nmp_cnt, usa300_nmp_sgd, col="blue", ylim = c(0,0.7), xlab="", ylab="")
par(new=TRUE)
plot(usa300_spt_cnt, usa300_spt_sgd, col="red", ylim = c(0,0.7), xlab="Epitope",
     ylab="Signal-Epitopdichte", main = expression(italic("S.aureus, USA300")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

# Epitopdichte~Signal-Epitopdichte
plot(usa300_nmp_dns, usa300_nmp_sgd, col="blue",
     ylim = c(0,0.4), xlim = c(0,0.1),xlab="", ylab="",
     main = expression(paste(italic("S.aureus, USA300"), " Sequences")))
par(new=TRUE)
plot(usa300_spt_dns, usa300_spt_sgd, xlab = "Density", ylab = "Signal-Epitopdichte",
     ylim = c(0,0.4), xlim = c(0,0.1), col = "red")
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch=1, col = c("blue", "red"))
par(new=FALSE)

par(mfrow=c(1,1))

# Plot Sequence Correlations
usa300_seq_cor_overview <- matrix(round(as.numeric(c(usa300_lenmas_cor, usa300_lensig_cor,
                                                   usa300_lencnt_nmp_cor, usa300_lencnt_spt_cor, 
                                                   usa300_dnssgd_nmp_cor, usa300_dnssgd_spt_cor,
                                                   usa300_lenmas_pvl, usa300_lensig_pvl,
                                                   usa300_lencnt_nmp_pvl, usa300_lencnt_spt_pvl, 
                                                   usa300_dnssgd_nmp_pvl, usa300_dnssgd_spt_pvl)),
                                      digits = 4),
                                nrow = 6, ncol = 2)
rownames(usa300_seq_cor_overview) <- c("Length ~ Mass", "Length ~ Signal Length", 
                                     "Length ~ Number of Epitopes (N)",
                                     "Length ~ Number of Epitopes (S)",
                                     "Density ~ Signal Density (N)",
                                     "Density ~ Signal Density (S)")
colnames(usa300_seq_cor_overview) <- c("Correlation [-1,1]", "P-Value")
usa300_seq_cor_table <- tableGrob(usa300_seq_cor_overview, theme = table_theme1)
usa300_seq_cor_title1 <- textGrob("Sequence Correlations")
usa300_seq_cor_title2 <- textGrob(expression(italic("S.aureus, USA300")))
grid.arrange(usa300_seq_cor_title1, usa300_seq_cor_title2, usa300_seq_cor_table, heights = unit(c(5,20,50), rep("mm", 3)))
