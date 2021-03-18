########################### Data collection Mu50  #############################
###############################################################################
# Before executing this file load the 'data_analysis_functions.R' and         #
# 'staphAureus_data_analysis.R' in this Environment.                          #
###############################################################################

# Read Mu50 / ATCC700699 Files
mu50_uniProt <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/mu50_uniProt.xlsx", sheet = "All")

# Length
mu50_len <- as.numeric(unlist(mu50_uniProt["LENGTH"]))

# Mass
mu50_mas <- as.numeric(gsub(",","",unlist(mu50_uniProt["MASS"], use.names = FALSE)))

# Signal regions
mu50_sig_na <- as.numeric( unlist( mu50_uniProt["SIGNAL"]))
mu50_sig <- na.omit(mu50_sig_na)

# Plot UniProt overview
mu50_UP_overview <- matrix(c(nrow(mu50_uniProt),
                               round(mean(mu50_len)), 
                               length(na.omit(mu50_sig)),
                               round(mean(na.omit(mu50_sig)))),
                             nrow = 1, ncol = 4)
colnames(mu50_UP_overview) <- c("Number of Sequences", "Length [AA]", "Number of Signals","Signal Length [AA]")
rownames(mu50_UP_overview) <- c("Mu50 / ATCC700699")
mu50_table_up <- tableGrob(mu50_UP_overview, theme = table_theme)
mu50_up_title <- textGrob("UniProt Overview")
mu50_up_subtt <- textGrob(expression(paste(italic("S.aureus"), "-Strain ", italic("Mu50"))))
grid.arrange(mu50_up_title, mu50_up_subtt, mu50_table_up, ncol=1,
             heights = unit(c(10,10,30),rep("mm",4)))




##### Epitopes #####

mu50_nmp <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/mu50_netMHCpan.xlsx", sheet = "All")
mu50_spt <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/mu50_SYFPEITHI.xlsx", sheet = "All")

# Densities
mu50_nmp_dns <- as.numeric(unlist(mu50_nmp["DENSITY"]))
mu50_spt_dns <- as.numeric(unlist(mu50_spt["DENSITY"]))

# Signal Densities
mu50_nmp_sgd <- as.numeric(unlist(mu50_nmp["SIG DENSITY"]))
mu50_spt_sgd <- as.numeric(unlist(mu50_spt["SIG DENSITY"]))

# Epitope counts
mu50_nmp_cnt <- count_epis(mu50_nmp)
mu50_spt_cnt <- count_epis(mu50_spt)

# Plot prediction overview
mu50_Pred_overview_nmp <- matrix(c(round( mean( mu50_nmp_cnt)),
                                     round( mean( mu50_nmp_dns),
                                            digits = 3),
                                     round( mean( na.omit( mu50_nmp_sgd)),
                                            digits = 3)),
                                   nrow = 1, ncol = 3)

mu50_Pred_overview_spt <- matrix(c(round( mean( mu50_spt_cnt)),
                                     round( mean( mu50_spt_dns),
                                            digits = 3),
                                     round( mean( na.omit( mu50_spt_sgd)),
                                            digits = 3)),
                                   nrow = 1, ncol = 3)
colnames(mu50_Pred_overview_nmp) <- c("Number of Epitopes", "Density", "Signal Density")
colnames(mu50_Pred_overview_spt) <- c("Number of Epitopes", "Density", "Signal Density")
rownames(mu50_Pred_overview_nmp) <- "Mu50 / ATCC700699"
rownames(mu50_Pred_overview_spt) <- "Mu50 / ATCC700699"
mu50_nmp_table <- tableGrob(mu50_Pred_overview_nmp, theme = table_theme)
mu50_spt_table <- tableGrob(mu50_Pred_overview_spt, theme = table_theme)
mu50_pred_title <- textGrob("Epitopes Overview")
mu50_nmp_subtt <- textGrob(expression(paste(italic("S.aureus"), "-Strain ",
                                            italic("Mu50"), ", NetMHCpan")))
mu50_spt_subtt <- textGrob(expression(paste(italic("S.aureus"), "-Strain ",
                                            italic("Mu50"), ", SYFPEITHI")))
grid.arrange(mu50_pred_title, mu50_nmp_subtt, mu50_nmp_table, mu50_spt_subtt, mu50_spt_table,
             ncol=1, heights = unit(c(20,1,30,1,30), rep("mm", 5)))



############################### Analysis ##################################

# Correlation Length - Mass
mu50_lenmas_cor <- cor.test(mu50_len, mu50_mas, method = "spearman")["estimate"]
mu50_lenmas_pvl <- cor.test(mu50_len, mu50_mas, method = "spearman", exact = FALSE)["p.value"]

# Correlation Length - Signal
mu50_lensig_cor <- cor.test(mu50_len, mu50_sig_na, method = "spearman")["estimate"]
mu50_lensig_pvl <- cor.test(mu50_len, mu50_sig_na, method = "spearman", exact = FALSE)["p.value"]

# Correlation Length - Epi count
mu50_lencnt_nmp_cor <- cor.test(mu50_len, mu50_nmp_cnt, method = "spearman")["estimate"]
mu50_lencnt_nmp_pvl <- cor.test(mu50_len, mu50_nmp_cnt, method = "spearman", exact = FALSE)["p.value"]
mu50_lencnt_spt_cor <- cor.test(mu50_len, mu50_spt_cnt, method = "spearman")["estimate"]
mu50_lencnt_spt_pvl <- cor.test(mu50_len, mu50_spt_cnt, method = "spearman", exact = FALSE)["p.value"]

# Correlation Length - Density
mu50_lendns_nmp_cor <- cor.test(mu50_len, mu50_nmp_dns, method = "spearman")["estimate"]
mu50_lendns_nmp_pvl <- cor.test(mu50_len, mu50_nmp_dns, method = "spearman", exact = FALSE)["p.value"]
mu50_lendns_spt_cor <- cor.test(mu50_len, mu50_spt_dns, method = "spearman")["estimate"]
mu50_lendns_spt_pvl <- cor.test(mu50_len, mu50_spt_dns, method = "spearman", exact = FALSE)["p.value"]

# Correlation Length - Signal Density
mu50_lensgd_nmp_cor <- cor.test(mu50_len, mu50_nmp_sgd, method = "spearman")["estimate"]
mu50_lensgd_nmp_pvl <- cor.test(mu50_len, mu50_nmp_sgd, method = "spearman", exact = FALSE)["p.value"]
mu50_lensgd_spt_cor <- cor.test(mu50_len, mu50_spt_sgd, method = "spearman")["estimate"]
mu50_lensgd_spt_pvl <- cor.test(mu50_len, mu50_spt_sgd, method = "spearman", exact = FALSE)["p.value"]

# Correlation Epi Count - Density
mu50_cntdns_nmp_cor <- cor.test(mu50_nmp_cnt, mu50_nmp_dns, method = "spearman")["estimate"]
mu50_cntdns_nmp_pvl <- cor.test(mu50_nmp_cnt, mu50_nmp_dns, method = "spearman", exact = FALSE)["p.value"]
mu50_cntdns_spt_cor <- cor.test(mu50_spt_cnt, mu50_spt_dns, method = "spearman")["estimate"]
mu50_cntdns_spt_pvl <- cor.test(mu50_spt_cnt, mu50_spt_dns, method = "spearman", exact = FALSE)["p.value"]

# Correlation Epi Count - Signal Density
mu50_cntsgd_nmp_cor <- cor.test(mu50_nmp_cnt, mu50_nmp_sgd, method = "spearman")["estimate"]
mu50_cntsgd_nmp_pvl <- cor.test(mu50_nmp_cnt, mu50_nmp_sgd, method = "spearman", exact = FALSE)["p.value"]
mu50_cntsgd_spt_cor <- cor.test(mu50_spt_cnt, mu50_spt_sgd, method = "spearman")["estimate"]
mu50_cntsgd_spt_pvl <- cor.test(mu50_spt_cnt, mu50_spt_sgd, method = "spearman", exact = FALSE)["p.value"]


# Correlation Density - Signal Density
mu50_dnssgd_nmp_cor <- cor.test(mu50_nmp_dns, mu50_nmp_sgd, method = "spearman")["estimate"]
mu50_dnssgd_nmp_pvl <- cor.test(mu50_nmp_dns, mu50_nmp_sgd, method = "spearman", exact = FALSE)["p.value"]
mu50_dnssgd_spt_cor <- cor.test(mu50_spt_dns, mu50_spt_sgd, method = "spearman")["estimate"]
mu50_dnssgd_spt_pvl <- cor.test(mu50_spt_dns, mu50_spt_sgd, method = "spearman", exact = FALSE)["p.value"]

# Correlation Type Plots
# Länge~Masse
par(mfrow=c(4,2))
plot(mu50_len, mu50_mas, xlab = "Länge", ylab = "Masse", 
     main = expression(italic("S.aureus, Mu50")))

# Länge~Signallänge
plot(mu50_len, mu50_sig_na, xlab = "Länge", ylab = "Signallänge", 
     main = expression(italic("S.aureus, Mu50")))

# Länge~Anz. Epitope
plot(mu50_len, mu50_nmp_cnt, col="blue", ylim = c(0,90), xlab="", ylab="")
par(new=TRUE)
plot(mu50_len, mu50_spt_cnt, col="red", ylim = c(0,90), xlab="Länge",
     ylab="Epitope", main = expression(italic("S.aureus, Mu50")))
legend("right", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

# Länge~Epitopdichte
plot(mu50_len, mu50_nmp_dns, col="blue", ylim = c(0,0.25), xlab="", ylab="")
par(new=TRUE)
plot(mu50_len, mu50_spt_dns, col="red", ylim = c(0,0.25), xlab="Länge",
     ylab="Epitopdichte", main = expression(italic("S.aureus, Mu50")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

# Länge~Signal-Epitopdichte
plot(mu50_len, mu50_nmp_sgd, col="blue", ylim = c(0,0.7), xlab="", ylab="")
par(new=TRUE)
plot(mu50_len, mu50_spt_sgd, col="red", ylim = c(0,0.7), xlab="Länge",
     ylab="Signal-Epitopdichte", main = expression(italic("S.aureus, Mu50")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

# Anz.Epitope~Epitopdichte
plot(mu50_nmp_cnt, mu50_nmp_dns, col="blue", ylim = c(0,0.25), xlab="", ylab="")
par(new=TRUE)
plot(mu50_spt_cnt, mu50_spt_dns, col="red", ylim = c(0,0.25), xlab="Länge",
     ylab="Epitopdichte", main = expression(italic("S.aureus, Mu50")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

# Anz.Epitope~Signal-Epitopdichte
plot(mu50_nmp_cnt, mu50_nmp_sgd, col="blue", ylim = c(0,0.7), xlab="", ylab="")
par(new=TRUE)
plot(mu50_spt_cnt, mu50_spt_sgd, col="red", ylim = c(0,0.7), xlab="Epitope",
     ylab="Signal-Epitopdichte", main = expression(italic("S.aureus, Mu50")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

# Epitopdichte~Signal-Epitopdichte
plot(mu50_nmp_dns, mu50_nmp_sgd, col="blue",
     ylim = c(0,0.4), xlim = c(0,0.1),xlab="", ylab="",
     main = expression(paste(italic("S.aureus, Mu50"), " Sequences")))
par(new=TRUE)
plot(mu50_spt_dns, mu50_spt_sgd, xlab = "Density", ylab = "Signal-Epitopdichte",
     ylim = c(0,0.4), xlim = c(0,0.1), col = "red")
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch=1, col = c("blue", "red"))
par(new=FALSE)

par(mfrow=c(1,1))

# Plot Sequence Correlations
mu50_seq_cor_overview <- matrix(round(as.numeric(c(mu50_lenmas_cor, mu50_lensig_cor,
                                                  mu50_lencnt_nmp_cor, mu50_lencnt_spt_cor, 
                                                  mu50_dnssgd_nmp_cor, mu50_dnssgd_spt_cor,
                                                  mu50_lenmas_pvl, mu50_lensig_pvl,
                                                  mu50_lencnt_nmp_pvl, mu50_lencnt_spt_pvl, 
                                                  mu50_dnssgd_nmp_pvl, mu50_dnssgd_spt_pvl)),
                                     digits = 4),
                               nrow = 6, ncol = 2)
rownames(mu50_seq_cor_overview) <- c("Length ~ Mass", "Length ~ Signal Length", 
                                    "Length ~ Number of Epitopes (N)",
                                    "Length ~ Number of Epitopes (S)",
                                    "Density ~ Signal Density (N)",
                                    "Density ~ Signal Density (S)")
colnames(mu50_seq_cor_overview) <- c("Correlation [-1,1]", "P-Value")
mu50_seq_cor_table <- tableGrob(mu50_seq_cor_overview, theme = table_theme1)
mu50_seq_cor_title1 <- textGrob("Sequence Correlations")
mu50_seq_cor_title2 <- textGrob(expression(italic("S.aureus, Mu50 / ATCC700699")))
grid.arrange(mu50_seq_cor_title1, mu50_seq_cor_title2, mu50_seq_cor_table, heights = unit(c(5,20,50), rep("mm", 3)))


