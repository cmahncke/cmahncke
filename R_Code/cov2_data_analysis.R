######################### Data collection SARS-CoV-2 ############################
###############################################################################
# Before executing this file load the 'data_analysis_functions.R in this      #
# Environment.                                                                #
###############################################################################

##### Proteins #####

# Read UniProt Files
cov_uniProt <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/cov2_uniProt.xlsx", sheet = "All")

# Length
cov_len <- as.numeric(unlist(cov_uniProt["LENGTH"]))

# Mass
cov_mas <- as.numeric(gsub(",","",unlist(cov_uniProt["MASS"], use.names = FALSE)))

# Signal regions
cov_sig_na <- as.numeric( unlist( cov_uniProt["SIGNAL"]))
cov_sig <- na.omit(cov_sig_na)

# Lineages
cov_lineage_na <- find_strain_na( cov_uniProt["LINEAGE"])
cov_lineage <- find_strain( cov_uniProt["LINEAGE"])

# Plot UniProt overview
cov_UP_overview <- matrix(c(nrow(cov_uniProt),
                        round(mean(cov_len)), 
                        length(na.omit(cov_sig)),
                        round(mean(na.omit(cov_sig)))),
                      nrow = 1, ncol = 4)
colnames(cov_UP_overview) <- c("Number of Sequences", "Length [AA]", "Number of Signals","Signal Length [AA]")
rownames(cov_UP_overview) <- c("SARS-CoV-2")
cov_up_title <- textGrob("UniProt Overview")
cov_up_table <- tableGrob(cov_UP_overview, theme = table_theme)
grid.arrange(cov_up_title, cov_up_table, ncol=1, heights = unit(c(10,30),rep("mm",2)))


##### Epitopes #####

# Read NetMHCpan files
cov_nmp <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/cov2_netMHCPan.xlsx", sheet = "All")
# Read SYFPEITHI files
cov_spt <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/Spreadsheets/cov2_SYFPEITHI.xlsx", sheet = "All")


# Densities
cov_nmp_dns <- as.numeric(unlist(cov_nmp["DENSITY"]))
cov_spt_dns <- as.numeric(unlist(cov_spt["DENSITY"]))

# Signal Densities
cov_nmp_sgd <- as.numeric(unlist(cov_nmp["SIG DENSITY"]))
cov_spt_sgd <- as.numeric(unlist(cov_spt["SIG DENSITY"]))

# Epitope counts
cov_nmp_cnt <- count_epis(cov_nmp)
cov_spt_cnt <- count_epis(cov_spt)

# Plot prediction overview
cov_Pred_overview_nmp <- matrix(c(round( mean( cov_nmp_cnt)),
                              round( mean(  cov_nmp_dns),
                                     digits = 3),
                              round( mean( na.omit( cov_nmp_sgd)),
                                     digits = 3)),
                            nrow = 1, ncol = 3)

cov_Pred_overview_spt <- matrix(c(round( mean( cov_spt_cnt)),
                              round( mean( cov_spt_dns),
                                     digits = 3),
                              round( mean( na.omit( cov_spt_sgd)),
                                     digits = 3)),
                              nrow = 1, ncol = 3)
colnames(cov_Pred_overview_nmp) <- c("Number of Epitopes", "Density", "Signal Density")
colnames(cov_Pred_overview_spt) <- c("Number of Epitopes", "Density", "Signal Density")
rownames(cov_Pred_overview_nmp) <- "SARS-CoV-2"
rownames(cov_Pred_overview_spt) <- "SARS-CoV-2"
cov_nmp_table <- tableGrob(cov_Pred_overview_nmp, theme = table_theme)
cov_spt_table <- tableGrob(cov_Pred_overview_spt, theme = table_theme)
cov_pred_title <- textGrob("Epitopes Overview")
cov_nmp_subtt <- textGrob("NetMHCpan")
cov_spt_subtt <- textGrob("SYFPEITHI")
grid.arrange(cov_pred_title, cov_nmp_subtt, cov_nmp_table, cov_spt_subtt, cov_spt_table,
             ncol=1, heights = unit(c(15,2,25,1,25), rep("mm", 5)))




############################### Analysis ##################################

# Correlation Length - Mass
cov_lenmas_cor <- cor.test(cov_len, cov_mas, method = "spearman")["estimate"]
cov_lenmas_pvl <- cor.test(cov_len, cov_mas, method = "spearman", exact = FALSE)["p.value"]

# Correlation Length - Signal
cov_lensig_cor <- cor.test(cov_len, cov_sig_na, method = "spearman")["estimate"]
cov_lensig_pvl <- cor.test(cov_len, cov_sig_na, method = "spearman", exact = FALSE)["p.value"]

# Correlation Length - Epi count
cov_lencnt_nmp_cor <- cor.test(cov_len, cov_nmp_cnt, method = "spearman")["estimate"]
cov_lencnt_nmp_pvl <- cor.test(cov_len, cov_nmp_cnt, method = "spearman", exact = FALSE)["p.value"]
cov_lencnt_spt_cor <- cor.test(cov_len, cov_spt_cnt, method = "spearman")["estimate"]
cov_lencnt_spt_pvl <- cor.test(cov_len, cov_spt_cnt, method = "spearman", exact = FALSE)["p.value"]

# Correlation Length - Density
cov_lendns_nmp_cor <- cor.test(cov_len, cov_nmp_dns, method = "spearman")["estimate"]
cov_lendns_nmp_pvl <- cor.test(cov_len, cov_nmp_dns, method = "spearman", exact = FALSE)["p.value"]
cov_lendns_spt_cor <- cor.test(cov_len, cov_spt_dns, method = "spearman")["estimate"]
cov_lendns_spt_pvl <- cor.test(cov_len, cov_spt_dns, method = "spearman", exact = FALSE)["p.value"]

# Correlation Length - Signal Density
cov_lensgd_nmp_cor <- cor.test(cov_len, cov_nmp_sgd, method = "spearman")["estimate"]
cov_lensgd_nmp_pvl <- cor.test(cov_len, cov_nmp_sgd, method = "spearman", exact = FALSE)["p.value"]
cov_lensgd_spt_cor <- cor.test(cov_len, cov_spt_sgd, method = "spearman")["estimate"]
cov_lensgd_spt_pvl <- cor.test(cov_len, cov_spt_sgd, method = "spearman", exact = FALSE)["p.value"]

# Correlation Epi Count - Density
cov_cntdns_nmp_cor <- cor.test(cov_nmp_cnt, cov_nmp_dns, method = "spearman")["estimate"]
cov_cntdns_nmp_pvl <- cor.test(cov_nmp_cnt, cov_nmp_dns, method = "spearman", exact = FALSE)["p.value"]
cov_cntdns_spt_cor <- cor.test(cov_spt_cnt, cov_spt_dns, method = "spearman")["estimate"]
cov_cntdns_spt_pvl <- cor.test(cov_spt_cnt, cov_spt_dns, method = "spearman", exact = FALSE)["p.value"]

# Correlation Epi Count - Signal Density
cov_cntsgd_nmp_cor <- cor.test(cov_nmp_cnt, cov_nmp_sgd, method = "spearman")["estimate"]
cov_cntsgd_nmp_pvl <- cor.test(cov_nmp_cnt, cov_nmp_sgd, method = "spearman", exact = FALSE)["p.value"]
cov_cntsgd_spt_cor <- cor.test(cov_spt_cnt, cov_spt_sgd, method = "spearman")["estimate"]
cov_cntsgd_spt_pvl <- cor.test(cov_spt_cnt, cov_spt_sgd, method = "spearman", exact = FALSE)["p.value"]


# Correlation Density - Signal Density
cov_dnssgd_nmp_cor <- cor.test(cov_nmp_dns, cov_nmp_sgd, method = "spearman")["estimate"]
cov_dnssgd_nmp_pvl <- cor.test(cov_nmp_dns, cov_nmp_sgd, method = "spearman", exact = FALSE)["p.value"]
cov_dnssgd_spt_cor <- cor.test(cov_spt_dns, cov_spt_sgd, method = "spearman")["estimate"]
cov_dnssgd_spt_pvl <- cor.test(cov_spt_dns, cov_spt_sgd, method = "spearman", exact = FALSE)["p.value"]

# Correlation Type Plots
# Länge~Masse
par(mfrow=c(4,2))
plot(cov_len, cov_mas, xlab = "Länge", ylab = "Masse", 
     main = expression(paste(italic("SARS-CoV-2"), " Sequenzen")))

# Länge~Signallänge
plot(cov_len, cov_sig_na, xlab = "Länge", ylab = "Signallänge", 
     main = expression(italic("SARS-CoV-2")))

# Länge~Anz. Epitope
plot(cov_len, cov_nmp_cnt, col="blue", ylim = c(0,90), xlab="", ylab="")
par(new=TRUE)
plot(cov_len, cov_spt_cnt, col="red", ylim = c(0,90), xlab="Länge",
     ylab="Epitope", main = expression(italic("SARS-CoV-2")))
legend("right", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

# Länge~Epitopdichte
plot(cov_len, cov_nmp_dns, col="blue", ylim = c(0,0.25), xlab="", ylab="")
par(new=TRUE)
plot(cov_len, cov_spt_dns, col="red", ylim = c(0,0.25), xlab="Länge",
     ylab="Epitopdichte", main = expression(italic("SARS-CoV-2")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

# Länge~Signal-Epitopdichte
plot(cov_len, cov_nmp_sgd, col="blue", ylim = c(0,0.7), xlab="", ylab="")
par(new=TRUE)
plot(cov_len, cov_spt_sgd, col="red", ylim = c(0,0.7), xlab="Länge",
     ylab="Signal-Epitopdichte", main = expression(italic("SARS-CoV-2")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

# Anz.Epitope~Epitopdichte
plot(cov_nmp_cnt, cov_nmp_dns, col="blue", ylim = c(0,0.25), xlab="", ylab="")
par(new=TRUE)
plot(cov_spt_cnt, cov_spt_dns, col="red", ylim = c(0,0.25), xlab="Länge",
     ylab="Epitopdichte", main = expression(italic("SARS-CoV-2")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

# Anz.Epitope~Signal-Epitopdichte
plot(cov_nmp_cnt, cov_nmp_sgd, col="blue", ylim = c(0,0.7), xlab="", ylab="")
par(new=TRUE)
plot(cov_spt_cnt, cov_spt_sgd, col="red", ylim = c(0,0.7), xlab="Epitope",
     ylab="Signal-Epitopdichte", main = expression(italic("SARS-CoV-2")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

# Epitopdichte~Signal-Epitopdichte
plot(cov_nmp_dns, cov_nmp_sgd, col="blue",
     ylim = c(0,0.4), xlim = c(0,0.1),xlab="", ylab="",
     main = expression(paste(italic("SARS-CoV-2"), " Sequences")))
par(new=TRUE)
plot(cov_spt_dns, cov_spt_sgd, xlab = "Density", ylab = "Signal-Epitopdichte",
     ylim = c(0,0.4), xlim = c(0,0.1), col = "red")
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch=1, col = c("blue", "red"))
par(new=FALSE)

par(mfrow=c(1,1))


# Plot Sequence Correlations
cov_seq_cor_overview <- matrix(round(as.numeric(c(cov_lenmas_cor, cov_lensig_cor,
                                                 cov_lencnt_nmp_cor, cov_lencnt_spt_cor, 
                                                 cov_dnssgd_nmp_cor, cov_dnssgd_spt_cor,
                                                 cov_lenmas_pvl, cov_lensig_pvl,
                                                 cov_lencnt_nmp_pvl, cov_lencnt_spt_pvl, 
                                                 cov_dnssgd_nmp_pvl, cov_dnssgd_spt_pvl)),
                                    digits = 4),
                              nrow = 6, ncol = 2)
rownames(cov_seq_cor_overview) <- c("Length ~ Mass", "Length ~ Signal Length", 
                                   "Length ~ Number of Epitopes (N)",
                                   "Length ~ Number of Epitopes (S)",
                                   "Density ~ Signal Density (N)",
                                   "Density ~ Signal Density (S)")
colnames(cov_seq_cor_overview) <- c("Correlation [-1,1]", "P-Value")
cov_seq_cor_table <- tableGrob(cov_seq_cor_overview, theme = table_theme1)
cov_seq_cor_title1 <- textGrob("Sequence Correlations")
cov_seq_cor_title2 <- textGrob(expression(italic("SARS-CoV-2")))
grid.arrange(cov_seq_cor_title1, cov_seq_cor_title2, cov_seq_cor_table, heights = unit(c(5,20,50), rep("mm", 3)))

