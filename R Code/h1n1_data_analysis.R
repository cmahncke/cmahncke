########################### Data collection H1N1 ##############################
###############################################################################
# Before executing this file load the 'data_analysis_functions.R in this      #
# Environment.                                                                #
###############################################################################

##### Proteins #####

# Read UniProt Files
pur_h1n1_uniProt <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/pur_uniProt.xlsx", sheet = "All")
scb_h1n1_uniProt <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/scb_uniProt.xlsx", sheet = "All")
stp_h1n1_uniProt <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/stp_uniProt.xlsx", sheet = "All")
tex_h1n1_uniProt <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/tex_uniProt.xlsx", sheet = "All")

# Length
pur_h1n1_len <- as.numeric(unlist(pur_h1n1_uniProt["LENGTH"]))
scb_h1n1_len <- as.numeric(unlist(scb_h1n1_uniProt["LENGTH"]))
stp_h1n1_len <- as.numeric(unlist(stp_h1n1_uniProt["LENGTH"]))
tex_h1n1_len <- as.numeric(unlist(tex_h1n1_uniProt["LENGTH"]))

# Mass
pur_h1n1_mas <- as.numeric(gsub(",","",unlist(pur_h1n1_uniProt["MASS"], use.names = FALSE)))
scb_h1n1_mas <- as.numeric(gsub(",","",unlist(scb_h1n1_uniProt["MASS"], use.names = FALSE)))
stp_h1n1_mas <- as.numeric(gsub(",","",unlist(stp_h1n1_uniProt["MASS"], use.names = FALSE)))
tex_h1n1_mas <- as.numeric(gsub(",","",unlist(tex_h1n1_uniProt["MASS"], use.names = FALSE)))

# Signal regions
pur_h1n1_sig_na <- as.numeric( unlist( pur_h1n1_uniProt["SIGNAL"]))
scb_h1n1_sig_na <- as.numeric( unlist( scb_h1n1_uniProt["SIGNAL"]))
stp_h1n1_sig_na <- as.numeric( unlist( stp_h1n1_uniProt["SIGNAL"]))
tex_h1n1_sig_na <- as.numeric( unlist( tex_h1n1_uniProt["SIGNAL"]))

pur_h1n1_sig <- na.omit(pur_h1n1_sig_na)
scb_h1n1_sig <- na.omit(scb_h1n1_sig_na)
stp_h1n1_sig <- na.omit(stp_h1n1_sig_na)
tex_h1n1_sig <- na.omit(tex_h1n1_sig_na)

# Plot UniProt overview
h1n1_UP_overview <- t(matrix(c(nrow(pur_h1n1_uniProt),
                               round(mean(pur_h1n1_len)),
                               length(na.omit(pur_h1n1_sig)),
                               round(mean(na.omit(pur_h1n1_sig))),
                               nrow(scb_h1n1_uniProt),
                               round(mean(scb_h1n1_len)),
                               length(na.omit(scb_h1n1_sig)),
                               round(mean(na.omit(scb_h1n1_sig))),
                               nrow(stp_h1n1_uniProt),
                               round(mean(stp_h1n1_len)),
                               length(na.omit(stp_h1n1_sig)),
                               round(mean(na.omit(stp_h1n1_sig))),
                               nrow(tex_h1n1_uniProt),
                               round(mean(tex_h1n1_len)),
                               length(na.omit(tex_h1n1_sig)),
                               round(mean(na.omit(tex_h1n1_sig)))),
                           nrow = 4, ncol = 4))
colnames(h1n1_UP_overview) <- c("Anz. Sequenzen", "Länge [AS]",
                                "Anz. Signalpeptide","Signallänge [AS]")
rownames(h1n1_UP_overview) <- c("Puerto Rico 1934", "South Cantebury 2000",
                                "St.Petersburg 2006", "Texas 2007")
h1n1_up_table <- tableGrob(h1n1_UP_overview, theme = table_theme1)
h1n1_up_title <- textGrob("UniProt Overview")
h1n1_up_subtt <- textGrob(expression(paste(italic("Influenza A"), ", H1N1-Strains")))
grid.arrange(h1n1_up_title, h1n1_up_subtt, h1n1_table_up, ncol=1,
             heights = unit(c(10,10,55),rep("mm",3)))


##### Epitopes #####

# Read NetMHCpan files
pur_h1n1_nmp <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/pur_netMHCPan.xlsx", sheet = "All")
scb_h1n1_nmp <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/scb_netMHCPan.xlsx", sheet = "All")
stp_h1n1_nmp <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/stp_netMHCPan.xlsx", sheet = "All")
tex_h1n1_nmp <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/tex_netMHCPan.xlsx", sheet = "All")
# Read SYFPEITHI files
pur_h1n1_spt <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/pur_SYFPEITHI.xlsx", sheet = "All")
scb_h1n1_spt <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/scb_SYFPEITHI.xlsx", sheet = "All")
stp_h1n1_spt <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/stp_SYFPEITHI.xlsx", sheet = "All")
tex_h1n1_spt <- read_excel("C:/Users/cedri/OneDrive/Dokumente/BA/Zwischenergebnisse/results/tex_SYFPEITHI.xlsx", sheet = "All")


# Densities
pur_h1n1_nmp_dns <- as.numeric(unlist(pur_h1n1_nmp["DENSITY"]))
scb_h1n1_nmp_dns <- as.numeric(unlist(scb_h1n1_nmp["DENSITY"]))
stp_h1n1_nmp_dns <- as.numeric(unlist(stp_h1n1_nmp["DENSITY"]))
tex_h1n1_nmp_dns <- as.numeric(unlist(tex_h1n1_nmp["DENSITY"]))

pur_h1n1_spt_dns <- as.numeric(unlist(pur_h1n1_spt["DENSITY"]))
scb_h1n1_spt_dns <- as.numeric(unlist(scb_h1n1_spt["DENSITY"]))
stp_h1n1_spt_dns <- as.numeric(unlist(stp_h1n1_spt["DENSITY"]))
tex_h1n1_spt_dns <- as.numeric(unlist(tex_h1n1_spt["DENSITY"]))

# Signal Densities
pur_h1n1_nmp_sgd <- as.numeric(unlist(pur_h1n1_nmp["SIG DENSITY"]))
scb_h1n1_nmp_sgd <- as.numeric(unlist(scb_h1n1_nmp["SIG DENSITY"]))
stp_h1n1_nmp_sgd <- as.numeric(unlist(stp_h1n1_nmp["SIG DENSITY"]))
tex_h1n1_nmp_sgd <- as.numeric(unlist(tex_h1n1_nmp["SIG DENSITY"]))

pur_h1n1_spt_sgd <- as.numeric(unlist(pur_h1n1_spt["SIG DENSITY"]))
scb_h1n1_spt_sgd <- as.numeric(unlist(scb_h1n1_spt["SIG DENSITY"]))
stp_h1n1_spt_sgd <- as.numeric(unlist(stp_h1n1_spt["SIG DENSITY"]))
tex_h1n1_spt_sgd <- as.numeric(unlist(tex_h1n1_spt["SIG DENSITY"]))

# Epitope counts
pur_h1n1_nmp_cnt <- count_epis(pur_h1n1_nmp)
scb_h1n1_nmp_cnt <- count_epis(scb_h1n1_nmp)
stp_h1n1_nmp_cnt <- count_epis(stp_h1n1_nmp)
tex_h1n1_nmp_cnt <- count_epis(tex_h1n1_nmp)

pur_h1n1_spt_cnt <- count_epis(pur_h1n1_spt)
scb_h1n1_spt_cnt <- count_epis(scb_h1n1_spt)
stp_h1n1_spt_cnt <- count_epis(stp_h1n1_spt)
tex_h1n1_spt_cnt <- count_epis(tex_h1n1_spt)

# Plot prediction overview
h1n1_Pred_overview_nmp <- matrix(c(round(c(mean(pur_h1n1_nmp_cnt),
                                           mean(scb_h1n1_nmp_cnt),
                                           mean(stp_h1n1_nmp_cnt),
                                           mean(tex_h1n1_nmp_cnt))),
                                   round(c(mean(as.numeric(unlist(pur_h1n1_nmp_dns))),
                                           mean(as.numeric(unlist(scb_h1n1_nmp_dns))),
                                           mean(as.numeric(unlist(stp_h1n1_nmp_dns))),
                                           mean(as.numeric(unlist(tex_h1n1_nmp_dns)))),
                                         digits = 3),
                                   round(c(mean(na.omit(as.numeric(unlist(pur_h1n1_nmp_sgd)))),
                                           mean(na.omit(as.numeric(unlist(scb_h1n1_nmp_sgd)))),
                                           mean(na.omit(as.numeric(unlist(stp_h1n1_nmp_sgd)))),
                                           mean(na.omit(as.numeric(unlist(tex_h1n1_nmp_sgd))))),
                                         digits = 3)),
                                 nrow = 4, ncol = 3)

h1n1_Pred_overview_spt <- matrix(c(round(c(mean(pur_h1n1_spt_cnt),
                                           mean(scb_h1n1_spt_cnt),
                                           mean(stp_h1n1_spt_cnt),
                                           mean(tex_h1n1_spt_cnt))),
                                   round(c(mean(as.numeric(unlist(pur_h1n1_spt_dns))),
                                           mean(as.numeric(unlist(scb_h1n1_spt_dns))),
                                           mean(as.numeric(unlist(stp_h1n1_spt_dns))),
                                           mean(as.numeric(unlist(tex_h1n1_spt_dns)))),
                                         digits = 3),
                                   round(c(mean(na.omit(as.numeric(unlist(pur_h1n1_spt_sgd)))),
                                           mean(na.omit(as.numeric(unlist(scb_h1n1_spt_sgd)))),
                                           mean(na.omit(as.numeric(unlist(stp_h1n1_spt_sgd)))),
                                           mean(na.omit(as.numeric(unlist(tex_h1n1_spt_sgd))))),
                                         digits = 3)),
                                 nrow = 4, ncol = 3)
                                 
colnames(h1n1_Pred_overview_nmp) <- c("Number of Epitopes", "Density", "Signal Density")
colnames(h1n1_Pred_overview_spt) <- c("Number of Epitopes", "Density", "Signal Density")
rownames(h1n1_Pred_overview_nmp) <- c("Puerto Rico / 1934", "South Cantebury / 2000",
                                      "St.Petersburg / 2006", "Texas / 2007")
rownames(h1n1_Pred_overview_spt) <- c("Puerto Rico / 1934", "South Cantebury / 2000",
                                      "St.Petersburg / 2006", "Texas / 2007")
h1n1_nmp_table <- tableGrob(h1n1_Pred_overview_nmp, theme = table_theme)
h1n1_spt_table <- tableGrob(h1n1_Pred_overview_spt, theme = table_theme)
h1n1_pred_title <- textGrob("Epitopes Overview")
h1n1_nmp_subtt <- textGrob(expression(paste(italic("Influenza A - H1N1"), ", NetMHCpan")))
h1n1_spt_subtt <- textGrob(expression(paste(italic("Influenza A - H1N1"), ", SYFPEITHI")))
grid.arrange(h1n1_pred_title, h1n1_nmp_subtt, h1n1_nmp_table, h1n1_spt_subtt,
             h1n1_spt_table, ncol=1, heights = unit(c(10,10,40,10,40), rep("mm", 4)))




############################### Analysis ##################################

# Correlation Length - Mass
pur_h1n1_lenmas_cor <- cor.test(pur_h1n1_len, pur_h1n1_mas, method = "spearman")["estimate"]
pur_h1n1_lenmas_pvl <- cor.test(pur_h1n1_len, pur_h1n1_mas, method = "spearman", exact = FALSE)["p.value"]
scb_h1n1_lenmas_cor <- cor.test(scb_h1n1_len, scb_h1n1_mas, method = "spearman")["estimate"]
scb_h1n1_lenmas_pvl <- cor.test(scb_h1n1_len, scb_h1n1_mas, method = "spearman", exact = FALSE)["p.value"]
stp_h1n1_lenmas_cor <- cor.test(stp_h1n1_len, stp_h1n1_mas, method = "spearman")["estimate"]
stp_h1n1_lenmas_pvl <- cor.test(stp_h1n1_len, stp_h1n1_mas, method = "spearman", exact = FALSE)["p.value"]
tex_h1n1_lenmas_cor <- cor.test(tex_h1n1_len, tex_h1n1_mas, method = "spearman")["estimate"]
tex_h1n1_lenmas_pvl <- cor.test(tex_h1n1_len, tex_h1n1_mas, method = "spearman", exact = FALSE)["p.value"]

par(mfrow=c(2,2))
plot(pur_h1n1_len, pur_h1n1_mas, xlab = "Length", ylab = "Mass", 
     main = expression(paste(italic("H1N1, Puerto Rico/1934"), " Sequences")))
plot(scb_h1n1_len, scb_h1n1_mas, xlab = "Length", ylab = "Mass", 
     main = expression(paste(italic("H1N1, S.Canterury/2000"), " Sequences")))
plot(stp_h1n1_len, stp_h1n1_mas, xlab = "Length", ylab = "Mass", 
     main = expression(paste(italic("H1N1, St.Petersburg/2006"), " Sequences")))
plot(tex_h1n1_len, tex_h1n1_mas, xlab = "Length", ylab = "Mass", 
     main = expression(paste(italic("H1N1, Texas/2007"), " Sequences")))

# Correlation Length - Epitope Count
pur_h1n1_lencnt_nmp_cor <- cor.test(pur_h1n1_len, pur_h1n1_nmp_cnt, method = "spearman")["estimate"]
pur_h1n1_lencnt_nmp_pvl <- cor.test(pur_h1n1_len, pur_h1n1_nmp_cnt, method = "spearman", exact = FALSE)["p.value"]
scb_h1n1_lencnt_nmp_cor <- cor.test(scb_h1n1_len, scb_h1n1_nmp_cnt, method = "spearman")["estimate"]
scb_h1n1_lencnt_nmp_pvl <- cor.test(scb_h1n1_len, scb_h1n1_nmp_cnt, method = "spearman", exact = FALSE)["p.value"]
stp_h1n1_lencnt_nmp_cor <- cor.test(stp_h1n1_len, stp_h1n1_nmp_cnt, method = "spearman")["estimate"]
stp_h1n1_lencnt_nmp_pvl <- cor.test(stp_h1n1_len, stp_h1n1_nmp_cnt, method = "spearman", exact = FALSE)["p.value"]
tex_h1n1_lencnt_nmp_cor <- cor.test(tex_h1n1_len, tex_h1n1_nmp_cnt, method = "spearman")["estimate"]
tex_h1n1_lencnt_nmp_pvl <- cor.test(tex_h1n1_len, tex_h1n1_nmp_cnt, method = "spearman", exact = FALSE)["p.value"]

pur_h1n1_lencnt_spt_cor <- cor.test(pur_h1n1_len, pur_h1n1_spt_cnt, method = "spearman")["estimate"]
pur_h1n1_lencnt_spt_pvl <- cor.test(pur_h1n1_len, pur_h1n1_spt_cnt, method = "spearman", exact = FALSE)["p.value"]
scb_h1n1_lencnt_spt_cor <- cor.test(scb_h1n1_len, scb_h1n1_spt_cnt, method = "spearman")["estimate"]
scb_h1n1_lencnt_spt_pvl <- cor.test(scb_h1n1_len, scb_h1n1_spt_cnt, method = "spearman", exact = FALSE)["p.value"]
stp_h1n1_lencnt_spt_cor <- cor.test(stp_h1n1_len, stp_h1n1_spt_cnt, method = "spearman")["estimate"]
stp_h1n1_lencnt_spt_pvl <- cor.test(stp_h1n1_len, stp_h1n1_spt_cnt, method = "spearman", exact = FALSE)["p.value"]
tex_h1n1_lencnt_spt_cor <- cor.test(tex_h1n1_len, tex_h1n1_spt_cnt, method = "spearman")["estimate"]
tex_h1n1_lencnt_spt_pvl <- cor.test(tex_h1n1_len, tex_h1n1_spt_cnt, method = "spearman", exact = FALSE)["p.value"]

par(mfrow=c(2,2))
plot(pur_h1n1_len, pur_h1n1_nmp_cnt, col="blue", ylim = c(0,60), xlab="", ylab="")
par(new=TRUE)
plot(pur_h1n1_len, pur_h1n1_spt_cnt, col="red", ylim = c(0,60), xlab="Sequence Length",
     ylab="Epitopes", main = expression(italic("H1N1, Puerto Rico/1934")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

plot(scb_h1n1_len, scb_h1n1_nmp_cnt, col="blue", ylim = c(0,60), xlab="", ylab="")
par(new=TRUE)
plot(scb_h1n1_len, scb_h1n1_spt_cnt, col="red", ylim = c(0,60), xlab="Sequence Length",
     ylab="Epitopes", main = expression(italic("H1N1, S.Canterbury/2000")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

plot(stp_h1n1_len, stp_h1n1_nmp_cnt, col="blue", ylim = c(0,60), xlab="", ylab="")
par(new=TRUE)
plot(stp_h1n1_len, stp_h1n1_spt_cnt, col="red", ylim = c(0,60), xlab="Sequence Length",
     ylab="Epitopes", main = expression(italic("H1N1, St.Petersburg/2006")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

plot(tex_h1n1_len, tex_h1n1_nmp_cnt, col="blue", ylim = c(0,60), xlab="", ylab="")
par(new=TRUE)
plot(tex_h1n1_len, tex_h1n1_spt_cnt, col="red", ylim = c(0,60), xlab="Sequence Length",
     ylab="Epitopes", main = expression(italic("H1N1, Texas/2007")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

# Correlation Length - Density
pur_h1n1_lendns_nmp_cor <- cor.test(pur_h1n1_len, pur_h1n1_nmp_dns, method = "spearman")["estimate"]
pur_h1n1_lendns_nmp_pvl <- cor.test(pur_h1n1_len, pur_h1n1_nmp_dns, method = "spearman", exact = FALSE)["p.value"]
scb_h1n1_lendns_nmp_cor <- cor.test(scb_h1n1_len, scb_h1n1_nmp_dns, method = "spearman")["estimate"]
scb_h1n1_lendns_nmp_pvl <- cor.test(scb_h1n1_len, scb_h1n1_nmp_dns, method = "spearman", exact = FALSE)["p.value"]
stp_h1n1_lendns_nmp_cor <- cor.test(stp_h1n1_len, stp_h1n1_nmp_dns, method = "spearman")["estimate"]
stp_h1n1_lendns_nmp_pvl <- cor.test(stp_h1n1_len, stp_h1n1_nmp_dns, method = "spearman", exact = FALSE)["p.value"]
tex_h1n1_lendns_nmp_cor <- cor.test(tex_h1n1_len, tex_h1n1_nmp_dns, method = "spearman")["estimate"]
tex_h1n1_lendns_nmp_pvl <- cor.test(tex_h1n1_len, tex_h1n1_nmp_dns, method = "spearman", exact = FALSE)["p.value"]

pur_h1n1_lendns_spt_cor <- cor.test(pur_h1n1_len, pur_h1n1_spt_dns, method = "spearman")["estimate"]
pur_h1n1_lendns_spt_pvl <- cor.test(pur_h1n1_len, pur_h1n1_spt_dns, method = "spearman", exact = FALSE)["p.value"]
scb_h1n1_lendns_spt_cor <- cor.test(scb_h1n1_len, scb_h1n1_spt_dns, method = "spearman")["estimate"]
scb_h1n1_lendns_spt_pvl <- cor.test(scb_h1n1_len, scb_h1n1_spt_dns, method = "spearman", exact = FALSE)["p.value"]
stp_h1n1_lendns_spt_cor <- cor.test(stp_h1n1_len, stp_h1n1_spt_dns, method = "spearman")["estimate"]
stp_h1n1_lendns_spt_pvl <- cor.test(stp_h1n1_len, stp_h1n1_spt_dns, method = "spearman", exact = FALSE)["p.value"]
tex_h1n1_lendns_spt_cor <- cor.test(tex_h1n1_len, tex_h1n1_spt_dns, method = "spearman")["estimate"]
tex_h1n1_lendns_spt_pvl <- cor.test(tex_h1n1_len, tex_h1n1_spt_dns, method = "spearman", exact = FALSE)["p.value"]

par(mfrow=c(2,2))
plot(pur_h1n1_len, pur_h1n1_nmp_dns, col="blue", ylim = c(0,0.1), xlab="", ylab="")
par(new=TRUE)
plot(pur_h1n1_len, pur_h1n1_spt_dns, col="red", ylim = c(0,0.1), xlab="Länge",
     ylab="Epitopdichte", main = expression(italic("H1N1, Puerto Rico/1934")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

plot(scb_h1n1_len, scb_h1n1_nmp_dns, col="blue", ylim = c(0,0.1), xlab="", ylab="")
par(new=TRUE)
plot(scb_h1n1_len, scb_h1n1_spt_dns, col="red", ylim = c(0,0.1), xlab="Länge",
     ylab="Epitopdichte", main = expression(italic("H1N1, S.Canterbury/2000")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

plot(stp_h1n1_len, stp_h1n1_nmp_dns, col="blue", ylim = c(0,0.1), xlab="", ylab="")
par(new=TRUE)
plot(stp_h1n1_len, stp_h1n1_spt_dns, col="red", ylim = c(0,0.1), xlab="Länge",
     ylab="Epitopdichte", main = expression(italic("H1N1, St.Petersburg/2006")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

plot(tex_h1n1_len, tex_h1n1_nmp_dns, col="blue", ylim = c(0,0.1), xlab="", ylab="")
par(new=TRUE)
plot(tex_h1n1_len, tex_h1n1_spt_dns, col="red", ylim = c(0,0.1), xlab="Länge",
     ylab="Epitopdichte", main = expression(italic("H1N1, Texas/2007")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

# Correlation Epi Count - Density
pur_h1n1_cntdns_nmp_cor <- cor.test(pur_h1n1_nmp_cnt, pur_h1n1_nmp_dns, method = "spearman")["estimate"]
pur_h1n1_cntdns_nmp_pvl <- cor.test(pur_h1n1_nmp_cnt, pur_h1n1_nmp_dns, method = "spearman", exact = FALSE)["p.value"]
scb_h1n1_cntdns_nmp_cor <- cor.test(scb_h1n1_nmp_cnt, scb_h1n1_nmp_dns, method = "spearman")["estimate"]
scb_h1n1_cntdns_nmp_pvl <- cor.test(scb_h1n1_nmp_cnt, scb_h1n1_nmp_dns, method = "spearman", exact = FALSE)["p.value"]
stp_h1n1_cntdns_nmp_cor <- cor.test(stp_h1n1_nmp_cnt, stp_h1n1_nmp_dns, method = "spearman")["estimate"]
stp_h1n1_cntdns_nmp_pvl <- cor.test(stp_h1n1_nmp_cnt, stp_h1n1_nmp_dns, method = "spearman", exact = FALSE)["p.value"]
tex_h1n1_cntdns_nmp_cor <- cor.test(tex_h1n1_nmp_cnt, tex_h1n1_nmp_dns, method = "spearman")["estimate"]
tex_h1n1_cntdns_nmp_pvl <- cor.test(tex_h1n1_nmp_cnt, tex_h1n1_nmp_dns, method = "spearman", exact = FALSE)["p.value"]

pur_h1n1_cntdns_spt_cor <- cor.test(pur_h1n1_spt_cnt, pur_h1n1_spt_dns, method = "spearman")["estimate"]
pur_h1n1_cntdns_spt_pvl <- cor.test(pur_h1n1_spt_cnt, pur_h1n1_spt_dns, method = "spearman", exact = FALSE)["p.value"]
scb_h1n1_cntdns_spt_cor <- cor.test(scb_h1n1_spt_cnt, scb_h1n1_spt_dns, method = "spearman")["estimate"]
scb_h1n1_cntdns_spt_pvl <- cor.test(scb_h1n1_spt_cnt, scb_h1n1_spt_dns, method = "spearman", exact = FALSE)["p.value"]
stp_h1n1_cntdns_spt_cor <- cor.test(stp_h1n1_spt_cnt, stp_h1n1_spt_dns, method = "spearman")["estimate"]
stp_h1n1_cntdns_spt_pvl <- cor.test(stp_h1n1_spt_cnt, stp_h1n1_spt_dns, method = "spearman", exact = FALSE)["p.value"]
tex_h1n1_cntdns_spt_cor <- cor.test(tex_h1n1_spt_cnt, tex_h1n1_spt_dns, method = "spearman")["estimate"]
tex_h1n1_cntdns_spt_pvl <- cor.test(tex_h1n1_spt_cnt, tex_h1n1_spt_dns, method = "spearman", exact = FALSE)["p.value"]

par(mfrow=c(2,2))
plot(pur_h1n1_nmp_cnt, pur_h1n1_nmp_dns, col="blue", ylim = c(0,0.1), xlim = c(0,30), xlab="", ylab="")
par(new=TRUE)
plot(pur_h1n1_spt_cnt, pur_h1n1_spt_dns, col="red", ylim = c(0,0.1), xlim = c(0,30), xlab="Anz. Epitope",
     ylab="Epitopdichte", main = expression(italic("H1N1, Puerto Rico/1934")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

plot(scb_h1n1_nmp_cnt, scb_h1n1_nmp_dns, col="blue", ylim = c(0,0.1), xlim = c(0,30), xlab="", ylab="")
par(new=TRUE)
plot(scb_h1n1_spt_cnt, scb_h1n1_spt_dns, col="red", ylim = c(0,0.1), xlim = c(0,30), xlab="Anz. Epitope",
     ylab="Epitopdichte", main = expression(italic("H1N1, S.Canterbury/2000")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

plot(stp_h1n1_nmp_cnt, stp_h1n1_nmp_dns, col="blue", ylim = c(0,0.1), xlim = c(0,30), xlab="", ylab="")
par(new=TRUE)
plot(stp_h1n1_spt_cnt, stp_h1n1_spt_dns, col="red", ylim = c(0,0.1), xlim = c(0,30), xlab="Anz. Epitope",
     ylab="Epitopdichte", main = expression(italic("H1N1, St.Petersburg/2006")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

plot(tex_h1n1_nmp_cnt, tex_h1n1_nmp_dns, col="blue", ylim = c(0,0.1), xlim = c(0,30), xlab="", ylab="")
par(new=TRUE)
plot(tex_h1n1_spt_cnt, tex_h1n1_spt_dns, col="red", ylim = c(0,0.1), xlim = c(0,30), xlab="Anz. Epitope",
     ylab="Epitopdichte", main = expression(italic("H1N1, Texas/2007")))
legend("topright", legend = c("NetMHCpan", "SYFPEITHI"), pch = 1, col = c("blue", "red"))
par(new=FALSE)

# Plot Sequence Correlations
h1n1_pur_up_cor_overview <- matrix(round(c(as.numeric(pur_h1n1_lenmas_cor),
                                         as.numeric(pur_h1n1_lenmas_pvl),
                                         as.numeric(pur_h1n1_lencnt_nmp_cor),
                                         as.numeric(pur_h1n1_lencnt_nmp_pvl),
                                         as.numeric(pur_h1n1_lencnt_spt_cor),
                                         as.numeric(pur_h1n1_lencnt_spt_pvl)),digits = 4),
                                   ncol = 3)
h1n1_scb_up_cor_overview <- matrix(round(c(as.numeric(scb_h1n1_lenmas_cor),
                                         as.numeric(scb_h1n1_lenmas_pvl),
                                         as.numeric(scb_h1n1_lencnt_nmp_cor),
                                         as.numeric(scb_h1n1_lencnt_nmp_pvl),
                                         as.numeric(scb_h1n1_lencnt_spt_cor),
                                         as.numeric(scb_h1n1_lencnt_spt_pvl)),digits = 4),
                                   ncol = 3)
h1n1_stp_up_cor_overview <- matrix(round(c(as.numeric(stp_h1n1_lenmas_cor),
                                         as.numeric(stp_h1n1_lenmas_pvl),
                                         as.numeric(stp_h1n1_lencnt_nmp_cor),
                                         as.numeric(stp_h1n1_lencnt_nmp_pvl),
                                         as.numeric(stp_h1n1_lencnt_spt_cor),
                                         as.numeric(stp_h1n1_lencnt_spt_pvl)),digits = 4),
                                   ncol = 3)
h1n1_tex_up_cor_overview <- matrix(round(c(as.numeric(tex_h1n1_lenmas_cor),
                                         as.numeric(tex_h1n1_lenmas_pvl),
                                         as.numeric(tex_h1n1_lencnt_nmp_cor),
                                         as.numeric(tex_h1n1_lencnt_nmp_pvl),
                                         as.numeric(tex_h1n1_lencnt_spt_cor),
                                         as.numeric(tex_h1n1_lencnt_spt_pvl)),digits = 4),
                                   ncol = 3)
colnames(h1n1_pur_up_cor_overview) <- c("Length ~ Mass",
                                   "Length ~ Epitopes (N)",
                                   "Length ~ Epitopes (S)")
colnames(h1n1_scb_up_cor_overview) <- c("Length ~ Mass",
                                        "Length ~ Epitopes (N)",
                                        "Length ~ Epitopes (S)")
colnames(h1n1_stp_up_cor_overview) <- c("Length ~ Mass",
                                        "Length ~ Epitopes (N)",
                                        "Length ~ Epitopes (S)")
colnames(h1n1_tex_up_cor_overview) <- c("Length ~ Mass",
                                        "Length ~ Epitopes (N)",
                                        "Length ~ Epitopes (S)")
rownames(h1n1_pur_up_cor_overview) <- c("Correlation", "P-Val")
rownames(h1n1_scb_up_cor_overview) <- c("Correlation", "P-Val")
rownames(h1n1_stp_up_cor_overview) <- c("Correlation", "P-Val")
rownames(h1n1_tex_up_cor_overview) <- c("Correlation", "P-Val")

h1n1_up_cor_title1.1 <- textGrob(" ")
h1n1_up_cor_title1.2 <- textGrob("Sequence Correlations")
h1n1_up_cor_subtt1.1 <- textGrob(" ")
h1n1_up_cor_subtt1.2 <- textGrob(expression(paste(italic("Influenza A"), "-Strains")))
h1n1_pur_up_cor_title2 <- textGrob(expression(italic("Puerto Rico 1934")))
h1n1_scb_up_cor_title2 <- textGrob(expression(italic("S.Canterbury 2000")))
h1n1_stp_up_cor_title2 <- textGrob(expression(italic("St.Petersburg 2006")))
h1n1_tex_up_cor_title2 <- textGrob(expression(italic("Texas 2007")))

h1n1_pur_up_cor_table <- tableGrob(h1n1_pur_up_cor_overview, theme = table_theme1)
h1n1_scb_up_cor_table <- tableGrob(h1n1_scb_up_cor_overview, theme = table_theme2)
h1n1_stp_up_cor_table <- tableGrob(h1n1_stp_up_cor_overview, theme = table_theme2)
h1n1_tex_up_cor_table <- tableGrob(h1n1_tex_up_cor_overview, theme = table_theme2)

grid.arrange(h1n1_up_cor_title1.1, h1n1_up_cor_title1.2,
             h1n1_up_cor_subtt1.1, h1n1_up_cor_subtt1.2,
        h1n1_pur_up_cor_title2, h1n1_pur_up_cor_table,
             h1n1_scb_up_cor_title2, h1n1_scb_up_cor_table,
             h1n1_stp_up_cor_title2, h1n1_stp_up_cor_table,
             h1n1_tex_up_cor_title2, h1n1_tex_up_cor_table,
             ncol=2, widths = unit(c(100,100), rep("mm", 2)),
             heights = unit(c(1,rep(25,5)), rep("mm",5)))

