###################### Data Analysis Functions ################################

# Libraries
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(plyr)
library(purrr)
library(readxl)
library(rlist)
library(stringr)

# Lineages
find_strain_na <- function(lineages) {
  i = 1
  strains = vector()
  for (lineage in t(lineages)) {
    strain <- unlist( regmatches(lineage, gregexpr("\\(.*?\\)", lineage)))[1]
    na_strain <- lapply( strain, function(x) if(identical(x, NA_character_)) "N/A" else x  )
    clean_strain <- gsub("[\\(\\)]|\\bstrain\\b|\\s", "", unlist(na_strain))
    strains[i] = clean_strain
    i = i+1
  }
  return(strains)
}
find_strain <- function(lineages) {
  i = 1
  strains = vector()
  for (lineage in t(lineages)) {
    strain <- unlist( regmatches(lineage, gregexpr("\\(.*?\\)", lineage)))[1]
    clean_strain <- gsub("[\\(\\)]|\\bstrain\\b|\\s", "", unlist(strain))
    strains[i] = clean_strain
    i = i+1
  }
  return(strains)
}

# tableTheme
table_theme <- ttheme_minimal(colhead = list(bg_params = list(fill="white")),
                              core = list(bg_params = list(fill=c(blues9[1],blues9[2]))))
table_theme1 <- ttheme_minimal(colhead = list(bg_params = list(fill="white"),
                                             fg_params = list(fontface=2)),
                              rowhead = list(fg_params = list(fontface=2)),
                              core = list(bg_params = list(fill=c(blues9[1],blues9[2]))))
table_theme2 <- ttheme_minimal(colhead = list(bg_params = list(fill="white"),
                                              fg_params = list(fontface=2, col="white")),
                               rowhead = list(fg_params = list(fontface=2)),
                               core = list(bg_params = list(fill=c(blues9[1],blues9[2]))))
table_theme3 <- ttheme_minimal(colhead = list(bg_params = list(fill="white"),
                                              fg_params = list(fontface=4)),
                               rowhead = list(fg_params = list(fontface=2)),
                               core = list(bg_params = list(fill=c(blues9[1],blues9[2]))))
table_theme4 <- ttheme_minimal(colhead = list(bg_params = list(fill="white"),
                                              fg_params = list(fontface=3)),
                               rowhead = list(fg_params = list(fontface=2)),
                               core = list(bg_params = list(fill=c(blues9[1],blues9[1],blues9[2],blues9[2]))))

# Plot strains
rotate_x <- function(data, labels_vec, rot_angle, hgt, title, n, m) {
  par(mar = c(12,5,hgt,1))
  plt <- barplot(data, col='steelblue', las=2, xaxt="n", main = title, ylim = range(pretty(c(0, data))),
                 space = 0.5)
  text(plt, par("usr")[3], labels = labels_vec, srt = rot_angle, adj = c(1.1,1.1), xpd = TRUE, cex=1)
  legend("topright", legend = c(paste("N", n, sep = " = "),
                                paste("Methode", m, sep = ": ")),
         inset = c(0,-0.375), xpd = TRUE, cex = 0.925)
  par(mar = c(5.1, 4.1, 4.1, 2.1))
}

# Epitope counts
count_epis <- function(eps_data) {
  count = c()
  for (row in 1:nrow(eps_data)) {
    x <- c(eps_data[row,])
    count[row] <- length(x[!is.na(x)])-6
  }
  return(count)
}

tuk_plot <- function (x, xlab, ylab, ylabels = NULL, cexaxis, las, ...) {
  for (i in seq_along(x)) {
    xi <- x[[i]][, -4L, drop = FALSE]
    pi <- x[[i]][, 4, drop = FALSE]
    yvals <- nrow(xi):1L
    dev.hold()
    on.exit(dev.flush())
    plot(c(xi[, "lwr"], xi[, "upr"]), rep.int(yvals, 2L), 
         type = "n", axes = FALSE, xlab = "", ylab = "", main = NULL,
         ...)
    axis(1, ...)
    # change for custom axis labels
    if (is.null(ylabels)) ylabels <- dimnames(xi)[[1L]]
    
    if (nrow(xi)>1) axis(2, at = rev(seq(2,nrow(xi),2)), labels = ylabels[rev(seq(2,nrow(xi),2))]
                         ,las=las,srt = 0, cex.axis=cexaxis, ...)
    axis(2, at = rev(seq(1,nrow(xi),2)), labels = ylabels[rev(seq(1,nrow(xi),2))],
         las=las, srt = 0, cex.axis=cexaxis, ...)
    abline(h = yvals, lty = 1, lwd = 0.5, col = "lightgray")
    abline(v = 0, lty = 2, lwd = 0.5, ...)
    
    color <- rep("black",nrow(pi))
    color[as.vector(pi)<0.05] <- "red"
    segments(xi[, "lwr"], yvals, xi[, "upr"], yvals, col = color, ...)
    segments(as.vector(xi), rep.int(yvals - 0.1, 3L), as.vector(xi), 
             rep.int(yvals + 0.1, 3L), col = color, ...)
    title(# change for custom axis titles
          xlab = xlab, ylab = ylab)
    
    box()
    dev.flush()
    on.exit()
  }
}
