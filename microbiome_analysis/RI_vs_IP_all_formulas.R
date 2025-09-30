#restart/refresh
rm(list = ls())

setwd("C:/Users/kirchon/Documents/Projects/Thomas_vaccine/R")

source("string_splitter_metaphlan.R")

require(ggplot2)
require(vegan)
require(coin)
require(qvalue)
require(dunn.test)
require(gridExtra)
require(cartography)
require(Maaslin2)

#kreport <- read.table("combined_kreports.tsv", fill = TRUE, sep = "\t", row.names = NULL, header = TRUE)
meta <- read.table("metadata.csv", sep = ",",header = TRUE)
#colnames(meta) <- c("sampleid","groupid", "group", "fourTU", "host" )
phlan <- read.table("merged_abundance_table.txt", header = TRUE)
tax <- as.matrix(phlan[,1])
tax.split <- as.data.frame(string_splitter(tax))

pt <- cbind(tax.split, phlan)

#pt2 <- pt[-which(pt$Phylum == "Bacteria_unclassified"),]
pt2 <- pt[,-c(1,2,3,4,5,8)]


#get genera df
pt3 <- pt2[!is.na(pt2$Genus),]
pt4 <- pt3[is.na(pt3$Species),]  
pt5 <- pt4[,-2]
rownames(pt5) <-pt5[,1]
pt6 <- pt5[,-c(1,2)]

names <- meta$Mouse_ID

colnames(pt6) <- (c(names))

#confirm the ids are in the same order
match(meta$Mouse_ID, colnames(pt6))

#remove the GGBs from the data frame
pt6 <- pt6[-grep("GGB", rownames(pt6)),]

#remove the None administration routes
ids <- meta$Mouse_ID[which(meta$Route == "None")]
pt6 <- pt6[,-which(colnames(pt6) %in% ids)]
meta <- meta[-which(meta$Mouse_ID %in% ids),]
###CHANGE ME
meta <- meta
abund <- pt6

classes <- meta$Route
df <- abund
dist.meth <- "bray"   
dist <- vegdist(t(df), method=dist.meth, na.rm = TRUE)
groups <- classes
names(groups) <- colnames(df)

write.csv(abund, "results/RI_vs_IP_all_formulas/csvs/RI_vs_IP_all_formulas_abundance.csv")

pdf("results/RI_vs_IP_all_formulas/taxonomy/genus/plots/pcoa/RI_vs_IP_all_formuals.pdf", onefile = TRUE, width = 10, height = 8)
# PCoA
bray_curtis_pcoa <- ecodist::pco(dist)

# All components could be found here: 
# bray_curtis_pcoa$vectors
# But we only need the first two to demonstrate what we can do:
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])

bray_curtis_pcoa_df$Group <- factor(meta$Route, levels = c("RI", "IP"))

write.csv(bray_curtis_pcoa_df, "results/RI_vs_IP_all_formulas/csvs/RI_vs_IP_all_formulas_pcoa.csv")

pc1 <- bray_curtis_pcoa$values[1]/sum(bray_curtis_pcoa$values)
pc2 <- bray_curtis_pcoa$values[2]/sum(bray_curtis_pcoa$values)

pc1.round <- format(round((pc1*100), 2), nsmall = 2)
pc2.round <- format(round((pc2*100), 2), nsmall = 2)


colors <- c(carto.pal(pal1 = "multi.pal", n1 = 4))

# Creates a plot
plot <- ggplot(data = bray_curtis_pcoa_df, aes_string(x = "pcoa1", y = "pcoa2", color = "Group")) +
  geom_point() +
  labs(x = paste0("PC1 (", pc1.round, "%)"),
       y = paste0("PC2 (", pc2.round, "%)"), 
       title = "Bray-Curtis PCoA") +
  theme_classic(26) +
  geom_point(size=5) +
  scale_color_manual(values = colors) +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    legend.title=element_blank())

plot
dev.off()

classes.df <- as.data.frame(classes)
colnames(classes.df) <- "Combination"

#comparing beta diversity values with adonis
set.seed(234)
adon <- adonis2(dist ~ Combination, classes.df)
per.pval <- adon$`Pr(>F)`[1]

#comparing beta dispersion values with adonis
bd <- betadisper(dist, groups, type = "centroid")
bd.adon <- adonis2(bd$distances ~ bd$group)
bd.pval <- bd.adon$`Pr(>F)`[1]

#saving adonis pvalues
pval.df <- as.data.frame(matrix(c(per.pval,bd.pval), nrow = 1, ncol = 2))
colnames(pval.df) <- c("beta-diversity", "beta-dispersion")
rownames(pval.df) <- "pvalue"
write.csv(pval.df, "results/RI_vs_IP_all_formulas/taxonomy/genus/pvalues/beta_diversity/allgroups.csv")

####diff comparisons kruskal

ctu.txt <- t(abund)

#set kw function
get_kruskal_pvalue<- function( dataframe.row, classes){
  df = data.frame(abundance = dataframe.row, feature = classes )
  p <- pvalue( 
    kruskal_test(
      abundance~factor(feature), 
      data = df,
      distribution= approximate(nresample = 5000)
    )
  )
  return(p)
}


#kruskal
kdf <- as.data.frame(matrix(data = 0, nrow = dim(ctu.txt)[2], ncol = 2))
set.seed(431)
for (i in 1:dim(ctu.txt)[2]){
  print(i)
  kr <- get_kruskal_pvalue(ctu.txt[,i], classes)[1]
  c <- colnames(ctu.txt)[i]
  rownames(kdf)[i] <- c
  kdf[i,1] <- kr
}

pa <- p.adjust(kdf[,1], method = "fdr")
kdf[,2] <- pa
colnames(kdf) <- c("pval", "fdr")
kdf.ord <- kdf[order(kdf[,2], decreasing = FALSE),]
write.csv(kdf.ord, "results/RI_vs_IP_all_formulas/taxonomy/genus/pvalues/diff_comparisons/RI_vs_IP_all_formulas.csv")

#Dunn's test

pdf("results/RI_vs_IP_all_formulas/taxonomy/genus/pvalues/diff_comparisons/RI_vs_IP_all_formulas_dunn.pdf", onefile = TRUE, height = 10, width = 7)
new.kdf.ord <- kdf.ord[kdf.ord[,2] < 0.05,]
for(i in 1:dim(new.kdf.ord)[1]){
  print(c("this is the i",i))
  dun1 <- dunn.test(as.numeric(abund[i,]), groups, method = "bh")
  col1 <- dun1$comparisons
  col2 <- dun1$P.adjusted
  dunn.df <- cbind(col1, col2)
  colnames(dunn.df) <- c(rownames(kdf.ord)[i], "Adjusted P-value")
  dunn.df <- dunn.df[order(dunn.df[,2], decreasing = FALSE),]
  dunn.df <- as.data.frame(dunn.df[dunn.df[,2] < 0.05,])
  if (dim(dunn.df)[1] == 0) {
    dunn.df <- data.frame(col1 = "None",col2 = "None")
    colnames(dunn.df) <- c(rownames(kdf.ord)[i], "Adjusted P-value")
  }
  if (dim(dunn.df)[2] == 1){
    dunn.df <- t(dunn.df)
  }
  rownames(dunn.df) <- NULL
  #which(dunn.df$Adjust)
  #table <- tableGrob(dunn.df)
  #grid.arrange(table, heights = unit(c(10, 0), "npc"))
  grid.table(dunn.df)
  plot.new()
}

dev.off()

############################
# Boxplots #
###########################

taxa_in_order <- rownames(kdf.ord)

var1 <- "RI"
var2 <- "IP"

vars <- c(var1, var2)

#set up ids
for (i in 1:length(vars)){
  assign(paste0("ids.", i), meta$Mouse_ID[which(meta$Route == vars[i])])
}

#set up groups in required dataframe format
for (i in 1:length(vars)){
  vex <- eval(sym(paste0("ids.", i)))
  assign(paste0("group", i), abund[,vex])
}

#transpose
for (i in 1:length(vars)){
  hey <- eval(sym(paste0("group", i)))
  assign(paste0("group.", i), t(hey))
}

#combine
combo <- NULL
for (i in 1:length(vars)){
  com <- eval(sym(paste0("group.", i)))
  combo <- as.data.frame(rbind(combo, com))
}

fact <- NULL
for (i in 1:length(vars)){
  grp <- eval(sym(paste0("group.", i)))
  reps <- rep(vars[i], nrow(grp))
  fact <- c(fact, reps)
}

combo$Group <- factor(fact, levels = vars)

combo.ord <- combo[,match(rownames(kdf.ord), colnames(combo) )]

#all groups boxplots
pdf("results/RI_vs_IP_all_formulas/taxonomy/genus/plots/box_plots/RI_vs_IP_all_formulas.pdf", onefile = TRUE, width = 10, height = 7)

for(i in 1:dim(kdf.ord)[1]){
  print(i)
  ind <- i
  
  taxa_abund <- combo.ord[,ind]
  tax_title <- (colnames(combo.ord)[ind])
  pv <- kdf.ord[ind,2]
  
  
  
  p <- ggplot(combo, aes(x=Group, y=taxa_abund, fill=Group)) + 
    geom_boxplot(outlier.shape = NA)
  
  print(p + theme_classic(26) + 
          scale_fill_manual(values = colors) +
          labs(title = tax_title,
               subtitle = paste("p-value = ",pv)) +
          geom_jitter(size = 2) +
          theme(axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                legend.title=element_blank(),
                axis.title.x=element_blank(),
          )+
          ylab("Relative Abundance"))
}
dev.off()

#PAM
library(cluster)
library(factoextra)

scaleddata = scale(t(pt6))
fviz_nbclust(scaleddata, pam, method ="silhouette")+theme_minimal()

pamResult <-pam(scaleddata, k = 2)
pamResult

#Masslin2
require(Maaslin2)

rownames(meta) <- meta$Sample_ID


fit_data2 = Maaslin2(
  input_data = pt6, 
  input_metadata = meta, 
  output = "maaslin2_output",
  reference = c("Group,High_Cal_No_Pain"),
  max_significance = 0.05,
  min_prevalence = 0.1,
  correction = "BH",
  analysis_method = "LM",
  #normalization = "CSS",
  #transform = "LOG",
  fixed_effects = c("Group", "Age"))

#
alphadiv <- as.data.frame(diversity(t(abund), index = "shannon"))
match(names(alphadiv), meta$Mouse_ID)

alphadiv[,2] <- meta$Route

colnames(alphadiv) <- c("alphadiv", "Group")

set.seed(123)
pv <- pvalue(kruskal_test(alphadiv~factor(classes), data = alphadiv, distribution= approximate(nresample = 5000)))[1]
combo <- alphadiv

write.csv(combo, "results/RI_vs_IP_all_formulas/csvs/RI_vs_IP_all_formulas_alpha_div.csv")

pdf("results/RI_vs_IP_all_formulas/taxonomy/genus/plots/box_plots/RI_vs_IP_all_formulas_alpha_diversity.pdf", onefile = TRUE, width = 10, height = 7)

p <- ggplot(combo, aes(x=Group, y=alphadiv, fill=Group)) + 
  geom_boxplot(outlier.shape = NA)

print(p + theme_classic(26) + 
        scale_fill_manual(values = colors) +
        labs(title = "Alpha Diversity",
             subtitle = paste("p-value = ",pv)) +
        geom_jitter(size = 2) +
        theme(axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              legend.title=element_blank(),
              axis.title.x=element_blank(),
        )+
        ylab("Shannon Diversity Index"))


dev.off()


