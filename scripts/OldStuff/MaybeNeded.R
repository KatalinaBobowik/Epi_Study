# compare vivax vs falciparum in LFC values
allENSEMBL=data.frame(c(topTable.FalcvsHealthy$ENSEMBL,topTable.VivaxvsHealthy$ENSEMBL))
colnames(allENSEMBL)=c("ENSEMBL")
vivax_df=data.frame(topTable.VivaxvsHealthy$ENSEMBL, topTable.VivaxvsHealthy$AveExpr)
colnames(vivax_df)=c("ENSEMBL","Vivax")
combined_vvx=merge(allENSEMBL, vivax_df, all = TRUE)
falc_df=data.frame(topTable.FalcvsHealthy$ENSEMBL, topTable.FalcvsHealthy$AveExpr)
colnames(falc_df)=c("ENSEMBL","Falciparum")
combined_all=merge(combined_vvx, falc_df, all = TRUE)

# plot density of LFC values
meltedLFC=melt(combined_all)
colnames(meltedLFC)=c("ENSEMBL","Species","LogFC")
pdf(paste0(outputdir,"FalciparumAndVivax_LFC.pdf"))
ggplot(meltedLFC, aes(x=AveExpr, fill=Species)) + geom_density(alpha=0.4) + theme_bw(base_size = 14) + scale_fill_manual(values=c(species.col[3], species.col[1])) +
  labs(title = "LFC Distribution of Significant Genes", y="Density" ) +
  annotate("text", x = -3.8, y = 0.65, label = paste0("Vivax Variance = ",round(var(topTable.VivaxvsHealthy$logFC),2)), size = 4) +
  annotate("text", x = -3.5, y = 0.62, label = paste0("Falciparum Variance = ",round(var(topTable.FalcvsHealthy$logFC),2)), size = 4)
dev.off()




--
# get row variance and plot
VVX=data.frame(rownames(topTable.VivaxvsHealthy),rowVars(vDup$E[rownames(topTable.VivaxvsHealthy),which(new$samples$species=="vivax")]))
colnames(VVX)=c("ENSEMBL","Vivax")
FALC=data.frame(rownames(vDup$E[rownames(topTable.VivaxvsHealthy),which(new$samples$species=="falciparum")]),rowVars(vDup$E[rownames(topTable.VivaxvsHealthy),which(new$samples$species=="falciparum")]))
colnames(FALC)=c("ENSEMBL","Falciparum")
combined_all_var=merge(VVX, FALC, all = TRUE)
meltedVar=melt(combined_all_var)
colnames(meltedVar)=c("ENSEMBL","Species","Variance")
p=ggplot(meltedVar, aes(x=Variance, fill=Species)) + geom_density(alpha=0.4) + theme_bw(base_size = 14) +
labs(title = "Variance Distribution of Significant Genes", y="Density" )
pdf(paste0(outputdir,"FalciparumAndVivax_Variance_DensityPlot.pdf"), width=10)
p+scale_fill_manual(values=c(species.col[3], species.col[1]))
dev.off()


---

VVX=data.frame(rownames(vDup$E[rownames(topTable.FalcvsHealthy),which(new$samples$species=="vivax")]),rowVars(vDup$E[rownames(topTable.FalcvsHealthy),which(new$samples$species=="vivax")]))
colnames(VVX)=c("ENSEMBL","Vivax")
FALC=data.frame(rownames(vDup$E[rownames(topTable.FalcvsHealthy),which(new$samples$species=="falciparum")]),rowVars(vDup$E[rownames(topTable.FalcvsHealthy),which(new$samples$species=="falciparum")]))
colnames(FALC)=c("ENSEMBL","Falciparum")
combined_all_var=merge(VVX, FALC, all = TRUE)
meltedVar=melt(combined_all_var)
colnames(meltedVar)=c("ENSEMBL","Species","Variance")
p=ggplot(meltedVar, aes(x=Variance, fill=Species)) + geom_density(alpha=0.4) + theme_bw(base_size = 14) +
labs(title = "Variance Distribution of Significant Genes", y="Density" )
pdf(paste0(outputdir,"FalciparumAndVivax_Variance_DensityPlot.pdf"), width=10)
p+scale_fill_manual(values=c(species.col[3], species.col[1]))
dev.off()




---

new=merged[,c(which(merged$samples$species=="control"),which(merged$samples$species=="vivax"),which(merged$samples$species=="falciparum"))]

# Set up design matrix
design <- model.matrix(~0 + new$samples$species + new$samples$Island + new$samples$sex + new$samples$EEF_falciparum + new$samples$Merozoite_falciparum + new$samples$oocyst_falciparum + new$samples$Ring_falciparum + new$samples$ookoo_falciparum + new$samples$Trophozoite_falciparum + new$samples$bbSpz_vivax + new$samples$Female_vivax + new$samples$Male_vivax + new$samples$Merozoite_vivax + new$samples$oocyst_vivax + new$samples$ook_vivax  + new$samples$Ring_vivax + new$samples$ookoo_vivax + new$samples$Schizont_vivax + new$samples$sgSpz_vivax + new$samples$Trophozoite_vivax + new$samples$Gran + new$samples$Bcell + new$samples$CD4T + new$samples$CD8T + new$samples$NK + new$samples$Mono)


design <- model.matrix(~0 + new$samples$species + new$samples$Island + new$samples$sex + new$samples$EEF_falciparum + new$samples$Merozoite_falciparum + new$samples$oocyst_falciparum + new$samples$Ring_falciparum + new$samples$ookoo_falciparum + new$samples$Trophozoite_falciparum + new$samples$Gran + new$samples$Bcell + new$samples$CD4T + new$samples$CD8T + new$samples$NK + new$samples$Mono)


# rename columns to exclude spaces and unrecognised characters
colnames(design)=gsub("\\$", "", colnames(design)) %>% gsub("newsamplesspecies", "", .) %>% gsub("newsamplessex", "", .)
colnames(design)=gsub("\\$", "", colnames(design)) %>% gsub("newsamplesIsland", "", .)
colnames(design)=gsub("West Papua", "Mappi",  colnames(design))
colnames(design)=gsub("newsamples", "",  colnames(design))

# set up contrast matrix
contr.matrix <- makeContrasts(SickvsHealthy=(vivax + falciparum)/2 - control, FalcvsHealthy=falciparum - control, VivaxvsHealthy=vivax - control, VivaxvsFalciparum=vivax - falciparum, levels=colnames(design))

# Using duplicate correlation and blocking -----------------------------------------------------

# First, we need to perform voom normalisation
v <- voom(new, design, plot=T)


# With duplicate correction and blocking:
# the inter-subject correlation is input into the linear model fit
voomDupVfit <- lmFit(v, design)
voomDupVfit <- contrasts.fit(voomDupVfit, contrasts=contr.matrix)
voomDupEfit <- eBayes(voomDupVfit, robust=T)
dt <- decideTests(voomDupEfit, p.value=0.05, lfc=0)

# explore different pvalue thresholds
summary(dt)

#        SickvsHealthy FalcvsHealthy VivaxvsHealthy VivaxvsFalciparum
# Down             578           169            338                 0
# NotSig         10013         11184          10657             11729
# Up              1144           382            740                 6

dt <- decideTests(voomDupEfit, p.value=0.01, lfc=0)
summary(dt)





VVX=data.frame(rownames(merged),rowVars(vDup$E[,which(merged$samples$species=="vivax")]))
colnames(VVX)=c("ENSEMBL","Vivax")
FALC=data.frame(rownames(merged),rowVars(vDup$E[,which(merged$samples$species=="falciparum")]))
colnames(FALC)=c("ENSEMBL","Falciparum")
combined_all_var=merge(VVX, FALC, all = TRUE)
meltedVar=melt(combined_all_var)
colnames(meltedVar)=c("ENSEMBL","Species","Variance")

# do Levene's test
vivax_df$Species="Vivax"
colnames(vivax_df)[2]="LogFC"
falc_df$Species="Falciparum"
colnames(falc_df)[2]="LogFC"
levene=rbind(vivax_df,falc_df)
levene.test(levene[,"LogFC"], location="mean", levene[,"Species"], correction.method=c("zero.correction"))

topTable.FalcvsHealthy <- topTable(voomDupEfit, coef=2, p.value=0.05, n=Inf, lfc=0, sort.by="p")
topTable.VivaxvsHealthy <- topTable(voomDupEfit, coef=3, p.value=0.05, n=Inf, lfc=0, sort.by="p")

x=data.frame(topTable.FalcvsHealthy[rownames(dt[which(dt[,2]!=0 & dt[,3]!=0),]),"logFC"],topTable.VivaxvsHealthy[rownames(dt[which(dt[,2]!=0 & dt[,3]!=0),]),"logFC"])
colnames(x) = c("Falc","Vivax")
rownames(x) = merged[rownames(dt[which(dt[,1]!=0 & dt[,3]!=0 & dt[,4]!=0),]),]$genes$SYMBOL
pdf(paste0(outputdir,"HeatmapFalciparum_AllComparisonPops.pdf"))
Heatmap(as.matrix(x))
