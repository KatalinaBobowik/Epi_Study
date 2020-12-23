# library
library(tidyverse)

# this is what I want, but I need a gap----
data=as.matrix(as.data.frame(otu_table(relative_phyloCounts)))
melt=melt(data)
melt$Island='NA'
melt$Island=sapply(strsplit(as.character(melt$Var2), "[-.]"), `[`, 1)


# Get the name and the y position of each label
base_data <- melt %>% 
  group_by(Island) %>% 
  dplyr::summarize(start=min(as.numeric(as.factor(Var2))), end=max(as.numeric(as.factor(Var2)))) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

ggplot(melt, aes(x=as.factor(Var2), y=value, fill=Var1)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(aes(x=as.factor(Var2), y=value, fill=Var1, width=0.8), stat="identity", alpha=0.8) + 
  theme_minimal() + coord_polar() +
  theme(axis.ticks =element_blank(), legend.position = "none",  axis.title = element_blank(), plot.margin = unit(rep(-1,4), "cm"), axis.text.x=element_blank(), axis.text.y=element_blank(),panel.grid = element_blank(),) + 
# Add base line information
  geom_segment(data=base_data, aes(x = start, y = 1.1, xend = end, yend = 1.1), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = 1.2, label=Island), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE) +
  scale_fill_manual(values=Merged_Palette)


relative_phyloCounts <- microbiome::transform(aggregated_phyloCounts, "compositional")
p=plot_bar(relative_phyloCounts, fill = "Superkingdom")

# set colour palette
families=levels(p$data$Superkingdom)
# get number of families in each kingdom
table(sapply(strsplit(families, "[_.]"), `[`, 1))

PaletteBacteria = colorRampPalette(c("#023858","#74a9cf"))(13)
PaletteEukaryote = colorRampPalette(c("#4c0000","#b20000","#ff4c4c"))(4)
PaletteOther = colorRampPalette(c("black"))(1)
PaletteUnk = colorRampPalette(c("#5a5a5a"))(1)
PaletteVirus = colorRampPalette(c("#78c679","#006837"))(2)

Merged_Palette <- c(PaletteBacteria,PaletteEukaryote,PaletteOther,PaletteUnk,PaletteVirus)


----

melt$Island = as.factor(melt$Island)
# Get the name and the y position of each label
base_data <- melt %>% 
  group_by(Island) %>% 
  dplyr::summarize(start=min(as.numeric(as.factor(Var2))), end=max(as.numeric(as.factor(Var2)))) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 1
to_add <- data.frame( matrix(NA, empty_bar*nlevels(as.factor(melt$Var2)), ncol(melt)) )
colnames(to_add) <- colnames(melt)
to_add$Var2 <- rep(levels(as.factor(melt$Var2)), each=empty_bar)
to_add$Island <- sapply(strsplit(as.character(to_add$Var2), "[-.]"), `[`, 1)

melt <- rbind(melt, to_add)
melt <- melt %>% arrange(Island)

ggplot(melt, aes(x=as.factor(Var2), y=value, fill=Var1)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(aes(x=as.factor(Var2), y=value, fill=Var1, width=0.5), stat="identity", alpha=0.5) + 
  theme_minimal() + coord_polar() +
  theme(axis.ticks =element_blank(), legend.position = "none",  axis.title = element_blank(), plot.margin = unit(rep(-1,4), "cm"), axis.text.x=element_blank(), axis.text.y=element_blank(),panel.grid = element_blank(),) + 
# Add base line information
  geom_segment(data=base_data, aes(x = start, y = 1.1, xend = end, yend = 1.1), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = 1.2, label=Island), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE) 

---

# this is what I want, but I need a gap----
data=as.matrix(as.data.frame(otu_table(relative_phyloCounts)))
data=t(data)
data=as.data.frame(data)
data$group=sapply(strsplit(rownames(data), "[-.]"), `[`, 1)
data$individual=rownames(data)
data$individual=as.factor(data$individual)
data$group=as.factor(data$group)
 
# Transform data in a tidy format (long format)
data = data %>% gather(key = "observation", value="value", -c(22,23)) 

# Set a number of 'empty bar' to add at the end of each group
empty_bar=3
nObsType=nlevels(as.factor(data$observation))
to_add = data.frame( matrix(NA, empty_bar*nlevels(data$group)*nObsType, ncol(data)) )
colnames(to_add) = colnames(data)
to_add$group=rep(levels(data$group), each=empty_bar*nObsType )
data=rbind(data, to_add)
data=data %>% arrange(group, individual)
data$id=rep( seq(1, nrow(data)/nObsType) , each=nObsType)
 
 
# prepare a data frame for base lines
base_data=data %>% 
  group_by(group) %>% 
  dplyr::summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
#grid_data=grid_data[-1,]
 
# Make the plot
p = ggplot(data) +      
  
  # Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.8) +
  # Add a valu=100/75/50/25 lines. I do it at the beginning to make sure barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.25, xend = start, yend = 0.25), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.50, xend = start, yend = 0.50), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.75, xend = start, yend = 0.75), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 1, xend = start, yend = 1), colour = "white", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  #annotate("text", x = rep(max(data$id),5), y = c(0, 0.25, 0.50, 0.75, 1), label = c("0", "0.25", "0.50", "0.75", "1") , color="black", size=2 , angle=0, fontface="bold", hjust=1) +

  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=base_data, aes(x = title, y = 1.2, label=group), hjust=c(0,0,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE) +
  scale_fill_manual(values=Merged_Palette)

pdf(paste0(outputdir,"relativeTaxa_Compositional_CircularBarplot.pdf"), width=15)
p
dev.off()

