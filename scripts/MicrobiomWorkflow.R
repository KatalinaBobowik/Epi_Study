# created by KSB, 20.01.19
# Script created to inspect where reads are mapping to. 

# load packages
require(ggplot2)
require(RColorBrewer)
library(dplyr)
library(plyr)
library(phyloseq)
library(reshape2)
library(ggpubr)
library(vegan)

# set ggplot colour theme to white
theme_set(theme_bw())

# set up directories
inputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/ReferenceFiles/EpiStudy/"
outputdir = "/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/Epi_Study/"

# Import data and convert to a phyloseq object
# raw_CCMetagen_data <- read.csv(paste0(inputdir,"IndoUnmapped_species_table_noRepeats_RPM.csv"),check.names=FALSE)
raw_CCMetagen_data <- read.csv(paste0(inputdir,"IndoUnmapped_species_table_noRepeats_RPM.csv"),check.names=FALSE)

# Data preprocessing ------------------------------------

# remove plants from the analysis
raw_CCMetagen_data=raw_CCMetagen_data[-which(raw_CCMetagen_data$Kingdom %in% c("Viridiplantae","Fungi")),]
raw_CCMetagen_data=raw_CCMetagen_data[-which(raw_CCMetagen_data$Phylum %in% c("Bacillariophyta","Cnidaria","Echinodermata","Mollusca")),]
raw_CCMetagen_data=raw_CCMetagen_data[-which(raw_CCMetagen_data$Order %in% "Albuginales"),]

# Convert to PhyloSeq object ------------------------------------------

# Separate species' abundances and taxonomy columns
taxa_raw <- as.matrix(raw_CCMetagen_data[,c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species")])
#rownames(taxa_raw) <- taxa_raw[,"Family"]
abund_raw <- as.matrix(raw_CCMetagen_data[,-which(colnames(raw_CCMetagen_data) %in% c("Superkingdom","Kingdom","Phylum", "Class", "Order","Family","Genus","Species"))])
#rownames(abund_raw) <- CCMetagen_data[,which(names(CCMetagen_data)=="Family")]

# convert to Phyloseq object
tax = tax_table(taxa_raw)
taxa = otu_table(abund_raw, taxa_are_rows = TRUE)

CCMeta_physeq = phyloseq(taxa, tax)
# TopNOTUs <- names(sort(taxa_sums(CCMeta_physeq), TRUE))
# TopFamilies <- prune_taxa(TopNOTUs, CCMeta_physeq)
p=plot_bar(CCMeta_physeq, fill = "Family")

samplenames <- colnames(otu_table(CCMeta_physeq))
island <- sapply(strsplit(samplenames, "[-.]"), `[`, 1)

samples_df=data.frame(SampleName=colnames(otu_table(CCMeta_physeq)), SamplePop=island)
samples = sample_data(samples_df)
rownames(samples)=samples$SampleName
sample_data(CCMeta_physeq) <- samples

#######

library(BiocStyle)
library(knitr)
library(rmarkdown)
library(ggplot2)
library(gridExtra)
library(dada2)
library(phyloseq)
library(DECIPHER)
library(phangorn)


# Show available ranks in the dataset
rank_names(ps)
# Create table, number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)



Step 1: Look at taxa prevalence
ggplot(prevdf1, aes(TotalAbundance, Prevalence, color = Phylum)) +
  geom_hline(yintercept = prevalenceThreshold, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
  scale_y_log10() + scale_x_log10() +
  xlab("Total Abundance") +
  facet_wrap(~Phylum)

Step 2: Rarefaction curves

Step 3: Filter data if necessary

Step 5: Convert to RPM
It is usually necessary to transform microbiome count data to account for differences in library size, variance, scale, etc.

Step 6: Plot abundance of reads (do this in counts and RPM)

Step 7: look at principal coordinates analysis (PCoA) for outliers
keepSamples=names(which(colSums(otu_table(CCMeta_physeq))>0))
ps=prune_samples(keepSamples, CCMeta_physeq)

sample_data(pslog)$age_binned <- cut(sample_data(pslog)$age,
  				          breaks = c(0, 100, 200, 400))
out.wuf.log <- ordinate(ps, method = "PCoA", distance = "euclidean")
plot_ordination(ps, out.wuf.log)

evals <- out.wuf.log$values$Eigenvalues
plot_ordination(ps, out.wuf.log, color = "SamplePop") +
  labs(col = "SamplePop") +
  coord_fixed(sqrt(evals[2] / evals[1]))

pca <- prcomp(b[1:228], scale. = F, center=T)
autoplot(pca, data = b, colour = 'group', x=1, y=2, label = TRUE)
autoplot(pca, data = b, colour = 'group', x=2, y=3, label = TRUE)
autoplot(pca, data = b, colour = 'group', x=3, y=4, label = TRUE)


This shows a few phyla for which only one feature was observed. Those
may be worth filtering, and we’ll check that next.



The following ensures that features with ambiguous phylum annotation are
also removed. Note the flexibility in defining strings that should be
considered ambiguous annotation.

```{r removeNAphyla}
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
```


A useful next step is to explore feature *prevalence* in the dataset,
which we will define here as the number of samples in which a taxon
appears at least once.

```{r prevfilter0}
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
```

Are there phyla that are comprised of mostly low-prevalence features? Compute the total and average prevalences of the features in each phylum.

```{r lowprev}
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
```

_Deinococcus-Thermus_ appeared in just over one percent of samples, and Fusobacteria appeared in just 2 samples total. In some cases it might be worthwhile to explore these two phyla in more detail despite this (though probably not Fusobacteria’s two samples). For the purposes of this example, though, they will be filtered from the dataset.

```{r taxfilter}
# Define phyla to filter
filterPhyla = c("Fusobacteria", "Deinococcus-Thermus")
# Filter entries with unidentified Phylum.
ps1 = subset_taxa(ps, !Phylum %in% filterPhyla)
ps1
```

### Prevalence Filtering {#prevalence-filtering .unnumbered}

First, explore the relationship of prevalence and total read count for
each feature. Sometimes this reveals outliers that should probably be
removed, and also provides insight into the ranges of either feature
that might be useful. This aspect depends quite a lot on the
experimental design and goals of the downstream inference, so keep these
in mind. It may even be the case that different types of downstream
inference require different choices here. There is no reason to expect
ahead of time that a single filtering workflow is appropriate for all
analysis.




```{r plotprevalence, fig.width=9, fig.height=5, fig.cap="Taxa prevalence versus total counts."}
# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
```

Each point in Figure \@ref(fig:plotprevalence) is a different taxa.
Exploration of the data in this way is often useful for selecting
filtering parameters, like the minimum prevalence criteria we will used
to filter the data
above.


Sometimes a natural separation in the dataset reveals itself, or at
least, a conservative choice that is in a stable region for which small
changes to the choice would have minor or no effect on the biological
interpreation (stability). Here no natural separation is immediately
evident, but it looks like we might reasonably define a prevalence
threshold in a range of zero to ten percent or so. Take care that this
choice does not introduce bias into a downstream analysis of association
of differential abundance.

The following uses five percent of all samples as the prevalence
threshold.



# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)



Abundance value transformation {#abundance-value-transformation .unnumbered}
------------------------------

It is usually necessary to transform microbiome count data to account
for differences in library size, variance, scale, etc. The `phyloseq`
package provides a flexible interface for defining new functions to
accomplish these transformations of the abundance values via the
function `transform_sample_counts()`.
 The first argument to this function is the phyloseq object you
want to transform, and the second argument is an R function that defines
the transformation. The R function is applied sample-wise, expecting
that the first unnamed argument is a vector of taxa counts in the same
order as the phyloseq object. Additional arguments are passed on to the
function specified in the second argument, providing an explicit means
to include pre-computed values, previously defined
parameters/thresholds, or any other object that might be appropriate for
computing the transformed values of interest.

This example begins by defining a custom plot function, 
`plot_abundance()`, that uses
phyloseq’s function to define a relative abundance graphic. 
We will use
this to compare more easily differences in scale and distribution of the
abundance values in our phyloseq object before and after transformation.


```{r abundancetransformation}
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "sex",y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}
```


library("shiny")
library("miniUI")
library("caret")
library("pls")
library("e1071")
library("ggplot2")
library("randomForest")
library("dplyr")
library("ggrepel")
library("nlme")
library("devtools")
library("reshape2")
library("PMA")
library("structSSI")
library("ade4")
library("ggnetwork")
library("intergraph")
library("scales")
library("jfukuyama/phyloseqGraphTest")
library("genefilter")
library("impute")

We back these claims by illustrating several analyses on the mouse
data prepared above. We experiment with several flavors of exploratory
ordination before shifting to more formal testing and modeling,
explaining the settings in which the different points of view are most
appropriate. Finally, we provide example analyses of multitable data,
using a study in which both metabolomic and microbial abundance
measurements were collected on the same samples, to demonstrate that
the general workflow presented here can be adapted to the multitable
setting.



### Preprocessing {#preprocessing .unnumbered}

Before doing the multivariate projections, we will add a few columns to
our sample data, which can then be used to annotate plots. From Figure
\@ref(fig:preprocessing-setup), we see that the ages of the mice come in a
couple of groups, and so we make a categorical variable corresponding to
young, middle-aged, and old mice. We also record the total number of
counts seen in each sample and log-transform the data as an approximate
variance stabilizing transformation.

```{r preprocessing-setup, fig.cap="Histogram of island groupings",fig.show="hold"}
qplot(sample_data(ps)$SamplePop, geom = "histogram",binwidth=10) + xlab("island")
```  

Figure \@ref(fig:preprocessing-setup) shows that the age covariate
belongs to three separate clusters.



```{r preprocessing2, fig.cap="Histograms comparing raw and log transformed read depths",fig.show="hold"}
qplot(log10(rowSums(otu_table(ps))),binwidth=0.2) +
  xlab("Logged counts-per-sample")
```  

These preliminary plots suggest certain preprocessing steps. 
The
  histogram in Figure \@ref(fig:preprocessing-setup)
  motivates the creation of a new categorical
  variable, binning age into one of the three peaks. 
  
  The histogram in Figure \@ref(fig:preprocessing2)
   suggests that a $\log\left(1 + x\right)$ transformation 
  might be
  sufficient for `normalizing` the abundance data for the exploratory analyses.
  
  In fact this transformation is not sufficient for 
  testing purposes and when performing differential abundances
  we recommend the variance stabilizing transformations
  available in DESeq2 through the `phyloseq_to_deseq2` function,
  see the [phyloseq_to_deseq2 tutorial here](https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html).

As our first step, we look at principal coordinates analysis (PCoA) with
either the Bray-Curtis dissimilarity on the weighted Unifrac
distance. 

```{r outlier-detect, fig.cap="Exploratory ordination analysis with log abundances.",fig.wide= TRUE}
sample_data(ps)$age_binned <- cut(sample_data(ps)$age,
  				          breaks = c(0, 100, 200, 400))
levels(sample_data(ps)$age_binned) <- list(Young100="(0,100]", Mid100to200="(100,200]", Old200="(200,400]")
sample_data(ps)$family_relationship=gsub(" ","",sample_data(ps)$family_relationship)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
out.wuf.log <- ordinate(ps, method = "MDS", distance = "wunifrac")
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned") +
  labs(col = "Binned Age") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

Figure \@ref(fig:outlier-detect) showing the ordination on the logged abundance data reveals a few outliers.

These turn
out to be the samples from females 5 and 6 on day 165 and the samples
from males 3, 4, 5, and 6 on day 175. We will take them out, since we
are mainly interested in the relationships between the non-outlier
points.

Before we continue, we should check the two female
outliers -- they have been taken over by the same OTU/ASV, which has a
relative abundance of over 90\% in each of them. This is the only time
in the entire data set that this ASV has such a high relative
abundance -- the rest of the time it is below 20\%. In particular, its
diversity is by far the lowest of all the samples.


```{r outlier-analyze, fig.width=9, fig.height=5, fig.cap="The outlier samples are dominated by a single ASV."}
rel_abund <- t(apply(otu_table(ps), 1, function(x) x / sum(x)))
qplot(rel_abund[, 12], geom = "histogram",binwidth=0.05) +
  xlab("Relative abundance")
```

# Different Ordination Projections {#different-ordination-projections .unnumbered}

As we have seen, an important first step in analyzing microbiome data is
to do unsupervised, exploratory analysis. This is simple to do in 
`phyloseq`,
which provides many distances and ordination methods.

After documenting the outliers, we are going to compute ordinations with
these outliers removed and more carefully study the output.

```{r remove-outliers}
outliers <- c("F5D165", "F6D165", "M3D175", "M4D175", "M5D175", "M6D175")
ps <- prune_samples(!(sample_names(ps) %in% outliers), ps)
```

We are also going to remove samples with fewer than 1000 reads:

```{r removelowreads}
which(!rowSums(otu_table(ps)) > 1000)
ps <- prune_samples(rowSums(otu_table(ps)) > 1000, ps)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
```

We’ll first
perform a PCoA using Bray-Curtis dissimilarity.

```{r ordinations-bray,fig.cap="A PCoA plot using Bray-Curtis between  samples."}
out.pcoa.log <- ordinate(pslog,  method = "MDS", distance = "bray")
evals <- out.pcoa.log$values[,1]
plot_ordination(pslog, out.pcoa.log, color = "age_binned",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

We see that
there is a fairly substantial age effect that is consistent between all the mice, male and female, and from different litters. 

Next we look at double principal coordinates analysis (DPCoA)
[@Pavoine:2004; @Purdom2010; @Fukuyama:2012], which is a phylogenetic
ordination method and that provides a biplot representation of both
samples and taxonomic categories. We see again that the second axis
corresponds to young vs. old mice, and the biplot suggests an
interpretation of the second axis: samples that have larger scores on
the second axis have more taxa from Bacteroidetes and one subset of
Firmicutes.

```{r ordinations-dpcoa, fig.wide=TRUE,fig.cap="A DPCoA plot incorporates phylogenetic information, but is  dominated by the first axis."}
out.dpcoa.log <- ordinate(pslog, method = "DPCoA")
evals <- out.dpcoa.log$eig
plot_ordination(pslog, out.dpcoa.log, color = "age_binned", label= "SampleID",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

In Figure \@ref(fig:ordinations-dpcoa) we 
have the first axis explains 
`r round(evals[1]/sum(evals)*100)` \% of the variability, about
`r round(evals[1]/evals[2])` times that of the second axis;
this translates into the elongated form of the ordination plot.


```{r dpcoabiplot,fig.wide=TRUE,fig.cap="Taxa responsible for Axis 1 and 2"}
plot_ordination(pslog, out.dpcoa.log, type = "species", color = "Phylum") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

Finally, we can look at the results of PCoA with weighted Unifrac. As
before, we find that the second axis is associated with an age effect,
which is fairly similar to DPCoA. This is not surprising, because both
are phylogenetic ordination methods taking abundance into account.
However, when we compare biplots, we see that the DPCoA gave a much
cleaner interpretation of the second axis, compared to weighted Unifrac.


```{r ordinations-wuf,fig.wide =TRUE, fig.cap="The sample positions produced by a PCoA using weighted Unifrac."}
out.wuf.log <- ordinate(pslog, method = "PCoA", distance ="wunifrac")
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned",
                  shape = "family_relationship") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Binned Age", shape = "Litter")
```



## Why are the ordination plots so far from square? {.unnumbered}
###Aspect ratio of ordination plots {.unnumbered}

In the ordination plots in Figure 8–Figure 14, you may have noticed as did the reviewers of the first version of the paper, that the maps are not presented as square representations as is often the case in standard PCoA and PCA plots in the literature.

The reason for this is that as we are trying to represent the distances between samples as faithfully as possible; we have to take into account that the second eigenvalue is always smaller than the first, sometimes considerably so, thus we normalize the axis norm ratios to the relevant eigenvalue ratios.
This ensures that the variability represented in the plots is done
so faithfully.


### PCA on ranks {#pca-on-ranks .unnumbered}

Microbial abundance data is often heavy-tailed, and sometimes it can be
hard to identify a transformation that brings the data to normality. In
these cases, it can be safer to ignore the raw abundances altogether,
and work instead with ranks. We demonstrate this idea using a
rank-transformed version of the data to perform PCA. First, we create a
new matrix, representing the abundances by their ranks, where the
microbe with the smallest in a sample gets mapped to rank 1, second
smallest rank 2, etc.

```{r rankab}
abund <- otu_table(pslog)
abund_ranks <- t(apply(abund, 1, rank))
```

Naively using these ranks could make differences between pairs of low
and high abundance microbes comparable. In the case where many bacteria
are absent or present at trace amounts, an artificially large difference
in rank could occur[@holmes2011] for minimally abundant taxa. To avoid
this, all those microbes with rank below some threshold are set to be
tied at 1. The ranks for the other microbes are shifted down, so there
is no large gap between ranks. 

```{r rankthreshold}
abund_ranks <- abund_ranks - 329
abund_ranks[abund_ranks < 1] <- 1
```

```{r pca-rank-visualize-procedure, fig.cap="Rank threshold transformation"}
library(dplyr)
library(reshape2)
abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
colnames(abund_df) <- c("sample", "seq", "abund", "rank")
abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
colnames(abund_df) <- c("sample", "seq", "abund", "rank")
sample_ix <- sample(1:nrow(abund_df), 8)
ggplot(abund_df %>%
         filter(sample %in% abund_df$sample[sample_ix])) +
  geom_point(aes(x = abund, y = rank, col = sample),
             position = position_jitter(width = 0.2), size = 1.5) +
  labs(x = "Abundance", y = "Thresholded rank") +
  scale_color_brewer(palette = "Set2")
```  

This transformation is illustrated in
Figure \@ref(fig:pca-rank-visualize-procedure).

The association between abundance and rank, for a few randomly
selected samples. The numbers of the $y$-axis are those supplied to
PCA.

We can now perform PCA and study the resulting biplot, given in 
the Figure
below. To produce annotation for this figure, we
used the following block.

```{r pca-rank-pca-setup}
library(ade4)
ranks_pca <- dudi.pca(abund_ranks, scannf = F, nf = 3)
row_scores <- data.frame(li = ranks_pca$li,
                         SampleID = rownames(abund_ranks))
col_scores <- data.frame(co = ranks_pca$co,
                         seq = colnames(abund_ranks))
tax <- tax_table(ps) %>%
  data.frame(stringsAsFactors = FALSE)
tax$seq <- rownames(tax)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
tax$otu_id <- seq_len(ncol(otu_table(ps)))
row_scores <- row_scores %>%
  left_join(sample_data(pslog))
col_scores <- col_scores %>%
  left_join(tax)
```


```{r pca-rank-pca-plot, fig.wide=TRUE,fig.height=8,fig.cap="The biplot resulting from the PCA after the truncated-ranking transformation."}
evals_prop <- 100 * (ranks_pca$eig / sum(ranks_pca$eig))
ggplot() +
  geom_point(data = row_scores, aes(x = li.Axis1, y = li.Axis2), shape = 2) +
  geom_point(data = col_scores, aes(x = 25 * co.Comp1, y = 25 * co.Comp2, col = Order),
             size = .3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  facet_grid(~ age_binned) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
       y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  coord_fixed(sqrt(ranks_pca$eig[2] / ranks_pca$eig[1])) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```
The results are similar to the PCoA analyses computed without applying a
truncated-ranking transformation, reinforcing our confidence in the
analysis on the original data.

### Canonical correspondence {#canonical-correspondence .unnumbered}

Canonical Correspondence Analysis (CCpnA) is an approach to ordination
of a species by sample table that incorporates supplemental information
about the samples. As before, the purpose of creating biplots is to
determine which types of bacterial communities are most prominent in
different mouse sample types. It can be easier to interpret these
biplots when the ordering between samples reflects sample
characteristics – variations in age or litter status in the mouse data,
for example – and this central to the design of CCpnA.

The function allows us to create biplots where the positions of samples
are determined by similarity in both species signatures and
environmental characteristics; in contrast, principal components
analysis or correspondence analysis only look at species signatures.
More formally, it ensures that the resulting CCpnA directions lie in the
span of the environmental variables; thorough treatments are available
in [@terBraak:1985; @greenacre2007correspondence].

Like PCoA and DPCoA, this method can be run using `ordinate` from the `phyloseq` package . In order to use
supplemental sample data, it is necessary to provide an extra argument,
specifying which of the features to consider – otherwise, defaults to
using all measurements when producing the ordination.


```{r ccpna-correspondence-analysis}
ps_ccpna <- ordinate(pslog, "CCA", formula = pslog ~ age_binned + family_relationship)
```

To access the positions for the biplot,
we can use the function `ordinate` in `phyloseq`.
Further, to facilitate figure annotation, we also join the site scores
with the environmental data in the slot. Of the 23 total taxonomic
orders, we only explicitly annotate the four most abundant – this makes
the biplot easier to read.

```{r ccpna-join-data, fig.cap="The mouse and bacteria scores generated by CCpnA.", fig.wide=TRUE, fig.height=10}
library(ggrepel)
ps_scores <- vegan::scores(ps_ccpna)
sites <- data.frame(ps_scores$sites)
sites$SampleID <- rownames(sites)
sites <- sites %>%
  left_join(sample_data(ps))
species <- data.frame(ps_scores$species)
species$otu_id <- seq_along(colnames(otu_table(ps)))
species <- species %>%
  left_join(tax)
evals_prop <- 100 * ps_ccpna$CCA$eig[1:2] / sum(ps_ccpna$CA$eig)
ggplot() +
  geom_point(data = sites, aes(x = CCA1, y = CCA2), shape = 2, alpha = 0.5) +
  geom_point(data = species, aes(x = CCA1, y = CCA2, col = Order), size = 0.5) +
  geom_text_repel(data = species %>% filter(CCA2 < -2),
                    aes(x = CCA1, y = CCA2, label = otu_id),
		    size = 1.5, segment.size = 0.1) +
  facet_grid(. ~ family_relationship) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
        y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  scale_color_brewer(palette = "Set2") +
  coord_fixed(sqrt(ps_ccpna$CCA$eig[2] / ps_ccpna$CCA$eig[1])*0.45   ) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```  

# Supervised learning {#sec:supervised-learning .unnumbered}

Here we illustrate some supervised learning methods that can be easily
run in R. The package wraps many prediction algorithms available in R
and performs parameter tuning automatically. Since we saw that
microbiome signatures change with age, we’ll apply supervised techniques
to try to predict age from microbiome composition.

We’ll first look at Partial Least Squares (PLS)[@PLSwold]. The first
step is to divide the data into training and test sets, with assignments
done by mouse, rather than by sample, to ensure that the test set
realistically simulates the collection of new data. Once we split the
data, we can use the function `train` to fit the PLS model.

```{r caret-pls}
library(caret)
sample_data(pslog)$age2 <- cut(sample_data(pslog)$age, c(0, 100, 400))
dataMatrix <- data.frame(age = sample_data(pslog)$age2, otu_table(pslog))
# take 8 mice at random to be the training set, and the remaining 4 the test set
trainingMice <- sample(unique(sample_data(pslog)$host_subject_id), size = 8)
inTrain <- which(sample_data(pslog)$host_subject_id %in% trainingMice)
training <- dataMatrix[inTrain,]
testing <- dataMatrix[-inTrain,]
plsFit <- train(age ~ ., data = training,
                method = "pls", preProc = "center")
```



Next we can predict class labels on the test set using the function 
`predict` and
compare to the truth. We see that the method does an excellent job of
predicting age.

```{r caret-pls-confusion}
plsClasses <- predict(plsFit, newdata = testing)
table(plsClasses, testing$age)
```





As another example, we can try out random forests. This is run in
exactly the same way as PLS, by switching the argument from to . Random
forests also perform well at the prediction task on this test set,
though there are more old mice misclassified as young.

```{r caret-rf, eval =T}
library(randomForest)
rfFit <- train(age ~ ., data = training, method = "rf",
               preProc = "center", proximity = TRUE)
rfClasses <- predict(rfFit, newdata = testing)
table(rfClasses, testing$age)
```


To interpret these PLS and random forest results, it is standard to
produce biplots and proximity plots, respectively. The code below
extracts coordinates and supplies annotation for points to include on
the PLS biplot.

```{r caret-pls-scores-loadings}
pls_biplot <- list("loadings" = loadings(plsFit$finalModel),
                   "scores" = scores(plsFit$finalModel))
class(pls_biplot$scores) <- "matrix"
pls_biplot$scores <- data.frame(sample_data(pslog)[inTrain, ],
                                pls_biplot$scores)
tax <- tax_table(ps)@.Data %>%
  data.frame(stringsAsFactors = FALSE)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
class(pls_biplot$loadings) <- "matrix"
pls_biplot$loadings <- data.frame(tax, pls_biplot$loadings)
```

```{r caret-pls-scores-plot, fig.height = 7, fig.width=8, fig.wide=TRUE, fig.cap="PLS produces a biplot representation designed to separate  samples by a response variable."}
ggplot() +
  geom_point(data = pls_biplot$scores,
             aes(x = Comp.1, y = Comp.2), shape = 2) +
  geom_point(data = pls_biplot$loadings,
             aes(x = 25 * Comp.1, y = 25 * Comp.2, col = Order),
             size = 0.3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Axis1", y = "Axis2", col = "Binned Age") +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  facet_grid( ~ age2) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

The resulting biplot is displayed in Figure
\@ref(fig:caret-pls-scores-plot); it can be interpreted similarly to
earlier ordination diagrams, with the exception that the projection is
chosen with an explicit reference to the binned age variable.
Specifically, PLS identifies a subspace to maximize discrimination
between classes, and the biplot displays sample projections and ASV
coefficients with respect to this subspace.


```{r caret-proximity, fig.wide=TRUE, fig.cap="The random forest model determines a distance between samples, which can be input into PCoA to produce a proximity plot."}
rf_prox <- cmdscale(1 - rfFit$finalModel$proximity) %>%
  data.frame(sample_data(pslog)[inTrain, ])
ggplot(rf_prox) +
  geom_point(aes(x = X1, y = X2, col = age_binned),
             size = 1, alpha = 0.7) +
  scale_color_manual(values = c("#A66EB8", "#238DB5", "#748B4F")) +
  guides(col = guide_legend(override.aes = list(size = 4))) +
  labs(col = "Binned Age", x = "Axis1", y = "Axis2")
```

A random forest proximity plot is displayed in Figure
\@ref(fig:caret-proximity). To generate this representation, a distance
is calculated between samples based on how frequently sample occur in
the same tree partition in the random forest's bootstrapping
procedure. If a pair of samples frequently occur in the same
partition, the pair is assigned a low distance. The resulting
distances are then input to PCoA, giving a glimpse into the random
forests' otherwise complex classification mechanism. The separation
between classes is clear, and manually inspecting points would reveal
what types of samples are easier or harder to classify.


```{r caret-rf-importance, fig.height=6, fig.width=9, fig.cap="Random forest aids to interpretation: importance scores."}
as.vector(tax_table(ps)[which.max(importance(rfFit$finalModel)), c("Family", "Genus")])
impOtu <- as.vector(otu_table(pslog)[,which.max(importance(rfFit$finalModel))])
maxImpDF <- data.frame(sample_data(pslog), abund = impOtu)
ggplot(maxImpDF) +   geom_histogram(aes(x = abund)) +
  facet_grid(age2 ~ .) +
  labs(x = "Abundance of discriminative bacteria", y = "Number of samples")
```

To further understand the fitted random forest model, we identify the
microbe with the most influence in the random forest prediction. This
turns out to be a microbe in family _Lachnospiraceae_ and genus
_Roseburia_. Figure \@ref(fig:caret-rf-importance) plots its abundance
across samples; we see that it is uniformly very low from age 0 to 100
days and much higher from age 100 to 400 days.




Graph-based analyses {#graph-based-visualization-and-testing .unnumbered}
=====================================

Creating and plotting graphs {#creating-and-plotting-graphs .unnumbered}
----------------------------

Phyloseq has functionality for creating graphs based on thresholding a
distance matrix, and the resulting networks can be plotting using the 
`r CRANpkg("ggnetwork")` package.
This package overlays onto the `ggplot` syntax, so you can use the function
`ggplot` on an `igraph` object and add and geoms to plot the network. To be
able to color the nodes or edges a certain way, we need to add these
attributes to the igraph object. Below we create a network by
thresholding the Jaccard dissimilarity (the default distance for the
function `make_network`) at .35, and then we add an attribute to the vertices
indicating which mouse the sample came from and which litter the mouse
was in. Then we can plot the network with the coloring by mouse and
shape by litter. 

      

```{r xtrasetupknitr, echo=FALSE}
theme_set(theme_bw())
min_theme <- theme_update(panel.border = element_blank(),
                          panel.grid = element_blank(),
                          axis.ticks = element_blank(),
                          legend.title = element_text(size = 8),
                          legend.text = element_text(size = 6),
                          axis.text = element_text(size = 6),
                          axis.title = element_text(size = 8),
                          strip.background = element_blank(),
                          
                          strip.text = element_text(size = 8),
                          legend.key = element_blank())
```


```{r setup-ggnetwork}
library("phyloseqGraphTest")
library("igraph")
library("ggnetwork")
net <- make_network(ps, max.dist=0.35)
sampledata <- data.frame(sample_data(ps))
V(net)$id <- sampledata[names(V(net)), "host_subject_id"]
V(net)$litter <- sampledata[names(V(net)), "family_relationship"]
```


```{r ggnetwork-plot,fig.wide=TRUE, fig.cap="A network created by thresholding the Jaccard dissimilarity  matrix."}
ggplot(net, aes(x = x, y = y, xend = xend, yend = yend), layout = "fruchtermanreingold") +
  geom_edges(color = "darkgray") +
  geom_nodes(aes(color = id, shape = litter),  size = 3 ) +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        legend.key.height = unit(0.5,"line")) +
  guides(col = guide_legend(override.aes = list(size = .5)))
```

We see the resulting network in Figure
\@ref(fig:ggnetwork-plot). The colors in the Figure represent which mouse the sample came from and  the shape represents which litter the mouse was in.
We can see that there is grouping of the samples
by both mouse and litter.



## Graph-based two-sample tests {#graph-based-two-sample-tests .unnumbered}
Graph-based two-sample tests were introduced by Friedman and Rafsky
[@friedman1979multivariate] as a generalization of the Wald-Wolfowitz
runs test. They proposed the use of a minimum spanning tree (MST) based
on the distances between the samples, and then counting the number of
edges on the tree that were between samples in different groups. It is
not necessary to use a minimum spanning tree (MST), the graph made by
linking nearest neighbors [@schilling1986multivariate] or distance
thresholding can also be used as the input graph. No matter what graph
we build between the samples, we can approximate a null distribution by
permuting the labels of the nodes of the graph.

### Minimum Spanning Tree (MST) {#minimum-spanning-tree-mst .unnumbered}

We first perform a test using an MST with Jaccard dissimilarity. We want
to know whether the two litters (`family_relationship`) come from the same distribution.
Since there is a grouping in the data by individual (`host_subject_id`), we can’t simply
permute all the labels, we need to maintain this nested structure – this
is what the `grouping` argument does. Here we permute the labels but keep the
structure intact.

```{r mst}
gt <- graph_perm_test(ps, "family_relationship", grouping = "host_subject_id",
                      distance = "jaccard", type = "mst")
gt$pval
```

```{r mst-plot, fig.width=8.5, fig.height=5, fig.cap="The graph and permutation histogram obtained from the minimal  spanning tree with Jaccard similarity."}
plotNet1=plot_test_network(gt) + theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))
plotPerm1=plot_permutations(gt)
grid.arrange(ncol = 2,  plotNet1, plotPerm1)
```


This test has a small $p$-value, and we reject the null hypothesis that
the two samples come from the same distribution. From the plot of the
minimum spanning tree in Figure \@ref(fig:mst-plot), we see by eye that the
samples group by litter more than we would expect by chance.

### Nearest neighbors {#nearest-neighbors .unnumbered}

The $k$-nearest neighbors graph is obtained by putting an edge between
two samples whenever one of them is in the set of $k$-nearest neighbors
of the other. 

```{r knn-1}
gt <- graph_perm_test(ps, "family_relationship", grouping = "host_subject_id",
                      distance = "jaccard", type = "knn", knn = 1)
```

```{r knn-1-plot,fig.cap="k=1 nearest-neighbor network and permutation histogram", fig.wide=TRUE, fig.height=5}
plotNet2=plot_test_network(gt) + theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))
plotPerm2=plot_permutations(gt)
grid.arrange(ncol = 2,  plotNet2, plotPerm2)
```



We see from Figure \@ref(fig:knn-1-plot) that if a pair of
samples has an edge between them in the nearest neighbor graph, they are
overwhelmingly likely to be in the same litter.


Linear modeling {#linear-modeling .unnumbered}
---------------

It is often of interest to evaluate the degree to which microbial
community diversity reflects characteristics of the environment from
which it was sampled. Unlike ordination, the purpose of this analysis is
not to develop a representation of many bacteria with respect to sample
characteristics; rather, it is to describe how a single measure of
overall community structure[^1] is associated with sample
characteristics. This is a somewhat simpler statistical goal, and can be
addressed through linear modeling, for which there are a range of
approaches in R. As an example, we will used a mixed-effects model to
study the relationship between mouse microbial community diversity and
the age and litter variables that have been our focus so far. This
choice was motivated by the observation that younger mice have
noticeably lower Shannon diversities, but that different mice have
different baseline diversities. The mixed-effects model is a starting
point for formalizing this observation.

We first compute the Shannon diversity associated with each sample and
join it with sample annotation.

```{r lm-get-alpha-diversity}
library("nlme")
library("reshape2")
ps_alpha_div <- estimate_richness(ps, split = TRUE, measure = "Shannon")
ps_alpha_div$SampleID <- rownames(ps_alpha_div) %>%
  as.factor()
ps_samp <- sample_data(ps) %>%
  unclass() %>%
  data.frame() %>%
  left_join(ps_alpha_div, by = "SampleID") %>%
  melt(measure.vars = "Shannon",
       variable.name = "diversity_measure",
       value.name = "alpha_diversity")
# reorder's facet from lowest to highest diversity
diversity_means <- ps_samp %>%
  group_by(host_subject_id) %>%
  summarise(mean_div = mean(alpha_diversity)) %>%
  arrange(mean_div)
ps_samp$host_subject_id <- factor(ps_samp$host_subject_id)
#                                  diversity_means$host_subject_id)
```

```{r lm-age}
alpha_div_model <- lme(fixed = alpha_diversity ~ age_binned, data = ps_samp,
                       random = ~ 1 | host_subject_id)
```

```{r lm-prediction-intervals}
new_data <- expand.grid(host_subject_id = levels(ps_samp$host_subject_id),
                        age_binned = levels(ps_samp$age_binned))
new_data$pred <- predict(alpha_div_model, newdata = new_data)
X <- model.matrix(eval(eval(alpha_div_model$call$fixed)[-2]),
                  new_data[-ncol(new_data)])
pred_var_fixed <- diag(X %*% alpha_div_model$varFix %*% t(X))
new_data$pred_var <- pred_var_fixed + alpha_div_model$sigma ^ 2
```


```{r lm-fitted-plot}
# fitted values, with error bars
ggplot(ps_samp %>% left_join(new_data)) +
  geom_errorbar(aes(x = age_binned, ymin = pred - 2 * sqrt(pred_var),
                    ymax = pred + 2 * sqrt(pred_var)),
                col = "#858585", size = .1) +
  geom_point(aes(x = age_binned, y = alpha_diversity,
                 col = family_relationship), size = 0.8) +
  facet_wrap(~host_subject_id) +
  scale_y_continuous(limits = c(2.4, 4.6), breaks = seq(0, 5, .5)) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Binned Age", y = "Shannon Diversity", color = "Litter") +
  guides(col = guide_legend(override.aes = list(size = 4))) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)),
        axis.text.x = element_text(angle = -90, size = 6),
        axis.text.y = element_text(size = 6))
```


Hierarchical multiple testing {#hierarchical-multiple-testing .unnumbered}
-----------------------------

Hypothesis testing can be used to identify individual microbes whose
abundance relates to sample variables of interest. A standard approach
is to compute a test statistic for each bacteria individually, measuring
its association with sample characteristics, and then jointly adjust
$p$-values to ensure a False Discovery Rate upper bound. This can be
accomplished through the Benjamini-Hochberg procedure, for example
[@BH:1995]. However, this procedure does not exploit any structure among
the tested hypotheses – for example, it is likely that if one
Ruminococcus species is strongly associated with age, then others are as
well. To integrate this information,
[@benjamini2003hierarchical; @benjamini2014selective] proposed a
hierarchical testing procedure, where taxonomic groups are only tested
if higher levels are found to be be associated. In the case where many
related species have a slight signal, this pooling of information can
increase power.

We apply this method to test the association between microbial abundance
and age. This provides a complementary view of the earlier analyses,
identifying individual bacteria that are responsible for the differences
between young and old mice.

We digress briefly from hierarchical testing to describe an alternative
form of count normalization. Rather than working with the logged data as
in our earlier analysis, we consider a variance stabilizing
transformation introduced by [@LoveDESeq2] for RNA-seq data and in
[@mcmurdie2014] for 16S rRNA generated count data and available in the
`r Biocpkg("DESeq2")`
package. The two transformations yield similar sets of significant
microbes. One difference is that, after accounting for size factors, the
histogram of row sums for DESeq is more spread out in the lower values,
refer to Figure \@ref(fig:deseq-vis). This is the motivation of using such
a transformation, although for high abundance counts, it is equivalent
to the log, for lower and mid range abundances it does not crush the
data and yields more powerful results. The code below illustrates the
mechanics of computing ’s variance stabilizing transformation on a
object.

```{r deseq-transform}
library("reshape2")
library("DESeq2")
#New version of DESeq2 needs special levels
sample_data(ps)$age_binned <- cut(sample_data(ps)$age,
  				          breaks = c(0, 100, 200, 400))
levels(sample_data(ps)$age_binned) <- list(Young100="(0,100]", Mid100to200="(100,200]", Old200="(200,400]")
sample_data(ps)$family_relationship = gsub(" ", "", sample_data(ps)$family_relationship)
ps_dds <- phyloseq_to_deseq2(ps, design = ~ age_binned + family_relationship)
# geometric mean, set to zero when all coordinates are zero
geo_mean_protected <- function(x) {
  if (all(x == 0)) {
    return (0)
  }
  exp(mean(log(x[x != 0])))
}
geoMeans <- apply(counts(ps_dds), 1, geo_mean_protected)
ps_dds <- estimateSizeFactors(ps_dds, geoMeans = geoMeans)
ps_dds <- estimateDispersions(ps_dds)
abund <- getVarianceStabilizedData(ps_dds)
```


We use the `r CRANpkg("structSSI")` package to perform the hierarchical
testing [@sankaran2014structssi], we first shorten the names for each taxa/ASV.


```{r structssi-shorten-names}
short_names <- substr(rownames(abund), 1, 5)%>%
  make.names(unique = TRUE)
rownames(abund) <- short_names
```

```{r deseq-vis, fig.cap="DEseq transformation abundance"}
abund_sums <- rbind(data.frame(sum = colSums(abund),
                               sample = colnames(abund),
                               type = "DESeq2"),
                    data.frame(sum = rowSums(otu_table(pslog)),
                               sample = rownames(otu_table(pslog)),
                               type = "log(1 + x)"))
ggplot(abund_sums) +
  geom_histogram(aes(x = sum), binwidth = 20) +
  facet_grid(type ~ .) +
  xlab("Total abundance within sample")
```  


The histogram on the top gives the total DESeq2 transformed
  abundance within each sample. The bottom histogram is the same as
  that in Figure \ref@(fig:preprocessing-setup), and is included to
  facilitate comparison.

Unlike standard multiple hypothesis testing, the hierarchical testing
procedure needs univariate tests for each higher-level taxonomic
group, not just every species. A helper function,
`treePValues`, is available for this; it expects an edgelist
encoding parent-child relationships, with the first row specifying the
root node.

```{r structssi-unadjp }
library("structSSI")
el <- phy_tree(pslog)$edge
el0 <- el
el0 <- el0[nrow(el):1, ]
el_names <- c(short_names, seq_len(phy_tree(pslog)$Nnode))
el[, 1] <- el_names[el0[, 1]]
el[, 2] <- el_names[as.numeric(el0[, 2])]
unadj_p <- treePValues(el, abund, sample_data(pslog)$age_binned)
```

We can now correct $p$-value using the hierarchical testing procedure.
The test results are guaranteed to control several variants of FDR
control, but at different levels; we defer details to
[@benjamini2003hierarchical; @benjamini2014selective; @sankaran2014structssi].
```{r structssi-test}
hfdr_res <- hFDR.adjust(unadj_p, el, .75)
summary(hfdr_res)
```

```{r  structssi-test-plotres, eval=FALSE}
#interactive part: not run
plot(hfdr_res, height = 5000) # opens in a browser
```

![Screen shot](figure/structssi-screenshot.png)
Figure :A screenshot of a subtree with many differentially abundant
    bacteria, as determined by the hierarchical testing
    procedure. Currently the user is hovering over the node associated
    with bacteria GCGAG.33; this causes the adjusted $p$-value (0.0295)
    to appear.
    
The plot opens in a new browser -- a static screenshot of a subtree is
displayed above. Nodes are shaded
according to $p$-values, from blue to orange, representing the
strongest to weakest associations. Grey nodes were never
tested, to focus power on more promising subtrees. Scanning the full
tree, it becomes clear that the association between age group and
bacterial abundance is present in only a few isolated taxonomic
groups, but that it is quite strong in those groups. To give context
to these results, we can retrieve the taxonomic identity of the
rejected hypotheses.


```{r structssi-tax}
tax <- tax_table(pslog)[, c("Family", "Genus")] %>%
  data.frame()
tax$seq <- short_names
```

```{r structssi-test-res}
options(digits=3)
hfdr_res@p.vals$seq <- rownames(hfdr_res@p.vals)
tax %>%
  left_join(hfdr_res@p.vals) %>%
  arrange(adjp) %>% head(10)
```


It seems that the most strongly associated bacteria all belong to
family _Lachnospiraceae_, which is consistent with the random forest
results.


# Multitable techniques {#multitable-techniques .unnumbered}

Many microbiome studies attempt to quantify variation in the microbial,
genomic, and metabolic measurements across different experimental
conditions. As a result, it is common to perform multiple assays on the
same biological samples and ask what features – bacteria, genes, or
metabolites, for example – are associated with different sample
conditions. There are many ways to approach these questions, which to
apply depends on the study’s focus.

Here, we will focus on one specific workflow that uses sparse Canonical
Correlation Analysis (sparse CCA), a method well-suited to both
exploratory comparisons between samples and the identification of
features with interesting variation. We will use an implementation from
the `r CRANpkg("PMA")` [@witten2009pma].

Since the mouse data used above included only a single table, we use a
new data set, collected by the study [@kashyap2013genetically]. There
are two tables here, one for bacteria and another with metabolites. 12
samples were obtained, each with measurements at 637 m/z values and
20,609 OTUs; however, about 96% of the entries of the microbial
abundance table are exactly zero. The code below retrieves 
and filters these data.


```{r multitable-setup}
metab <- read.csv("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/metabolites.csv",row.names = 1)
microbe_connect <-url("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/microbe.rda")
load(microbe_connect)
microbe
```

We see that `microbe` is a `phyloseq` object.
Our preprocessing mirrors that done for the mouse data. We first filter
down to microbes and metabolites of interest, removing those that are
zero across many samples. Then, we transform them to weaken the heavy
tails.
We also take the log of the metabolites.

```{r multitable-filtering}
library("genefilter")
keep_ix <- rowSums(metab == 0) <= 3
metab <- metab[keep_ix, ]
microbe <- prune_taxa(taxa_sums(microbe) > 4, microbe)
microbe <- filter_taxa(microbe, filterfun(kOverA(3, 2)), TRUE)
metab <- log(1 + metab, base = 10)
X <- otu_table(microbe)
X[X > 50] <- 50
dim(X)
dim(metab)
```
We see that both X and metab have 12 columns, these are actually
the samples and we will transpose them.

We can now apply sparse CCA. This method compares sets of features
across high-dimensional data tables, where there may be more measured
features than samples. In the process, it chooses a subset of available
features that capture the most covariance – these are the features that
reflect signals present across multiple tables. We then apply PCA to
this selected subset of features. In this sense, we use sparse CCA as a
screening procedure, rather than as an ordination method.

Our implementation is below. The parameters `penaltyx` and `penaltyz` are sparsity penalties. Smaller values of penaltyx will result in fewer selected microbes, similarly penaltyz modulates the number of selected metabolites.
See the documentation of the `CCA` in the `r CRANpkg("PMA")` package.

We tune them manually to
facilitate subsequent interpretation – we generally prefer more sparsity
than the default parameters would provide.


```{r multitable-sparse-cca}
library(PMA)
cca_res <- CCA(t(X),  t(metab), penaltyx = .15, penaltyz = .15)
cca_res
```

With these parameters, 5 microbes and 15 metabolites have been selected,
based on their ability to explain covariation between tables. Further,
these 20 features result in a correlation of 0.974 between the two
tables. We interpret this to mean that the microbial and metabolomic
data reflect similar underlying signals, and that these signals can be
approximated well by the 20 selected features. Be wary of the
correlation value, however, since the scores are far from the usual
bivariate normal cloud. Further, note that it is possible that other
subsets of features could explain the data just as well – sparse CCA has
minimized redundancy across features, but makes no guarantee that these
are the “true” features in any sense.

Nonetheless, we can still use these 20 features to compress information
from the two tables without much loss. To relate the recovered
metabolites and OTUs to characteristics of the samples on which they
were measured, we use them as input to an ordinary PCA.

```{r loadade4,echo=FALSE}
library(ade4)
```

```{r multitable-plug-in-pca }
combined <- cbind(t(X[cca_res$u != 0, ]),
                  t(metab[cca_res$v != 0, ]))
pca_res <- dudi.pca(combined, scannf = F, nf = 3)
```

```{r annotation}
genotype <- substr(rownames(pca_res$li), 1, 2)
sample_type <- substr(rownames(pca_res$l1), 3, 4)
feature_type <- grepl("\\.", colnames(combined))
feature_type <- ifelse(feature_type, "Metabolite", "OTU")
sample_info <- data.frame(pca_res$li, genotype, sample_type)
feature_info <- data.frame(pca_res$c1,
                           feature = substr(colnames(combined), 1, 6))
```

```{r multitable-interpret-pca, fig.wide=TRUE, fig.cap="A PCA triplot produced from the CCA selected features in from muliple data types;metabolites and OTUs."}
ggplot() +  geom_point(data = sample_info,
            aes(x = Axis1, y = Axis2, col = sample_type, shape = genotype), size = 3) + 
  geom_label_repel(data = feature_info,
                   aes(x = 5.5 * CS1, y = 5.5 * CS2, label = feature, fill = feature_type),
                   size = 2, segment.size = 0.3,
                   label.padding = unit(0.1, "lines"), label.size = 0) +
  geom_point(data = feature_info,
             aes(x = 5.5 * CS1, y = 5.5 * CS2, fill = feature_type),
             size = 1, shape = 23, col = "#383838") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#a6d854", "#e78ac3")) +
  guides(fill = guide_legend(override.aes = list(shape = 32, size = 0))) +
  coord_fixed(sqrt(pca_res$eig[2] / pca_res$eig[2])) +
  labs(x = sprintf("Axis1 [%s%% Variance]",
                   100 * round(pca_res$eig[1] / sum(pca_res$eig), 2)),
       y = sprintf("Axis2 [%s%% Variance]",
                   100 * round(pca_res$eig[2] / sum(pca_res$eig), 2)),
       fill = "Feature Type", col = "Sample Type")
```       

Note that we have departed from our convention of fixing the aspect ratio 
here as the second axis
represents very little of the variability and the plot would actually become unreadable.

Figure \@ref(fig:multitable-interpret-pca) displays a PCA *triplot*,
where we show different types of samples and the multidomain features
(Metabolites and OTUs). This allows comparison across the measured
samples – triangles for Knockout and circles for wild type – and
characterizes the influence the different features – diamonds with text
labels. For example, we see that the main variation in the data is
across PD and ST samples,
which correspond to the different diets.

Further, large values of 15 of the features are associated with ST
status, while small values for 5 of them indicate PD status. The
advantage of the sparse CCA screening is now clear – we can display most
of the variation across samples using a relatively simple plot, and can
avoid plotting the hundreds of additional points that would be needed to
display all of the features.


Conclusions {#conclusions .unnumbered}
===========

We have shown how a complete workflow in R is now available to
denoise, identify and normalize next generation amplicon sequencing
reads using probabilistic models with parameters fit using the data at
hand.

We have provided a brief overview of all the analyses that become
possible once the data has been imported into the _R_ environment.
Multivariate projections using the phylogenetic tree as the relevant
distance between OTUs/ASVs can be done using weighted unifrac or double
principal coordinate analyses using the
[*phyloseq*](http://www.bioconductor.org/packages/release/bioc/html/phyloseq.html)
package. Biplots provide the user with an interpretation key. These
biplots have been extended to triplots in the case of multidomain data
incorporating genetic, metabolic and taxa abundances. We illustrate the
use of network based analyses, whether the community graph is provided
from other sources or from a taxa co-occurrence computation using a
Jaccard distance.

We have briefly covered a small example of using supervised
learning functions (random forests and partial least squares) to
predict a response variable,

The main challenges in tackling microbiome data come from the many
different levels of heterogeneity both at the input and output levels.
These are easily accommodated through R’s capacity to combine data into
S4 classes. We are able to include layers of information, trees, sample
data description matrices, contingency table in the phyloseq data
structures. The plotting facilities of and allow for the layering of
information in the output into plots that combine graphs, multivariate
information and maps of the relationships between covariates and taxa
abundances. The layering concept allows the user to provide reproducible
publication level figures with multiple heterogeneous sources of
information. Our main goal in providing these tools has been to enhance
the statistical power of the analyses by enabling the user to combine
frequencies, quality scores and covariate information into complete and
testable projections.

Summary {#summary .unnumbered}
=======

This illustration of possible workflows for microbiome data combining
trees, networks, normalized read counts and sample information showcases
the capabilities and reproducibility of an R based system for analysing
bacterial communities. We have implemented key components in C
wrapped within the Bioconductor package
[*dada2*](http://www.bioconductor.org/packages/release/bioc/html/dada2.html)
to enable the different steps to be undertaken on a laptop.

Once the sequences have been filtered and tagged they can be assembled
into a phylogenetic tree directly in R using the maximum likelihood tree
estimation available in `r `. The sequences are then assembled into a
phyloseq object containing all the sample covariates, the phylogenetic
tree and the sample-taxa contingency table.

These data can then be visualized interactively with Shiny-phyloseq,
plotted with one line wrappers in phyloseq and filtered or transformed
very easily.

The last part of the paper shows more complex analyses that require
direct plotting and advanced statistical analyses.

Multivariate ordination methods allow useful lower dimensional
projections in the presence of phylogenetic information or multi-domain
data as shown in an example combining metabolites, OTU abundances and sample covariate information.


Supervised learning methods provide lists of the most relevant taxa in
discriminating between groups.

Bacterial communities can be represented as co-occurrence graphs using
network based plotting procedures available in R. We have also provided
examples where these graphs can be used to test community structure
through non parametric permutation resampling. This provides
implementations of the Friedman Rafsky[@friedman1979multivariate] tests
for microbiome data which have not been published previously.

Data availability {#data-availability .unnumbered}
=================

Intermediary data for the analyses are made available both on GitHub at
<https://github.com/spholmes/F1000_workflow> and at the Stanford digital
repository permanent url for this paper:\
<http://purl.stanford.edu/wh250nn9648>. All other data have been
previously published and the links are included in the paper.

Operation  and Session Information {#operation .unnumbered}
-----------------------------------

The programs and source for this article can be run using version
[3.1] or above of [R](https://cran.r-project.org) and version [3.3] or above of
[Bioconductor](https://www.bioconductor.org/install/).

