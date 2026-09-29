library(WGCNA)
library(rrcov)
library(ggplot2) 
library(ComplexHeatmap)
library(clusterProfiler)
source("Analysis/data.r")
options(stringsAsFactors = FALSE)
plot_sft <- function(sft, filename = 'Network-Topology-Analysis.png') {
  library(magrittr)
  library(patchwork)
  p1 <- sft$fitIndices %>% 
    ggplot(aes(x = Power, y = -sign(slope)*SFT.R.sq, label = Power)) +
    geom_text() +
    geom_abline(aes(color = 'red', slope = 0, intercept = 0.85)) +
    ylim(c(-0.75, 1)) + 
    xlab('Soft Threshold (power)') + 
    ylab('Scale Free Topology Model Fit, signed R^2') +
    labs(title = paste('Scale independence. Beta =', sft$powerEstimate)) +
    theme_minimal() +
    theme(legend.position = 'none')

  p2 <- sft$fitIndices %>% 
    ggplot(aes(x = Power, y =mean.k., label = Power)) +
    geom_text() +
    xlab('Soft Threshold (power)') + 
    ylab('Mean Connectivity') +
    labs(title = 'Mean Connectivity') +
    theme_minimal()

  if (!is.null(filename)) {
    cowplot::save_plot(plot = p1 + p2, filename = filename, base_asp=2)
  }
  return(p1 + p2)
}

## wrapper for plotDendroAndColors
plotDendroAndColors.wrapper <- function(net, path=NULL, prefix='Dendrogram', dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05, ...) {
  if (!is.null(path)) { pdf(paste0(path, '/', prefix, '.pdf'), wi=8, he=6) }
  mergedColors <- labels2colors(net$colors)
  WGCNA::plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors", dendroLabels = dendroLabels, hang = hang, addGuide = addGuide, guideHang=guideHang, ...)
  if (!is.null(path)) { dev.off() }
}

## Wrapper for module-trait correlations.
# net  : network
# meta : data frame of traits you want to correlate with. Ensure that rownames can re-arrange net$MEs
#      into the correct order.

plotModuleTraitCorrelation.wrapper <- function(datExpr, datTraits, net, name, path='./', cor_method = c("pearson", "kendall", "spearman"), ...) {
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  moduleColors = labels2colors(net$colors)
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  write.table(MEs, paste0(path, '/', 'MEs - ', name, '.txt'), sep='\t', row.names=T, col.names=T, quote=F)
  moduleTraitCor = cor(MEs, as.matrix(datTraits), use="p", method=cor_method)
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  pdf(paste0(path, '/', 'MTC - ', name, '.pdf'), wi=9, he=11+0.1*(ncol(net$MEs)-20)*(ncol(net$MEs) > 20))
  textMatrix=paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep="")
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar=c(6,8.5+0.1*(ncol(net$MEs)-20)*(ncol(net$MEs) > 20),3,3))
#  par(mar=c(1,1,1,1))
  labeledHeatmap(Matrix=moduleTraitCor,
                 xLabels=names(datTraits),
                 yLabels=names(MEs),
                 ySymbols=names(MEs),
                 colorLabels=FALSE,
                 colors=blueWhiteRed(50),
                 textMatrix=textMatrix,
                 setStdMargins=FALSE,
                 cex.text=0.5,
                 zlim=c(-1,1),
                 main=paste("Module-trait relationships -", name))
  dev.off()
  
  
  return(list(datExpr=datExpr,
              datTraits=datTraits,
              MEs=MEs,
              moduleTraitCor=moduleTraitCor,
              moduleTraitPvalue=moduleTraitPvalue,
              textMatrix=textMatrix))
}
xl <- load_lipid(demo=FALSE)
xp <- load_protein(demo=FALSE)

pca <- PcaClassic(t(assay(xl, 'quantile')), k=4)
pca <- PcaClassic(t(assay(xp, 'log10')[rowSums(is.na(assay(xp))) == 0,]), k=4)
pdf(); plot(pca); dev.off()

pdf()
pca$scores %>%
  as.data.frame %>%
  tibble::rownames_to_column('Sample') %>%
  ggplot(aes(x=PC1, y=PC2)) + 
  geom_point()
dev.off()

# ------------ Lipidomics WGCNA ----------------
fit <- list()
resids <- matrix(NA, nrow=nrow(xl), ncol=ncol(xl))
rownames(resids) = rownames(xl)
colnames(resids) = colnames(xl)
for (gg in rownames(xl)) {
  y <- data.frame(gg = assay(xl, 'log10')[gg,])
  fit[[gg]] <- lm(t(as.matrix(y)) ~ lipid.sum, data = data.frame(lipid.sum = colSums(assay(xl, 'raw'), na.rm=TRUE)), na.action = na.omit)
  to_add = resid(fit[[gg]])
  resids[gg, names(to_add)] = to_add
}

assay(xl, 'resids') <- resids
assay='resids'
powers <- c(c(1:10), seq(from=12, to=20, by=2))
sft_lipid <- WGCNA::pickSoftThreshold(t(assay(xl, assay)), powerVector=powers, verbose=5)
plot_sft(sft_lipid, filename = file.path('Results-WGCNA/lipid', assay, 'lipid Network-Topology-Analysis.png'))

enableWGCNAThreads()
net_lipid = blockwiseModules(
  t(assay(xl, assay)),
  power = 5,
  TOMType = 'signed',
  networkType = 'signed hybrid',
  minModuleSize = 2, # The minimum module size should be that of the smallest lipid class -- minus the free cholesterol.
  reassignThreshold = 0,
  mergeCutHeight = 0.6, # To reduce the total number of modules we should increase height. My protocol is to just increase this number in 0.1 increments until I get a number of modules less than the number of lipid classes we started with.
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  maxBlockSize = 3000,
  corType = 'bicor',
  nThreads = 20,
  deepSplit = 2,
  verbose = 3
)

pdf(file.path('Results-WGCNA/lipid/', assay, 'lipid dendro.pdf'), width=10, height=7)
#plotDendroAndColors.wrapper(net_lipid, path=NULL, prefix=NA)
#dev.off()

datTraits = WGCNA::binarizeCategoricalColumns.pairwise(data = as.data.frame(colData(xl)), convertColumns = c('Genotype', 'Sex'))[, c('Genotype.E4.vs.E3', 'Sex.M.vs.F')]
datTraits <- cbind(datTraits, lipid.sum = colSums(log10(assay(xl, 'raw'))))
lipid_result <- plotModuleTraitCorrelation.wrapper(t(assay(xl, assay)), datTraits = datTraits, net=net_lipid, name = 'Lipid', path=file.path('Results-WGCNA/lipid/', assay), cor_method = 'spearman')
writexl::write_xlsx(
  purrr::map(
    lipid_result, 
    ~as.data.frame(.x, row.names=rownames(.x)) %>% 
      tibble::rownames_to_column(var="rowname") %>%
      mutate(rowname = stringr::str_remove(rowname, "^ME"))
    ), paste0("Results-WGCNA/lipid/", assay, "/lipid-ME-correlation-results.xlsx"))

moduleColors = labels2colors(net_lipid$colors)
MEs0 <- moduleEigengenes(t(assay(xl, assay)), moduleColors)$eigengenes
MEs = orderMEs(MEs0)
geneModuleMembership = as.data.frame(cor(t(assay(xl, assay)), MEs, use = "p", method='spearman'))
geneModuleMembership$modules = moduleColors
write.csv(geneModuleMembership, file.path('Results-WGCNA/lipid', assay, 'lipid module membership.csv'))

geneModuleMembership$modules <- NULL
module_colors <- as.vector(gsub(pattern="ME", x=colnames(MEs), replacement=""))
names(module_colors) <- colnames(MEs)
p <- bind_cols(
  Sex = xl$Sex,
  Genotype = xl$Genotype,
  MEs
) %>%
  tidyr::pivot_longer(-c(Sex, Genotype), names_to = 'module', values_to='ME') %>%
  ggplot(aes(x = Genotype, y=ME)) +
  geom_boxplot(aes(fill=module), outlier.size = 0) + 
  geom_jitter(width=0.1, shape=21, fill='white') + 
  facet_wrap(~module) + 
  scale_fill_manual(values=module_colors)
cowplot::save_plot(
  plot = p,
  filename = file.path('Results-WGCNA/lipid/', assay, 'lipid ME boxplots by genotype.pdf'),
  base_height=8.5,
  base_width=11
)

module_colors <- as.vector(gsub(pattern="ME", x=colnames(geneModuleMembership), replacement=""))
names(module_colors) <- colnames(geneModuleMembership)
lipids.keep <- rowSums(is.na(geneModuleMembership)) == 0
hc <- columnAnnotation(
    name = 'Module',
    df = data.frame(module = colnames(geneModuleMembership)),
    col = list(module=module_colors)
  )

source('Analysis/colors.r')
hr <- rowAnnotation(
    name = 'Lipid Class',
    df = data.frame(LipidClass = rowData(xl)[lipids.keep,]$lipid.class),
    col = list(LipidClass=colors$lipid.class)
  )

h <- Heatmap(
  as.matrix(geneModuleMembership[lipids.keep,]),
  top_annotation = hc,
  left_annotation = hr, 
  row_names_gp  = gpar(fontsize=8)
  )
pdf(file.path('Results-WGCNA/lipid/', assay, 'lipid membership.pdf'), width=10, height=25); h; dev.off()

# lipid class gsea
gsea_result <- list()
TERM2GENE = as.data.frame(rowData(xl)) %>%
  dplyr::select(lipid.class, lipid.species) %>%
  dplyr::rename(term=lipid.class, gene=lipid.species)
TERM2NAME = as.data.frame(rowData(xl)) %>%
  dplyr::select(lipid.class, lipid.class) %>%
  dplyr::rename(term=lipid.class)
for (module in colnames(geneModuleMembership)) {
  geneList <- geneModuleMembership[, module]
  names(geneList) <- rownames(geneModuleMembership)
  gsea_result[['lipid']][[module]] <- GSEA(
    sort(geneList, decreasing = TRUE),
    TERM2GENE = TERM2GENE, 
    eps=0)
}
lipid_gsea <- purrr::map(gsea_result$lipid, ~as.data.frame(.x@result))
writexl::write_xlsx(lipid_gsea, file.path('Results-WGCNA/lipid', assay, 'lipid-module-enrichment-results.xlsx'))

# ----- Proteomics WGCNA ---------
fit <- list()
resids <- matrix(NA, nrow=nrow(xp), ncol=ncol(xp))
rownames(resids) = rownames(xp)
colnames(resids) = colnames(xp)
for (gg in rownames(xp)) {
  y <- data.frame(gg = assay(xp, 'log10')[gg,])
  fit[[gg]] <- lm(t(as.matrix(y)) ~ protein.sum, data = data.frame(protein.sum = colSums(assay(xp, 'raw'), na.rm=TRUE)), na.action = na.omit)
  to_add = resid(fit[[gg]])
  resids[gg, names(to_add)] = to_add
}

assay(xp, 'resids') <- resids
assay='resids'
powers <- c(c(1:10), seq(from=12, to=20, by=2))
sft_protein <- WGCNA::pickSoftThreshold(t(assay(xp, assay)), powerVector=powers, verbose=5)
plot_sft(sft_protein, filename = paste0('Results-WGCNA/protein/', assay, '/protein Network-Topology-Analysis.png'))
net_protein = blockwiseModules(
  t(assay(xp, assay)),
  power = 10,
  TOMType = 'signed',
  networkType = 'signed hybrid',
  minModuleSize = 20, # The minimum module size should be that of the smallest protein class.
  reassignThreshold = 0,
  mergeCutHeight = 0.4, # It's hard to say how many modules is a good amount. Papers typically arrive at a number of modules in the dozens.
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  maxBlockSize = 5000,
  corType = 'bicor',
  nThreads = 20,
  deepSplit = 2,
  verbose = 3
)

pdf(paste0('Results-WGCNA/protein/', assay, '/protein dendro.pdf'), width=10, height=7)
#plotDendroAndColors.wrapper(net_protein, path=NULL, prefix=NA)
#dev.off()

datTraits = WGCNA::binarizeCategoricalColumns.pairwise(data = as.data.frame(colData(xp)), convertColumns = c('Genotype', 'Sex'))[, c('Genotype.E4.vs.E3', 'Sex.M.vs.F')]
datTraits <- cbind(datTraits, protein.sum = log10(colSums(assay(xp, 'raw'), na.rm=TRUE)))
protein_result = plotModuleTraitCorrelation.wrapper(t(assay(xp, assay)), datTraits = datTraits, net=net_protein, name = 'Protein', path=paste0('Results-WGCNA/protein/', assay), cor_method = 'spearman')
writexl::write_xlsx(purrr::map(protein_result, ~as.data.frame(.x, row.names = rownames(.x))), paste0("Results-WGCNA/protein/", assay, "/protein-ME-correlation-results.xlsx"))
writexl::write_xlsx(
  purrr::map(
    protein_result, 
    ~as.data.frame(.x, row.names=rownames(.x)) %>% 
      tibble::rownames_to_column(var="rowname") %>%
      mutate(rowname = stringr::str_remove(rowname, "^ME"))
    ), paste0("Results-WGCNA/protein/", assay, "/protein-ME-correlation-results.xlsx"))

moduleColors = labels2colors(net_protein$colors)
MEs0 <- moduleEigengenes(t(assay(xp, assay)), moduleColors)$eigengenes
MEs = orderMEs(MEs0)
geneModuleMembership = as.data.frame(cor(t(assay(xp, assay)), MEs, use = "p", method='spearman'))
geneModuleMembership$modules = moduleColors
write.csv(geneModuleMembership, file.path('Results-WGCNA/protein', assay, 'protein module membership.csv'))

module_colors <- as.vector(gsub(pattern="ME", x=colnames(MEs), replacement=""))
names(module_colors) <- colnames(MEs)
p <- bind_cols(
  Sex = xp$Sex,
  Genotype = xp$Genotype,
  MEs
) %>%
  tidyr::pivot_longer(-c(Sex, Genotype), names_to = 'module', values_to='ME') %>%
  ggplot(aes(x = Genotype, y=ME)) +
  geom_boxplot(aes(fill=module), outlier.size = 0) + 
  geom_jitter(width=0.1, shape=21, fill='white') + 
  facet_wrap(~module) + 
  scale_fill_manual(values=module_colors)
cowplot::save_plot(
  plot = p,
  filename = file.path('Results-WGCNA/protein/', assay, 'protein ME boxplots by genotype.pdf'),
  base_height=8.5,
  base_width=11
)

geneModuleMembership$modules = NULL
module_colors <- as.vector(gsub(pattern="ME", x=colnames(geneModuleMembership), replacement=""))
names(module_colors) <- colnames(geneModuleMembership)
ha <- columnAnnotation(
    name = 'Module',
    df = data.frame(module = colnames(geneModuleMembership)),
    col = list(module=module_colors)
  )

h <- Heatmap(
  as.matrix(geneModuleMembership),
  top_annotation = ha,
  row_names_gp  = gpar(fontsize=8)
  )
pdf(paste0('Results-WGCNA/protein/', assay, '/protein membership.pdf'), width=10, height=30); h; dev.off()

library(org.Hs.eg.db)
all_proteins = select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns='UNIPROT')
for (module in colnames(geneModuleMembership)) {
  geneList <- geneModuleMembership[, module]
  names(geneList) <- rownames(geneModuleMembership)
  geneList = geneList[names(geneList) %in% all_proteins$UNIPROT]
  gsea_result[['protein']][[module]] <- tryCatch(
    {
      gseGO(
        sort(geneList, decreasing = TRUE),
        keyType = 'UNIPROT',
        ont = 'ALL',
        OrgDb = org.Hs.eg.db,
        eps=0)
    }, error = function(cond) {
      gseGO(
        sort(geneList, decreasing = TRUE),
        keyType = 'UNIPROT',
        ont = 'BP',
        OrgDb = org.Hs.eg.db,
        eps=0)
    })
}
prot_gsea <- purrr::map(gsea_result$protein, ~as.data.frame(.x@result))
writexl::write_xlsx(prot_gsea, file.path('Results-WGCNA/protein', assay, 'protein-module-enrichment-results.xlsx'))

# -------- Plot correlation of protein and lipid modules -------
MEs_lipid = read.table("Results-WGCNA/lipid/resids/MEs - Lipid.txt")
MEs_protein = read.table("Results-WGCNA/protein/resids/MEs - Protein.txt")
ME_cor <- cor(MEs_lipid, MEs_protein, method='spearman')
pME_cor <- corPvalueStudent(ME_cor, nSamples=40)

module_colors <- as.vector(gsub(pattern="ME", x=colnames(ME_cor), replacement=""))
names(module_colors) <- colnames(ME_cor)
hc <- columnAnnotation(
  name = 'Module',
  df = data.frame(`Protein MEs` = colnames(ME_cor), check.names = FALSE),
  col = list(`Protein MEs` = module_colors)
)

module_colors <- as.vector(gsub(pattern="ME", x=rownames(ME_cor), replacement=""))
names(module_colors) <- rownames(ME_cor)
hr <- rowAnnotation(
  name = 'Module',
  df = data.frame(`Lipid MEs` = rownames(ME_cor), check.names=FALSE),
  col = list(`Lipid MEs` = module_colors)
)

h <- Heatmap(
  as.matrix(ME_cor),
  top_annotation = hc,
  left_annotation = hr,
  row_title = 'Lipid modules',
  column_title = 'Protein modules',
  heatmap_legend_param = list(title = "pearson's cor"),
  #row_names_gp  = gpar(fontsize=8)
)
pdf(file.path('Results-WGCNA/lipid-protein', assay, 'ME_cor.pdf'), width=10, height=7); h; dev.off()

write.table(ME_cor, file.path('Results-WGCNA/lipid-protein', assay, 'ME_cor.txt'), sep=',', row.names = TRUE, col.names = TRUE)
write.table(pME_cor, file.path('Results-WGCNA/lipid-protein', assay, 'ME_cor p-values.txt'), sep=',', row.names = TRUE, col.names = TRUE)

# ------ Genotype comparisons -----
E3 <- xp$Genotype == 'E3'
E4 <- xp$Genotype == 'E4'

#pdf(width=6, height=5)
#bind_cols(
#  Sex = xp$Sex,
#  Genotype = xp$Genotype,
#  MEs_lipid %>% dplyr::rename_with(~paste0('lipid_', .x)),
#  MEs_protein %>% dplyr::rename_with(~paste0('protein_', .x))
#) %>%
#  ggplot(aes(x = lipid_MEtan, y=protein_MEgreen)) + 
#  geom_point(aes(color=Sex)) + 
#  #facet_wrap(~Genotype) + 
#  geom_smooth(method = 'lm')
#dev.off()


