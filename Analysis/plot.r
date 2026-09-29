source("Analysis/colors.r")

plot_density <- function(x, type, assay, file.dir=NULL) {
  p <- assay(x[[type]], assay) %>%
    tibble::rownames_to_column('UniqueIdentifier') %>%
    tidyr::pivot_longer(-UniqueIdentifier, names_to = 'ID Number', values_to = assay) %>%
    ggplot(aes(x = .data[[assay]], group=`ID Number`)) +
    geom_density() + 
    labs(title = paste(type, assay, 'distribution by sample'))
  
  if (!is.null(file.dir)) {
    cowplot::save_plot(
      plot = p,
      filename = file.path(file.dir, paste(type, assay, 'density plot.pdf')),
      base_asp = 1.3,
      base_height=5
    )
  } 
  return(p)
}
plot_pca_biplot <- function(x, type, assay, npcs=4, file.dir=NULL) {
  remove <- rowSums(is.na(assay(x[[type]], assay))) > 0 | rowSums(apply(assay(x[[type]], assay), MARGIN = c(1,2), is.infinite)) > 0
  keep = !remove
  pca <- prcomp(t(na.omit(assay(x[[type]], assay)[keep, ])), center = TRUE, scale = TRUE)
  #pca <- PcaClassic(t(na.omit(assay(x[[type]], assay)[keep, ])), k=npcs, scale=TRUE, signflip = FALSE)
  pca_props = pca$sdev^2 / sum(pca$sdev^2)
  p <- bind_cols(
    pca$x,
    colData(x[[type]]) %>% as_tibble,
    total = colSums(na.omit(assay(x[[type]], assay)[keep, ]))
  ) %>%
    ggplot()
  
  pca_biplots <- list()
  for (pc in 1:(npcs-1)) {
    pc_x = paste0('PC', pc)
    pc_y = paste0('PC', pc+1)
    x_lab = paste0(pc_x, ' (', round(pca_props[pc]*100, digits=1), '%)')
    y_lab = paste0(pc_y, ' (', round(pca_props[pc+1]*100, digits=1), '%)')
    for (color in c('Group', 'Genotype', 'Sex')) {
      pca_biplots[[paste0(pc, '_', color)]] <- p + geom_point(aes(x = .data[[pc_x]], y=.data[[pc_y]], fill=.data[[color]]), shape=21, color='black', size=3) + 
        xlab(x_lab) + 
        ylab(y_lab) + 
        scale_fill_manual(values = colors[[color]])
    }
    color = 'total'
    pca_biplots[[paste0(pc, '_', color)]] <- p + geom_point(aes(x = .data[[pc_x]], y=.data[[pc_y]], fill=.data[[color]]), shape=21, color='black', size=3) + 
      xlab(x_lab) + 
      ylab(y_lab) + 
      scale_fill_viridis_c(option='C')
  }
  p <- wrap_plots(pca_biplots, nrow=3)
  if (!is.null(file.dir)) {
    cowplot::save_plot(
      plot = p,
      filename = file.path(file.dir, paste(type, assay, 'PCA biplots.pdf')),
      base_height=10,
      base_asp=1.6
    )
  } 
  return(p)
}
