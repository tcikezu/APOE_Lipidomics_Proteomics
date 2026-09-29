# With this script we can load lipidomics and proteomics data.
# Let's do this as follows.
# - load_lipid()
# - load_protein()
# Both will return a summarized experiment
# Sanity check passed 10/13/2023: `all(colData(load_lipid()) == colData(load_protein())) == TRUE`

load_lipid <- function(demo=TRUE) {
  library(SummarizedExperiment)
  library(magrittr)
  library(dplyr)

  if (demo) {
    input.file <- "Data/lipids_demo.xlsx"
  } else {
    input.file <- 'Data/lipidomics-173-female and male.xlsx'
  }
  sample_info <- t(readxl::read_excel(input.file, n_max = 1)) %>%
    as.data.frame %>%
    tibble::rownames_to_column('rowname')
  colnames(sample_info) <- sample_info[1,]
  sample_info <- sample_info[2:nrow(sample_info),]
  sample_info <- sample_info %>%
    mutate( across(everything(), ~stringr::str_replace_all(.x, '-', '_')) )
  rownames(sample_info) <- sample_info$`ID number`
  sample_info <- sample_info %>%
    tidyr::separate(Group, sep = '_', into=c('Genotype', 'Sex'), remove = FALSE)
  
  data <- readxl::read_excel(input.file, skip=1, .name_repair = 'minimal') %>%
    tibble::column_to_rownames('Group')
  colnames(data) = rownames(sample_info)
  
  row_info <- data.frame(
    lipid.species = rownames(data)
  ) %>%
    mutate(
      lipid.class = ifelse(
        lipid.species != 'free',
        stringr::str_extract(lipid.species, '^[A-Z]+'),
        'free')
    )
  
  warn_data(data)
  
  are_zero <- rowSums(data <= 0) > 0
  data <- data[!are_zero, ]
  row_info <- row_info[!are_zero, ]
  
  x <- SummarizedExperiment(
    assays = list(raw = data,
                  log10 = log10(data),
                  log10_composition = log10(data / colSums(data, na.rm=TRUE)),
                  quantile = limma::normalizeBetweenArrays(data, method='quantile')
                  ),
    colData = sample_info,
    rowData = row_info
  )
  return(x)
}

load_protein <- function(resid = FALSE, demo = TRUE) {
  library(SummarizedExperiment)
  library(magrittr)
  library(dplyr)
  library(limma)

  if (demo) {
    input.file <- 'Data/proteins_demo.xlsx'
  } else {
    input.file <- 'Data/Proteomics-APOE-EV-female and male.xlsx'
  }

  sample_info <- readxl::read_excel(input.file, n_max = 2, .name_repair = 'minimal', col_names = FALSE)
  sample_info <- t(sample_info[,3:ncol(sample_info)])
  colnames(sample_info) <- sample_info[1,]
  colnames(sample_info)[2] = 'ID number'
  sample_info <- sample_info[2:nrow(sample_info),]
  sample_info <- sample_info %>%
    as.data.frame(row.names = FALSE) %>%
    mutate(Group = stringr::str_replace(Group, '-', '_')) %>%
    tidyr::separate(Group, sep = '_', into=c('Genotype', 'Sex'), remove = FALSE) %>%
    mutate(`ID number` = stringr::str_remove(`ID number`, '^APO')) %>%
    dplyr::relocate(`ID number`)
  rownames(sample_info) <- sample_info[['ID number']]
  
  row_info <- readxl::read_excel(input.file, skip=1) %>%
    dplyr::select(Accession.Number, Genes, Protein.Description) %>%
    tidyr::separate(Accession.Number, sep='\\|', into = c('db', 'UniqueIdentifier', 'EntryName'), remove = FALSE) %>%
    data.frame(row.names = .[['UniqueIdentifier']])
  
  data <- readxl::read_excel(input.file, skip=1) %>%
    dplyr::select(!(Accession.Number : Protein.Description)) %>%
    data.frame(row.names=row_info$UniqueIdentifier)
  colnames(data) <- stringr::str_remove(colnames(data), '^APO')
  
  warn_data(data)
  
  are_zero <- rowSums(data <= 0, na.rm = TRUE) > 0
  data <- data[!are_zero, ]
  row_info <- row_info[!are_zero, ]
  
  x <- SummarizedExperiment(
   assays = list(raw = data,
                 log10 = log10(data),
                 log10_composition = log10(data / colSums(data, na.rm=TRUE)),
                 quantile = limma::normalizeBetweenArrays(data, method='quantile')
                 ),
    colData = sample_info,
    rowData = row_info
  )
  
  # "limma" normalized
  #assay(x, 'limmaVSN') <- limma::normalizeVSN(assay(x, 'raw'))
  #meanSdPlot(as.matrix(assay(x, 'raw')))
  
  #limma::plotMD(assay(x, 'log10_composition'))
  #limma::plotSA(assay(x, 'raw'))
  #
  #assay(x, 'log10') %>%
  #  as.data.frame %>%
  #  tibble::rownames_to_column('Gene') %>%
  #  tidyr::pivot_longer(-Gene) %>%
  #  mutate(genotype = stringr::str_split_i(name, '_', 1)) %>%
  #  group_by(name) %>%
  #  ggplot(aes(x = genotype, y=value)) + 
  #  geom_boxplot()
  #
  #data.frame(
  #  sd = rowSds(as.matrix(assay(x, 'log10')), na.rm = TRUE),
  #  mean = rowMeans(assay(x, 'log10'), na.rm=TRUE)
  #) %>%
  #  mutate(
  #    mean_rank = rank(mean)
  #  ) %>%
  #  ggplot(aes(x=mean_rank, y=sd)) + 
  #  geom_point() #+ 
  #  #scale_x_continuous(trans='log') +
  #  #scale_y_continuous(trans='log')
  #
  return(x)
}

warn_data <- function(data) {
  if (any(rowSums(data < 0, na.rm = TRUE) > 0)) {
    rlog::log_warn(paste(rownames(data)[rowSums(data < 0, na.rm=TRUE) > 0], 'has negative values.'))
    rlog::log_warn(paste(sum(rowSums(data < 0) > 0), 'features have negative values.'))
  }
  if (any(rowSums(data == 0, na.rm=TRUE) > 0)) {
    rlog::log_warn(paste(paste0(rownames(data)[rowSums(data == 0, na.rm=TRUE) > 0], collapse=','), 'has zero values.'))
    rlog::log_warn(paste(sum(rowSums(data == 0, na.rm=TRUE) > 0), 'features have zero values.'))
  }
  if (any(rowSums(is.na(data)) > 0)) {
    rlog::log_warn(paste(paste0(rownames(data)[rowSums(is.na(data)) > 0], collapse=','), 'has NA values.'))
    rlog::log_warn(paste(sum(rowSums(is.na(data)) > 0), 'features have NA values.'))
  }
}

load_corr <- function(x, assay='log10', method='spearman') {
  corr = rcorr(x = t(assay(x$lipid, assay)),
               y = t(assay(x$protein, assay)),
               type = method)
  
  corr$n <- corr$n[rownames(x$lipid), rownames(x$protein)]
  corr$r <- corr$r[rownames(x$lipid), rownames(x$protein)]
  corr$P <- corr$P[rownames(x$lipid), rownames(x$protein)]
  xn = corr$n
  
  corr$n <- corr$n[rowSums(xn == 0) < ncol(xn), colSums(xn == 0) < nrow(xn)]
  corr$r <- corr$r[rowSums(xn == 0) < ncol(xn), colSums(xn == 0) < nrow(xn)]
  corr$P <- corr$P[rowSums(xn == 0) < ncol(xn), colSums(xn == 0) < nrow(xn)]
  return(corr)
}

compute_resid <- function(x) {
  resids <- matrix(0, nrow=nrow(x), ncol=ncol(x))
  rownames(resids) = rownames(x)
  colnames(resids) = colnames(x)
  for (p in rownames(x)) {
    fit <- lm(t(as.matrix(assay(x, 'log10')[p,])) ~ Group,
            data=as.data.frame(colData(x)), na.action = na.exclude)
    resids[p, names(resid(fit))] = resid(fit)
  }
  assay(x, 'resid') = resids
  return(x)
}
