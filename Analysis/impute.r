# Imputation strategies
# 1. Pairwise correlation followed by rank or scale matching
# 2. Nearest-neighbors followed by rank or scale matching

# 1. Pairwise correlation followed by quantile matching
imputation_strategy_1 <- function(data, k=5, cor.threshold=0.5) {
  cor_data <- cor(data, method = 'spearman', use = 'pairwise.complete.obs')
  na_features = colnames(data)[colSums(is.na(data)) > 0]
  imputed_data = data
  
  for (ff in na_features) {
    scores_topk = sort(cor_data[ff,], decreasing=TRUE)[2:ncol(data)][1:k]
    scores_topk = scores_topk[scores_topk > cor.threshold]
    
    if (length(scores_topk) == 0) {
      next
    }
    
    impute_features = names(scores_topk)[1:length(scores_topk)]
    ranks = apply(data[,c(ff, impute_features)], MARGIN=2, FUN=rank, na.last = 'keep')
    scales = scale(data[, c(ff, impute_features)])
    
    #data[,c(ff, impute_features)] %>%
    ranks %>%  
      as.data.frame %>%
      tibble::rownames_to_column('ID Number') %>%
      tidyr::pivot_longer(-`ID Number`, names_to='Unique.Identifier', values_to='Protein') %>%
      group_by(`ID Number`) %>%
      mutate(has_NA = sum(is.na(Protein)) > 0) %>%
      mutate(Unique.Identifier = factor(Unique.Identifier, levels=c(ff, impute_features))) %>%
      ggplot(aes(x = Unique.Identifier, y=Protein, group=`ID Number`, label = `ID Number`, color=has_NA)) + 
      geom_point(na.rm=FALSE) + 
      geom_path() + 
      geom_text() + 
      facet_wrap(~`ID Number` == 'E3_22')
      #facet_wrap(~Unique.Identifier)
    
    
    mean_rank = stats::weighted.mean(ranks[is.na(imputed_data[,ff]),-1], w=scores_topk^2 / sum(scores_topk^2), na.rm=TRUE)
    median_rank = median(ranks[is.na(imputed_data[,ff]),-1], na.rm=TRUE)
    mean_scale = stats::weighted.mean(scales[is.na(imputed_data[,ff]),-1], w=scores_topk^2 / sum(scores_topk^2), na.rm=TRUE)
    
    mean_rank_impute = stats::weighted.mean(c(data[ranks[,ff] == floor(mean_rank) & !is.na(ranks[,ff]),ff], data[ranks[,ff] == ceiling(mean_rank) & !is.na(ranks[,ff]),ff]),
                         w = c(1 - (mean_rank - floor(mean_rank)),
                               mean_rank - floor(mean_rank))
    )
    
    median_imputed = data[ranks[,ff] == median_rank & !is.na(ranks[,ff]),ff]
   
    slow =  max(scales[scales[,ff] < mean_scale,ff], na.rm=TRUE)
    shigh = min(scales[scales[,ff] >= mean_scale, ff], na.rm=TRUE)
    xlow =  data[scales[,ff] == slow & !is.na(scales[,ff]), ff]
    xhigh = data[scales[,ff] == shigh & !is.na(scales[,ff]), ff]
    scale_imputed = stats::weighted.mean(
      x = c(xlow, xhigh),
      w = 1 - c((mean_scale - slow)/(shigh - slow), (shigh - mean_scale)/(shigh-slow))
    )
    
    for (gg in c()) {
      if (sum(is.na(imputed_data[,ff])) == 0) {
        break
      }
      else {
        rank(data[,gg], na.last = 'keep')[is.na(imputed_data[,ff])]
        
      }
    }
  }
  
  
}