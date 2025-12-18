# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we run MOFAcell models for each of the available datasets:
library(MOFAcellulaR)
library(tidyverse)
library(cowplot)

tmm_trns <- function(pb_dat_list, scale_factor = 1000000) {
  
  pb_dat_red <- purrr::map(pb_dat_list, function(x) {
    if (dim(x)[2] != 0){
      all_nf <- edgeR::calcNormFactors(x, method = "TMM")
      sfs <- all_nf$samples$lib.size * all_nf$samples$norm.factors
      pb <- base::sweep(assay(x, "counts"), MARGIN = 2, sfs, FUN = "/")
      SummarizedExperiment::assay(x, "logcounts") <- base::log1p(pb * scale_factor)
      
      return(x)}
    else{
      return(NULL)
    }
  })
  # Remove NULL elements from the list
  pb_dat_red <- purrr::compact(pb_dat_red)
  return(pb_dat_red)
  
}


files = snakemake@input$pbulk
coldata_files = snakemake@input$coldata
metadata_files = snakemake@input$meta
marker_files = snakemake@input$markers

source(snakemake@params$aesthestics_R)
mofa_dir = snakemake@params$mofa_dir
model_stats_plt_file = snakemake@output$plot
model_stats_file = snakemake@output$file

studies <- sapply(strsplit(files, "/"), function(x) {
  parts <- unlist(x)
  parts[length(parts)]
}, USE.NAMES = FALSE) %>% gsub(".csv","",.)

mofa_pdf = paste0(mofa_dir, studies, '_mofa.pdf')
mofa_out = paste0(mofa_dir, studies, '_mofa.hdf5')
umap_pdf = paste0(mofa_dir, studies, '_umap.pdf')
factor_pdf = paste0(mofa_dir, studies, '_factors.pdf')


input_df <- data.frame(
  file = I(files),
  study = I(studies),
  coldata_file = I(coldata_files),
  metadata_file = I(metadata_files),
  marker_csv = I(marker_files),
  mofa_out = mofa_out,
  mofa_pdf = mofa_pdf,
  umap_pdf = umap_pdf,
  factor_pdf = factor_pdf,
  check.names = FALSE
)

model_outs <- pmap(input_df, function(study, file, 
                                      coldata_file, metadata_file, 
                                      marker_csv, mofa_out, 
                                      mofa_pdf, umap_pdf, factor_pdf) {
  
  print(study)
  print(file)
  print(coldata_file)
  print(marker_csv)
  print(metadata_file)
  print(mofa_pdf)
  print(factor_pdf)
  
  # Importing pb data
  pb_data <- read_csv(file,
                      show_col_types = FALSE)
  
  colnames(pb_data)[1] <- "sample"
  
  pb_data <- pb_data %>%
    column_to_rownames("sample") %>%
    as.matrix() %>%
    t()
  
  # Importing coldata of the matrices - ignoring not annotated cell types
  coldat <- read_csv(coldata_file,
                     show_col_types = FALSE)[,-1]  %>%
    column_to_rownames("colname") %>%
    dplyr::rename(ncells = "counts") %>%
    dplyr::filter(annotation_MOFA != "none")
  
  pb_data <- pb_data[,rownames(coldat)]
  
  # Importing meta_data
  meta_data <- read_csv(metadata_file,
                        show_col_types = FALSE)[,-1]
  # depending on no. of samples in study, diff. number of factors will be used
  n_samples <- nrow(meta_data)
  if (n_samples <= 6) { 
    factor_nr <- 4
  } else if ((n_samples > 6) & (n_samples <= 10)) {
    factor_nr <- 6
  } else if ((n_samples > 10) & (n_samples <= 14)) {
    factor_nr <- 8
  } else if (n_samples > 14) {
    factor_nr <- 10
  } 
  
  
  # Importing markers
  mrkrs <- read_csv(marker_csv, show_col_types = FALSE) %>%
    dplyr::filter(FDR < 0.01, logFC > 1) %>%
    dplyr::rename("annotation_MOFA" = name) %>%
    group_by(annotation_MOFA) %>%
    nest() %>%
    dplyr::mutate(data = map(data, ~.x[[1]])) %>%
    deframe()
  
  # Defining cts
  cts <- coldat$annotation_MOFA %>% 
    unique() %>%
    set_names()
  print(cts)
  
  # Defining parameters
  min_samples <- (n_samples * 0.4) %>% floor()
  
  # Creating MOFAcell object to make manipulations
  MOFAcell_obj <- MOFAcellulaR::create_init_exp(pb_data, coldat) %>%
    MOFAcellulaR::filt_profiles(pb_dat = .,
                                cts = cts,
                                ncells = 25,
                                counts_col = "ncells",
                                ct_col = "annotation_MOFA") %>%
    MOFAcellulaR::filt_views_bysamples(pb_dat_list = .,
                                       nsamples = min_samples) %>%
    MOFAcellulaR::filt_gex_byexpr(pb_dat_list = .,
                                  min.count = 100,
                                  min.prop = 0.4) %>%
    MOFAcellulaR::filt_views_bygenes(pb_dat_list = .,
                                     ngenes = 50) %>%
    MOFAcellulaR::filt_samples_bycov(pb_dat_list = ., # Filtering of low quality samples
                                     prop_coverage = 0.97) %>%
    tmm_trns(pb_dat_list = .,
              scale_factor = 1000000) %>% 
    MOFAcellulaR::filt_gex_bybckgrnd(pb_dat_list = .,
                                     prior_mrks = mrkrs) %>%
    MOFAcellulaR::filt_views_bygenes(pb_dat_list = .,
                                     ngenes = 50) %>%
    MOFAcellulaR::pb_dat2MOFA(pb_dat_list = .,
                              sample_column = "sample")
  
  # Generating QC of feature and sample coverage
  view_completion <- MOFAcell_obj %>%
    group_by(view) %>%
    summarise(n_genes = length(unique(feature)),
              n_ind = length(unique(sample))) %>%
    dplyr::mutate(perc_samples = (n_ind/n_samples) * 100)
  
  
  # Building MOFA model
  
  MOFAobject <- create_mofa(MOFAcell_obj)
  
  data_opts <- get_default_data_options(MOFAobject)
  
  data_opts$center_groups <- TRUE
  
  model_opts <- get_default_model_options(MOFAobject)
  
  model_opts$num_factors <- factor_nr
  
  model_opts$spikeslab_weights <- FALSE
  
  train_opts <- get_default_training_options(MOFAobject)
  
  print(MOFAobject)  
  # Prepare MOFA model:
  MOFAobject <- prepare_mofa(
    object = MOFAobject,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
  
  

  # Train model:
  model <- run_mofa(MOFAobject, mofa_out)
  
  # Summarize model
  
  # R2 - of model
  var_explained <- model@cache$variance_explained$r2_total[[1]] %>%
    enframe("view","R2")
  
  
  print('model:')
  print(model)
  # Association with HF
  expl_var_HF <- MOFAcellulaR::get_associations(model = model,
                                  metadata = meta_data,
                                  sample_id_column = "sample",
                                  test_variable = "cond_test",
                                  test_type = "categorical",
                                  categorical_type = "parametric",
                                  group = FALSE)
  
  # Complete summary table
  useful_factors <- expl_var_HF %>%
    dplyr::filter(adj_pvalue < 0.01) %>%
    pull(Factor)
  
  if(rlang::is_empty(useful_factors)) {
    var_explained <- var_explained %>%
      dplyr::mutate(R2_HF = 0)
    
  } else {
    R2_HF <- model@cache$variance_explained$r2_per_factor[[1]][useful_factors,, drop = F] %>%
      colSums() %>%
      enframe("view","R2_HF")
    
    var_explained <- var_explained %>%
      left_join(R2_HF, by = "view")
    
  }
  
  model_summ <- left_join(view_completion, 
                          var_explained, 
                          by = "view") %>%
    dplyr::mutate(study = study)
  
  # Visualize
  assoc_list <- list(disease = expl_var_HF)
  
  
  col_list <- list(cond_test = etiology_colors[meta_data$cond_test %>% unique],
                   study = study_colors[meta_data$study %>% unique],
                   grouping = etiology_colors[meta_data$grouping %>% unique])
  scores_hmap <- MOFAcellulaR::plot_MOFA_hmap(model = model,
                                              group = FALSE,
                                              metadata = meta_data,
                                              sample_id_column = "sample",
                                              sample_anns = c("cond_test","study", "grouping"),
                                              assoc_list = assoc_list,
                                              col_rows = col_list)
  
  
  pdf(mofa_pdf, height = 8, width = 4.5)
  draw(scores_hmap)
  dev.off()
  
  # Generate UMAP
  UMAP_dat <- MOFAcellulaR::plot_sample_2D(model, 
                                           metadata = meta_data, 
                                           color_by = "cond_test",
                                           sample_id_column = "sample",
                                           n_neighbors = factor_nr, seed = 1)
  
  UMAP_plt <- last_plot()
  
  
  # Generate UMAP
  UMAP_dat <- MOFAcellulaR::plot_sample_2D(model, 
                                           metadata = meta_data, 
                                           color_by = "grouping",
                                           sample_id_column = "sample",
                                           n_neighbors = factor_nr, seed = 1)
  
  UMAP_plt2 <- last_plot()
  
  
  pdf(umap_pdf, height = 3, width = 3.9)
  plot(UMAP_plt +
         scale_color_manual(values = etiology_colors[meta_data$cond_test %>% unique])  +
         ggtitle(study) )
  plot(UMAP_plt2+
         scale_color_manual(values = etiology_colors[meta_data$grouping %>% unique])  +
         ggtitle(study))
  dev.off()
  
  
  
  
  # extract which Factors are the three most significant ones
  smallest_indices <- order(expl_var_HF$adj_pvalue)[1:3]
  sig_factors <- expl_var_HF$Factor[smallest_indices]
  
  factors <- MOFA2::get_factors(model, factors = "all") %>%
    base::do.call(base::rbind, .) %>% 
    as_tibble(rownames = NA) %>% 
    rownames_to_column()
  
  factors_meta <- factors %>%
    dplyr::rename("sample" = rowname) %>%
    dplyr::left_join(meta_data, by = "sample")
  
  
  print(factor_pdf)
  # plot the most significant factor against the other two
  pdf(factor_pdf, height = 3, width = 4.9)
  
  for (fac in sig_factors[2:3]){
    plot1 <- ggplot(factors_meta, aes(x = .data[[sig_factors[1]]], y = .data[[fac]], color = cond_test))+
      geom_point()+
      labs(x = sig_factors[1], y = fac, color = "condition") +
      theme_minimal()+
      scale_color_manual(values = etiology_colors[meta_data$cond_test %>% unique])
    
    plot4 <- ggplot(factors_meta, aes(x = .data[[sig_factors[1]]], y = .data[[fac]], color = grouping))+
      geom_point()+
      labs(x = sig_factors[1], y = fac, color = "disease") +
      theme_minimal()+
      scale_color_manual(values = etiology_colors[meta_data$grouping %>% unique])
    
    pl <- align_plots(plot1, plot4, align="v")
    plot(ggdraw(pl[[1]]))
    plot(ggdraw(pl[[2]]))
  }
  dev.off()
  
  
  
  
  
  return(model_summ)
  
})

# Generate final map 

model_stats <- enframe(model_outs) %>%
  unnest(c(value)) %>%
  dplyr::select(-name) %>%
  dplyr::mutate(R2 = as.numeric(R2))

write_csv(model_stats, model_stats_file)

R2_qc <- model_stats %>%
  ggplot(aes(x = study, y = R2)) +
  geom_bar(stat = "identity") +
  facet_wrap(. ~ view, nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("") +
  ylab("Model R2")

genes_qc <- model_stats %>%
  ggplot(aes(x = study, y = n_genes)) +
  geom_bar(stat = "identity") +
  facet_wrap(. ~ view, nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("") +
  ylab("Number of genes")

sample_qc <- model_stats %>%
  ggplot(aes(x = study, y = perc_samples)) +
  geom_bar(stat = "identity") +
  facet_wrap(. ~ view, nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("") +
  ylab("% profiled samples")

R2_hf_qc <- model_stats %>%
  ggplot(aes(x = study, y = R2_HF)) +
  geom_bar(stat = "identity") +
  facet_wrap(. ~ view, nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("") +
  ylab("Model R2 disease")


model_stats_plt <- cowplot::plot_grid(genes_qc, sample_qc, R2_qc, R2_hf_qc, ncol = 2, nrow = 2, align = "hv")

pdf(model_stats_plt_file, height = 8, width = 10)

plot(model_stats_plt)

dev.off()





