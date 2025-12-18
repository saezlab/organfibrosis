# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

library(MOFAcellulaR)
library(tidyverse)
library(cowplot)
# Main ---------------------------------------------------------------
files = snakemake@input$pbulk
coldata_files = snakemake@input$coldata
metadata_files = snakemake@input$meta
marker_files = snakemake@input$markers

source(snakemake@params$aesthestics_R)
mofa_dir = snakemake@params$mofa_dir
mofa_out = snakemake@output$file
umap_pdf = snakemake@output$umap_pdf
mofa_pdf = snakemake@output$mofa_pdf
stats_pdf = snakemake@output$stats_pdf
factor_pdf = snakemake@output$factor_pdf

studies <- sapply(strsplit(files, "/"), function(x) {
  parts <- unlist(x)
  parts[length(parts)]
}, USE.NAMES = FALSE) %>% gsub(".csv","",.)


quality_plots = paste0(mofa_dir, studies, '_model_filters.pdf')

input_df <- data.frame(
  file = I(files),
  study = I(studies),
  coldata_file = I(coldata_files),
  metadata_file = I(metadata_files),
  marker_csv = I(marker_files),
  quality_pdf = I(quality_plots),
  check.names = FALSE
)


# We get all meta_data
all_meta <- tibble(file = metadata_files) %>%
  dplyr::mutate(metadata = map(file, ~ read_csv(.x, show_col_types = FALSE)[,-1])) %>%
  unnest(metadata) %>%
  dplyr::select(-file) %>%
  mutate(condition = cond_test)

all_samples <- all_meta$sample %>%
  unique() %>%
  length()
print(all_samples)

study_samples <- all_meta %>%
  dplyr::select(study, sample) %>%
  unique() %>%
  group_by(study) %>%
  summarize(n_samples = n())

plot_genes_samples <- function(sum_experiment, title){
  vectorgenes <- c()
  vectorsamples <- c()
  vectorctype <- c() 
  for (i in sum_experiment) {
    genes = dim(i)[1]
    samples = dim(i)[2]
    ctype = i$annotation_MOFA[1]
    vectorgenes <- c(vectorgenes, genes)
    vectorsamples <- c(vectorsamples, samples)
    vectorctype <- c(vectorctype, ctype)
  }
  data <- data.frame(genes = vectorgenes, samples = vectorsamples, ctype = vectorctype)
  
  plot <- ggplot(data, aes(x = vectorgenes, y = vectorsamples, color = ctype)) +
    geom_point() +
    theme_minimal()+
    scale_x_continuous(limits = c(0, NA)) +  # Set x-axis limits
    scale_y_continuous(limits = c(0, NA)) + 
    xlab("genes") +
    ylab("samples") +
    ggtitle(title)
  
  return(plot)
}


get_ct_list <- function (pb_dat, cts, ct_col = "cell_type") 
{
  cts <- purrr::set_names(cts)
  pb_dat_list <- purrr::map(cts, function(ctype) {
    ix <- base::which(SummarizedExperiment::colData(pb_dat)[,ct_col] == ctype)
    return(pb_dat[, ix])
  })
  return(pb_dat_list)
}

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


# We get all meta_data
all_meta <- tibble(file = metadata_files) %>%
  dplyr::mutate(metadata = map(file, ~ read_csv(.x, show_col_types = FALSE)[,-1])) %>%
  unnest(metadata) %>%
  dplyr::select(-file) %>%
  mutate(heart_failure = cond_test)

all_samples <- all_meta$sample %>%
  unique() %>%
  length()

study_samples <- all_meta %>%
  dplyr::select(study, sample) %>%
  unique() %>%
  group_by(study) %>%
  summarize(n_samples = n())

# First we need to get all markers together [Union]
all_mrks <- input_df %>%
  dplyr::select(study, marker_csv) %>%
  dplyr::mutate(markers = map(marker_csv, read_csv,
                              show_col_types = FALSE)) %>%
  dplyr::select(-marker_csv) %>%
  unnest(c(markers)) %>%
  dplyr::filter(FDR < 0.01, logFC > 2) %>%
  dplyr::rename("annotation_MOFA" = name) %>%
  dplyr::select(annotation_MOFA, gene) %>%
  unique() %>%
  group_by(annotation_MOFA) %>%
  nest() %>%
  dplyr::mutate(data = map(data, ~.x[[1]])) %>%
  deframe()

# We need all pbulks together
# Processing will be handled for each study individually

input_df_pb <- input_df %>%
  dplyr::select(-marker_csv)


all_pbs <- pmap(input_df_pb, function(file, study, coldata_file, metadata_file, quality_pdf) {
  
  print(study)
  
  # Importing pb data
  pb_data <- read_csv(file,
                      show_col_types = FALSE)
  
  colnames(pb_data)[1] <- "sample_id"
  
  pb_data <- pb_data %>%
    column_to_rownames("sample_id") %>%
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
  
  # Defining cts
  cts <- coldat$annotation_MOFA %>%
    unique() %>%
    set_names()
  
  # Defining parameters
  n_samples <- nrow(meta_data)
  min_samples <- (n_samples * 0.4) %>% floor()
  
  pb_obj <- MOFAcellulaR::create_init_exp(pb_data, coldat) 
  pdf(quality_pdf, height = 4, width = 5)
  
  ct_list <- get_ct_list(pb_dat = pb_obj,cts = cts,ct_col = "annotation_MOFA")
  plot <-plot_genes_samples(ct_list, "no filter")
  plot(plot)
  
  ct_list <- MOFAcellulaR::filt_profiles(pb_dat = pb_obj,
                                         cts = cts,
                                         ncells = 25,
                                         counts_col = "ncells",
                                         ct_col = "annotation_MOFA") # This refers to the column name in testcoldata where the cell-type label was stored
  
  
  plot <- plot_genes_samples(ct_list, "filt_profiles\nviews with low number of cells")
  plot(plot)
  
  ct_list <- MOFAcellulaR::filt_views_bysamples(pb_dat_list = ct_list,
                                                nsamples = min_samples)
  plot <- plot_genes_samples(ct_list,"filt_views_bysamples\nviews with low number of samples")
  plot(plot)
  
  ct_list <- MOFAcellulaR::filt_gex_byexpr(pb_dat_list = ct_list,
                                           min.count = 100,
                                           min.prop = 0.4) 
  plot <- plot_genes_samples(ct_list,"filt_gex_byexpr\nlowly expressed genes")
  plot(plot)
  
  ct_list <- MOFAcellulaR::filt_views_bygenes(pb_dat_list = ct_list,
                                              ngenes = 50)
  plot <- plot_genes_samples(ct_list,"filt_views_bygenes\nviews with low number of genes")
  plot(plot)
  
  ct_list <- MOFAcellulaR::filt_samples_bycov(pb_dat_list = ct_list,
                                              prop_coverage = 0.97)
  plot <- plot_genes_samples(ct_list,"filt_samples_bycov\nsamples from views with a low coverage")
  plot(plot)
  
  ct_list <- tmm_trns(pb_dat_list = ct_list,
                        scale_factor = 1000000)
  # only normalization step, no need for new plot
  
  ct_list <- MOFAcellulaR::filt_gex_bybckgrnd(pb_dat_list = ct_list,
                                              prior_mrks = all_mrks)
  plot <- plot_genes_samples(ct_list,"filt_gex_bybckgrnd\nmarker genes")
  plot(plot)
  
  ct_list <- MOFAcellulaR::filt_views_bygenes(pb_dat_list = ct_list,
                                              ngenes = 50)
  
  plot <- plot_genes_samples(ct_list,"filt_views_bygenes \nviews with low number of genes")
  plot(plot)
  
  
  dev.off()
  MOFAcell_obj <- pb_dat2MOFA(pb_dat_list = ct_list, 
                              sample_column = "sample")%>%
    plyr::mutate(group = study)
  
})

all_pbs <- bind_rows(all_pbs)

# Model completion
# Generating QC of feature and sample coverage
view_completion <- all_pbs %>%
  group_by(view) %>%
  summarise(n_genes = length(unique(feature)),
            n_ind = length(unique(sample))) %>%
  dplyr::mutate(perc_samples = (n_ind/all_samples) * 100)

view_completion_study <- all_pbs %>%
  group_by(group, view) %>%
  summarise(n_genes = length(unique(feature)),
            n_ind = length(unique(sample)))  %>%
  left_join(study_samples, by = c("group" = "study")) %>%
  dplyr::mutate(perc_samples = (n_ind/n_samples) * 100)

bar_plt0 <- ggplot(view_completion_study, aes(x = view, y = perc_samples)) +
  geom_bar(stat = "identity") +
  facet_wrap(.~group, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("") +
  ylab("Percentage of samples")

bar_plt1 <- ggplot(view_completion_study, aes(x = view, y = n_genes)) +
  geom_bar(stat = "identity") +
  facet_wrap(.~group, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("") +
  ylab("Number of genes")


pdf(stats_pdf, height = 3, width = 3.9)
plot(bar_plt0)
plot(bar_plt1)
dev.off()




# Fitting the MOFA model

# Building MOFA model

MOFAobject <- create_mofa(all_pbs)

data_opts <- get_default_data_options(MOFAobject)

data_opts$center_groups <- TRUE

model_opts <- get_default_model_options(MOFAobject)

model_opts$num_factors <- 15

model_opts$spikeslab_weights <- FALSE

train_opts <- get_default_training_options(MOFAobject)

# Prepare MOFA model:
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# Train model:
model <- run_mofa(MOFAobject, mofa_out)


n_samples <- sum(model@dimensions$N)
# Variance Explained
var_explained <- model@cache$variance_explained$r2_total %>%
  map(., ~ enframe(.x, "view","R2")) %>%
  enframe("study") %>%
  unnest(c(value))

# Association with HF (all studies)
# Association with HF
expl_var_HF <- MOFAcellulaR::get_associations(model = model,
                                              metadata = all_meta,
                                              sample_id_column = "sample",
                                              #test_variable = "heart_failure",
                                              test_variable = "cond_test",
                                              test_type = "categorical",
                                              categorical_type = "parametric",
                                              group = TRUE)

# Association per study

studies <- set_names(all_meta$study %>% unique)

assoc_list <- map(studies, function(s) {
  
  MOFAcellulaR::get_associations(model = model,
                                 metadata = all_meta%>%
                                   dplyr::filter(study == s),
                                 sample_id_column = "sample",
                                 test_variable = "cond_test",
                                 test_type = "categorical",
                                 categorical_type = "parametric",
                                 group = TRUE)%>%
    replace(is.na(.), 1)
  
  
  
})

assoc_list[["All"]] <- expl_var_HF

# Visualize
col_list <- list(heart_failure = c("control" = "darkgrey",
                                   "fibrosis" = "black"),
                 cond_test = etiology_colors[all_meta$cond_test %>% unique],
                 study = study_colors[all_meta$study %>% unique],
                 grouping = etiology_colors[all_meta$grouping %>% unique],
                 organ = organ_colors[all_meta$organ %>% unique])

scores_hmap <- MOFAcellulaR::plot_MOFA_hmap(model = model,
                                            group = TRUE,
                                            metadata = all_meta,
                                            sample_id_column = "sample",
                                            sample_anns = c("cond_test","study",'grouping', 'organ'),
                                            assoc_list = assoc_list,
                                            col_rows = col_list)
pdf(mofa_pdf, height = (length(files) *3), width = 7)

draw(scores_hmap)

dev.off()

# UMAP


# Generate UMAP
UMAP_dat <- MOFAcellulaR::plot_sample_2D(model,
                                         metadata = all_meta,
                                         color_by = "grouping",
                                         sample_id_column = "sample",
                                         group = T,
                                         n_neighbors = ceiling(n_samples * 0.25),
                                         seed = 1)
UMAP_plt0 <- last_plot()

UMAP_dat <- MOFAcellulaR::plot_sample_2D(model,
                                         metadata = all_meta,
                                         color_by = "study",
                                         sample_id_column = "sample",
                                         group = T,
                                         n_neighbors = ceiling(n_samples * 0.25),
                                         seed = 1)
UMAP_plt1 <- last_plot()
UMAP_plt1 <-UMAP_plt1 +
  scale_color_manual(values = study_colors[all_meta$study %>% unique])

# Generate UMAP
UMAP_dat <- MOFAcellulaR::plot_sample_2D(model,
                                         metadata = all_meta,
                                         color_by = "cond_test",
                                         sample_id_column = "sample",
                                         group = T,
                                         n_neighbors = ceiling(n_samples * 0.25),
                                         seed = 1)

UMAP_plt2 <- last_plot()

UMAP_plt2 <- UMAP_plt2 +
  scale_color_manual(values = etiology_colors[all_meta$cond_test %>% unique])

# Generate UMAP
UMAP_dat <- MOFAcellulaR::plot_sample_2D(model,
                                         metadata = all_meta,
                                         color_by = "organ",
                                         sample_id_column = "sample",
                                         group = T,
                                         n_neighbors = ceiling(n_samples * 0.25),
                                         seed = 1)

UMAP_plt3 <- last_plot()

UMAP_plt3 <- UMAP_plt3 +
  scale_color_manual(values = organ_colors[all_meta$organ %>% unique])




pl <- align_plots(UMAP_plt0, UMAP_plt1, UMAP_plt2, UMAP_plt3, align="v")

pdf(umap_pdf, height = 3, width = 4.9)



plot(ggdraw(pl[[1]]))
plot(ggdraw(pl[[2]]))
plot(ggdraw(pl[[3]]))
plot(ggdraw(pl[[4]]))


dev.off()



factors <- MOFA2::get_factors(model, factors = "all") %>%
  base::do.call(base::rbind, .) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column()

factors_meta <- factors %>%
  dplyr::rename("sample" = rowname) %>%
  dplyr::left_join(all_meta, by = "sample") %>%
  dplyr::rename("condition" = heart_failure)



# extract which Factors are the three most significant ones
smallest_indices <- order(assoc_list$All$adj_pvalue)[1:3]
sig_factors <- assoc_list$All$Factor[smallest_indices]


print(factor_pdf)
# plot the most significant factor against the other two
pdf(factor_pdf, height = 3, width = 4.9)

for (fac in sig_factors[2:3]){
  plot1 <- ggplot(factors_meta, aes(x = .data[[sig_factors[1]]], y = .data[[fac]], color = condition))+
    geom_point()+
    labs(x = sig_factors[1], y = fac, color = "condition") +
    theme_minimal()+
    scale_color_manual(values = etiology_colors[all_meta$cond_test %>% unique])
  

  plot2 <- ggplot(factors_meta, aes(x = .data[[sig_factors[1]]], y = .data[[fac]], color = study))+
    geom_point()+
    labs(x = sig_factors[1], y = fac, color = "study") +
    theme_minimal()+
    scale_color_manual(values = study_colors[all_meta$study %>% unique])
  
  plot3 <- ggplot(factors_meta, aes(x = .data[[sig_factors[1]]], y = .data[[fac]], color = organ))+
    geom_point()+
    labs(x = sig_factors[1], y = fac, color = "organ") +
    theme_minimal()+
    scale_color_manual(values = organ_colors[all_meta$organ %>% unique])
  
  plot4 <- ggplot(factors_meta, aes(x = .data[[sig_factors[1]]], y = .data[[fac]], color = grouping))+
    geom_point()+
    labs(x = sig_factors[1], y = fac, color = "organ") +
    theme_minimal()+
    scale_color_manual(values = etiology_colors[all_meta$grouping %>% unique])
  
  pl <- align_plots(plot1, plot2, plot3, plot4, align="v")
  plot(ggdraw(pl[[1]]))
  plot(ggdraw(pl[[2]]))
  plot(ggdraw(pl[[3]]))
  plot(ggdraw(pl[[4]]))
}
dev.off()


