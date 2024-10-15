# library(DoddLabDatabase)
# annotate_metabolite(ms1_file = 'Peak_Table.csv',
#                     ms2_type = 'mzML',
#                     path = '~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow',
#                     polarity = 'positive',
#                     lib = 'dodd',
#                     column = 'hilic',
#                     ce = '20',
#                     adduct_list = '[M+H]+',
#                     mz_tol = 15,
#                     mz_ppm_thr = 150,
#                     pf_rt_range = 10,
#                     tolerance_rt_range = 20,
#                     is_rt_score = TRUE,
#                     is_ms2_score = TRUE)
#
# load('~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow/05_object_hilic_pos_outlier_removal.RData')
# annotate_metabolite(object = object_hilic_pos,
#                     ms2_type = 'mzML',
#                     path = '~/Project/00_IBD_project/Data/20230227_develop_metabolite_ID_workflow',
#                     polarity = 'positive',
#                     lib = 'dodd',
#                     column = 'hilic',
#                     ce = '20',
#                     adduct_list = '[M+H]+',
#                     mz_tol = 15,
#                     mz_ppm_thr = 150,
#                     pf_rt_range = 10,
#                     tolerance_rt_range = 20,
#                     is_rt_score = TRUE,
#                     is_ms2_score = TRUE)
#
#
#
# annotate_metabolite(ms1_file = 'data.csv',
#                     ms2_type = 'mzML',
#                     path = '~/Project/00_IBD_project/Data/20230327_raw_data_processing_test/DemoData/',
#                     polarity = 'positive',
#                     lib = 'dodd',
#                     column = 'hilic',
#                     ce = '20',
#                     adduct_list = '[M+H]+',
#                     mz_tol = 15,
#                     mz_ppm_thr = 150,
#                     pf_rt_range = 10,
#                     tolerance_rt_range = 20,
#                     is_rt_score = TRUE,
#                     is_ms2_score = TRUE)
#
#
# load(file.path(path_output, "00_intermediate_data", 'ms2'))
# load(file.path(path_output, "00_intermediate_data", 'ms1_data'))
# load(file.path(path_output, "00_intermediate_data", 'ms1_result'))
# load(file.path(path_output, "00_intermediate_data", 'lib_meta'))
# load(file.path(path_output, "00_intermediate_data", 'lib_spec'))
# load(file.path(path_output, "00_intermediate_data", 'ms2_result'))
# load(file.path(path_output, "00_intermediate_data", 'exp_spec'))
#
# path <- '~/Project/00_IBD_project/Data/20230327_raw_data_processing_test/DemoData/'
# merge_one_modes(path = path,
#                 column = 'hilic',
#                 polarity = 'positive')


################################################################################
# library(DoddLabDatabase)
# parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20230327_raw_data_processing_test/DemoData/',
#                                                                   lib = 'dodd',
#                                                                   ce = '20',
#                                                                   column = 'hilic',
#                                                                   polarity = 'positive',
#                                                                   is_rt_score = TRUE,
#                                                                   is_ms2_score = TRUE)
# annotate_metabolite(parameter_set_annotation = parameter_set_annotation)

# load('~/Project/00_IBD_project/Data/20230327_raw_data_processing_test/DemoData/01_metabolite_annotation/00_intermediate_data/result_annotation')
# load('~/Project/00_IBD_project/Data/20230327_raw_data_processing_test/DemoData/01_metabolite_annotation/00_intermediate_data/ms1_data')


################################################################################
# library(DoddLabDatabase)

# parameter_set_annotation@
# parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/05_coworker/Aedan/230727_untargeted_data_processing/HILIC_pos/',
#                                                                   lib = 'peptide',
#                                                                   ce = '20',
#                                                                   column = 'hilic',
#                                                                   polarity = 'positive',
#                                                                   is_rt_score = FALSE,
#                                                                   is_ms2_score = TRUE)
# parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.5
# parameter_set_annotation@para_ms2_match$matched_frag_cutoff <- 2
#
#
# annotate_metabolite(parameter_set_annotation = parameter_set_annotation)
#
#
# parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/05_coworker/Aedan/230727_untargeted_data_processing/HILIC_neg/',
#                                                                   lib = 'peptide',
#                                                                   ce = '20',
#                                                                   column = 'hilic',
#                                                                   polarity = 'negative',
#                                                                   is_rt_score = FALSE,
#                                                                   is_ms2_score = TRUE)
# parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.5
# parameter_set_annotation@para_ms2_match$matched_frag_cutoff <- 2
#
#
# annotate_metabolite(parameter_set_annotation = parameter_set_annotation)


###############################################################################
# library(DoddLabDatabase)
#
# parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20231205_IBD_enorllment/c18_pos/',
#                                                                   lib = 'gnps_bile_acid',
#                                                                   column = 'hilic',
#                                                                   polarity = 'positive',
#                                                                   is_rt_score = FALSE,
#                                                                   is_ms2_score = TRUE)
# parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.7
# parameter_set_annotation@para_ms2_match$matched_frag_cutoff <- 2
#
#
# annotate_metabolite(parameter_set_annotation = parameter_set_annotation)

#
# parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20231211_IBD_ID_test/c18_pos/',
#                                                                   lib = 'gnps_bile_acid',
#                                                                   column = 'hilic',
#                                                                   polarity = 'positive',
#                                                                   is_rt_score = FALSE,
#                                                                   is_ms2_score = TRUE)
# parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.7
# parameter_set_annotation@para_ms2_match$matched_frag_cutoff <- 2
# parameter_set_annotation@para_ms2_match$direction <- 'forward'
#
# annotate_metabolite(parameter_set_annotation = parameter_set_annotation)


# parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20231211_IBD_ID_test/c18_pos/',
#                                                                   lib = 'dodd',
#                                                                   column = 'hilic',
#                                                                   polarity = 'positive',
#                                                                   is_rt_score = TRUE,
#                                                                   is_ms2_score = TRUE)
# parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.8
# parameter_set_annotation@para_ms2_match$direction <- 'reverse'
# annotate_metabolite(parameter_set_annotation = parameter_set_annotation)
#
# load('~/Project/00_IBD_project/Data/20231211_IBD_ID_test/c18_pos/01_metabolite_annotation/00_intermediate_data/lib_meta')
# load('~/Project/00_IBD_project/Data/20231211_IBD_ID_test/c18_pos/01_metabolite_annotation/00_intermediate_data/ms1_result')
#
#
# get_variable_name <- function(x) {
#   deparse(substitute(x))
# }
#
# # Example
# my_variable <- 42
#
# deparse(substitute(my_variable))
# get_variable_name(my_variable)


# load('~/Project/00_IBD_project/Data/20231211_IBD_ID_test/c18_pos/01_input_data_cleaning/05_object_c18_pos_outlier_removal.RData')
#
# # mz + RT + ms2
# parameter_set_annotation <- initialize_annotation_parameter_class(object = object_c18_pos,
#                                                                   path = '~/Project/00_IBD_project/Data/20231205_IBD_enorllment/c18_pos/',
#                                                                   lib = 'dodd',
#                                                                   column = 'c18',
#                                                                   polarity = 'positive',
#                                                                   is_rt_score = TRUE,
#                                                                   is_ms2_score = TRUE)
# parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.8
# parameter_set_annotation@para_ms2_match$direction <- 'reverse'
# annotate_metabolite(parameter_set_annotation = parameter_set_annotation)

# file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_mz_rt_ms2')

# path_output <- '~/Project/00_IBD_project/Data/20231205_IBD_enorllment/c18_pos/01_metabolite_annotation_dodd_mz_rt_ms2'
# load(file.path(path_output, "00_intermediate_data", 'ms1_result'))

#
################################################################################
# Metabolite identification ----------------------------------------------------
# C18 pos --------------------------------------------------------------------
# setwd('~/Project/00_IBD_project/Data/20240329_update_annotation_table_for_discussion/')
#
#
# # library(DoddLabPackages)
# library(DoddLabMetID)
# library(DoddLabDatabase)
#
# load('~/Project/00_IBD_project/Data/20240329_update_annotation_table_for_discussion//c18_pos/05_object_c18_pos_outlier_removal.RData')
#
# # mz + RT + ms2
# parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240329_update_annotation_table_for_discussion/c18_pos/',
#                                                                   lib = 'dodd',
#                                                                   column = 'c18',
#                                                                   polarity = 'positive',
#                                                                   is_rt_score = TRUE,
#                                                                   is_ms2_score = TRUE)
# parameter_set_annotation@para_ms1_match$mz_tol <- 10
# parameter_set_annotation@para_ms1_match$mz_ppm_thr <- 200
# parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.8
# parameter_set_annotation@para_ms2_match$direction <- 'reverse'
# annotate_metabolite(object = object_c18_pos, parameter_set_annotation)
#
# file.rename(from = '~/Project/00_IBD_project/Data/20240329_update_annotation_table_for_discussion/c18_pos/01_metabolite_annotation', to = './01_metabolite_annotation_dodd_mz_rt_ms2')



# ################################################################################
# library(DoddLabMetID)
# library(DoddLabDatabase)
#
# load('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/c18_pos/01_input_data_cleaning/05_object_c18_pos_outlier_removal.RData')
#
# # mz + RT + ms2
# parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/c18_pos/',
#                                                                   lib = 'dodd',
#                                                                   column = 'c18',
#                                                                   polarity = 'positive',
#                                                                   is_rt_score = TRUE,
#                                                                   is_ms2_score = TRUE)
# parameter_set_annotation@para_ms1_match$mz_tol <- 10
# parameter_set_annotation@para_ms1_match$mz_ppm_thr <- 200
# parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.8
# parameter_set_annotation@para_ms2_match$direction <- 'reverse'
# annotate_metabolite(object = object_c18_pos, parameter_set_annotation)
#
# file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_dodd_mz_rt_ms2')
#
# merge_one_modes(path = '~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/c18_pos/',column = 'c18', polarity = 'positive')



# load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_neg/01_input_data_cleaning/03_object_c18_neg_serrf.RData')
#
# # mz + RT + ms2
# parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_neg/',
#                                                                   lib = 'dodd',
#                                                                   column = 'c18',
#                                                                   polarity = 'negative',
#                                                                   is_rt_score = TRUE,
#                                                                   is_ms2_score = TRUE)
# parameter_set_annotation@para_ms1_match$mz_tol <- 10
# parameter_set_annotation@para_ms1_match$mz_ppm_thr <- 200
# parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.8
# parameter_set_annotation@para_ms2_match$direction <- 'reverse'
# annotate_metabolite(object = object_c18_neg_serrf, parameter_set_annotation)
#
# file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_dodd_mz_rt_ms2')
