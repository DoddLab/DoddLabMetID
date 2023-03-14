################################################################################

merge_one_modes <- function(path = '.',
                            column = c('hilic', 'c18'),
                            polarity = c('positive', 'negative')) {
  list_dir <- list.files(path)
  list_result <- vector(mode='list', length=3)

  if ('01_metabolite_annotation_mz_rt_ms2' %in% list_dir) {
    message(crayon::blue('Read annotation result from dodd lib (mz + rt + ms2) ...\n'))
    annot_table_mz_rt_ms2 <- readxl::read_xlsx(file.path(path, '01_metabolite_annotation_mz_rt_ms2/annotation_summary.xlsx'))
    load(file.path(path, '01_metabolite_annotation_mz_rt_ms2', '00_intermediate_data', 'ms2'))

    # add inchikey1 & confidence level
    annot_table_mz_rt_ms2 <- annot_table_mz_rt_ms2 %>%
      dplyr::mutate(inchikey1 = annot_table_mz_rt_ms2$inchikey %>%
                      ImmsTools::SplitInchiKey() %>%
                      dplyr::pull(inchikey1)
      ) %>%
      dplyr::select(feature_name:inchikey, inchikey1, everything()) %>%
      dplyr::mutate(confidence_level = 'Level1')

    # remove replication: keep the highest ms2 score candidate
    annot_table_mz_rt_ms2 <- pbapply::pblapply(unique(annot_table_mz_rt_ms2$feature_name), function(x){
      annot_table_mz_rt_ms2 %>%
        dplyr::filter(feature_name == x) %>%
        dplyr::arrange(desc(msms_score_reverse)) %>%
        dplyr::distinct(inchikey1, .keep_all = TRUE)
    }) %>%
      dplyr::bind_rows()

    annot_table_mz_rt_ms2 <- annot_table_mz_rt_ms2 %>%
      dplyr::filter((mz <= 150) | (mz_error <= 10)) %>%
      dplyr::filter(rt_error <= 15 | msms_matched_frag >= 3)

    list_result[[1]] <- annot_table_mz_rt_ms2
  }

  if ('01_metabolite_annotation_mz_rt' %in% list_dir) {
    message(crayon::blue('Read annotation result from dodd lib (mz + rt) ...\n'))
    annot_table_mz_rt <- readxl::read_xlsx(file.path(path, '01_metabolite_annotation_mz_rt/annotation_summary.xlsx'))

    # add inchikey1 & confidence level
    annot_table_mz_rt <- annot_table_mz_rt %>%
      dplyr::mutate(inchikey1 = annot_table_mz_rt$inchikey %>%
                      ImmsTools::SplitInchiKey() %>%
                      dplyr::pull(inchikey1)
      ) %>%
      dplyr::select(feature_name:inchikey, inchikey1, everything()) %>%
      dplyr::mutate(confidence_level = 'Level2.1')

    # remove replication: keep the candidate follow RT score & mz_error
    annot_table_mz_rt <- pbapply::pblapply(unique(annot_table_mz_rt$feature_name), function(x){
      annot_table_mz_rt %>%
        dplyr::filter(feature_name == x) %>%
        dplyr::arrange(desc(rt_score),
                       mz_error) %>%
        dplyr::distinct(inchikey1, .keep_all = TRUE)
    }) %>%
      dplyr::bind_rows()

    # additional filtering
    annot_table_mz_rt <- annot_table_mz_rt %>%
      dplyr::filter((mz <= 150) | (mz_error <= 10)) %>%
      dplyr::filter(rt_error <= 15)

    list_result[[2]] <- annot_table_mz_rt
  }

  if ('01_metabolite_annotation_msdial' %in% list_dir) {
    message(crayon::blue('Read annotation result from MS-DIAL lib (mz + ms2) ...\n'))
    annot_table_msdial <- readxl::read_xlsx(file.path(path, '01_metabolite_annotation_msdial/annotation_summary.xlsx'))

    if (any(ls() == 'ms2')) {
      load(file.path(path, '01_metabolite_annotation_msdial', '00_intermediate_data', 'ms2'))
    }

    # modify msdial annotation
    annot_table_msdial <- annot_table_msdial %>%
      dplyr::mutate(inchikey1 = annot_table_msdial$inchikey %>%
                      ImmsTools::SplitInchiKey() %>%
                      dplyr::pull(inchikey1)
      ) %>%
      dplyr::select(feature_name:inchikey, inchikey1, everything()) %>%
      dplyr::mutate(msms_score_reverse = msms_score_forward,
                    confidence_level = 'Level2.2')

    # dereplication: keep the highest ms2 score candidate
    annot_table_msdial <- pbapply::pblapply(unique(annot_table_msdial$feature_name), function(x){
      annot_table_msdial %>%
        dplyr::filter(feature_name == x) %>%
        dplyr::arrange(desc(msms_score_forward)) %>%
        dplyr::distinct(inchikey1, .keep_all = TRUE)
    }) %>%
      dplyr::bind_rows()

    # additional filtering for msdial results
    annot_table_msdial <- annot_table_msdial %>%
      dplyr::filter(msms_matched_frag >= 2) %>%
      dplyr::filter((mz <= 150) | (mz_error <= 10))

    list_result[[3]] <- annot_table_msdial
  }


  # merge annotation table from all modes
  annot_table_merge <- list_result %>%
    dplyr::bind_rows()

  annot_table_merge <- pbapply::pblapply(unique(annot_table_merge$feature_name), function(x){
    result <- annot_table_merge %>%
      dplyr::filter(feature_name == x) %>%
      dplyr::arrange(confidence_level) %>%
      dplyr::distinct(inchikey1, .keep_all = TRUE)

    highest_confidence_level <- result$confidence_level[1]

    result <- result %>%
      dplyr::filter(confidence_level == highest_confidence_level)

    return(result)
  }) %>%
    dplyr::bind_rows()


  annot_table_merge <- pbapply::pblapply(unique(annot_table_merge$inchikey1), function(y){
    result <- annot_table_merge %>%
      dplyr::filter(inchikey1 == y)

    highest_confidence_level <- result$confidence_level[1]

    result <- result %>%
      dplyr::filter(confidence_level == highest_confidence_level)

    return(result)
  }) %>%
    dplyr::bind_rows()

  # add with ms2 label
  if (any(ls() == 'ms2')) {
    feature_with_ms2 <- names(ms2)
    temp <- rep(-1, length = nrow(annot_table_merge))
    temp[annot_table_merge$feature_name %in% feature_with_ms2] <- 1

    annot_table_merge <- annot_table_merge %>%
      dplyr::mutate(with_ms2 = temp) %>%
      dplyr::select(feature_name:rt, with_ms2, dplyr::everything())
  } else {
    annot_table_merge <- annot_table_merge %>%
      dplyr::mutate(with_ms2 = -1) %>%
      dplyr::select(feature_name:rt, with_ms2, dplyr::everything())
  }

  # add column and polarity label
  annot_table_merge <- annot_table_merge %>%
    dplyr::mutate(column = column,
                  polarity = polarity) %>%
    dplyr::select(feature_name:with_ms2, column, polarity, dplyr::everything())

  writexl::write_xlsx(annot_table_merge,
                      path = file.path(path,
                                       paste0('annotation_table_merge_', column, '_', polarity, '.xlsx')),
                      format_headers = FALSE)
}
