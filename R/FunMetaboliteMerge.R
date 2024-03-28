################################################################################
  # merge_one_modes -----------------------------------------------------------
#' @title merge_one_modes
#' @description merge multiple annotation using different database and parameter. Note: The input table should be same
#' @author Zhiwei Zhou
#' @param path '.'
#' @param column 'hilic', 'c18'
#' @param polarity 'positive', 'negative'
#' @return
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export


merge_one_modes <- function(path = '.',
                            column = c('hilic', 'c18'),
                            polarity = c('positive', 'negative')) {
  list_dir <- list.files(path)
  list_result <- vector(mode='list', length=3)

  if ('01_metabolite_annotation_dodd_mz_rt_ms2' %in% list_dir) {
    message(crayon::blue('Read annotation result from dodd lib (mz + rt + ms2) ...\n'))
    annot_table_mz_rt_ms2 <- readxl::read_xlsx(file.path(path, '01_metabolite_annotation_dodd_mz_rt_ms2/annotation_summary.xlsx'))
    load(file.path(path, '01_metabolite_annotation_dodd_mz_rt_ms2', '00_intermediate_data', 'ms2'))

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
      dplyr::filter(rt_error <= 15)

    list_result[[1]] <- annot_table_mz_rt_ms2
  }

  if ('01_metabolite_annotation_dodd_mz_rt' %in% list_dir) {
    message(crayon::blue('Read annotation result from dodd lib (mz + rt) ...\n'))
    annot_table_mz_rt <- readxl::read_xlsx(file.path(path, '01_metabolite_annotation_dodd_mz_rt/annotation_summary.xlsx'))

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

  if ('01_metabolite_annotation_mz_ms2_public_db' %in% list_dir) {
    message(crayon::blue('Read annotation result from MS-DIAL lib (mz + ms2) ...\n'))
    annot_table_publicdb <- readxl::read_xlsx(file.path(path, '01_metabolite_annotation_mz_ms2_public_db/annotation_summary.xlsx'))

    if (any(ls() == 'ms2')) {
      load(file.path(path, '01_metabolite_annotation_mz_ms2_public_db', '00_intermediate_data', 'ms2'))
    }

    # modify publicdb annotation
    annot_table_publicdb <- annot_table_publicdb %>%
      dplyr::mutate(inchikey1 = annot_table_publicdb$inchikey %>%
                      ImmsTools::SplitInchiKey() %>%
                      dplyr::pull(inchikey1)
      ) %>%
      dplyr::select(feature_name:inchikey, inchikey1, everything()) %>%
      dplyr::mutate(msms_score_reverse = msms_score_forward,
                    confidence_level = 'Level2.2')

    # dereplication: keep the highest ms2 score candidate
    annot_table_publicdb <- pbapply::pblapply(unique(annot_table_publicdb$feature_name), function(x){
      annot_table_publicdb %>%
        dplyr::filter(feature_name == x) %>%
        dplyr::arrange(desc(msms_score_forward)) %>%
        dplyr::distinct(inchikey1, .keep_all = TRUE)
    }) %>%
      dplyr::bind_rows()

    # additional filtering for publicdb results
    annot_table_publicdb <- annot_table_publicdb %>%
      dplyr::filter(msms_matched_frag >= 2) %>%
      dplyr::filter((mz <= 150) | (mz_error <= 10))

    list_result[[3]] <- annot_table_publicdb
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

  message(crayon::green('Done!'))
}



################################################################################
# export ms2 for each merged candidates for quickly inspection -----------------

# load_ms2_obj ---------------------------------------------------------------
load_ms2_obj <- function(path = '.'){
  list_mode <- c('c18_pos', 'c18_neg', 'hilic_pos', 'hilic_neg')

  ms2_dodd_list <- lapply(list_mode, function(x){
    temp_path <- file.path(path, x, '01_metabolite_annotation_dodd_mz_rt_ms2/00_intermediate_data/ms2_result')
    load(temp_path)
    ms2_result
  })
  names(ms2_dodd_list) <- list_mode

  ms2_public_list <- lapply(list_mode, function(x){
    temp_path <- file.path(path, x, '01_metabolite_annotation_mz_ms2_public_db/00_intermediate_data/ms2_result')
    load(temp_path)
    ms2_result
  })
  names(ms2_public_list) <- list_mode

  result <- list(dodd_lib = ms2_dodd_list,
                 public_lib = ms2_public_list)

  return(result)
}



################################################################################
  # plot_merge_ms2 -------------------------------------------------------------

#' @title plot_merge_ms2
#' @description plot ms2 of merged table for manual check
#' @author Zhiwei Zhou
#' @param path '.' the dirname of 4 folders that contains metabolite identification result
#' @param annot_merge_file "merge_annotation_pos.xlsx"
#' @return
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export


plot_merge_ms2 <- function(path = '.',
                           annot_merge_file = 'merge_annotation_pos.xlsx') {

  cat('Read annotation merge file ...\n')
  annot_merge <- readxl::read_xlsx(file.path(path, annot_merge_file))

  # only keep the ms2 annotated records
  annot_merge <- annot_merge %>%
    dplyr::filter(confidence_level %in% c('Level1', 'Level2.2'))

  cat('Load ms2 object ...\n')
  ms2_list_obj <- load_ms2_obj(path)

  ms2_hilic_pos_dodd <- ms2_list_obj$dodd_lib$hilic_pos
  ms2_hilic_neg_dodd <- ms2_list_obj$dodd_lib$hilic_neg
  ms2_c18_pos_dodd <- ms2_list_obj$dodd_lib$c18_pos
  ms2_c18_neg_dodd <- ms2_list_obj$dodd_lib$c18_neg

  ms2_hilic_pos_public <- ms2_list_obj$public_lib$hilic_pos
  ms2_hilic_neg_public <- ms2_list_obj$public_lib$hilic_neg
  ms2_c18_pos_public <- ms2_list_obj$public_lib$c18_pos
  ms2_c18_neg_public <- ms2_list_obj$public_lib$c18_neg

  cat('Plot ms2 for manually check ...\n')
  plot_list <- pbapply::pblapply(seq_along(annot_merge$feature_name), function(i){
    temp_result <- annot_merge %>%
      dplyr::slice(i)

    temp_confidence_level <- temp_result$confidence_level
    temp_column <- temp_result$column
    temp_polarity <- temp_result$polarity
    temp_feature_name <- temp_result$feature_name
    temp_feature_mz <- temp_result$mz
    temp_feature_rt <- temp_result$rt
    temp_cpd_id <- temp_result$id
    temp_cpd_name <- temp_result$name
    temp_adduct <- temp_result$adduct
    temp_mz_error <- temp_result$mz_error
    temp_rt_error <- temp_result$rt_error
    temp_ms2_score_forward <- temp_result$msms_score_forward
    temp_ms2_score_reverse <- temp_result$msms_score_reverse
    temp_ms2_match_fragment <- temp_result$msms_matched_frag
    temp_mz_lib <- temp_result$mz_lib
    temp_rt_lib <- temp_result$rt_lib
    temp_smiles <- temp_result$smiles
    temp_inchikey <- temp_result$inchikey

    if (temp_confidence_level == 'Level1'){

      if (temp_column == 'hilic') {
        if (temp_polarity == 'positive') {
          temp_ms2_obj <- which(names(ms2_hilic_pos_dodd) == temp_feature_name) %>%
            ms2_hilic_pos_dodd[[.]]
          temp_ms2_obj <- which(temp_ms2_obj@info$name == temp_cpd_id) %>%
            temp_ms2_obj@matchedFragments[[.]]
        } else {
          temp_ms2_obj <- which(names(ms2_hilic_neg_dodd) == temp_feature_name) %>%
            ms2_hilic_neg_dodd[[.]]
          temp_ms2_obj <- which(temp_ms2_obj@info$name == temp_cpd_id) %>%
            temp_ms2_obj@matchedFragments[[.]]
        }
      } else {
        if (temp_polarity == 'positive') {
          temp_ms2_obj <- which(names(ms2_c18_pos_dodd) == temp_feature_name) %>%
            ms2_c18_pos_dodd[[.]]
          temp_ms2_obj <- which(temp_ms2_obj@info$name == temp_cpd_id) %>%
            temp_ms2_obj@matchedFragments[[.]]
        } else {
          temp_ms2_obj <- which(names(ms2_c18_neg_dodd) == temp_feature_name) %>%
            ms2_c18_neg_dodd[[.]]
          temp_ms2_obj <- which(temp_ms2_obj@info$name == temp_cpd_id) %>%
            temp_ms2_obj@matchedFragments[[.]]
        }
      }
    }

    if (temp_confidence_level == 'Level2.1'){
      return(NULL)
    }

    if (temp_confidence_level == 'Level2.2'){
      if (temp_column == 'hilic') {
        if (temp_polarity == 'positive') {
          temp_ms2_obj <- which(names(ms2_hilic_pos_public) == temp_feature_name) %>%
            ms2_hilic_pos_public[[.]]
          temp_ms2_obj <- which(temp_ms2_obj@info$name == temp_cpd_id) %>%
            temp_ms2_obj@matchedFragments[[.]]
        } else {
          temp_ms2_obj <- which(names(ms2_hilic_neg_public) == temp_feature_name) %>%
            ms2_hilic_neg_public[[.]]
          temp_ms2_obj <- which(temp_ms2_obj@info$name == temp_cpd_id) %>%
            temp_ms2_obj@matchedFragments[[.]]
        }
      } else {
        if (temp_polarity == 'positive') {
          temp_ms2_obj <- which(names(ms2_c18_pos_public) == temp_feature_name) %>%
            ms2_c18_pos_public[[.]]
          temp_ms2_obj <- which(temp_ms2_obj@info$name == temp_cpd_id) %>%
            temp_ms2_obj@matchedFragments[[.]]
        } else {
          temp_ms2_obj <- which(names(ms2_c18_neg_public) == temp_feature_name) %>%
            ms2_c18_neg_public[[.]]
          temp_ms2_obj <- which(temp_ms2_obj@info$name == temp_cpd_id) %>%
            temp_ms2_obj@matchedFragments[[.]]
        }
      }
    }


    text <- paste(c(
      paste0('Feature: ', temp_feature_name),
      paste0('Feature m/z: ', round(temp_feature_mz, 4)),
      paste0('Feature RT: ', round(temp_feature_rt, 1)),
      paste0('Column: ', temp_column),
      paste0('Polarity: ', temp_polarity),
      paste0('Compound ID: ', temp_cpd_id),
      paste0('Compound name: ', temp_cpd_name),
      paste0('Confidence level: ', temp_confidence_level),
      paste0('Adduct: ', temp_adduct),
      paste0('SMILES: ', temp_smiles),
      paste0('InChIKey: ', temp_inchikey),
      paste0('m/z lib: ', round(temp_mz_lib, 4)),
      paste0('m/z error: ', temp_mz_error),
      paste0('RT lib: ', temp_rt_lib),
      paste0('RT error: ', temp_rt_error),
      paste0('MS2 score forward: ', temp_ms2_score_forward),
      paste0('MS2 score reverse: ', temp_ms2_score_reverse),
      paste0('MS2 matched frag.: ', temp_ms2_match_fragment)
    ),
    collapse = '\n')

    suppressMessages(
      temp_plot <- DoddLabMetID::plot_id_ms2(obj_spec = temp_ms2_obj) +
        ggplot2::scale_colour_manual(
          name = 'Attribute',
          labels= c(paste0('Experiment'),
                    'Unmatched fragments',
                    paste0('Library ')),
          values = c(
            'experiment' = 'black',
            'library' = 'red',
            'frag_unmatch' = 'gray'
          )
        ) +
        ggplot2::scale_shape_manual(
          name = 'Label',
          labels= c('matched' = "Matched",
                    'unmatched' = "Unmatched"),
          values = c(
            'matched' = 16,
            'unmatched' = 4
          )
        ) +
        ggplot2::ggtitle(label = paste0(temp_feature_name,
                                        ': ', temp_cpd_name,
                                        ' (ID: ', temp_cpd_id,
                                        ', DP: ', temp_ms2_score_forward,
                                        ')')) +
        ZZWtool::ZZW_annotate_text2(label = text, x = 0, y = 1) +
        ggplot2::theme(legend.position = c(0.85, 0.85),
                       title = ggplot2::element_text(vjust = 0.5))
    )

    return(temp_plot)

  })


  dir.create(file.path(path, '00_manual_check_merge'), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(
    filename = file.path(path, '00_manual_check_merge',
                         'ms2_plot_merge_manual_check.pdf'),
    plot = gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1),
    width = 15, height = 9
  )

  cat('Done\n')

}
