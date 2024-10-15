#' @title generate_input_from_msp
#' @description generate input files from msp file for annotate_metabolite
#' @param msp_file path to msp file
#' @param dir_path path to directory where input files will be saved
#' @return input files for annotate_metabolite
#' @export
#' @examples
#' generate_input_from_msp(msp_file = '~/Project/05_coworker/Kazuma/241014/ms2_spec.msp',
#'                        dir_path = '~/Project/05_coworker/Kazuma/241014/')


# generate_input_from_msp(msp_file = '~/Project/05_coworker/Kazuma/241014/ms2_spec.msp',
#                         dir_path = '~/Project/05_coworker/Kazuma/241014/')
# tem_spec <- readMSP('~/Project/05_coworker/Kazuma/241014/ms2_spec.msp')

# dir_path <- '~/Project/05_coworker/Kazuma/241014/'

generate_input_from_msp <- function(msp_file = '~/Project/05_coworker/Kazuma/241014/ms2_spec.msp',
                                    dir_path = '.'){

  message(crayon::blue("Start generate input files from msp file\n"))
  # generate ms1 table
  raw_spec <- readMSP(msp_file)
  temp_table <- lapply(raw_spec, function(x){
    x$info
  }) %>%
    dplyr::bind_rows()

  temp_table <- temp_table %>%
    dplyr::select(name, mz, RETENTIONTIME, COLLISIONENERGY) %>%
    dplyr::rename(mz = mz, rt = RETENTIONTIME, collision_energy = COLLISIONENERGY)

  temp_variable_id <- featureReName(temp_table$name)
  temp_table$name <- temp_variable_id

  readr::write_csv(temp_table,
                   file = file.path(dir_path, 'input_ms1_table.csv'))

  # modify ms2 msp file
  raw_spec2 <- lapply(seq_along(raw_spec), function(i){
    x <- raw_spec[[i]]
    x$info$name <- temp_variable_id[i]
    x
  })

  purrr::walk(seq_along(raw_spec2), function(i){
    # cat(i, ' ')
    DoddLabDatabase::generate_msp(file_name = file.path(dir_path, 'input.msp'),
                                  cmp_name = raw_spec2[[i]]$info$name,
                                  precusormz = raw_spec2[[i]]$info$mz,
                                  rt = as.numeric(raw_spec2[[i]]$info$RETENTIONTIME),
                                  instrument = raw_spec2[[i]]$info$INSTRUMENT,
                                  instrument_type = raw_spec2[[i]]$info$INSTRUMENTTYPE,
                                  ce = raw_spec2[[i]]$info$COLLISIONENERGY,
                                  polarity = raw_spec2[[i]]$info$POLARITY,
                                  spec = raw_spec2[[i]]$spec)
  })

  message(crayon::blue("Done!\n"))

}
