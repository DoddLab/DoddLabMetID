################################################################################
#   plot_id_ms2 ----------------------------------------------------------------

#' @title plot_id_ms2
#' @description generate ms2 mirror plot
#' @author Zhiwei Zhou
#' @param obj_ms2 object of SpectraData
#' @param obj_spec
#' @export
#' @examples

plot_id_ms2 <- function(
    obj_ms2 = NULL,
    obj_spec = NULL
) {
  # mode <- match(mode)

  if (length(obj_ms2) > 0) {
    obj_spec <- obj_ms2@matchedFragments[[1]]
  }

  if (length(obj_spec) == 0) {
    stop('Please input obj_spec')
  }

  temp_spec <- obj_spec %>%
    # dplyr::mutate(int = tidyr::replace_na(int, 0)) %>%
    dplyr::mutate(int_lib = intensity/max(intensity),
                  int_exp = intensityExp/max(intensityExp)) %>%
    dplyr::rename(mz_lib = mz,
                  mz_exp = mzExp) %>%
    dplyr::select(mz_lib, int_lib, mz_exp, int_exp, fragPrecursor) %>%
    dplyr::mutate(label = dplyr::case_when(
      int_lib > 0 & int_exp > 0 ~ 'matched',
      !(int_lib > 0 & int_exp > 0) ~ 'unmatched'
    ))

  # switch (mode,
  #         'SpecLibMatch' = {
  #
  #         },
  #         'NeighborMatch' = {
  #
  #         }
  # )

  temp_spec1 <- temp_spec %>%
    dplyr::select(mz_lib, int_lib, fragPrecursor, label) %>%
    dplyr::rename(mz = mz_lib, int = int_lib) %>%
    # dplyr::mutate(int = 0-int,
    #               tag = 'library')
    dplyr::mutate(int = 0-int,
                  tag = dplyr::case_when(
                    label == 'unmatched' ~ 'frag_unmatch',
                    label == 'matched' ~ 'library'
                  ))

  temp_spec2 <- temp_spec %>%
    dplyr::select(mz_exp, int_exp, fragPrecursor, label) %>%
    dplyr::rename(mz = mz_exp, int = int_exp) %>%
    # dplyr::mutate(tag = 'experiment')
    dplyr::mutate(tag = dplyr::case_when(
      label == 'unmatched' ~ 'frag_unmatch',
      label == 'matched' ~ 'experiment'
    ))


  temp_data <- temp_spec1 %>%
    dplyr::bind_rows(temp_spec2)


  temp_plot <- ggplot2::ggplot(temp_data) +
    ggplot2::geom_segment(ggplot2::aes(x = mz, xend = mz,
                                       y = 0, yend = int,
                                       colour = tag)) +
    ggplot2::geom_point(ggplot2::aes(x = mz,
                                     y = int,
                                     shape = label,
                                     colour = tag))+
    ggplot2::scale_colour_manual(values = c(
      'experiment' = 'dodgerblue',
      'library' = 'tomato',
      'frag_unmatch' = 'gray'
    )) +
    ggplot2::scale_shape_manual(values = c(
      'matched' = 16,
      'unmatched' = 4
    )) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::xlim(0.95*min(temp_data$mz),
                  1.05*max(temp_data$mz)) +
    ggplot2::xlab('m/z') +
    ggplot2::ylab('Relative intensity') +
    ZZWTheme() +
    # ggplot2::ggtitle(label = paste0(obj_ms2@info$name, '{',
    #                                 round(obj_ms2@info$scoreForward, 4), '}')) +
    ggplot2::theme(legend.position = c(0.8, 0.75),
                   title = ggplot2::element_text(vjust = 0.5))


  return(temp_plot)
}



#   plot_id_shift_ms2 ----------------------------------------------------------
#' @title plot_id_shift_ms2
#' @description generate ms2 mirror plot for shift match (bonanza, hybrid, gnps)
#' @author Zhiwei Zhou
#' @param obj_ms2 object of SpectraData
#' @export
#' @examples

# load('/home/zhouzw/Data_processing/20210224_metdna2_development_test/obj_ms2_cpd1_for_metdna2.RData')
# load('/home/zhouzw/Data_processing/20210224_metdna2_development_test/obj_ms2_cpd2_for_metdna2.RData')
# score_dp <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'dp')
# score_bonanza <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'bonanza')
# score_hybrid <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'hybrid')
# score_gnps <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'gnps')
#
# plotIdShiftMs2(obj_ms2 = score_bonanza)
# plotIdShiftMs2(obj_ms2 = score_hybrid)


plot_id_shift_ms2 <- function(
    obj_ms2 = NULL,
    obj_spec_frag_match,
    obj_spec_frag_nl
) {

  if (length(obj_ms2) > 0) {
    obj_spec_frag_match <- obj_ms2@matchedFragments[[1]]
    obj_spec_frag_nl <- obj_ms2@nlFragments[[1]]
  }

  if (length(obj_spec_frag_match) == 0 & length(obj_spec_frag_nl) == 0) {
    stop('Please input obj_spec_frag_match and obj_spec_frag_nl')
  }


  obj_spec_frag_match <- obj_spec_frag_match %>%
    dplyr::mutate(type = 'exact_match')
  obj_spec_frag_nl <- obj_spec_frag_nl %>%
    dplyr::mutate(type = 'nl_match')

  temp_spec_frag_match <- obj_spec_frag_match %>%
    dplyr::rename(mz_lib = mz,
                  mz_exp = mzExp,
                  int_lib = intensity,
                  int_exp = intensityExp) %>%
    dplyr::select(mz_lib, int_lib, mz_exp, int_exp, fragPrecursor, type) %>%
    dplyr::mutate(label = dplyr::case_when(
      int_lib > 0 & int_exp > 0 ~ 'matched',
      !(int_lib > 0 & int_exp > 0) ~ 'unmatched'
    ))

  temp_spec_frag_nl <- obj_spec_frag_nl %>%
    dplyr::rename(mz_lib = mz,
                  mz_exp = mzExp,
                  int_lib = intensity,
                  int_exp = intensityExp) %>%
    dplyr::select(mz_lib, int_lib, mz_exp, int_exp, fragPrecursor, type) %>%
    dplyr::filter(int_lib > 0 & int_exp > 0) %>%
    dplyr::mutate(mz_lib = dplyr::case_when(type == 'exact_match' ~ mz_lib,
                                            type == 'nl_match' ~ mz_exp)) %>%
    dplyr::mutate(label = dplyr::case_when(
      int_lib > 0 & int_exp > 0 ~ 'matched',
      !(int_lib > 0 & int_exp > 0) ~ 'unmatched'
    ))

  temp_spec <- temp_spec_frag_match %>%
    dplyr::bind_rows(temp_spec_frag_nl)

  temp_spec1 <- temp_spec %>%
    dplyr::select(mz_lib, int_lib, fragPrecursor, type, label) %>%
    dplyr::mutate(int_lib = int_lib/max(int_lib)) %>%
    dplyr::rename(mz = mz_lib, int = int_lib) %>%
    dplyr::mutate(int = 0-int,
                  tag = dplyr::case_when(
                    label == 'unmatched' ~ 'frag_unmatch',
                    label == 'matched' & type == 'exact_match' ~ 'library',
                    label == 'matched' & type == 'nl_match' ~ 'library_shift'
                  ))

  temp_spec2 <- temp_spec %>%
    dplyr::select(mz_exp, int_exp, fragPrecursor, type, label) %>%
    dplyr::mutate(int_exp = int_exp/max(int_exp)) %>%
    dplyr::rename(mz = mz_exp, int = int_exp) %>%
    dplyr::mutate(tag = dplyr::case_when(
      label == 'unmatched' ~ 'frag_unmatch',
      label == 'matched' ~ 'experiment'
    ))


  temp_data <- temp_spec1 %>%
    dplyr::bind_rows(temp_spec2)

  temp_plot <- ggplot2::ggplot(temp_data) +
    ggplot2::geom_segment(ggplot2::aes(x = mz, xend = mz,
                                       y = 0, yend = int,
                                       colour = tag)) +
    ggplot2::geom_point(ggplot2::aes(x = mz,
                                     y = int,
                                     shape = label,
                                     colour = tag))+
    ggplot2::scale_colour_manual(values = c(
      'experiment' = 'dodgerblue',
      'library' = 'tomato',
      'frag_unmatch' = 'gray',
      'library_shift' = 'orange'
    )) +
    ggplot2::scale_shape_manual(values = c(
      'matched' = 16,
      'unmatched' = 4
    )) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::xlim(0.95*min(temp_data$mz),
                  1.05*max(temp_data$mz)) +
    ggplot2::ylim(-1, 1) +
    ggplot2::xlab('m/z') +
    ggplot2::ylab('Relative intensity') +
    ZZWTheme() +
    # ggplot2::ggtitle(label = paste0(obj_ms2@info$name, '{',
    #                                 round(obj_ms2@info$scoreForward, 4), '}')) +
    ggplot2::theme(legend.position = c(0.8, 0.75),
                   title = ggplot2::element_text(vjust = 0.5))


  return(temp_plot)

}



#   annotate_text2 ---------------------------------------------------------------
annotate_text2 <- function(label, x, y, facets=NULL, hjust=0, vjust=0, color='black', alpha=NA,
                           family=thm$text$family, size=thm$text$size, fontface=1, lineheight=1.0,
                           box_just=ifelse(c(x,y)<0.5,0,1), margin=ggplot2::unit(size/2, 'pt'), thm=ggplot2::theme_get()) {
  x <- scales::squish_infinite(x)
  y <- scales::squish_infinite(y)
  data <- if (is.null(facets)) data.frame(x=NA) else data.frame(x=NA, facets)

  tg <- grid::textGrob(
    label, x=0, y=0, hjust=hjust, vjust=vjust,
    gp=grid::gpar(col=ggplot2::alpha(color, alpha), fontsize=size, fontfamily=family, fontface=fontface, lineheight=lineheight)
  )
  ts <- grid::unit.c(grid::grobWidth(tg), grid::grobHeight(tg))
  vp <- grid::viewport(x=x, y=y, width=ts[1], height=ts[2], just=box_just)
  tg <- grid::editGrob(tg, x=ts[1]*hjust, y=ts[2]*vjust, vp=vp)
  inner <- grid::grobTree(tg, vp=grid::viewport(width=ggplot2::unit(1, 'npc')-margin*2, height=ggplot2::unit(1, 'npc')-margin*2))

  ggplot2::layer(
    data = NULL,
    stat = ggplot2::StatIdentity,
    position = ggplot2::PositionIdentity,
    geom = ggplot2::GeomCustomAnn,
    inherit.aes = TRUE,
    params = list(
      grob=grid::grobTree(inner),
      xmin=-Inf,
      xmax=Inf,
      ymin=-Inf,
      ymax=Inf
    )
  )
}



setGeneric(name = 'ZZWTheme',
           def = function(
    type = c('common', 'classic')
           ){

             type <- match.arg(type)
             switch (type,
                     'common' = {
                       result <- ggplot2::theme_bw() +
                         ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.4,
                                                                           face = "bold",
                                                                           size = 14),
                                        panel.grid.major = ggplot2::element_blank(),
                                        panel.grid.minor = ggplot2::element_blank(),
                                        panel.background = ggplot2::element_blank(),
                                        panel.border = ggplot2::element_rect(fill = NA, colour = 'black'),
                                        legend.background = ggplot2::element_blank(),
                                        axis.text.x = ggplot2::element_text(size = 10, colour = 'black'),
                                        axis.text.y = ggplot2::element_text(size = 10, angle = 90,
                                                                            colour = 'black'),
                                        axis.title.x = ggplot2::element_text(size = 12, colour = 'black'),
                                        axis.title.y = ggplot2::element_text(size = 12, colour = 'black'),
                                        axis.ticks = ggplot2::element_line(colour = 'black'))
                     },
                     'classic' = {
                       result <- ggplot2::theme_classic() +
                         ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.4,
                                                                           face = "bold",
                                                                           size = 14),
                                        panel.grid.major = ggplot2::element_blank(),
                                        panel.grid.minor = ggplot2::element_blank(),
                                        panel.background = ggplot2::element_rect(fill = NA, colour = NA),
                                        # panel.border = ggplot2::element_rect(fill = NA, colour = 'black'),
                                        legend.background = ggplot2::element_blank(),
                                        axis.line = ggplot2::element_line(colour = 'black'),
                                        axis.text.x = ggplot2::element_text(size = 10, colour = 'black'),
                                        axis.text.y = ggplot2::element_text(size = 10, angle = 90,
                                                                            colour = 'black'),
                                        axis.title.x = ggplot2::element_text(size = 12, colour = 'black'),
                                        axis.title.y = ggplot2::element_text(size = 12, colour = 'black'),
                                        axis.ticks = ggplot2::element_line(colour = 'black'))
                     }
             )



             return(result)
           }
)
