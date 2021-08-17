#' Generate GGplot for methods comparison in simulation
#'
#' This function generates a figure for comparison in simulation and save it under the specified path
#' @param fig.dt a data frame of reaults for figure
#' @param fig.title a character as the title of figure
#' @param quants a string vector of quantities to plot
#' @param ref.quants a string vector of quantities to add refrence lines
#' @param ref.locs a real number vector of locations for reference lines
#' @param ref.labels a vector of refline lables
#' @param cols a vector of color codes for plots
#' @param facet.row a string vector of variable names in rows for facet
#' @param facet.col a string vector of variable names in columns for facet
#' @param fig.path path to save the figure
#' @param fig.width width of figure
#' @param fig.height height of figure
#' @param log_10 log10 transformaiton of y axis (used to enlarge difference of Type I error)
#' @export
#' @import dplyr ggplot2
#' @details Returns the ggplot subject for print and save it as .png under specified path
#' @author Xi "Ada" Wang
#' @author Ph.D. Student of Biostatistics
#' @author Penn State College of Medicine
#' @author xzw149@@psu.edu
#' @export

plot.simu <- function(fig.dt, fig.title, quants,
                      ref.quants, ref.locs, ref.labels,
                      cols,
                      facet.row, facet.col,
                      fig.path = "C:\\Users\\wangxi8\\OneDrive - Merck Sharp & Dohme, Corp\\AugRCT_new\\Simulation Results\\Figure\\",
                      fig.width = 12, fig.height = 9,
                      log_10 = FALSE){


  facet.fml <- as.formula(paste(paste(facet.row, collapse = "+"), "~", paste(facet.col, collapse="+"), sep=""))
  # reflines.dt: a dataframe with position of reflines
  reflines.dt <- data.frame(quantity = ref.quants,
                            refline = ref.locs)
  # reflines.text: a dataframe with text of reflines
  reflines.text <- data.frame(quantity = ref.quants,
                              label = ref.labels,
                              y = ref.locs,
                              Drift = "1. No CDD",
                              unm.confs.ind = "1. No",
                              set.num = 1)
  fig.title1 <- paste0(fig.title, " (",paste(quants, collapse = ", "),")")
  if (length(unique(fig.dt$unm.confs.ind))>1){
    p1 <- ggplot(fig.dt %>% dplyr::filter(quantity %in% quants),
                 aes(x=factor(label.new), y=value, group=set.num, col=Drift, shape=unm.confs.ind,
                     linetype=unm.confs.ind))
  }else {
    p1 <- ggplot(fig.dt %>% dplyr::filter(quantity %in% quants),
                 aes(x=factor(label.new), y=value, group=set.num, col=Drift))
  }


  p <-  p1 +
    geom_point() + #theme_bw()+ #remove background
    geom_line(size=0.7, alpha=0.7) +
    facet_grid(facet.fml, scales="free", space="free_x") +
    labs(title="", x="Method") +
    #labs(title=fig.title1, x="Method") +
    theme(axis.text.x = element_text(angle = 45,vjust=0.7),text = element_text(size=12)) +
    # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    geom_hline(data=reflines.dt, aes(yintercept=refline), colour="grey20", linetype="dashed") +
    geom_text(data=reflines.text,mapping=aes(x=-Inf, y=y,hjust   = 0.05, vjust = 0.05,label=label),size=3,col="grey20") +
    scale_color_manual(values = cols) #+ scale_colour_discrete(labels = parse_format())
  if (log_10 == TRUE){p <- p + scale_y_log10()}

  ggsave(filename = paste0(str_replace(fig.title,":",""),".png"), path = fig.path, plot=p,
         width = fig.width, height = fig.height, dpi = 300, units = "in", device='png')
  return(p)
}
