# Common helper functions
# by Matthew Galbraith

## Setting and modifying default theme for plots
theme_set(theme_gray(base_size=12, base_family="Arial") +
            theme(
              panel.border=element_rect(colour="black", fill="transparent"),
              plot.title=element_text(face="bold", hjust=0),
              axis.text=element_text(color="black", size=14),
              axis.text.x=element_text(angle=0, hjust=0.5),
              axis.ticks = element_line(color = "black"), # make sure tick marks are black
              panel.background=element_blank(),
              panel.grid=element_blank(),
              plot.background=element_blank(),
              strip.background = element_blank(), # facet label borders
              legend.key=element_blank(), legend.background=element_blank() # remove grey bg from legend
            )
)


## Density color function
getDenCols <- function(x, y, transform = TRUE) { # set to TRUE if using log2 transformation of data
  if(transform) {
    df <- data.frame(log2(x), log2(y))
  } else{
    df <- data.frame(x, y)
  }
  z <- grDevices::densCols(df, colramp = grDevices::colorRampPalette(c("black", "white")))
  df$dens <- grDevices::col2rgb(z)[1,] + 1L
  cols <-  grDevices::colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
  df$col <- cols[df$dens]
  return(df$dens)
} # End of function


## Excel export function
export_excel <- function(named_list, filename = "") {
  wb <- openxlsx::createWorkbook()
  ## Loop through the list of split tables as well as their names
  ## and add each one as a sheet to the workbook
  Map(function(data, name){
    openxlsx::addWorksheet(wb, name)
    openxlsx::writeData(wb, name, data)
  }, named_list, names(named_list))
  ## Save workbook to working directory
  openxlsx::saveWorkbook(wb, file = here("results", paste0(out_file_prefix, filename, ".xlsx")), overwrite = TRUE)
  cat("Saved as:", here("results", paste0(out_file_prefix, filename, ".xlsx")))
} # end of function



# get standard ggplot colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


# Volcano plot function for lms with labels -----
volcano_plot_lab_lm <- function(
    res, 
    y_lim = c(0, NA), 
    n_labels = 3,
    title = "Volcano plot needs title", 
    subtitle = "optional subtitle"
    ){
  res <- res %>% 
    mutate(
      color = if_else(BHadj_pval < 0.1, "q < 0.1", "All")
    )
  # get max for x-axis
  x_lim <- res %>% 
    summarize(max = max(log2(FoldChange), na.rm = TRUE), min = min(log2(FoldChange), na.rm = TRUE)) %>% 
    abs() %>% 
    max() %>% 
    plyr::round_any(accuracy = 1, f = ceiling) # CUSTOMIZE
  res %>% 
    ggplot(aes(log2(FoldChange), -log10(BHadj_pval), color = color)) + 
    geom_hline(yintercept = -log10(0.1), linetype = 2) + 
    geom_vline(xintercept = 0, linetype = 2) + 
    geom_point() + 
    scale_color_manual(values = c("q < 0.1" = "red", "All" = "black")) + 
    xlim(-x_lim, x_lim) +
    ylim(y_lim) +
    geom_text_repel(data = res %>% filter(!is.na(BHadj_pval) & FoldChange > 1) %>% slice_max(order_by = FoldChange, n = n_labels), aes(label = Cell_cluster_name), min.segment.length = 0, show.legend = FALSE, nudge_x = 0.25, nudge_y = 0.1) +
    geom_text_repel(data = res %>% filter(!is.na(BHadj_pval) & FoldChange < 1) %>% slice_min(order_by = FoldChange, n = n_labels), aes(label = Cell_cluster_name), min.segment.length = 0, show.legend = FALSE, nudge_x = -0.25, nudge_y = 0.1) +
    theme(aspect.ratio=1.2) +
    labs(
      title = title,
      subtitle = subtitle
    )
} # end of function
