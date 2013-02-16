library(grid)

theme_classic2 <- function (base_size = 12, base_family = "") 
{
    theme_grey(base_size = base_size, base_family = base_family) %+replace% 
        theme(axis.text = element_text(size = rel(0.8)), 
              axis.ticks = element_line(colour = "black"), 
              legend.key = element_rect(colour = "grey80"), 
              panel.background = element_rect(fill = "white", colour = NA), 
              panel.border = element_rect(fill = NA, colour = "black", size=1.0), 
              panel.grid.major = element_line(linetype="blank", size = 0.2), 
              panel.grid.minor = element_line(linetype = "blank", size = 0.5), 
              panel.margin     = unit(0.75, "lines"),
              strip.background = element_rect(fill = NA, colour = "black",linetype="blank"), 
              strip.text.x = element_text(face="bold"))
}

# pos - where to add new labels
# newpage, vp - see ?print.ggplot
facetAdjust <- function(x, pos = c("up", "down"), 
                        newpage = is.null(vp), vp = NULL)
{
  # part of print.ggplot
  ggplot2:::set_last_plot(x)
  if(newpage)
    grid.newpage()
  pos <- match.arg(pos)
  p <- ggplot_build(x)
  gtable <- ggplot_gtable(p)
  # finding dimensions
  dims <- apply(p$panel$layout[2:3], 2, max)
  nrow <- dims[1]
  ncol <- dims[2]
  # number of panels in the plot
  panels <- sum(grepl("panel", names(gtable$grobs)))
  space <- ncol * nrow
  # missing panels
  n <- space - panels
  # checking whether modifications are needed
  if(panels != space){
    # indices of panels to fix
    idx <- (space - ncol - n + 1):(space - ncol)
    # copying x-axis of the last existing panel to the chosen panels 
    # in the row above
    gtable$grobs[paste0("axis_b",idx)] <- list(gtable$grobs[[paste0("axis_b",panels)]])
    if(pos == "down"){
      # if pos == down then shifting labels down to the same level as 
      # the x-axis of last panel
      rows <- grep(paste0("axis_b\\-[", idx[1], "-", idx[n], "]"), 
                   gtable$layout$name)
      lastAxis <- grep(paste0("axis_b\\-", panels), gtable$layout$name)
      gtable$layout[rows, c("t","b")] <- gtable$layout[lastAxis, c("t")]
    }
  }
  # again part of print.ggplot, plotting adjusted version
  if(is.null(vp)){
    grid.draw(gtable)
  }
  else{
    if (is.character(vp)) 
      seekViewport(vp)
    else pushViewport(vp)
    grid.draw(gtable)
    upViewport()
  }
  invisible(p)
}
