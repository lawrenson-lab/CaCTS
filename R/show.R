# @title .onAttach
# @description  Load required data into gloval enviroment

.onAttach <- function (libname, pkgname){

  #if (!interactive() || stats::runif(1) > 0.1) return()
  welcome.message <- paste0(
    " ===================================================================\n",
    "         ___        ___ _____  ___             _                    \n",
    "        |          |      |   |            _  | |                   \n",
    "        |    ___   |      |   |___        | | | | _                 \n",
    "        |     __|  |      |       |       | |_| || |                \n",
    "        |___ |__|  |___   |    ___|       |___   __|                \n",
    "                                              | |                   \n",
    "                                              | |                   \n",
    " -------------------------------------------------------------------\n",
    "       Finding TF Cancer Cell-Type-Specificity (CaCTS)              \n",
    "       Version:",utils::packageVersion("CaCTS"),"\n",
    " ===================================================================\n"
  )
  packageStartupMessage(welcome.message)

}


#' @title Vizualize the overall CaCTS scores curve and the expression values of top N ranked CaCTS scores
#' @description the visualize_scores function helps an user to identify the expression profile of a given TF
#' across all cancer type.
#' @param tf.scores a numeric data.frame containing CaCTS scores.
#' @param rep.matrix a expression matrix of the representative samples
#' @param topn number of top scores TF to be shown (default = 5)
#' @param query.name name of query cancer type or subtype
#' @param filename Name for the plot outuput saved as pdf
#' @param h height value for plot
#' @param w width value for plot
#' @importFrom ggplot2 ggplot
#' @importFrom ggrepel geom_text_repel
#' @return  a prot_grid object containing a list of ggplot2 plots
#' @export

visualize_scores <- function(tf.scores,
                      rep.matrix,
                      topn=5,
                      query.name,
		      filename, 
                      h,
                      w) {


  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("cowplot package is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("ggplot2 package is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if(missing(tf.scores))   stop("Missing tf.scores argument")
  if(missing(rep.matrix))   stop("Missing rep.matrix argument")
  if(missing(query.name))   stop("Missing query.name argument")

  if(!is.null(filename)) {
     foldername = paste0("CaCTs_res")      
     pathname = paste0(getwd(), "/", foldername)
     pdf(paste0(pathname, "/", filename), width=w,height=h)
  }

  if(!(query.name %in% colnames(rep.matrix))){
    stop(paste0("column ", query.name," not found in the rep.matrix"))
  }

  if(!("LogValue" %in% colnames(tf.scores))){
    stop(paste0("column LogValue not found in the tf.scores"))
  }

  query.index <- which(colnames(rep.matrix) %in% as.character(query.name))

  ticknames <- colnames(rep.matrix)

  axiscolors = rep("black", length(ticknames))
  axiscolors[query.index] <- "red"

  topNplot = list()
  topNscores = tail(tf.scores, topn)

  message(paste0("** building top ", topn, " plots..."))

  for (j in 1:topn) {

    TF.name <- topNscores$Name[j]
    TF.index <- which(rownames(rep.matrix) %in% TF.name)

    tf.expr.profile <- data.frame(pos=1:ncol(rep.matrix), value=as.numeric(rep.matrix[TF.index,]))

    plot_single_2 <- ggplot(aes(y = value, x = pos), data = tf.expr.profile) +
      geom_point() +
      geom_line() +
      scale_x_continuous(breaks=1:length(ticknames), labels=ticknames) +
      labs(title =paste0(TF.name, " - Mean expression profile"), x = "Cell types", y = "Relative expresssion") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
	    axis.text.x = element_text(angle = 90, hjust = 1, color=axiscolors))

    topNplot[[j]] <- plot_single_2

  }


  ###########################

  show.names <- tail(tf.scores$Name)

  tf.scores$Rank = seq(nrow(tf.scores), 1)

  tf.scores$x = 1:nrow(tf.scores)

  plot_all <- ggplot(aes(y = LogValue, x = x), data = tf.scores) +
    geom_text_repel(
      data = tf.scores[which(tf.scores$Name %in% show.names),],
      aes(label = tf.scores$Name[which(tf.scores$Name %in% show.names)]),
      size = 5,
      box.padding = unit(0.5, "lines"),
      point.padding = unit(0.9, "lines")
    ) +
    geom_line() +
    theme_minimal() +
    labs(title =  paste0("JSD scores - ",   query.name, " - Representative sample by mean"), x = "Transcription Factors", y = "Significance scores (-log10)") +
    theme(plot.title = element_text(hjust = 0.5))

  library(cowplot)

  topNplot[[6]] <- plot_all

  p <- plot_grid(plotlist = topNplot, ncol = 3)
  print(p)

  if(!missing(filename)) {
      dev.off()
      message(paste0("File saved as: ", pathname, "/", filename ))
  }

  return (p)
}



