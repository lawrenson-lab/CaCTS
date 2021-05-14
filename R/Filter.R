#' @title Filter top N TF scores based on expression rank
#' @description
#'     This is .
#' @param rep.matrix is a expression matrix of the representative samples from prepare_representaive_samples
#' @param tf.scores a numeric data.frame containing CaCTS scores.
#' @param query.name name of query cancer type or subtype
#' @param pn percentage of top N scores to be considered. Default 0.05
#' @param pnE percentage of top N expression rank to be considered. Default 0.05
#' @return Return the filtered list of CaCTS scores considering expression rank.
#' @export

filter_by_expression_rank <- function (rep.matrix,
                           tf.scores, query.name, pn=0.05, pnE=0.05) {


  if(missing(rep.matrix))   stop("Missing rep.matrix argument")
  if(missing(tf.scores))   stop("Missing scores argument")
  if(missing(query.name))   stop("Missing query.name argument")

  if(!(query.name %in% colnames(rep.matrix))){
    stop(paste0("column ", query.name," not found in the rep.matrix"))
  }

  message("Filtering TFs by expression rank... ")
  ## select top N JSD scores
  topNp = round(pn*nrow(rep.matrix))
  query.index <- which(colnames(rep.matrix) %in% query.name)
  tf.scores <- tf.scores[order(tf.scores$value, decreasing = F),]
  tf.scores$Rank <- 1:nrow(tf.scores)
  names.topJSD <- as.character(tf.scores$Name[1:topNp])
 
  ## select top N Expression scores
  topNp = round(pnE*nrow(rep.matrix))
  info.Exp <- rep.matrix[which(rownames(rep.matrix) %in% tf.scores$Name), query.index, drop = FALSE]
  info.Exp <- info.Exp[order(info.Exp[,1], decreasing = T), , drop = FALSE]
  info.Exp$Rank <- seq(1:nrow(info.Exp))
  names.topEXP = as.character(rownames(info.Exp)[1:topNp])

  final <- names.topEXP[which(names.topEXP %in% names.topJSD)]

  filtered.TF <- tf.scores[which(tf.scores$Name %in% final),]
  filtered.TF <- merge(filtered.TF, info.Exp, by.x = "Name", by.y = "row.names")

  colnames(filtered.TF) <- c("Name", "value", "LogValue", "JSD.rank", "Expr.mean", "Expr.rank")
  filtered.TF <- filtered.TF[order(filtered.TF$JSD.rank),]

  filename = paste0(query.name, "-CaCTS-scores-expression-filtered.txt")
  foldername = paste0("CaCTs_res")
  pathname = paste0(getwd(), "/", foldername)
  write.table(filtered.TF, file = paste0(pathname, "/", filename), quote = F, row.names=F, sep = "\t")
  message(paste0("File saved as: ", pathname, "/", filename ))


  return(filtered.TF)
}



