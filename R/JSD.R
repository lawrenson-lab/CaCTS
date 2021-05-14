#' @title Calculate The Jensen-Shannon Divergence considering expression matrix of representative samples
#' @description
#'     This is a implementation of the Jensen-Shannon Divergence for Transcription Factors expression
#'     across diffrent cancer types. The expectation is that you have cancer cell-type expression specificity.
#'     That are considered as CaCTS scores.
#' @param rep.matrix is a expression matrix of the representative samples from prepare_representaive_samples
#' @param query.index numeric index of the cancer type to be evaluated as query
#' @return Return the list of CaCTS scores assigned to each TF.
#' @export

run_CaCTS_score <- function (rep.matrix,
                           query.name) {


  if(missing(rep.matrix))   stop("Missing rep.matrix argument")
  if(missing(query.name))   stop("Missing query.name argument")

  if(!(query.name %in% colnames(rep.matrix))){
    stop(paste0("column ", query.name," not found in the rep.matrix"))
  }

  query.index <- which(colnames(rep.matrix) %in% query.name)

  scores = NULL

  message("Calculating scores... ")
  pb <- txtProgressBar(min = 0, max = nrow(rep.matrix), style = 3)

  ### block TF
  for (i in 1:nrow(rep.matrix)) {


    tf.expr.profile = data.frame(pos=1:ncol(rep.matrix), value=as.numeric(rep.matrix[i,]))

    #For the observed pattern, the vector was formed by values from the expression profiles of the query cell type and the balanced background dataset.
    #The elements in this vector are divided by the sum so that the new normalized vector sums to 1.

    #For the idealized pattern, the vector was formed by a value of 1 at the position equivalent to that of the query cell type and zeroes at all other positions

    obs = tf.expr.profile$value/ sum(tf.expr.profile$value, na.rm = T)

    ideal = rep(0.00000000000000001, ncol(rep.matrix))
    ideal[query.index] = 1

    m <- 0.5 * (obs + ideal)
    res <- 0.5 * (sum(obs * log(obs / m)) + sum(ideal * log(ideal / m)))

    scores = rbind(scores, res)
    setTxtProgressBar(pb, i)

  }


  close(pb)


  tf.scores = suppressWarnings(data.frame(Name=rownames(rep.matrix), value = scores))
  tf.scores$LogValue = -log10(tf.scores$value)
  tf.scores = tf.scores[order(tf.scores$value, decreasing = T),]

  filename = paste0(query.name, "-CaCTS-scores.txt")
  foldername = paste0("CaCTs_res")
  pathname = paste0(getwd(), "/", foldername)
  save.scores <- tf.scores[order(tf.scores$value, decreasing = F),]
  save.scores$Rank <- 1:nrow(save.scores)
  write.table(save.scores, file = paste0(pathname, "/", filename), quote = F, row.names=F, sep = "\t")
  message(paste0("File saved as: ", pathname, "/", filename ))

  #tf.scores = cbind(tf.scores, x=1:length(scores))

  return(tf.scores)
}


KLD <- function (A, B)
{
  sum(A * log(A/B))
}

JSD <- function (P, Q)
{
  M = (P + Q)/2
  jsd = 0.5 * KLD(P, M) + 0.5 * KLD(Q, M)
  return(jsd)
}
