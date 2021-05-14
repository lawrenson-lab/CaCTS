
#' @title Prepare the expression matrix to contain representative samples
#' @description Representative samples are generated considering the expression mean across all samples
#' @param expr.matrix a gene expression matrix
#' @param sample.descr a sample annotation data.frame containing cancer type or cancer subtype information
#' @param save.file if True the matrix with representative samples with be saved.
#' @param rep.NA value to replace NAs. Default 1e-17
#' @param filename Name for the outuput expression matrix file.
#' @return Return a expression matrix of the representative samples.
#' @export

prepare_representaive_samples <- function(expr.matrix,
                                          sample.descr,
                                          save.file=T,
                                          rep.NA = 0.00000000000000001,
                                          filename = "expr.matrix.representative.sample.mean.txt") {

  if(missing(expr.matrix))   stop("Missing expr.matrix argument")
  if(missing(sample.descr))   stop("Missing sample.descr argument")


  if(!("sample.id" %in% colnames(sample.descr))){
    stop(paste0("column sample.id not found in the sample.descr"))
  }

  if (!all(sample.descr$sample.id %in% colnames(expr.matrix))){
    stop("Expressiom matrix columns should be the same at sample.descr sample.id")
  }

  if(!("group.name" %in% colnames(sample.descr))){
    stop(paste0("column group.name not found in the sample.descr"))
  }

  ordered.matrix <- data.frame(row.names =  rownames(expr.matrix))

  cancer.group.name <- as.data.frame(table(sample.descr$group.name))

  # library(dplyr)
  # group_by(my_data, group) %>%
  #   summarise(
  #     count = n(),
  #     mean = mean(weight, na.rm = TRUE),
  #     sd = sd(weight, na.rm = TRUE)
  # )


  for (name in cancer.group.name$Var1) {

    sample.names <- sample.descr$sample.id[which(sample.descr$group.name == name)]

    sub.expr.matrix = expr.matrix[,(which(colnames(expr.matrix) %in% sample.names))]

    if (length(sample.names) > 1) {
       sub.expr.matrix$mean = apply(sub.expr.matrix,1,mean,na.rm=T)
       ordered.matrix = cbind(ordered.matrix, sub.expr.matrix$mean)

    }else {
       ordered.matrix = cbind(ordered.matrix, sub.expr.matrix)
    }

  }

  colnames(ordered.matrix) <- cancer.group.name$Var1

  ordered.matrix[ is.na(ordered.matrix) ] <- rep.NA

  if (save.file) {
    if(!is.null(filename)) {
      foldername = paste0("CaCTs_res")
      pathname = paste0(getwd(), "/", foldername)
      dir.create(pathname)
      message(paste("file saved in: ", file.path(pathname, filename)))
      write.table(ordered.matrix, file = file.path(pathname, filename), quote = F, sep = "\t")
    }
  }

  return(ordered.matrix)

}


#            if(is.null(id)) {
#                message("=============== INNPUT ERROR =================")
#                message("I'm expecting one of these columns:")
#                message(" => barcode")
#                message("    Has the complete barcode (TCGA-AA-3833-01A-01D-0904-05)")
#                message(" => bcr_patient_barcode")
#                message("    Has the patient barcode (TCGA-AA-3833)")
#                message(" => patient")
#                message("    Has the patient barcode (TCGA-AA-3833)")
#                message(" => sample")
#                message("    Has the sample barcode (TCGA-AA-3833-01A)")
#                message("-----------------------------------------------")
#                message("Obs: The complete barcode is the recommended one, as the others might lead to errors")
#                return(NULL)
#            }
