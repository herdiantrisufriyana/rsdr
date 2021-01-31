#' Get table of estimated weights of transformation
#'
#' This function get table of estimated weights to transform the input to 
#' output dimensions by \code{rsdr()} function.
#'
#' @param rsdr_object RSDR object, a list of results and parameters. This is
#' an output from \code{rsdr()} function.
#' @param dimensions Dimension index, a vector of integers for the column 
#' indices of the table, indicating the output dimensions transformed by 
#' \code{rsdr()}. The maximum integer is the minimum number of output 
#' dimensions among many dimensional reduction models from differen resampled 
#' subsets.
#'
#' @return RSDR object, a list of results and parameters. Use \code{plot()} to 
#' visualize weights that are used to transformed all input dimension to each 
#' output dimension.
#'
#' @keywords weights, rotated matrix, dimensional reduction
#'
#' @export
#'
#' @examples
#'
#' ## Create input example
#' library(medhist)
#' data(medhistdata)
#' ps_remover=extract_nps_mh(medhistdata)
#' 
#' mh_bin_nps=
#'   medhistdata[ps_remover_train$key,] %>%
#'   `exprs<-`(
#'     exprs(.) %>%
#'       t() %>%
#'       as.data.frame() %>%
#'       rownames_to_column(var='id') %>%
#'       column_to_rownames(var='id') %>%
#'       t()
#'   ) %>%
#'   trans_binary(verbose=F)
#'  
#' input=
#'   mh_bin_nps %>%
#'   exprs() %>%
#'   t() %>%
#'   as.data.frame()
#'  
#' ## Fit dimensional reduction models with resampling
#' rsdr_bin_nps=rsdr(input,'CV',10,'PCA')
#' 
#' ## Get table of estimated weights of transformation
#' composition(rsdr_bin_nps)

composition=function(rsdr_object,dimensions=NULL){
  
  rsdr_object=rsdr_object$rsdr_object
  
  if(is.null(dimensions)) dimensions=seq(ncol(rsdr_object[[1]]$dmr$rotm))
  
  min_dim=
    rsdr_object %>%
    sapply(function(x)dim(x$dmr$rotm)[2]) %>%
    min()
  
  dimensions=
    dimensions %>%
    .[. %in% seq(min_dim)]
  
  all_rotm=
    rsdr_object %>%
    lapply(X=seq(length(.)),Y=.,function(X,Y){
      as.data.frame(Y[[X]]$dmr$rotm[,dimensions]) %>%
        rownames_to_column(var='old_dim')
    }) %>%
    do.call(rbind,.)
  
  avg_rotm=
    all_rotm %>%
    group_by(old_dim) %>%
    summarize_all(function(x)mean(x,na.rm=T)) %>%
    ungroup() %>%
    column_to_rownames(var='old_dim')
  
  avg_rotm
  
}
