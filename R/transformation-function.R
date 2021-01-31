#' Transform dimensions using an RSDR model
#'
#' This function transforms input dimensions using a re-sampled dimensional 
#' reduction model fitted by \code{rsdr()}.
#'
#' @param tidy_set A TidySet (i.e. ExpressionSet) containing the visits of 
#' subjects in outcome dataset, paid by any payment systems. This TidySet also 
#' accomodates outcome dataset. This is an output of 
#' \code{compile_mh_outcome()} from \code{medhist} package, or manually 
#' compiled using \code{ExpressionSet()} from \code{Biobase} package.
#' @param rsdr_object RSDR object, a list of results and parameters. This is
#' an output from \code{rsdr()} function.
#' @param dimensions Dimension index, a vector of integers for the column 
#' indices of the table, indicating the output dimensions transformed by 
#' \code{rsdr()}. The maximum integer is the minimum number of output 
#' dimensions among many dimensional reduction models from differen resampled 
#' subsets.
#' @param input_dim Input dimension, a vector of characters containing the 
#' names of input dimensions that are selected to transform into output 
#' dimensions.
#' @param output_dim Output dimension, a vector of characters containing the 
#' names of output dimensions that are selected after transformation.
#' @param verbose Verbosity, a logical indicating whether progress should be
#' shown.
#'
#' @return A TidySet (i.e. ExpressionSet) containing the transformed table 
#' accessed using \code{exprs()} function from Biobase package. Composition of 
#' input dimensions as a weight table can be accessed using \code{fData()} from 
#' the same package. RSDR models and proportion of variance explained (PVE) 
#' using \code{input_dim} and \code{output_dim} can be accessed using 
#' \code{preproc()}.
#'
#' @keywords RSDR transformer
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
#' ## Transform dimensions using an RSDR model
#' dr_bin_nps=
#'   mh_bin_nps %>%
#'   transformation(
#'     rsdr_object=rsdr_bin_nps
#'     ,input_dim=rownames(mh_bin_nps)
#'     ,top_n=6
#'     ,verbose=F
#'   )

transformation=function(tidy_set
                        ,rsdr_object
                        ,dimensions=NULL
                        ,input_dim=NULL
                        ,top_n=NULL
                        ,output_dim=NULL
                        ,verbose=T){
  
  if(verbose) cat(paste0('Transform dimensions\n'))
  
  if(verbose){
    pb=startpb(0,9)
    on.exit(closepb(pb))
    setpb(pb,0)
  }
  
  if(verbose) setpb(pb,1)
  avg_rotm=
    rsdr_object %>%
    composition(dimensions) %>%
    as.matrix()
  
  if(verbose) setpb(pb,2)
  scaler=
    rsdr_object$rsdr_object %>%
    lapply(X=seq(length(.)),Y=.,function(X,Y){
      data.frame(
        input_dim=names(Y[[X]]$avg)
        ,avg=Y[[X]]$avg
        ,std=Y[[X]]$std
      ) %>%
        `rownames<-`(NULL)
    }) %>%
    do.call(rbind,.)
  
  if(verbose) setpb(pb,3)
  scaler=
    scaler %>%
    group_by(input_dim) %>%
    summarize_all(function(x)mean(x,na.rm=T)) %>%
    ungroup() %>%
    column_to_rownames(var='input_dim')
  
  if(verbose) setpb(pb,4)
  scaled_data=
    tidy_set %>%
    exprs() %>%
    t() %>%
    as.data.frame() %>%
    cbind(
      pData(tidy_set) %>%
        mutate(outcome=as.numeric(outcome=='event'))
    ) %>%
    select(outcome,everything()) %>%
    .[,rownames(avg_rotm)] %>%
    sweep(2,scaler$avg,'-') %>%
    sweep(2,scaler$std,'/') %>%
    as.matrix()
  
  if(is.null(input_dim)) input_dim=colnames(scaled_data)
  
  if(verbose) setpb(pb,5)
  converted_data0=scaled_data %*% avg_rotm
  
  if(verbose) setpb(pb,6)
  converted_data=scaled_data[,input_dim,drop=F] %*% avg_rotm[input_dim,,drop=F]
  
  if(verbose) setpb(pb,7)
  converted_data=as.data.frame(converted_data)
  
  if(is.null(top_n)) top_n=ncol(converted_data)
  
  if(verbose) setpb(pb,7)
  dr_var0=apply(converted_data0,2,var)
  
  if(verbose) setpb(pb,8)
  dr_var=apply(converted_data,2,var)
  
  if(verbose) setpb(pb,9)
  pve0=dr_var/sum(dr_var0)
  pve=dr_var/sum(dr_var)
  
  if(!is.null(top_n)){
    top_n_pve_dimensions=
      pve %>%
      sort(decreasing=T) %>%
      names() %>%
      .[seq(top_n)]
  }
  
  if(!is.null(output_dim)){
    top_n_pve_dimensions=output_dim
  }
  
  if(!is.null(top_n) | !is.null(output_dim)){
    converted_data=converted_data[,top_n_pve_dimensions]
  }
  
  ExpressionSet(
    assayData=
      converted_data %>%
      t()
    ,phenoData=
      tidy_set %>%
      pData() %>%
      cbind(
        tidy_set %>%
          exprs() %>%
          t() %>%
          as.data.frame()
      ) %>%
      AnnotatedDataFrame(
        varMetadata=
          tidy_set %>%
          phenoData() %>%
          varMetadata() %>%
          rbind(
            tidy_set %>%
              rownames() %>%
              data.frame(rowname=.,labelDescription=NA) %>%
              column_to_rownames(var='rowname')
          )
      )
    ,featureData=
      rsdr_object %>%
      composition() %>%
      .[rownames(tidy_set)[rownames(tidy_set)%in%rownames(.)]
        ,colnames(converted_data)[colnames(converted_data)%in%colnames(.)]
        ,drop=F] %>%
      t() %>%
      as.data.frame() %>%
      AnnotatedDataFrame(
        varMetadata=
          tidy_set %>%
          fData()
      ) %>%
      `varMetadata<-`(
        varMetadata(.) %>%
          select(labelDescription,everything())
      )
    ,experimentData=
      tidy_set %>%
      `preproc<-`(
        preproc(.) %>%
          c(list(dmr=rsdr_object$rsdr_object,pve0=pve0,pve=pve))
      ) %>%
      experimentData()
    ,annotation=
      tidy_set %>%
      annotation()
    ,protocolData=
      tidy_set %>%
      protocolData()
  )
  
}
