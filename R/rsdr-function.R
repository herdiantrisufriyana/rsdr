#' Fit dimensional reduction models with resampling
#'
#' This function fits a dimensional reduction model with resampling. Dataset is
#' resampled a number of times. A model using the same dimensional reduction
#' algorithm is fitted each time. These models can be used to transform dataset
#' into new dimensions using an estimated weight per a pair of input and output
#' dimension. The output dimensions are sorted from the highest to the lowest
#' proportion of variance explained (PVE). Thus, one can choose dimensions with
#' top PVE as the feature candidates for developing a prediction model.
#'
#' @param data Input data, a data frame with rows of samples and columns of 
#' variables that will be transformed.
#' @param rs_method Resampling method, a character of \code{BS} for 
#' bootstrapping or \code{CV} for k-fold cross-validation.
#' @param rs_number Resampling time/fold, an integer of any number. A common 
#' number for bootstrapping and cross-validation are 30 and 10, respectively.
#' @param dr_method Dimensional reduction method, a character of \code{PCA} for 
#' principal componen analysis or \code{SVD} for singular value decomposition.
#' @param sd_cutoff Standard deviation cutoff, a non-negative numeric of which 
#' a variable is excluded if the standard deviation is equal to this number or 
#' lower. This number is conceivably 0 if all values in a variable are the 
#' same. This situation (i.e. zero variance) is not allowed for dimensional 
#' reduction.
#' @param state An integer to set random seed for reproducible results.
#' @param cl Parallel cluster, a non-negative integer of number of CPU cluster
#' used for computation in parallel. Set to 1 if no parallelism is expected.
#'
#' @return RSDR object, a list of results and parameters. Use \code{plot()} to 
#' visualize weights that are used to transformed all input dimension to each 
#' output dimension.
#'
#' @keywords dimensional reduction, principal component analysis, singular 
#' value decomposition, resampling, bootstrapping, cross validation
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
#' ## Show fitting results
#' rsdr_bin_nps
#' 
#' ## Plot weights to transform dimensions
#' plot(rsdr_bin_nps_train)

rsdr=function(data
              ,rs_method=c('BS','CV')
              ,rs_number=c(30,10)
              ,dr_method=c('PCA','SVD')
              ,sd_cutoff=0
              ,state=33
              ,cl=1){
  
  set.seed(state)
  if(rs_method=='BS'){
    idx_list=
      seq(rs_number) %>%
      lapply(FUN=seq,nrow(data)) %>%
      lapply(FUN=sample,nrow(data),T)
  }else if(rs_method=='CV'){
    idx_list=
      round(seq(1,nrow(data),len=rs_number)) %>%
      lapply(X=2:length(.),Y=.,function(X,Y){
        seq(ifelse(X==2,1,Y[X-1]+1),Y[X])
      })
    idx_list=
      idx_list %>%
      lapply(X=seq(rs_number),Y=.,function(X,Y){
        Y[-X] %>% unlist
      })
  }else{
    idx_list=NA
  }
  
  if(dr_method=='PCA'){
    dimr_func=prcomp
    rotm_name='rotation'
    prefix='PC'
  }else if(dr_method=='SVD'){
    dimr_func=svd
    rotm_name='v'
    prefix='SV'
  }else{
    dimr_func=NA
    rotm_name=NA
    prefix=NA
  }
  
  cat(paste0('Fit ',dr_method,' models with ',rs_method,' resampling\n'))
  
  dimr_exec=function(idx_list,data,dimr_func){
    rs_data=data %>% .[idx_list,,drop=F]
    rs_nzv=
      rs_data %>%
      gather() %>%
      group_by(key) %>%
      summarize(sd_value=sd(value,na.rm=T),.groups='drop') %>%
      filter(sd_value>sd_cutoff)
    rs_data=rs_data[,rs_nzv$key,drop=F]
    
    avg=
      summarize_all(rs_data,mean) %>%
      gather()
    avg=
      setNames(avg$value,avg$key)

    std=
      summarize_all(rs_data,sd) %>%
      gather()
    std=
      setNames(std$value,std$key)

    dmr=
      as.matrix(rs_data) %>%
      sweep(2,avg,'-') %>%
      sweep(2,std,'/') %>%
      dimr_func()

    dmr$rotm=
      dmr[[rotm_name]] %>%
      `dimnames<-`(list(names(avg),paste0(prefix,seq(length(avg)))))
    
    list(avg=avg,std=std,dmr=dmr['rotm'])
  }
  
  cat('Started:',as.character(now()),'\n')
  cl=makeCluster(cl)
  clusterEvalQ(cl,{
    library('tidyverse')
    library('pbapply')
  })
    
    rsdr_object=pblapply(idx_list,cl=cl,FUN=dimr_exec,data,dimr_func)
    
  stopCluster(cl)
  rm(cl)
  gc()
  cat('End:',as.character(now()))
  
  range_dim=
    rsdr_object %>%
    sapply(function(x)dim(x$dmr$rotm)[2]) %>%
    range() %>%
    .[!is.na(.)] %>%
    unique()
  
  output=
    list(
      rsdr_object=rsdr_object
      ,input_dim=colnames(data)
      ,rs_method=rs_method
      ,rs_number=rs_number
      ,dr_method=dr_method
      ,sd_cutoff=sd_cutoff
      ,state=state
      ,min_output_dim=range_dim[1]
      ,max_output_dim=range_dim %>% .[length(.)]
      ,dimr_func=dimr_func
      ,rotm_name=rotm_name
      ,prefix=prefix
    )
  
  class(output)='rsdr'
  
  assign(x='print.rsdr',envir=baseenv(),value=function(x){
    
    input_desc=
      paste0(
        length(x$input_dim)
        ,' dimension',ifelse(length(x$input_dim)>1,'s','')
      )
    
    output_desc=
      paste0(
        ifelse(
          length(unique(c(x$min_output_dim,x$max_output_dim)))==1 &
          x$min_output_dim==1
          ,''
          ,'Up to '
        )
        ,ifelse(
          length(unique(c(x$min_output_dim,x$max_output_dim)))==1
          ,paste0(
            x$min_output_dim
            ,' dimension',ifelse(x$min_output_dim>1,'s','')
           )
          ,paste0(
            x$min_output_dim,'~',x$max_output_dim
            ,' dimensions'
           )
        )
      )
    
    rs_method=case_when(
      x$rs_method=='BS'~paste0(x$rs_number,'x Bootstrapping')
      ,x$rs_method=='CV'~paste0(x$rs_number,'-Fold Cross Validation')
      ,TRUE~''
    )
    
    dr_method=case_when(
      x$dr_method=='PCA'~'Principle Component Analysis (PCA)'
      ,x$dr_method=='SVD'~'Singular Value Decomposition (SVD)'
      ,TRUE~''
    )
    
    cat(paste0(
      '
    Method: rsdr

    Re-Sampled Dimensional Reduction (RSDR)
    
    Input: ',input_desc,'
    
    Parameters:
      DR method: ',dr_method,'
      RS method: ',rs_method,'
      SD cutoff: ',x$sd_cutoff,' (to remove zero-/near-zero variance)
      Random state: ',x$state,'
    
    Output: ',output_desc,'
    Prefix: ',x$prefix,'
    
    Callable variables in this object:
      ',paste(names(x)[1:3],collapse=', '),',
      ',paste(names(x)[4:6],collapse=', '),',
      ',paste(names(x)[7:9],collapse=', '),',
      ',paste(names(x)[10:12],collapse=', '),'
    '
    ))
    
  })
  
  assign(x='plot.rsdr',envir=baseenv(),value=function(x
                                                      ,w_cutoff=0
                                                      ,label=F
                                                      ,label_digits=2
                                                      ,label_size=3
                                                      ,font_family='sans'
                                                      ,font_size=9
                                                      ,input_dim_order=NULL
                                                      ,output_dim_order=NULL
                                                      ,input_dim_desc=NULL
                                                      ,output_dim_desc=NULL){
    
    data=composition(x)
    
    if(!is.null(input_dim_order)){
      data=
        data %>%
        .[input_dim_order[input_dim_order%in%rownames(.)],,drop=F]
    }
    
    if(!is.null(output_dim_order)){
      data=
        data %>%
        .[,output_dim_order[output_dim_order%in%colnames(.)],drop=F]
    }
    
    data=
      data %>%
      rownames_to_column(var='input_dim') %>%
      gather(output_dim,weight,-input_dim)
    
    if(!is.null(input_dim_desc)){
      data=
        data %>%
        left_join(input_dim_desc,by='input_dim') %>%
        unite(input_dim,input_dim_desc,input_dim,sep=' - ')
    }
    
    if(!is.null(output_dim_desc)){
      data=
        data %>%
        left_join(output_dim_desc,by='output_dim') %>%
        unite(output_dim,output_dim_desc,output_dim,sep=' - ')
    }
    
    data=
      data %>%
      mutate(
        input_dim=input_dim %>% factor(unique(.))
        ,output_dim=output_dim %>% factor(unique(.))
      ) %>%
      mutate(weight=ifelse(abs(weight)<=w_cutoff,0,weight))
    
    p=data %>%
      qplot(output_dim,input_dim,fill=weight,data=.,geom='tile')
    
    if(label){
      p=p +
        geom_text(
          aes(
            label=ifelse(weight==0,NA,round(weight,label_digits))
            ,alpha=abs(weight)
          )
          ,size=label_size
          ,color='#FFFFFF'
          ,family=font_family
          ,na.rm=T
        )
    }
    
    p +
      scale_x_discrete('Output dimension') +
      scale_y_discrete('Input dimension') +
      scale_fill_gradient2(
        paste0('Cutoff\n<=|',w_cutoff,'|')
        ,low='red',mid='black',high='green',midpoint=0
      ) +
      scale_alpha(guide=F) +
      theme(
        axis.title=element_text(family=font_family,size=font_size)
        ,axis.text=element_text(family=font_family,size=font_size)
        ,axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)
        ,axis.text.y=element_text(angle=0,hjust=1,vjust=0.5)
        ,legend.title.align=0.5
        ,legend.title=element_text(family=font_family,size=font_size)
        ,legend.text=element_text(family=font_family,size=font_size)
      )
    
  })
  
  output
  
}
