# Install and load
  
## Install and load these packages
devtools::install_github('herdiantrisufriyana/rsdr')
library(rsdr)

devtools::install_github('herdiantrisufriyana/medhist')
library(medhist)

## Showing load packages
library(tidyverse)
dslabs::ds_theme_set()
library(pbapply)
library(parallel)
library(BiocGenerics)
library(lubridate)
library(Biobase)
library(caret)
library(MLeval)


# Load artificial data

## Load data
data(medhistdata)

## Split data for external validation
set.seed(33)
idx=createDataPartition(medhistdata$outcome,times=1,p=0.8)

## Construct train set
train_set=medhistdata[,idx$Resample1]

## Construct test set
test_set=medhistdata[,-idx$Resample1]

## Extract and transform medical history
ps_remover_train=extract_nps_mh(train_set)

mh_bin_nps_train=
  train_set[ps_remover_train$key,] %>%
  `exprs<-`(
    exprs(.) %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column(var='id') %>%
      column_to_rownames(var='id') %>%
      t()
  ) %>%
  trans_binary(verbose=F)

mh_bin_nps_test=
  test_set[ps_remover_train$key,] %>%
  `exprs<-`(
    exprs(.) %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column(var='id') %>%
      column_to_rownames(var='id') %>%
      t()
  ) %>%
  trans_binary(verbose=F)

## Describe dataset
cat('Input dimensions (excluding outcome):\n')
dim(mh_bin_nps_train)
table(`Outcome in the samples:`=mh_bin_nps_train$outcome)
floor(table(`Event/non-event per variable (EPV):`=mh_bin_nps_train$outcome)/28)
floor(table(`If 6 variables, then the EPVs are:`=mh_bin_nps_train$outcome)/6)


# Dimensional reduction

## Get table of features with outcome
input_train=
  mh_bin_nps_train %>%
  exprs() %>%
  t() %>%
  as.data.frame()

## Fit dimensional reduction models with resampling
rsdr_bin_nps_train=rsdr(input_train,'CV',10,'PCA')

## Show fitting results
rsdr_bin_nps_train

## Transform dimensions of train set
dr_bin_nps_train=
  mh_bin_nps_train %>%
  transformation(
    rsdr_object=rsdr_bin_nps_train
    ,top_n=6
    ,verbose=F
  )

## Transform dimensions of test set
dr_bin_nps_test=
  mh_bin_nps_test %>%
  transformation(
    rsdr_object=rsdr_bin_nps_train
    ,output_dim=rownames(dr_bin_nps_train)
    ,verbose=F
  )

## Get proportion of variance explained without outcome
pve_no_outcome=
  preproc(dr_bin_nps_train)$pve %>%
  data.frame(pve=.) %>%
  rownames_to_column(var='dim') %>%
  arrange(desc(pve)) %>%
  mutate(dim=factor(dim,unique(dim)))

## Proportion of variance explained without outcome

pve_no_outcome %>%
  head() %>%
  setNames(c('Dimension','PVE')) %>%
  knitr::kable(
    format='html'
    ,caption='Proportion of variance explained without outcome.'
    ,digits=2
  ) %>%
  kableExtra::kable_styling(full_width=T)

## Variance explained
pve_no_outcome %>%
  mutate(cpve=cumsum(pve)) %>%
  qplot(dim,cpve,data=.) +
  geom_vline(xintercept=6,lty=2) +
  geom_label(
    aes(label=ifelse(dim==pve_no_outcome$dim[6],round(cpve,2),NA))
    ,na.rm=T
  ) +
  scale_x_discrete('Output dimension') +
  scale_y_continuous(
    'Cumulative proportion of variance explained'
    ,breaks=seq(0,1,0.05)
  ) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

## Cluster of input dimensions
input_dim_order=
  colnames(input_train)[heatmap(t(input_train))$rowInd]

## Cluster of output dimensions
output_dim_order=
  rownames(dr_bin_nps_train)[heatmap(exprs(dr_bin_nps_train))$rowInd]

## Set up plotting parameters for variable description
input_dim_desc=
  mh_bin_nps_train %>%
  fData() %>%
  rename(input_dim_desc=desc) %>%
  rbind(
    data.frame(
      input_dim_desc='Preeclampsia/eclampsia'
      ,row.names='outcome'
    )
  ) %>%
  .[input_dim_order,,drop=F] %>%
  rownames_to_column(var='input_dim')

output_dim_desc=
  c('PC8','Problems with previous pregnancy'
    ,'PC3','Lacking/late antenatal screening'
    ,'PC4','Commond cold with lab check'
    ,'PC24','Various problems increasing visits'
    ,'PC9','High-risk pregnancy'
    ,'PC1','Unknown') %>%
  matrix(ncol=2,byrow=T) %>%
  data.frame(output_dim=.[,1],output_dim_desc=.[,2])


## Dimension transformer
rsdr_bin_nps_train %>%
  plot(
    w_cutoff=0.01
    ,label=T
    ,input_dim_order=input_dim_order
    ,output_dim_order=output_dim_order
    ,input_dim_desc=input_dim_desc
    ,output_dim_desc=output_dim_desc
  )


# Predictive performance

## Set up internal validation for training set
set.seed(66)
int_val=
  trainControl(
    method='boot'
    ,number=30
    ,savePredictions=T
    ,classProbs=T
    ,summaryFunction=twoClassSummary
    ,allowParallel=F
  )

## Training the models
set.seed(66)
rand_predictor_candidates=sample(rownames(mh_bin_nps_train),6,F)
mod_mh_bin_nps=
  suppressWarnings(train(
    outcome~.
    ,data=
      mh_bin_nps_train %>%
      exprs() %>%
      t() %>%
      as.data.frame() %>%
      select_at(rand_predictor_candidates) %>%
      cbind(data.frame(outcome=mh_bin_nps_train$outcome))
    ,method='glm'
    ,metric='ROC'
    ,trControl=int_val
    ,tuneLength=10
  ))

set.seed(66)
mod_dr_bin_nps=
  suppressWarnings(train(
    outcome~.
    ,data=
      dr_bin_nps_train %>%
      exprs() %>%
      t() %>%
      as.data.frame() %>%
      cbind(data.frame(outcome=dr_bin_nps_train$outcome))
    ,method='glm'
    ,metric='ROC'
    ,trControl=int_val
    ,tuneLength=10
  ))

## Evaluating the models using testing set
eval_mh_bin_nps=
  mod_mh_bin_nps %>%
  predict(
    newdata=
      mh_bin_nps_test %>%
      exprs() %>%
      t() %>%
      as.data.frame() %>%
      cbind(data.frame(outcome=mh_bin_nps_test$outcome))
    ,type='prob'
  ) %>%
  cbind(data.frame(obs=mh_bin_nps_test$outcome)) %>%
  evalm(showplots=F,silent=T) %>%
  .$optres %>%
  .$Group1 %>%
  .['AUC-ROC',]

eval_dr_bin_nps=
  mod_dr_bin_nps %>%
  predict(
    newdata=
      dr_bin_nps_test %>%
      exprs() %>%
      t() %>%
      as.data.frame() %>%
      cbind(data.frame(outcome=dr_bin_nps_test$outcome))
    ,type='prob'
  ) %>%
  cbind(data.frame(obs=dr_bin_nps_test$outcome)) %>%
  evalm(showplots=F,silent=T) %>%
  .$optres %>%
  .$Group1 %>%
  .['AUC-ROC',]

## Predictive performance using different feature representation

rbind(
  mutate(eval_mh_bin_nps,rsdr='No')
  ,mutate(eval_dr_bin_nps,rsdr='Yes')
) %>%
  select(rsdr,everything()) %>%
  rename(
    `Applying re-sampled dimensional reduction (RSDR)?`=rsdr
    ,AUROC=Score
  ) %>%
  knitr::kable(
    format='html'
    ,caption='Predictive performance using different feature representation.'
  ) %>%
  kableExtra::kable_styling(full_width=T)
