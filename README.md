# Re-sampled dimensional reduction (RSDR)

This package applies a resampling method to estimate rotated matrix
for dimensional reduction. By applying the procedure, a number of new
dimensions can be used as feature candidates for predictive modeling;
thus, the number of candidates can be optimized depending on the training
sample size. This helps to fulfill minimum events per variable (EPV) for a
machine learning algorithm while optimizing the proportion of variance
explained (PVE). Unlike other packages for dimensional reduction, this
package applies resampling methods to prevent overfitting. The PVE
optimization takes the predicted outcome into account without using it to
represent the features.

## Quick Start rsdr R

<a href="https://htmlpreview.github.io/?https://github.com/herdiantrisufriyana/rsdr/blob/main/vignettes/quick-start-R.html">
Read vignette for simple example in R</a>

<a href="https://github.com/herdiantrisufriyana/rsdr/blob/main/vignettes/quick-start.R">Download R script</a>

<a href="https://github.com/herdiantrisufriyana/medhist/blob/main/man/rsdr_0.1.0.pdf">Reference manual</a>