# Ivarima

`Ivarima` Description: This package provides an automated method for identifying optimal intervention impact patterns in interventional ARIMA models. 
The approach incorporates structure-guided matching rules to improve delay time identification by extending the Linear Transfer Function (LTF) method to non-ideal conditions.
 It further integrates a stepwise ARIMAX search strategy to automatically determine the noise component and evaluate multiple candidate models. 
The optimal impact pattern is selected based on a set of model adequacy and fit criteria.

## Installation

```r
install.packages("devtools")
devtools::install_github("XYanSerenity/Ivarima")
```

## Example

```r
library(Ivarima)

y <- c(
  4.9, 5.2, 5.7, 5.7, 6.2, 6.7, 6.9, 7.1, 6.6, 7.0,
  6.9, 6.4, 6.6, 6.4, 7.0, 7.3, 6.0, 6.3, 4.8, 5.3,
  5.4, 4.7, 4.9, 4.4, 5.1, 5.3, 6.0, 5.9, 5.9, 5.6,
  5.3, 4.5, 4.7, 4.6, 4.3, 5.0, 5.2, 6.2, 5.8, 6.7,
  5.7, 6.1, 7.2, 6.5, 6.1, 6.3, 6.4, 7.0, 7.6, 7.2,
  7.5, 7.8, 7.2, 7.5, 5.6, 5.7, 4.9, 5.1, 6.2, 6.0,
  6.1, 7.5, 7.8, 8.0, 8.0, 8.1, 7.6, 7.1, 6.6, 5.6,
  5.9, 6.6, 6.8, 7.8, 7.9, 8.7, 7.7, 7.3, 6.7, 7.5,
  6.4, 9.7, 7.5, 7.1, 6.4, 6.0, 5.7, 5.0, 4.2, 5.1,
  5.4, 5.1, 5.3, 5.0, 4.8, 4.7, 5.0, 5.4, 4.3, 3.5
)

result <- auto_identify_impact_pattern(
  y = y,
  intervention_point = 82,
  k = 11
)

result$selected_model$best_model 


```
