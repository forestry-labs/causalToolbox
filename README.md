`causalToolbox` provides functions for estimating heterogenous treatment effects.

## How to install

The latest development version can be installed directly from Github using [devtools](https://github.com/hadley/devtools):

```R
if (!require("devtools")) install.packages("devtools")
devtools::install_github("forestry-labs/causalToolbox")
```

The package contains compiled code, and you must have a development environment to install the development version. (Use `devtools::has_devel()` to check whether you do.) If no development environment exists, Windows users download and install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) and macOS users download and install [Xcode](https://itunes.apple.com/us/app/xcode/id497799835).


## Example

For details please read our paper: https://arxiv.org/pdf/1706.03461.pdf
This example will not execute because there is no data.

```R
library(causalToolbox)
packageVersion("causalToolbox")

# create example data set
simulated_experiment <- simulate_causal_experiment(
    ntrain = 1000,
    ntest = 1000,
    dim = 10)
    
feature_train <- simulated_experiment$feat_tr
w_train <- simulated_experiment$W_tr
yobs_train <- simulated_experiment$Yobs_tr

# create the hte object using honest Random Forests (RF)
xl_rf <- X_RF(feat = feature_train, tr = w_train, yobs = yobs_train)

# alternatively, use BART instead of honest Random Forests. If you are not going
# to be careful about hyperparemeter tuning, we suggest using BART.
xl_bart <- X_BART(feat = feature_train, tr = w_train, yobs = yobs_train)

# estimate the CATE
feature_test <- simulated_experiment$feat_te

cate_esti_rf <- EstimateCate(xl_rf, feature_test)
cate_esti_bart <- EstimateCate(xl_bart, feature_test)

# evaluate the performance
cate_true <- simulated_experiment$tau_te
mean((cate_esti_rf - cate_true) ^ 2)
mean((cate_esti_bart - cate_true) ^ 2)

# Create confidence intervals via bootstrapping. Note that honest Random
# Forests have theoretically valid CIs, BART does not. However, BART often
# performs well for prediction without extensive tuning. See our paper
# (https://arxiv.org/abs/1706.03461)

# do not run (takes a long time)
xl_ci_rf <- CateCI(xl_rf, feature_test, B = 500)
xl_ci_bart <- CateCI(xl_bart, feature_test, B = 500)
```

