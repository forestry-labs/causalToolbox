devtools::load_all()
data <- get(load('~/Dropbox/01d-data_files.RData'))
data$treatment[data$treatment == 2] <- 0
gotv <- data[data$treatment == 0 | data$treatment == 1, ]
set.seed(2345)
sample_size <- floor(0.8 * nrow(gotv))
train_index <-  sample(seq_len(nrow(gotv)), size = sample_size)
train <- gotv[train_index, ]
test <- gotv[-train_index, ]

feature_tr <- train[,1:7]
w_tr <- train$treatment
yobs_tr <- train$voted

feature_test <- test[,1:7]
cate_true <- test$voted

t_rf <- T_RF(feat = feature_tr, tr = w_tr, yobs = yobs_tr)
t_bart_db <- T_BART(feat = feature_tr, tr = w_tr, yobs = yobs_tr, 
                 tree_package = "dbarts",
                 nthread = 4)
t_bart_bayes <- T_BART(feat = feature_tr, tr = w_tr, yobs = yobs_tr, 
                    tree_package = "BayesTree",
                    nthread = 4)

cate_esti_rf <- EstimateCate(t_rf, feature_test)
cate_esti_bart_db <- EstimateCate(t_bart_db, feature_test)
cate_esti_bart_bayes <- EstimateCate(t_bart_bayes, feature_test)

mean((cate_esti_rf - cate_true) ^ 2)
mean((cate_esti_bart_db - cate_true) ^ 2)
mean((cate_esti_bart_bayes - cate_true) ^ 2)