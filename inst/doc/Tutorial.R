## ----style, echo = FALSE, results = 'asis'------------------------------------
  BiocStyle::markdown()

## ----setup, include=FALSE, echo=F, warning= F, message=F----------------------
knitr::opts_chunk$set(
  message = FALSE, #fig.height = 6, fig.width = 6,
  warning = FALSE, 
  error = FALSE, 
  tidy = FALSE,
  fig.align = "center", 
  #  dpi = 600, 
  cache = TRUE,
  progress = FALSE, 
  quite = TRUE
)

library(BiocStyle)
#BiocStyle::latex(relative.path = TRUE)

library(knitr)

## ---- message=FALSE-----------------------------------------------------------
# Install the released version from CRAN using
if (!requireNamespace("multiclassPairs", quietly = TRUE)) {
  install.packages("multiclassPairs")
}

# Or install the dev version from GitHub using
# if (!requireNamespace("multiclassPairs", quietly = TRUE)) {
#  if (!requireNamespace("devtools", quietly = TRUE)) {
#    install.packages("devtools")
#  }
#  library(devtools) # this package is needed to install from GitHub
#  install_github("NourMarzouka/multiclassPairs")
#}

# Install the dependencies from Bioconductor
# BiocManager, Biobase, and switchBox packages from Bioconductor are needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("Biobase", quietly = TRUE)) {
  BiocManager::install("Biobase")
}
if (!requireNamespace("switchBox", quietly = TRUE)) {
  BiocManager::install("switchBox")
}

# load multiclassPairs library
library(multiclassPairs)

## ---- echo=FALSE,fig.cap="Workflow in multiclassPairs R package: Functions are colored in green."----
knitr::include_graphics("images/workflow_v0_2_1.png")

## ---- message=FALSE-----------------------------------------------------------
library(multiclassPairs)

# example of creating data object from matrix
# we will generate fake data in this example
# matrix with random values
Data <- matrix(runif(100000), 
               nrow=100, 
               ncol=100, 
               dimnames = list(paste0("G",1:100), 
                               paste0("S",1:100)))
# class labels
L1 <- sample(x = c("A","B","C"), size = 100, replace = TRUE)

# platform/study labels
P1 <- sample(x = c("P1","P2"), size = 100, replace = TRUE)

# create the data object
object <- ReadData(Data = Data,
                   Labels = L1,
                   Platform = P1,
                   verbose = FALSE)
object

## ---- message=FALSE-----------------------------------------------------------
library(multiclassPairs, quietly = TRUE)

# install Leukemia cancers data
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!requireNamespace("leukemiasEset", quietly = TRUE)){
  BiocManager::install("leukemiasEset")
}

# load the data package
library(leukemiasEset, quietly = TRUE)
data(leukemiasEset)

## ---- message=FALSE-----------------------------------------------------------
# check the Expressionset
leukemiasEset

# explore the phenotypes data
knitr::kable(head(pData(leukemiasEset)))

# We are interested in LeukemiaType
knitr::kable(table(pData(leukemiasEset)[,"LeukemiaType"]))

# split the data
# 60% as training data and 40% as testing data
n <- ncol(leukemiasEset)
set.seed(1234)
training_samples <- sample(1:n,size = n*0.6)

train <- leukemiasEset[1:1000,training_samples]
test  <- leukemiasEset[1:1000,-training_samples]

# just to be sure there are no shared samples between the training and testing data
sum(sampleNames(test) %in% sampleNames(train)) == 0

# create the data object
# when we use Expressionset we can use the name of the phenotypes variable 
# ReadData will automatically extract the phenotype variable and use it as class labels
# the same can be used with the Platform/study labels
# in this example we are not using any platform labels, so leave it NULL
object <- ReadData(Data = train, 
                   Labels = "LeukemiaType", 
                   Platform = NULL, 
                   verbose = FALSE)
object

## ---- echo=FALSE, fig.cap="One-vs-rest scheme"--------------------------------
knitr::include_graphics("images/one_vs_rest_scheme.png")

## ---- echo=FALSE, fig.cap="Gene filtering options"----------------------------
knitr::include_graphics("images/gene_filtering_TSP.png")

## ---- echo=FALSE, fig.cap="Platform-wise gene filtering"----------------------
knitr::include_graphics("images/platform_wise_gene_filtering_TSP.png")

## -----------------------------------------------------------------------------
# let's go with gene filtering using one_vs_one option
# for featureNo argument, a sufficient number of returned features is 
# recommended if large number of rules is used in the downstream training steps.
filtered_genes <- filter_genes_TSP(data_object = object,
                                   filter = "one_vs_one",
                                   platform_wise = FALSE,
                                   featureNo = 1000,
                                   UpDown = TRUE,
                                   verbose = TRUE)
filtered_genes

## -----------------------------------------------------------------------------
# Let's train our model
classifier <- train_one_vs_rest_TSP(data_object = object,
                                    filtered_genes = filtered_genes,
                                    k_range = 5:50,
                                    include_pivot = FALSE,
                                    one_vs_one_scores = TRUE,
                                    platform_wise_scores = FALSE,
                                    seed = 1234,
                                    verbose = FALSE)
classifier

## -----------------------------------------------------------------------------
# apply on the training data
# To have the classes in output in specific order, we can use classes argument
results_train <- predict_one_vs_rest_TSP(classifier = classifier,
                                         Data = object,
                                         tolerate_missed_genes = TRUE,
                                         weighted_votes = TRUE,
                                         classes = c("ALL","AML","CLL","CML","NoL"),
                                         verbose = TRUE)

# apply on the testing data
results_test <- predict_one_vs_rest_TSP(classifier = classifier,
                                        Data = test,
                                        tolerate_missed_genes = TRUE,
                                        weighted_votes = TRUE,
                                        classes=c("ALL","AML","CLL","CML","NoL"),
                                        verbose = TRUE)
# get a look over the scores in the testing data
knitr::kable(head(results_test))

## -----------------------------------------------------------------------------
# Confusion Matrix and Statistics on training data
caret::confusionMatrix(data = factor(results_train$max_score, 
                                     levels = unique(object$data$Labels)),
                       reference = factor(object$data$Labels, 
                                          levels = unique(object$data$Labels)),
                       mode="everything")

# Confusion Matrix and Statistics on testing data
caret::confusionMatrix(data = factor(results_test$max_score, 
                                     levels = unique(object$data$Labels)),
                       reference = factor(pData(test)[,"LeukemiaType"], 
                                          levels = unique(object$data$Labels)),
                       mode="everything")

## ---- include=FALSE, results="hide"-------------------------------------------
# Confusion Matrix and Statistics on training data
x1 <- caret::confusionMatrix(data = factor(results_train$max_score, 
                                           levels = unique(object$data$Labels)),
                             reference = factor(object$data$Labels, 
                                                levels = unique(object$data$Labels)),
                             mode="everything")

# Confusion Matrix and Statistics on testing data
x2 <- caret::confusionMatrix(data = factor(results_test$max_score, 
                                           levels = unique(object$data$Labels)),
                             reference = factor(pData(test)[,"LeukemiaType"], 
                                                levels = unique(object$data$Labels)),
                             mode="everything")


## -----------------------------------------------------------------------------
# plot for the rules and scores in the training data
plot_binary_TSP(Data = object, # we are using the data object here
                classifier = classifier, 
                prediction = results_train, 
                classes = c("ALL","AML","CLL","CML","NoL"),
                #margin = c(0,5,0,10),
                title = "Training data")

# plot for the rules and scores in the testing data
plot_binary_TSP(Data = test, # ExpressionSet
                ref = "LeukemiaType", # ref label names in pData
                classifier = classifier, 
                prediction = results_test, 
                classes = c("ALL","AML","CLL","CML","NoL"),
                title = "Testing data"#, 
                #margin = c(0,5,0,10)
                )

## -----------------------------------------------------------------------------
# (500 trees here just for fast example)
genes_RF <- sort_genes_RF(data_object = object,
                          # featureNo_altogether, it is better not to specify a number here
                          # featureNo_one_vs_rest, it is better not to specify a number here
                          rank_data = TRUE,
                          platform_wise = FALSE,
                          num.trees = 500, # more features, more tress are recommended
                          seed=123456, # for reproducibility
                          verbose = TRUE)
genes_RF # sorted genes object

## -----------------------------------------------------------------------------
# to get an idea of how many genes we will use
# and how many rules will be generated
summary_genes <- summary_genes_RF(sorted_genes_RF = genes_RF,
                                  genes_altogether = c(10,20,50,100,150,200),
                                  genes_one_vs_rest = c(10,20,50,100,150,200))
knitr::kable(summary_genes)

# 50 genes_altogether and 50 genes_one_vs_rest seems 
# to give enough number of  rules and unique genes for our classes
# (500 trees here just for fast example)
# Now let's run sort_rules_RF to create the rules and sort them
rules_RF <- sort_rules_RF(data_object = object, 
                          sorted_genes_RF = genes_RF,
                          genes_altogether = 50,
                          genes_one_vs_rest = 50, 
                          num.trees = 500,# more rules, more tress are recommended 
                          seed=123456,
                          verbose = TRUE)
rules_RF # sorted rules object

## ---- results="hide"----------------------------------------------------------
# prepare the simple data.frame for the parameters I want to test
# names of arguments as column names
# this df has three sets (3 rows) of parameters
parameters <- data.frame(
  gene_repetition=c(3,2,1),
  rules_one_vs_rest=c(2,3,10),
  rules_altogether=c(2,3,10),
  run_boruta=c(FALSE,"make_error",TRUE), # I want to produce error in the 2nd trial
  plot_boruta = FALSE,
  num.trees=c(100,200,300),
  stringsAsFactors = FALSE)

# parameters
# for overall and byclass possible options, check the help files
para_opt <- optimize_RF(data_object = object,
                        sorted_rules_RF = rules_RF,
                        parameters = parameters,
                        test_object = NULL,
                        overall = c("Accuracy","Kappa"), # wanted overall measurements 
                        byclass = c("F1"), # wanted measurements per class
                        verbose = TRUE)

para_opt # results object
# para_opt$summary # the df of with summarized information
knitr::kable(para_opt$summary)

## -----------------------------------------------------------------------------
# train the final model
# it is preferred to increase the number of trees and rules in case you have
# large number of samples and features
# for quick example, we have small number of trees and rules here
# based on the optimize_RF results we will select the parameters
RF_classifier <- train_RF(data_object = object,
                          sorted_rules_RF = rules_RF,
                          gene_repetition = 1,
                          rules_altogether = 10,
                          rules_one_vs_rest = 10,
                          run_boruta = TRUE, 
                          plot_boruta = FALSE,
                          probability = TRUE,
                          num.trees = 300,
                          boruta_args = list(),
                          verbose = TRUE)

## -----------------------------------------------------------------------------
# plot proximity matrix of the out-of-bag samples
# Note: this takes a lot of time if the data is big
proximity_matrix_RF(object = object,
             classifier = RF_classifier, 
             plot = TRUE,
             return_matrix = FALSE, # if we need to get the matrix itself
             title = "Leukemias",
             cluster_cols = TRUE)

## -----------------------------------------------------------------------------
# training accuracy
# get the prediction labels from the trained model
# if the classifier trained using probability	= FALSE
training_pred <- RF_classifier$RF_scheme$RF_classifier$predictions
if (is.factor(training_pred)) {
  x <- as.character(training_pred)
}

# if the classifier trained using probability	= TRUE
if (is.matrix(training_pred)) {
  x <- colnames(training_pred)[max.col(training_pred)]
}

# training accuracy
caret::confusionMatrix(data =factor(x),
                       reference = factor(object$data$Labels),
                       mode = "everything")

## ---- include=FALSE, results="hide"-------------------------------------------
# training accuracy
x1 <- caret::confusionMatrix(data =factor(x),
                             reference = factor(object$data$Labels),
                             mode = "everything")

## -----------------------------------------------------------------------------
# apply on test data
results <- predict_RF(classifier = RF_classifier, 
                      Data = test,
                      impute = TRUE) # can handle missed genes by imputation

# get the prediction labels
# if the classifier trained using probability	= FALSE
test_pred <- results$predictions
if (is.factor(test_pred)) {
  x <- as.character(test_pred)
}

# if the classifier trained using probability	= TRUE
if (is.matrix(test_pred)) {
  x <- colnames(test_pred)[max.col(test_pred)]
}

# training accuracy
caret::confusionMatrix(data = factor(x),
                       reference = factor(pData(test)[,"LeukemiaType"]),
                       mode = "everything")

## ---- include=FALSE, results="hide"-------------------------------------------
x2 <- caret::confusionMatrix(data = factor(x),
                             reference = factor(pData(test)[,"LeukemiaType"]),
                             mode = "everything")

## -----------------------------------------------------------------------------
#visualize the binary rules in training dataset
plot_binary_RF(Data = object,
               classifier = RF_classifier,
               prediction = NULL, 
               as_training = TRUE, # to extract the scores from the model
               show_scores = TRUE,
               top_anno = "ref",
               show_predictions = TRUE, 
               #margin = c(0,5,0,8),
               title = "Training data")

# visualize the binary rules in testing dataset
plot_binary_RF(Data = test,
               ref = "LeukemiaType", # Get ref labels from the test ExpressionSet
               classifier = RF_classifier,
               prediction = results, 
               as_training = FALSE, 
               show_scores = TRUE,
               top_anno = "ref",
               show_predictions = TRUE,
               title = "Testing data")

## -----------------------------------------------------------------------------
sessionInfo()

