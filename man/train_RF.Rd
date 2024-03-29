\name{train_RF}
\alias{train_RF}
\title{
Train pair-based random forest model
}
\description{
\code{train_RF} trains random forest model based on binary gene rules (such as geneA<geneB).
Boruta package is used to remove the unimportant rules and ranger function from ranger package is used for the training.
}
\usage{
train_RF(data_object,
          sorted_rules_RF,
          gene_repetition = 1,
          rules_altogether = 50,
          rules_one_vs_rest = 50,
          run_boruta = FALSE,
          plot_boruta = FALSE,
          boruta_args = list(doTrace = 1),
          num.trees = 500,
          min.node.size = 1,
          importance = "impurity",
          write.forest = TRUE,
          keep.inbag = TRUE,
          probability	= TRUE,
          verbose = TRUE, ...)
}
\arguments{
  \item{data_object}{
data object generated by ReadData function. Object contains the data and labels.
}
  \item{sorted_rules_RF}{
RandomForest_sorted_rules object generated by sort_rules_RF function
}
  \item{gene_repetition}{
interger indicating how many times the gene is allowed to be repeated in the pairs/rules. Default is 1.
  }
  \item{rules_altogether}{
integer indicating how many unique rules to be used from altogether slot in the sorted rules object. Default is 200.
}
  \item{rules_one_vs_rest}{
integer indicating how many unique rules to be used from each one_vs_rest slot (class vs rest slots) in the sorted rules object. Default is 200.
}
  \item{run_boruta}{
  logical indicates if Boruta algorithm should be run before building the RF model. Boruta will be used to remove the unimportant rules. Default is FALSE.
}
  \item{plot_boruta}{
   logical indicates if Boruta is allowed to plot importance history plots. Default is FALSE.
}
  \item{boruta_args}{
  list of argument to be passed to Boruta algorithm. Default for doTrace argument in Boruta is 1.
}
  \item{num.trees}{
an integer. Number of trees. Default is 500. It is recommended to increase num.trees in case of having large number of features (ranger function argument).
}
  \item{min.node.size}{
an integer. Minimal node size. Default is 1. (ranger function argument)
}
  \item{importance}{
Variable importance mode, should be one of 'impurity', 'impurity_corrected', 'permutation'.
Defualt is 'impurity' (ranger function argument)
}
  \item{write.forest}{
Save ranger.forest object, required for prediction. Default is TRUE. (ranger function argument). Always should be true to return the trained RF model.
}
  \item{keep.inbag}{
Save how often observations are in-bag in each tree. Default is TRUE. (ranger function argument). Needed for co-clustering heatmaps.
}
  \item{probability}{
Grow a probability forest as in Malley et al. (2012). Default is TRUE. (ranger function argument). Needed to plot probability scores in the binary rules heatmaps. If TRUE, when the classifier is used to predict a sample class the user will get "ranger.prediction" object with a matrix with scores for each class. If FALSE, the classifier will give a "ranger.prediction" object with the predicted class without scores for each class.
}
  \item{verbose}{
a logical value indicating whether processing messages will be printed or not. Default is TRUE.
}
  \item{\dots}{
any additional arguments to be passed to ranger function (i.e. random forest function) in ranger package. For example, seed for reproducibility. Note, seed argument will be used also for Burota run.
}
}
\details{
train_RF function extracts the lists of the sorted rules from altogether and classes slots, then it keep the top rules those fit with the gene_repetition number, this step reduces the number of the rule dramatically. From the left rules, rules_altogether and rules_one_vs_rest determine how many rules will be used. In case rules_altogether and rules_one_vs_rest were larger than the left rules then all the rules will be used. After that these rules will be pooled in one list and fid to Boruta function to remove the unimportant rules. Then random forest will be trained on the important rules.
}
\value{
train_RF returns rule_based_RandomForest object which contains the final RF classifier and the used genes and rules in the final model.

Boruta results are also included in the object.

The object also contains TrainingMatrix which is a binary matrix for the rules in the training data, this is used for imputation purposes during the prediction if the sample misses some values.
}
\author{
Nour-al-dain Marzouka <nour-al-dain.marzouka at med.lu.se>
}

\examples{# generate random data
Data <- matrix(runif(800), nrow=100, ncol=80,
               dimnames = list(paste0("G",1:100), paste0("S",1:80)))

# generate random labels
L <- sample(x = c("A","B","C","D"), size = 80, replace = TRUE)

# generate random platform labels
P <- sample(c("P1","P2","P3"), size = 80, replace = TRUE)

# create data object
object <- ReadData(Data = Data,
                   Labels = L,
                   Platform = P,
                   verbose = FALSE)

# sort genes
genes_RF <- sort_genes_RF(data_object = object,
                          seed=123456, verbose = FALSE)

# to get an idea of how many genes we will use
# and how many rules will be generated
# summary_genes_RF(sorted_genes_RF = genes_RF,
#                  genes_altogether = c(10,20,50,100,150,200),
#                  genes_one_vs_rest = c(10,20,50,100,150,200))

# creat and sort rules
# rules_RF <- sort_rules_RF(data_object = object,
#                           sorted_genes_RF = genes_RF,
#                           genes_altogether = 100,
#                           genes_one_vs_rest = 100,
#                           seed=123456,
#                           verbose = FALSE)

# parameters <- data.frame(
#   gene_repetition=c(3,2,1),
#   rules_one_vs_rest=0,
#   rules_altogether=c(2,3,10),
#   run_boruta=c(FALSE,"produce_error",FALSE),
#   plot_boruta = FALSE,
#   num.trees=c(100,200,300),
#   stringsAsFactors = FALSE)

# parameters

# test <- optimize_RF(data_object = object,
#                     sorted_rules_RF = rules_RF,
#                     test_object = NULL,
#                     overall = c("Accuracy"),
#                     byclass = NULL, verbose = FALSE,
#                     parameters = parameters)
# test
# test$summary[which.max(test$summary$Accuracy),]
#
# # train the final model
# # it is preferred to increase the number of trees and rules in case you have
# # large number of samples and features
# # for quick example, we have small number of trees and rules here
# # based on the optimize_RF results we will select the parameters
# RF_classifier <- train_RF(data_object = object,
#                           gene_repetition = 1,
#                           rules_altogether = 0,
#                           rules_one_vs_rest = 10,
#                           run_boruta = FALSE,
#                           plot_boruta = FALSE,
#                           probability = TRUE,
#                           num.trees = 300,
#                           sorted_rules_RF = rules_RF,
#                           boruta_args = list(),
#                           verbose = TRUE)
#
# # training accuracy
# # get the prediction labels
# # if the classifier trained using probability	= FALSE
# training_pred <- RF_classifier$RF_scheme$RF_classifier$predictions
# if (is.factor(training_pred)) {
#   x <- as.character(training_pred)
# }
#
# # if the classifier trained using probability	= TRUE
# if (is.matrix(training_pred)) {
#   x <- colnames(training_pred)[max.col(training_pred)]
# }
#
# # training accuracy
# caret::confusionMatrix(data =factor(x),
#                 reference = factor(object$data$Labels),
#                 mode = "everything")

# not to run
# visualize the binary rules in training dataset
# plot_binary_RF(Data = object,
#                classifier = RF_classifier,
#                prediction = NULL, as_training = TRUE,
#                show_scores = TRUE,
#                top_anno = "ref",
#                show_predictions = TRUE,
#                title = "Training data")

# not to run
# predict
# test_object # any test data
# results <- predict_RF(classifier = RF_classifier, impute = TRUE,
#                       Data = test_object)
#
# # visualize the binary rules in training dataset
# plot_binary_RF(Data = test_object,
#                classifier = RF_classifier,
#                prediction = results, as_training = FALSE,
#                show_scores = TRUE,
#                top_anno = "ref",
#                show_predictions = TRUE,
#                title = "Test data")
}
