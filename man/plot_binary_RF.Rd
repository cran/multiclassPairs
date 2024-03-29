\name{plot_binary_RF}
\alias{plot_binary_RF}
\title{
  Plot binary rule-based heatmaps
}
\description{
  \code{plot_binary_RF} Plot binary heatmaps for datasets based on rule-based random forest classifier
}
\usage{
  plot_binary_RF(Data,
                  classifier,
                  ref              = NULL,
                  prediction       = NULL,
                  as_training      = FALSE,
                  platform         = NULL,
                  classes          = NULL,
                  platforms_ord    = NULL,
                  top_anno         = c("ref", "prediction", "platform")[1],
                  title            = "",
                  binary_col       = c("white", "black", "gray"),
                  ref_col          = NULL,
                  pred_col         = NULL,
                  platform_col     = NULL,
                  show_ref         = TRUE,
                  show_predictions = TRUE,
                  show_platform    = TRUE,
                  show_scores      = TRUE,
                  show_rule_name   = TRUE,
                  legend           = TRUE,
                  cluster_cols     = TRUE,
                  cluster_rows     = TRUE,
                  anno_height      = 0.03,
                  score_height     = 0.03,
                  margin           = c(0, 5, 0, 5))
}
\arguments{
  \item{Data}{
    a matrix or a dataframe, samples as columns and row as features/genes. Can also take ExpressionSet, or data_object generated by ReadData function.
  }
  \item{classifier}{
    Classifier as a rule_based_RandomForest object, generated by \code{train_RF} function
  }
  \item{ref}{
    Optional vector with the reference labels. Ref labels in data_object will be used if not ref input provided. For ExpressionSet, the name of the ref variable in the pheno data can be used.
  }
  \item{prediction}{
    Optional. "ranger.prediction" object for the class scores generated by \code{predict_RF} function.
  }
    \item{as_training}{
    Logical indicates if the plot is for the training data. It means the predictions will be extracted from the classifier itself and any prediction object will be ignored. If TRUE, then the training data/object should be used for Data argument.
  }

  \item{platform}{
    Optional vector with the platform/study labels or any additional annotation. Platform labels in data_object will be used if no platform input is provided. For ExpressionSet, the name of the variable in the pheno data can be used.
  }
  \item{classes}{
    Optional vector with class names. This will determine which classes will be plotted and in which order. It is not recommended to use both "classes" and "platforms_ord" arguments together to avoid sample order conflict and may result in an improper plotting for samples.
  }
  \item{platforms_ord}{
    Optional vector with the platform/study names. This will determine which platform/study will be plotted and in which order. This will be used when top_anno="platform". It is not recommended to use both "classes" and "platforms_ord" arguments together.
  }
  \item{top_anno}{
    Determine the top annotation level. Samples will be grouped based on the top_anno. Input can be one of three options: "ref", "prediction", "platform". Default is "ref".
  }
  \item{title}{
    Character input as a title for the whole heatmap. Default is "".
  }
  \item{binary_col}{
    Vector determines the colors of the binary heatmap. Default is c("white", "black", "gray"). First color for the color when the rule is false in the sample. Second color for the color when the rule is true. Third color is for NAs.
  }
  \item{ref_col}{
    Optional named vector determines the colors of classes for the reference labels. Default is NULL. Vector names should match with the ref labels.
  }
  \item{pred_col}{
    Optional named vector determines the colors of classes for the prediction labels. Default is NULL. Vector names should match with the prediction labels in the prediction labels.
  }
  \item{platform_col}{
    Optional named vector determines the colors of platforms/study labels. Default is NULL. Vector names should match with the platforms/study labels.
  }
  \item{show_ref}{
    Logical. Determines if the ref labels will be plotted or not. If the top_anno argument is "ref" then show_ref will be ignored and ref labels will be plotted.
  }
  \item{show_predictions}{
    Logical. Determines if the prediction labels will be plotted or not. If the top_anno argument is "prediction" then show_predictions will be ignored and predictions will be plotted.
  }
  \item{show_platform}{
    Logical. Determines if the platform/study labels will be plotted or not. If the top_anno argument is "platform" then show_platform will be ignored and platforms will be plotted.
  }
  \item{show_scores}{
    Logical. Determines if the prediction scores will be plotted or not. To visualize scores, the classifier should be trained with probability=TRUE otherwise show_scores will be turned FALSE automatically.
  }
  \item{show_rule_name}{
    Logical. Determines if the rule names will be plotted on the left side of the heatmapp or not.
  }
  \item{legend}{
    Logical. Determines if a legend will be plotted under the heatmap.
  }
  \item{cluster_cols}{
    Logical. Clustering the samples in each class (i.e. not all samples in the cohort) based on the binary rules for that class. If top_anno is "platform" then the rules from all classes are used to cluster the samples in each platform.
  }
  \item{cluster_rows}{
    Logical. Clustering the rules in each class.
  }
  \item{anno_height}{
    Determines the height of the annotations. It is recommended not to go out of this range 0.01<height<0.1. Default is 0.03.
  }
  \item{score_height}{
    Determines the height of the score bars. It is recommended not to go out of this range 0.01<height<0.1. Default is 0.03.
  }
  \item{margin}{
    Determines the margins of the heatmap. Default is c(0, 5, 0, 5).
  }
}

\value{
  returns a heatmap plot for the binary rule
}

\author{
  Nour-al-dain Marzouka <nour-al-dain.marzouka at med.lu.se>
}

\examples{# generate random data
Data <- matrix(runif(8000), nrow=100, ncol=80,
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

# Or you can use expand.grid to generate dataframe with all parameter combinations
# parameters <- expand.grid(
#   gene_repetition=c(3,2,1),
#   rules_one_vs_rest=0,
#   rules_altogether=c(2,3,10),
#   num.trees=c(100,500,1000),
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
# Extract and plot the proximity matrix from the classifier for the training data
# it takes long time for large data
# proximity_mat <- proximity_matrix_RF(object = object,
#                       classifier = RF_classifier,
#                       plot=TRUE,
#                       return_matrix=TRUE,
#                       title = "Test",
#                       cluster_cols = TRUE)

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
#                title = "Test data")}
