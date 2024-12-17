#####GECKO - Geographical Ecology and Conservation Knowledge Online
#####Version 1.0.0 (2022-05-28)
#####By Vasco Branco, Pedro Cardoso, Luís Correia
#####Maintainer: vasco.branco@helsinki.fi
#####Changed from v0.1.2:
#####Full release. Added new SPECTRE functionalities.
#####Renamed gecko.examples to gecko.data

################################################################################
##################################MAIN FUNCTIONS################################
################################################################################

#' Create a confusion matrix
#' @description Create a confusion matrix for any multiclass set of predicted vs observed labels
#' in a classification problem.
#' @param actual dataframe. Original labels.
#' @param predicted dataframe. Predicted labels.
#' @return data.frame. Predicted labels (rows) x Observed labels (cols).
#' @examples
#' x = c("FALSE", "TRUE", "FALSE", "TRUE", "TRUE")
#' y = c("TRUE", "TRUE", "TRUE", "TRUE", "TRUE")
#' confusion.matrix(x, y)
#' @export
confusion.matrix = function (actual, predicted){
  # predicted rows vs expected columns
  if ((is(actual, "factor") && is(predicted, "factor"))) {
    classes <- levels(actual)
  } else {
    classes <- unique(c(unique(actual), unique(predicted)))
  }
  error_frame <- data.frame(0, 0, 0, 0)
  for (i in seq(length(actual))) {
    if (actual[i] == classes[1] && actual[i] == predicted[i]) {
      error_frame[1] <- error_frame[1] + 1
    } else if (actual[i] == classes[2] && actual[i] == predicted[i]) {
      error_frame[4] <- error_frame[4] + 1
    } else if (actual[i] == classes[1] && predicted[i] == classes[2]) {
      error_frame[2] <- error_frame[2] + 1
    } else if (actual[i] == classes[2] && predicted[i] == classes[1]) {
      error_frame[3] <- error_frame[3] + 1
    }
  }
  output <- matrix(error_frame, ncol = 2)
  colnames(output) <- classes
  row.names(output) <- classes
  return(output)
}

#' Performance of model predictions
#' @description Calculate the performance of a model through a comparison 
#' between predicted and observed labels. Available metrics are \code{accuracy},
#' \code{F1} and \code{TSS}.
#' @param actual dataframe. Same formatting as \code{y}, containg some sort of classification data. 
#' @param predicted dataframe. Same formatting as \code{x}, containg the predicted classifications of a model trained over the data in \code{x}.
#' @param metric character. String specifying the metric used, one of \code{accuracy}, \code{F1} and \code{TSS}.
#' @details
#' \strong{The F-score or F-measure (F1)} is: \cr
#' \cr
#' \eqn{F1 = 2 \dfrac{Precision * Recall}{Precision + Recall}}, with \cr
#' \cr
#' \eqn{Precision = \dfrac{True Positive}{True Positive + False Positive}} \cr
#' \cr
#' \eqn{Recall = \dfrac{True Positive}{True Positive + False Negative}} \cr
#' \cr
#' \strong{Accuracy} is: \cr
#' \cr
#' \eqn{\dfrac{100 * (True Postives + True Negatives)}{True Postives + True Negatives + False Positives + False Negatives}}
#' \cr
#'  \cr
#' \strong{The Pierce's skill score (PSS),  Bookmaker's Informedness (BM) or True Skill Statistic (TSS)} is: \cr
#' \cr
#' \eqn{TSS = TPR + TNR - 1}, \cr
#' with \eqn{TPR} being the True Positive Rate, positives correctly labelled 
#' as such and \eqn{TNR}, the True Negative Rate, the rate of negatives correctly
#' labelled, such that:\cr
#' \cr
#' \eqn{TPR = \dfrac{True Positives}{True Positives + False Negatives}}
#' \cr
#' \eqn{TNR = \dfrac{True Negatives}{True Negatives + False Positives}}
#' \cr
#' Take in consideration the fact that the F1 score is not a robust metric in datasets with class imbalances.
#' @references
#' PSS:
#' Peirce, C. S. (1884). The numerical measure of the success of predictions. Science, 4, 453–454.
#' @return numeric.
#' @examples
#' observed = c("FALSE", "TRUE", "FALSE", "TRUE", "TRUE")
#' predicted = c("TRUE", "TRUE", "TRUE", "TRUE", "TRUE")
#' performance.metrics(observed, predicted, "TSS")
#' @export
performance.metrics = function (actual, predicted, metric){
  if (metric == "accuracy") {
   error_count <- 0
   for (i in seq(length(actual))) {
     if (actual[i] == predicted[i]) {} else {
       error_count <- error_count + 1
     }
   }
   accuracy <- round((length(actual) - error_count) / length(actual) * 100, 2)
   # cat(paste0("Predictions have an accuracy of ", accuracy, "% \n"))
   return(accuracy)
 } else if (metric == "F1") {
   # And Multiclass?
   cmat <- confusion.matrix(actual, predicted)
   cmat <- as.data.frame(cmat)

   precision <- cmat[[2, 2]] / (cmat[[2, 2]] + cmat[[2, 1]])
   recall <- cmat[[2, 2]] / (cmat[[2, 2]] + cmat[[1, 2]])
   F1 <- 2 * (precision * recall) / (precision + recall)
   print(paste("F1 score is ", F1, ".", sep = ""))
   return(F1)
 } else if (metric == "TSS") {
   cmat <- confusion.matrix(actual, predicted)
   cmat <- as.data.frame(cmat)

   recall <- cmat[[2, 2]] / (cmat[[2, 2]] + cmat[[1, 2]])
   TSS <- recall + (cmat[[1, 1]] / (cmat[[1, 1]] + cmat[[2, 1]])) - 1
   print(paste("True Skill Statistic score is ", TSS, ".", sep = ""))
   return(TSS)
 } else {
   warning("Unrecognized method. Please check documentation.")
   return(NULL)
 }
}


#' Split a dataset for model training
#' @description Split a dataset for model training while keeping class representativity.
#' @param data dataframe. Containg some sort of classification data. The last column
#' must contain the label data. 
#' @param proportion numeric. A value between 0 a 1 determining the proportion of the dataset 
#' split between training and testing.
#' @return list. First element is the train data, second element is the test data.
#' @examples
#' # Binary label case
#' my_data = data.frame(X = runif(20), Y = runif(20), Z = runif(20), Label =
#' c(rep("presence", 10), rep("outlier", 10)) )
#' splitDataset(my_data, 0.8)
#' 
#' # Multi label case
#' my_data = data.frame(X = runif(60), Y = runif(60), Z = runif(60), Label =
#' c(rep("A", 20), rep("B", 30), rep("C", 10)) )
#' splitDataset(my_data, 0.8)
#' @export
splitDataset = function(data, proportion){
  labels <- table(data[, ncol(data)])
  labels <- sort(labels, decreasing = TRUE)

  set <- matrix(nrow = 0, ncol = ncol(data), dimnames = list(c(), colnames(data)))
  set <- as.data.frame(set)

  presence_set_a <- set
  presence_set_b <- set

  for (i in 1:length(labels)) {
    presence_set <- data[data[ncol(data)] == names(labels[i]), ]
    presence_floor <- floor(nrow(presence_set) * proportion)
    if (presence_floor == 0) {
      presence_sample <- sample(1:nrow(presence_set), floor(nrow(presence_set) * proportion))
      presence_set_a <- rbind(presence_set_a, presence_set)
      warning(paste0(
        "Proportion invalid with current class count. All data of class ",
        labels[i], " placed under the training set."
      ))
    } else {
      presence_sample <- sample(1:nrow(presence_set), floor(nrow(presence_set) * proportion))
      presence_set_a <- rbind(presence_set_a, presence_set[presence_sample, ])
      presence_set_b <- rbind(presence_set_b, presence_set[-presence_sample, ])
    }
  }
  outcome_list <- list()
  outcome_list[[1]] <- presence_set_a
  outcome_list[[2]] <- presence_set_b
  return(outcome_list)
}


#' Detect outliers in a set of geographical coordinates
#' @description This function generates pseudo-abscences from an input data.frame
#' containing latitude and longitude coordinates by using environmental data and
#' then uses both presences and pseudo-absences to train a SVM model used to 
#' flag possible outliers for a given species.
#' @param longlat data.frame. With two columns containing latitude and longitude, describing
#' the locations of a species, which may contain outliers.
#' @param training data.frame. With the same formatting as \code{longlat}, indicating only known
#' locations where a target species occurs. Used exclusively as training data for
#' method 'svm'.
#' @param hi_res logical. Specifies if 1 KM resolution environmental data should be used. 
#' If \code{FALSE} 10 KM resolution data is used instead.
#' @param crop logical. Indicates whether environmental data should be cropped to
#'  an extent similar to what is given in \code{longlat} and \code{training}. Useful to avoid
#'  large processing times of higher resolutions.
#' @param threshold numeric. Value indicating the threshold for classifying 
#' outliers in methods \code{"geo"} and \code{"env"}. E.g.: under the default
#'  of 0.05, points that are at an average distance greater than the 95% quartil
#'  of the average distances of all points, will be classified as outliers.
#' @param method A string specifying the outlier detection method. \code{"geo"} 
#' calculates the euclidean distance between point coordinates and classifies as
#' outliers those outside the 0% and 95% interval in an assumed gaussian distribuition;
#' \code{"env"}
#' performs the same calculation but instead uses the environmental data extracted
#' from those points. \code{"svm"} will use the dataset given to \code{"longlat"} and it corresponding
#' extracted environmental data to train a support vector machine model that then
#' predicts outliers.
#' @details Environmental data used is WorldClim and requires a long download, see
#' \code{\link[gecko:gecko.setDir]{gecko::gecko.setDir()}}
#' This function is heavily based on the methods described in Liu et al. (2017). 
#' There the authors describe SVM_pdSDM, a pseudo-SDM method similar to a 
#' two-class presence only SVM that is capable of using pseudo-absence points, 
#' implemented with the ksvm function in the R package kernlab. 
#' It is suggested that, for each set of \code{"n"} occurence 
#' records, \code{"2 * n"} pseudo-absences points are generated.
#' Whilst using it keep in mind works highlighting limitations such as such as
#' Meynard et al. (2019). See References section.
#' @return list if \code{method = "all"}, containing whether or not a given point
#' was classified as \code{TRUE} or \code{FALSE} along with the confusion matrix
#' for the training data. If \code{method = "geo"} or 
#' \code{method = "env"} a data.frame is returned. 
#' @references Liu, C., White, M. and Newell, G. (2017) ‘Detecting outliers in species distribution data’, Journal of Biogeography, 45(1), pp. 164–176. doi:10.1111/jbi.13122. \cr
#' \cr
#' Meynard, C.N., Kaplan, D.M. and Leroy, B. (2019) ‘Detecting outliers in species distribution data: Some caveats and clarifications on a virtual species study’, Journal of Biogeography, 46(9), pp. 2141–2144. doi:10.1111/jbi.13626. \cr
#' @examples
#' \dontrun{
#' new_occurences = gecko.data("records")
#' old_occurences = data.frame(X = runif(10, -17.1, -17.05), Y = runif(10, 32.73, 32.76))
#' outliers.detect(new_occurences, old_occurences)
#' }
#' @export
outliers.detect = function(longlat, training = NULL, hi_res = TRUE, crop = FALSE, threshold = 0.05, method = "all"){
  if (hi_res){
    res = "1 km"
  } else {
    res = "10 km"
  }
  
  # needs to call for worldclim to be downloaded
  data_dir <- gecko.getDir()
  if (is.null(data_dir) && (method %in% c("env", "all"))) {
    cat("One or more methods you have selected require downloading large data.
    This is a large dataset that is prone to fail by timeout if downloaded
    through R. You can do this process now or later by running gecko.setDir() 
    and placing the files needed inside it after download them through your browser.
    See detailed instructions for this in the documention of gecko.worldclim().
    Would you like to continue?\n")
    answer <- readline("[y/n]?")
    if (tolower(answer) == "y") {
      gecko.setDir()
      gecko.worldclim(res) # download the data
    } else {
      return(NULL)
    }
  }
  
  if (!file.exists(paste0(gecko.getDir(), "worldclim/", res))){
    gecko.worldclim(res)
  }

  if (method %in% c("env", "all")) {
    worldclim_stack <- terra::rast(dir(paste0(gecko.getDir(), "worldclim/", res), full.names = T))
    if (crop && !is.null(training)){
      all_data = rbind(as.matrix(longlat),as.matrix(training))
      worldclim_stack = terra::crop(worldclim_stack, c(min(all_data[,1]), max(all_data[,1]), min(all_data[,2]), max(all_data[,2])))
    }
    
    
  }

  # error handling -------------------------------------------------------------
  if (!is(longlat, "data.frame") || (!is(training, "data.frame") && !is.null(training))) {
    warning("Both longlat and training (if not NULL) need to be of type data.frame.")
    return(NULL)
  }

  if (ncol(longlat) != 2) {
    warning("Both longlat and training (if not NULL) need to have only two columns with latitude and longitude.")
    return(NULL)
  }

  if (!is.null(training)) {
    if (ncol(training) != 2) {
      warning("Both longlat and training (if not NULL) need to have only two columns with latitude and longitude.")
      return(NULL)
    }
  }

  points <- longlat
  # points = as.data.frame(points)
  points[, 1] <- as.numeric(points[, 1])
  points[, 2] <- as.numeric(points[, 2])
  colnames(points) <- c("x_coords", "y_coords")
  out <- points

  # Euclid dist with environmental data ----------------------------------------
  if (method == "env" || method == "all") {
    env_data <- terra::extract(worldclim_stack,
      ID = FALSE,
      terra::vect(points, geom = c("x_coords", "y_coords"))
    )
    env_data <- scale(env_data)
    out <- cbind(out, "env" = dist.euclid(env_data, threshold))
  }

  # Euclid dist with just the coordinates --------------------------------------
  if (method == "geo" || method == "all") {
    out <- cbind(out, "geo" = dist.euclid(points, threshold))
  }

  # SVM ------------------------------------------------------------------------
  if (method == "svm" || method == "all") {
    if (is.null(training)) {
      warning("No training data has been supplied. Skipping 'svm'.")
      return(out)
    }
    ## Generating pseudo-absence points with biomod.
    ## Liu recommends that PA.nb.absences = 2 * number of points
    ## Does raster size affect the time needed?
    ## Should have an option for a different environmental stack

    myBiomodData <- biomod2::BIOMOD_FormatingData(
      resp.name = "species_presence",
      resp.var = c(rep(1, nrow(training))),
      expl.var = worldclim_stack, # used to be worldclim_stack
      resp.xy = training,
      PA.nb.rep = 1,
      PA.nb.absences = 2 * nrow(training),
      PA.strategy = "random",
      PA.dist.min = 200
    )

    ## from: https://rpubs.com/dgeorges/416446
    ## When you have extracted the PA table from data.formatted.Biomod.object you can
    ## easily select presences (filter(status == 1)), pseudo-absences (filter(is.na(status)))
    ## or absences (filter(status == 0)) even if no absences are defined in our case
    point_frame <- data.frame(
      myBiomodData@coord[, 1], myBiomodData@coord[, 2],
      myBiomodData@data.species
    )

    ## keep a copy of our pseudo-absence points for consulting later.
    pseudo_absence_points <- point_frame[is.na(point_frame$myBiomodData.data.species), 1:2]
    colnames(pseudo_absence_points) <- c("x_coords", "y_coords")
    ## reclassify all our classifications
    for (row in seq(nrow(point_frame))) {
      if (is.na(point_frame[row, 3])) {
        point_frame[row, 3] <- "outlier"
      } else {
        point_frame[row, 3] <- "presence"
      }
    }
    colnames(point_frame) <- c("x_coords", "y_coords", "class")

    # Get environ. data for both datasets ----------------------------------------
    l_frames <- list()
    for (i in 1:2) {
      if (i == 1) {
        points <- point_frame
      } else {
        # longlat = gpt_matrix[,4:5]
        # points = longlat
        #
        points <- cbind(longlat, rep("presence", nrow(longlat)))
      }
      colnames(points) <- c("x_coords", "y_coords", "class")
      check_numeric <- tryCatch(
        {
          points[, 1] <- as.numeric(points[, 1])
          points[, 2] <- as.numeric(points[, 2])
          check_numeric <- FALSE
        },
        error = function(e) {
          message("Test data.frame must have data coerceable to numeric.")
          check_numeric <- TRUE
        }
      )
      if (check_numeric) {
        return(NULL)
      }

      env_data <- terra::extract(worldclim_stack, ID = FALSE, terra::vect(points, geom = c("x_coords", "y_coords")))

      # can't create model prediction with NA, remove data points with NA.
      # if training, silently remove. if longlat, append as outliers in the end
      na_rows <- c(1:nrow(env_data))[rowSums(is.na(env_data[, 1:ncol(env_data)])) > 0]
      if (i == 1) {
        train_NA <- na_rows
      } else {
        GPT_NA <- na_rows
      }

      if (length(na_rows) != 0) {
        if (i == 1) {
        } else {
          env_data <- env_data[-na_rows, ]
          points <- points[-na_rows, ]
        }
      }

      points <- data.frame(env_data, points[, ncol(points)])
      colname_list <- list()
      for (item in seq(terra::nlyr(worldclim_stack))) {
        colname_list[[item]] <- paste("env_", item, sep = "")
      }
      # points = points[,-c(1)]
      colnames(points) <- c(colname_list, "species_presence")
      l_frames[[i]] <- points
    }

    split <- splitDataset(l_frames[[1]], 0.8)
    training_set <- split[[1]]
    validation_set <- split[[2]]


    ## SVM_pdSDM (Liu et al. 2017) with kernlab. ##
    training_set$species_presence <- as.factor(training_set$species_presence)
    validation_set$species_presence <- as.factor(validation_set$species_presence)

    pseudo_svm <- kernlab::ksvm(
      x = species_presence~.,
      data = training_set
    )
    
    ## Evaluating results

    model_predictions <- kernlab::predict(pseudo_svm, validation_set)
    # should be 
    cat("SVM TRAINING \n")
    val_accuracy <- performance.metrics(
      validation_set[, ncol(validation_set)],
      model_predictions, 
      "accuracy"
    )
    cat(paste0("Base model (training data):"))
    cat(paste0("Accuracy: ", val_accuracy, "% \n"))
    

    val_mat <- confusion.matrix(
      validation_set[, ncol(validation_set)],
      model_predictions
    )



    # Predict with model data --------------------------------------------------
    model_predictions <- kernlab::predict(pseudo_svm, l_frames[[2]])

    svm_out <- rep(NA, length(model_predictions))
    for (l in 1:length(model_predictions)) {
      temp <- row.names(l_frames[[2]])[l]
      # this 'presence' / 'absence' should be changed
      if (as.character(model_predictions[l]) == "presence") {
        svm_out[as.numeric(temp)] <- FALSE
      } else {
        svm_out[as.numeric(temp)] <- TRUE
      }
    }

    out <- cbind(out, "svm" = svm_out)
    sum <- 0
    for (r in 3:ncol(out)) {
      sum <- sum + out[, r]
    }
    out <- cbind(out, "possible.outliers" = sum)
    
    
    cat("SVM TESTING \n")
    performance.metrics(
      factor(l_frames[[2]][, ncol(l_frames[[2]])], levels = c("outlier", "presence")),
      model_predictions, "accuracy"
    )

    return(list(
      outlier_matrix = out,
      SVM_traindata_cmatrix = val_mat
    ))
  }

  return(out)
}
