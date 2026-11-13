#' Encode DNA sequences as categorical nucleotide features
#'
#' Converts equal-length DNA sequences into one factor column per
#' nucleotide position. Use uppercase A, T, C and G.
#'
#' @param dna_strings A non-empty character vector of DNA sequences
#'   of equal length. For the supplied prediction model, use five-base
#'   sequences with A at the central position.
#' @return A data frame with one row per sequence and columns named
#'   nt_pos1, nt_pos2, and so on. Each column is a factor with levels
#'   A, T, C and G. Unsupported characters become NA.
#' @keywords internal
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(
    unlist(strsplit(dna_strings, "")),
    ncol = nn,
    byrow = TRUE
  )
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(
    seq_df, factor, levels = c("A", "T", "C", "G")
  )
  return(seq_df)
}


#' Predict m6A sites for multiple observations
#'
#' Uses a fitted random forest model to predict m6A probabilities
#' and assign Positive or Negative status to each input row.
#'
#' @param ml_fit A fitted randomForest classification model supporting
#'   probability prediction with a class named Positive.
#' @param feature_df A data frame containing gc_content, RNA_type,
#'   RNA_region, exon_length, distance_to_junction,
#'   evolutionary_conservation and DNA_5mer.
#'   RNA_type must be mRNA, lincRNA, lncRNA or pseudogene.
#'   RNA_region must be CDS, intron, 3'UTR or 5'UTR.
#'   DNA_5mer must contain uppercase five-base DNA sequences with
#'   A at the central position. Numeric features must use the same
#'   preprocessing as the training data; the course exon length and
#'   junction distance inputs are already log2-transformed.
#' @param positive_threshold Numeric probability threshold between
#'   0 and 1. A score strictly greater than this threshold is Positive;
#'   a score equal to or below it is Negative. Defaults to 0.5.
#' @return A data frame containing the input columns, with RNA_type
#'   and RNA_region converted to factors, plus predicted_m6A_prob
#'   (numeric) and predicted_m6A_status (character).
#' @import randomForest
#' @importFrom stats predict
#' @examples
#' ml_fit <- readRDS(system.file(
#'   "extdata", "rf_fit.rds",
#'   package = "m6APrediction", mustWork = TRUE
#' ))
#' feature_df <- read.csv(system.file(
#'   "extdata", "m6A_input_example.csv",
#'   package = "m6APrediction", mustWork = TRUE
#' ))
#' prediction_multiple(
#'   ml_fit,
#'   feature_df[1:3, ],
#'   positive_threshold = 0.6
#' )
#' @export
prediction_multiple <- function(
    ml_fit, feature_df, positive_threshold = 0.5){

  stopifnot(all(c(
    "gc_content", "RNA_type", "RNA_region", "exon_length",
    "distance_to_junction", "evolutionary_conservation", "DNA_5mer"
  ) %in% colnames(feature_df)))

  feature_df$RNA_type <- factor(
    feature_df$RNA_type,
    levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")
  )

  feature_df$RNA_region <- factor(
    feature_df$RNA_region,
    levels = c("CDS", "intron", "3'UTR", "5'UTR")
  )

  infered_prob <- predict(
    ml_fit,
    newdata = cbind(
      feature_df,
      dna_encoding(feature_df$DNA_5mer)
    ),
    type = "prob"
  )

  feature_df$predicted_m6A_prob <- infered_prob[, "Positive"]
  feature_df$predicted_m6A_status <- "Negative"

  feature_df$predicted_m6A_status[
    feature_df$predicted_m6A_prob > positive_threshold
  ] <- "Positive"

  return(feature_df)
}


#' Predict m6A status for a single site
#'
#' Creates a feature table for one site and returns its predicted
#' m6A probability and classification.
#'
#' @param ml_fit A fitted randomForest classification model supporting
#'   probability prediction with a class named Positive.
#' @param gc_content Numeric GC content for one site, using the
#'   same scale as the training data.
#' @param RNA_type Character RNA type: mRNA, lincRNA, lncRNA
#'   or pseudogene.
#' @param RNA_region Character RNA region: CDS, intron, 3'UTR
#'   or 5'UTR.
#' @param exon_length Numeric exon length already log2-transformed
#'   as in the course training data.
#' @param distance_to_junction Numeric distance to a junction already
#'   log2-transformed as in the course training data.
#' @param evolutionary_conservation Numeric conservation score,
#'   using the same scale as the training data.
#' @param DNA_5mer A single uppercase five-base DNA sequence using
#'   A, T, C and G, with A at the central position, such as GGACA.
#' @param positive_threshold Numeric probability threshold between
#'   0 and 1. A score strictly greater than this threshold is Positive;
#'   otherwise it is Negative. Defaults to 0.5.
#' @return A named character vector of length two, with elements
#'   predicted_m6A_prob and predicted_m6A_status.
#'   The probability is converted to character because the vector
#'   also contains the classification text.
#' @export
#' @examples
#' ml_fit <- readRDS(system.file(
#'   "extdata", "rf_fit.rds",
#'   package = "m6APrediction", mustWork = TRUE
#' ))
#' prediction_single(
#'   ml_fit,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 12,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "GGACA",
#'   positive_threshold = 0.5
#' )
prediction_single <- function(
    ml_fit, gc_content, RNA_type, RNA_region, exon_length,
    distance_to_junction, evolutionary_conservation, DNA_5mer,
    positive_threshold = 0.5){

  pred_outcome <- prediction_multiple(
    ml_fit = ml_fit,
    data.frame(
      gc_content = gc_content,
      RNA_type = factor(
        RNA_type,
        levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")
      ),
      RNA_region = factor(
        RNA_region,
        levels = c("CDS", "intron", "3'UTR", "5'UTR")
      ),
      exon_length = exon_length,
      distance_to_junction = distance_to_junction,
      evolutionary_conservation = evolutionary_conservation,
      DNA_5mer = DNA_5mer
    ),
    positive_threshold = positive_threshold
  )

  returned_vector <- c(
    pred_outcome$predicted_m6A_prob,
    pred_outcome$predicted_m6A_status
  )

  names(returned_vector) <- c(
    "predicted_m6A_prob",
    "predicted_m6A_status"
  )

  return(returned_vector)
}

