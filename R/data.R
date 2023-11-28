#' infectionAgeHIV Example Data
#' 
#' Example data for use in the vignette for infectionAgeHIV
#' 
#' @format ## `df.single.example`
#' A data frame with 4 rows and 4 columns.
#' Contains values for which there is only measurement per patient
#' \describe{
#'   \item{patient_id}{A unique identifier for the patient}
#'   \item{last_neg_test_date}{The patient's last negative HIV test date, if available}
#'   \item{first_pos_test_date}{The patient's first positive HIV test date}
#'   \item{ART_start_date}{The date the patient started antiretroviral therapy}
#' }
"df.single.example"
#' 
#' @format ## `df.multiple.example`
#' A data frame with 11 rows and 4 columns
#' Contains values for which there may be multiple measurements per patient
#' \describe{
#'   \item{patient_id}{The patient_id as in df.single}
#'   \item{name}{The name of the biomarker measurement (BED, LAg, CD4, pol, pol2, or seq)}
#'   \item{date}{The date of the biomarker measurement}
#'   \item{value}{The value of the biomarker measurement}
#' }
#' "df.multiple.example"