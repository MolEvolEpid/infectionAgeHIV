#' @title Find infection age distributions
#' @description
#' Function to determine which biomarker values are usable in the multiple biomarker model and find distributions
#' for the amount of time between infection and diagnosis
#' @param df.single A data frame of last negative test date, first positive test date, and ART start date for each patient
#' @param df.multiple A data frame of biomarker values for each patient
#' @param n.adapt Number of iterations to adapt the model
#' @param n.burn Number of iterations to run before sampling starts
#' @param n.iter Number of iterations of sampling after burn-in is completed
#' @return A data frame with the inputs, quality controlled inputs (removing values too long after ART start),
#' and inferred distributions
#' @export
find.infection.ages <- function(df.single, df.multiple,
                                n.adapt = 1e4, n.burn = 1e4, n.iter = 1e5){
  #make sure patients in df.single are unique
  if(length(unique(df.single$patient_id)) != dim(df.single)[1]){
    stop("Patient IDs in df.single must be unique")
  } else{
    patient_ids <- df.single$patient_id
  }
  
  #number of individuals
  nInd <- length(patient_ids)
  
  #format dates
  df.single$last_neg_test_date <- as.Date(df.single$last_neg_test_date, format = "%Y-%m-%d")
  df.single$first_pos_test_date <- as.Date(df.single$first_pos_test_date, format = "%Y-%m-%d")
  df.single$ART_start_date <- as.Date(df.single$ART_start_date, format = "%Y-%m-%d")
  df.multiple$date <- as.Date(df.multiple$date, format = "%Y-%m-%d")
  
  #types of biomarker measurements
  biomarker_names <- c("BED", "LAg", "CD4", "pol", "pol2", "seq")
  
  #create list to store biomarker measurements
  biomarkers <- rep(list(vector(mode = "list", length = nInd)), length(biomarker_names))
  names(biomarkers) <- biomarker_names
  #list for biomarker measurement dates
  biomarker_dates <- biomarkers
  
  #extract data from df.multiple and put it into lists
  for(i in seq_along(patient_ids)){
    #find indices of df.multiple for this patient
    pat_indices <- (df.multiple$patient_id == patient_ids[i])
    for(j in seq_along(biomarker_names)){
      #find indices in df.multiple for these measurements
      indices <- which(pat_indices & df.multiple$name == biomarker_names[j])
      #find order of dates for these measurements
      order <- order(df.multiple$date[indices])
      #put biomarker measurement dates into list
      biomarker_dates[[biomarker_names[j]]][[i]] <- df.multiple$date[indices[order]]
      #put biomarker measurement values into list
      biomarkers[[biomarker_names[j]]][[i]] <- df.multiple$value[indices[order]]
    }
  }
  
  #check to make sure dates are logically consistent 
  #(e.g., last negative before first positive, ART start and measurements after diagnosis)
  if(!all(df.single$last_neg_test_date < df.single$first_pos_test_date, na.rm = TRUE)){
    stop("The last negative test date must be before the first positive test date")
  }
  if(any(df.single$ART_start_date < df.single$first_pos_test_date)){
    stop("The ART start date must not be before the first positive test date")
  }
  if(any(sapply(biomarker_dates$BED, FUN = function(x){if(length(x) > 0) min(x) else NA}) < df.single$first_pos_test_date, na.rm = TRUE)){
    stop("BED measurements must not be before the first positive test date")
  }
  if(any(sapply(biomarker_dates$LAg, FUN = function(x){if(length(x) > 0) min(x) else NA}) < df.single$first_pos_test_date, na.rm = TRUE)){
    stop("LAg measurements must not be before the first positive test date")
  }
  if(any(sapply(biomarker_dates$CD4, FUN = function(x){if(length(x) > 0) min(x) else NA}) < df.single$first_pos_test_date, na.rm = TRUE)){
    stop("CD4 measurements must not be before the first positive test date")
  }
  if(any(sapply(biomarker_dates$pol, FUN = function(x){if(length(x) > 0) min(x) else NA}) < df.single$first_pos_test_date, na.rm = TRUE)){
    stop("pol measurements must not be before the first positive test date")
  }
  if(any(sapply(biomarker_dates$pol2, FUN = function(x){if(length(x) > 0) min(x) else NA}) < df.single$first_pos_test_date, na.rm = TRUE)){
    stop("pol2 measurements must not be before the first positive test date")
  }
  if(any(sapply(biomarker_dates$seq, FUN = function(x){if(length(x) > 0) min(x) else NA}) < df.single$first_pos_test_date, na.rm = TRUE)){
    stop("seq measurements must not be before the first positive test date")
  }
  
  df <- prepare.HIV.data(patient_ID = df.single$patient_id, 
                         last_neg_test_date = df.single$last_neg_test_date,
                         first_pos_test_date = df.single$first_pos_test_date,
                         ART_start_date = df.single$ART_start_date,
                         BED_dates = biomarker_dates$BED,
                         BED = lapply(biomarkers$BED, FUN = as.numeric),
                         LAg_dates = biomarker_dates$LAg,
                         LAg = lapply(biomarkers$LAg, FUN = as.numeric),
                         CD4_dates = biomarker_dates$CD4,
                         CD4 = lapply(biomarkers$CD4, FUN = as.numeric),
                         seq_dates = biomarker_dates$seq,
                         seqs = biomarkers$seq,
                         pol2_dates = biomarker_dates$pol2,
                         pol2 = lapply(biomarkers$pol2, FUN = as.numeric),
                         find_infection_age_distributions = TRUE,
                         cluster_ID = rep(1,nInd),
                         n.adapt = n.adapt,
                         n.burn = n.burn,
                         n.iter = n.iter)
  return(df)
}