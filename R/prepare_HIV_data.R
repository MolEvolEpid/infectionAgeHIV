#'@title Prepare data for/use multiple biomarker model
#'@description Function to put biomarker and other date into a dataframe, calculate the polypmorphism count,
#' and (possibly) run the multiple biomarker model on it
#'@param patient_ID An identifier for each patient
#'@param last_neg_test_date The last negative test date, if available
#'@param first_pos_test_date The first positive test date. It is overriden by the earliest biomarker sample date if that is earlier
#'or if the first positive test date is unavailable.
#'@param ART_start_date The start date of antiretroviral therapy
#'@param BED_dates A list of vectors for the sample dates of the BED values
#'@param BED A list of vectors for values for the BED test
#'@param LAg_dates A list of vectors for the sample dates of the LAg values
#'@param LAg A list of vectors for values for the LAg test
#'@param CD4_dates A list of vectors for the sample dates of the CD4 counts
#'@param CD4 A list of vectors for the CD4+ T-cell counts
#'@param seq_dates A list of vectors for the sample dates of the HIV sequences
#'@param seq A list of HIV sequences (or a list of lists of sequences) for each patient in character or DNAbin format
#'@param seq_names The names of the sequences
#'@param pol_override A list of vectors for values for the polymorphism count if they have already been calculated externally. 
#'These values will replace the pol values calculated by this function.
#'@param seq_length_override A list of vectors for the sequence lengths 
#'(to be used along with pol_override if they have already been calculated externally)
#'@param pol2_dates A list of vectors for the sample dates of the pol2 values
#'@param pol A list of vectors for values for the pol2 value
#'@param VL_dates A list of vectors for the sample dates of the viral load
#'@param VL A list of vectors for values for the viral load
#'@param ART_cutoff_delay The amount of time (in days) after the start of antiretroviral therapy 
#'after which biomarker values should be disregarded
#'@param cluster_ID An identifier for the transmission cluster that the patient belongs to
#'@param gender A vector for the genders of each patient
#'@param age_at_sampling A list of vectors for the age of each patient at the time of each sequence
#'@param birth_location A vector for the birth location of each patient
#'@param suspected_infection_location A vector for the suspected infection location of each patient
#'@param risk_group A vector for the transmission risk group of each patient
#'@param aids_diagnosis_date A vector for the AIDS diagnosis date of each patient if they have been diagnosed with AIDS
#'@param death_date The date of death of each patient if the patient has died
#'@param date_format A character string specifying the format of the dates (to be passed to as.Date())
#'@param find_infection_age_distributions Whether or not to run the multiple biomarker model (MBM) on the data
#'@param prior.type The type of prior used for the time between infection and diagnosis in the multiple biomarker model.
#'1 indicates a gamma distribution with mean 2 years and standard deviation 1.5 years.
#'2 indicates a continuous uniform distribution with minimum 0 and maximum 12 years.
#'3 indicates a distribution for an individual that is HIV positive but has not developed AIDS symptoms, 
#'assuming that AIDS symptoms develop after a length of time according to a gamma distribution with shape 3.349 and rate 0.327
#'4 indicates a user-supplied distribution from the parameter user.prior.pdf
#'@param user.prior.pdf A pdf to use as the prior distribution for the time between infection and diagnosis in the multiple biomarker model.
#'It must be a list with x and y components corresponding to the time before diagnosis and probability density.
#'@param n.adapt Number of adaptation iterations if the MBM is run
#'@param n.burn Number of burn-in iterations if the MBM is run
#'@param n.iter Number of sampling iterations if the MBM is run
#'@param seed The RNG seed
#'@param ... Additional lists of patient data to put into the dataframe
#'@return A dataframe containing the input data, quality controlled versions of the data, possibly infection age distributions
#'@export
prepare.HIV.data <- function(patient_ID,
                             last_neg_test_date = NULL,
                             first_pos_test_date = NULL,
                             ART_start_date = NULL,
                             BED_dates = NULL,
                             BED = NULL,
                             LAg_dates = NULL,
                             LAg = NULL,
                             CD4_dates = NULL,
                             CD4 = NULL,
                             seq_dates,
                             seqs,
                             seq_names = NULL,
                             pol_override = NULL,
                             seq_length_override = NULL,
                             pol2_dates = NULL,
                             pol2 = NULL, 
                             VL_dates = NULL,
                             VL = NULL,
                             ART_cutoff_delay = 3, #in days
                             cluster_ID = NULL,
                             gender = NULL,
                             age_at_sampling = NULL,
                             birth_location = NULL,
                             suspected_infection_location = NULL,
                             risk_group = NULL,
                             aids_diagnosis_date = NULL,
                             death_date = NULL,
                             date_format = "%Y-%m-%d",
                             find_infection_age_distributions = FALSE,
                             prior.type = 1,
                             user.prior.pdf = list(x = c(0, 10), y = c(1/10, 1/10)),
                             n.adapt = 1e4, 
                             n.burn = 1e5, 
                             n.iter = 1e6,
                             seed = sample(2^31-1, 1),
                             ...){
  
  #find total number of individuals
  total_inds <- length(patient_ID)
  
  #check inputs
  
  #convert patient IDs to characters if they are not already
  patient_ID <- as.character(patient_ID)
  
  #convert dates to years (decimal valued, so July 1st 2000 would be approximately 2000.5 (2000.497))
  if(!is.null(last_neg_test_date)){
    last_neg_test_date <- as.double(as.Date(last_neg_test_date, date_format))/365.25 + 1970
  } else{
    last_neg_test_date <- rep(NA, total_inds)
  }
  if(!is.null(first_pos_test_date)){
    first_pos_test_date <- as.double(as.Date(first_pos_test_date, date_format))/365.25 + 1970
  } else{
    first_pos_test_date <- rep(NA, total_inds)
  }
  if(!is.null(ART_start_date)){
    ART_start_date <- as.double(as.Date(ART_start_date, date_format))/365.25 + 1970
  } else{
    ART_start_date <- rep(NA, total_inds)
  }
  if(!is.null(aids_diagnosis_date)){
    aids_diagnosis_date <- as.double(as.Date(aids_diagnosis_date, date_format))/365.25 + 1970
  } else{
    aids_diagnosis_date <- rep(NA, total_inds)
  }
  if(!is.null(death_date)){
    death_date <- as.double(as.Date(death_date, date_format))/365.25 + 1970
  } else{
    death_date <- rep(NA, total_inds)
  }
  
  #function to parse dates/replace empty values with NA
  parse_dates <- function(x, date_format){
    if(length(x) >=1) as.double(as.Date(x, date_format))/365.25+1970
    else NA
  }
  
  #BED sample dates and values
  if(!is.null(BED_dates)){
    BED_dates <- lapply(as.list(BED_dates), 
                        FUN = parse_dates, 
                        date_format = date_format)
  } else{
    BED_dates <- as.list(rep(NA, total_inds))
  }
  if(!is.null(BED)){
    BED <- lapply(BED, FUN = function(x) if(length(x) >= 1) as.numeric(x) else NA)
  } else{
    BED <- as.list(rep(NA, total_inds))
  }
  
  #LAg sample dates and values
  if(!is.null(LAg_dates)){
    LAg_dates <- lapply(as.list(LAg_dates), 
                        FUN = parse_dates, 
                        date_format = date_format)
  } else{
    LAg_dates <- as.list(rep(NA, total_inds))
  }
  if(!is.null(LAg)){
    LAg <- lapply(LAg, FUN = function(x) if(length(x) >= 1) as.numeric(x) else NA)
  } else{
    LAg <- as.list(rep(NA, total_inds))
  }
  
  #CD4 sample dates and values
  if(!is.null(CD4_dates)){
    CD4_dates <- lapply(as.list(CD4_dates), 
                        FUN = parse_dates, 
                        date_format = date_format)
  } else{
    CD4_dates <- as.list(rep(NA, total_inds))
  }
  if(!is.null(CD4)){
    CD4 <- lapply(CD4, FUN = function(x) if(length(x) >= 1) as.numeric(x) else NA)
  } else{
    CD4 <- as.list(rep(NA, total_inds))
  }
  
  #sequence sample dates and values
  if(!is.null(seq_dates)){
    seq_dates <- lapply(as.list(seq_dates), 
                        FUN = parse_dates, 
                        date_format = date_format)
  } else{
    stop("seq_dates must be provided")
  }
  if(!is.null(seqs)){
    seq <- lapply(BED, FUN = function(x) if(length(x) >= 1) x else NA)
  } else{
    seqs <- as.list(rep("n", total_inds))
    warning("No sequences found. Using blank sequences.")
  }
  
  if(!is.null(seq_names)){
    seq_names <- as.list(seq_names)
  } else{
    n_seqs <- lapply(seqs, FUN = length)
    #print(n_seqs)
    seq_names <- mapply(FUN = function(a,b) if(length(b) > 0) paste(a,b,sep = "_") else character(0), 
                        a = patient_ID, b = lapply(n_seqs, FUN = seq_len))
  }
  
  #pol2 sample dates and values
  if(!is.null(pol2_dates)){
    pol2_dates <- lapply(as.list(pol2_dates), 
                         FUN = parse_dates, 
                         date_format = date_format)
  } else{
    pol2_dates <- as.list(rep(NA, total_inds))
  }
  if(!is.null(pol2)){
    pol2 <- lapply(pol2, FUN = function(x) if(length(x) >= 1) as.numeric(x) else NA)
  } else{
    pol2 <- as.list(rep(NA, total_inds))
  }
  
  #viral load sample dates and values
  if(!is.null(VL_dates)){
    VL_dates <- lapply(as.list(VL_dates), 
                       FUN = parse_dates, 
                       date_format = date_format)
  } else{
    VL_dates <- as.list(rep(NA, total_inds))
  }
  if(!is.null(VL)){
    VL <- lapply(VL, FUN = function(x) if(length(x) >= 1) as.numeric(x) else NA)
  } else{
    VL <- as.list(rep(NA, total_inds))
  }
  
  #cluster_ID
  if(!is.null(cluster_ID)){
    cluster_ID <- as.character(cluster_ID)
  } else{
    stop("cluster_ID must be provided")
  }
  
  #gender
  if(!is.null(gender)){
    gender <- as.character(gender)
  } else{
    gender <- as.character(rep(NA, total_inds))
  }
  
  #age at the time of sequence samples
  if(!is.null(age_at_sampling)){
    age_at_sampling <- as.list(age_at_sampling)
  } else{
    age_at_sampling <- as.list(rep(NA, total_inds))
  }
  
  #birth location
  if(!is.null(birth_location)){
    birth_location <- as.character(birth_location)
  } else{
    birth_location <- as.character(rep(NA, total_inds))
  }
  
  #suspected infection location
  if(!is.null(suspected_infection_location)){
    suspected_infection_location <- as.character(suspected_infection_location)
  } else{
    suspected_infection_location <- as.character(rep(NA, total_inds))
  }
  
  #risk group
  if(!is.null(risk_group)){
    risk_group <- as.character(risk_group)
  } else{
    risk_group <- as.character(rep(NA, total_inds))
  }
  
  #use the earliest sampled biomarker as the first positive date if it is blank or after biomarker measurements (add VL dates?)
  first_pos_test_date_adj <- mapply(FUN = min, 
                                    first_pos_test_date, 
                                    ART_start_date,
                                    BED_dates, 
                                    LAg_dates,
                                    CD4_dates,
                                    seq_dates, 
                                    pol2_dates, 
                                    aids_diagnosis_date, 
                                    na.rm = TRUE)
  
  #calculate polymorphism count from sequences if not provided
  if(is.null(pol_override)){
    #collapse inner lists in seqs 
    #TODO check more about what format seqs is in
    #seqs <- lapply(seqs, FUN = strsplit, split = "")
    #function to calculate polymorphism counts
    calculate_pol <- function(sequences){
      nseq <- length(sequences)
      length_seq <- integer(length = nseq)
      pol <- double(length = nseq)
      for(i in seq_len(nseq)){
        seq.split <- strsplit(tolower(as.character(sequences[[i]])), split = "")
        frequencies <- table(seq.split)
        polymorphic_count <- sum(frequencies[names(frequencies) != "-" & names(frequencies) != "n" & names(frequencies) != "a" & 
                                               names(frequencies) != "c" & names(frequencies) != "g" & names(frequencies) != "t"])
        length_seq[i] <- sum(frequencies[names(frequencies) != "-" & names(frequencies) != "n"]) #length of sequence excluding gaps and Ns
        pol[i] <- polymorphic_count/length_seq[i]
        #change NaN to NA
        if(is.na(pol[i])) pol[i] <- NA
      }
      return(list(length = length_seq, pol = pol))
    }
    length_and_pol <- lapply(seqs, FUN = calculate_pol)
    seq_length <- lapply(length_and_pol, FUN = function(x) x$length)
    pol <- lapply(length_and_pol, FUN = function(x) x$pol)
  } else{ #override pol values if provided
    pol <- as.list(pol_override)
    seq_length = as.list(seq_length_override)
  }
  
  #find number of sequences per individual
  nSeq <- sapply(pol, FUN = length)
  
  #convert ART cutoff delay from days to years
  ART_cutoff_delay_years <- ART_cutoff_delay/365.25
  #remove biomarker values that are taken too late after the start of ART
  
  #function to remove values that are unusable due to being too long after ART start
  remove_late_samples <- function(samples, sample_dates, ART_dates, ART_cutoff_delay_years){
    samples[sample_dates > ART_dates+ART_cutoff_delay_years] <- NA
    return(samples)
  }
  
  BED_qc <- mapply(remove_late_samples,
                   samples = BED, sample_dates = BED_dates, ART_dates = ART_start_date, ART_cutoff_delay_years = ART_cutoff_delay_years, 
                   SIMPLIFY = FALSE)
  LAg_qc <- mapply(remove_late_samples, 
                   samples = LAg, sample_dates = LAg_dates, ART_dates = ART_start_date, ART_cutoff_delay_years = ART_cutoff_delay_years, 
                   SIMPLIFY = FALSE)
  CD4_qc <- mapply(remove_late_samples, 
                   samples = CD4, sample_dates = CD4_dates, ART_dates = ART_start_date, ART_cutoff_delay_years = ART_cutoff_delay_years, 
                   SIMPLIFY = FALSE)
  BED_qc <- mapply(remove_late_samples, 
                   samples = BED, sample_dates = BED_dates, ART_dates = ART_start_date, ART_cutoff_delay_years = ART_cutoff_delay_years, 
                   SIMPLIFY = FALSE)
  pol_qc <- mapply(remove_late_samples, 
                   samples = pol, sample_dates = seq_dates, ART_dates = ART_start_date, ART_cutoff_delay_years = ART_cutoff_delay_years, 
                   SIMPLIFY = FALSE)
  pol2_qc <- mapply(remove_late_samples, 
                    samples = pol2, sample_dates = pol2_dates, ART_dates = ART_start_date, ART_cutoff_delay_years = ART_cutoff_delay_years, 
                    SIMPLIFY = FALSE)
  VL_qc <- mapply(remove_late_samples, 
                  samples = VL, sample_dates = VL_dates, ART_dates = ART_start_date, ART_cutoff_delay_years = ART_cutoff_delay_years, 
                  SIMPLIFY = FALSE)
  
  #is there at least one usable pol value
  usable_pol <- sapply(pol_qc, FUN = function(x) !all(is.na(x)))
  
  #make sure other things are character vectors
  cluster_ID <- as.character(cluster_ID)
  gender <- as.character(gender)
  birth_location <- as.character(birth_location)
  suspected_infection_location <- as.character(suspected_infection_location)
  risk_group <- as.character(risk_group)
  
  #put everything into a data.frame
  df <- data.frame(patient_ID = patient_ID,
                   last_neg_test_date = last_neg_test_date,
                   first_pos_test_date = first_pos_test_date,
                   first_pos_test_date_adj = first_pos_test_date_adj,
                   ART_start_date = ART_start_date,
                   BED_dates = I(BED_dates),
                   BED = I(BED),
                   BED_qc = I(BED_qc),
                   LAg_dates = I(LAg_dates),
                   LAg = I(LAg),
                   LAg_qc = I(LAg_qc),
                   CD4_dates = I(CD4_dates),
                   CD4 = I(CD4),
                   CD4_qc = I(CD4_qc),
                   seq_dates = I(seq_dates),
                   seqs = I(seqs),
                   seq_names = I(seq_names),
                   seq_length = I(seq_length),
                   nSeq = nSeq,
                   pol = I(pol),
                   pol_qc = I(pol_qc),
                   usable_pol = usable_pol,
                   pol2_dates = I(pol2_dates),
                   pol2 = I(pol2), 
                   pol2_qc = I(pol2_qc),
                   VL_dates = I(VL_dates),
                   VL = I(VL),
                   VL_qc = I(VL_qc),
                   ART_cutoff_delay_years = rep(ART_cutoff_delay_years, length = total_inds), #in days
                   cluster_ID = cluster_ID,
                   gender = gender,
                   age_at_sampling = I(age_at_sampling),
                   birth_location = birth_location,
                   suspected_infection_location = suspected_infection_location,
                   risk_group = risk_group,
                   aids_diagnosis_date = aids_diagnosis_date,
                   death_date = death_date,
                   stringsAsFactors = FALSE)
  
  #add additional user supplied columns
  additional_input <- list(...)
  if("temp_col_name" %in% names(additional_input)) stop("Please use a name other than 'temp_col_name' for additional input columns")
  #return(additional_input)
  for(i in seq_along(additional_input)){
    df$temp_col_name <- I(additional_input[[i]])
    names(df)[names(df) == "temp_col_name"] <- names(additional_input)[i]
  }
  
  if(find_infection_age_distributions){
    #find probability distributions for the time between infection and diagnosis
    infection_age_dists <- run.mbm(df, n.adapt = n.adapt, n.burn = n.burn, n.iter = n.iter, 
                                   prior.type = prior.type, user.prior.pdf = user.prior.pdf,
                                   overall.seed = seed)
    #return(infection_age_dists$infection_age_dists_diag)
    #find cdfs
    cdfs <- lapply(infection_age_dists$infection_age_dists_diag, FUN = find.cdf)
    
    #find icdfs
    icdfs <- lapply(cdfs, FUN = find.icdf)
    
    #find pdf and cdf in terms of real time
    #cdfs_rt <- lapply(infection_age_dists$infection_age_dists_diag, FUN = function(dist){find_cdf(-rev(dist$x), rev(dist$y))})
    real_times <- mapply(FUN = function(diag, dist){rev(diag - dist$x)}, 
                         diag = df$first_pos_test_date_adj, dist = infection_age_dists$infection_age_dists_diag, SIMPLIFY = FALSE)
    #put times and (forward time) cdfs together
    cdf <- vector(mode = "list", length = length(cdfs))
    for(i in seq_along(cdfs)){
      cdf[[i]]$x <- real_times[[i]]
      cdf[[i]]$y <- rev(cdfs[[i]]$y)
    }
    #put times and (forward time) pdfs together
    pdf <- vector(mode = "list", length = length(cdfs))
    for(i in seq_along(cdfs)){
      pdf[[i]]$x <- real_times[[i]]
      pdf[[i]]$y <- rev(infection_age_dists$infection_age_dists_diag[[i]]$y)
    }
    
    #put distributions into dataframe
    #distributions from diagnosis
    df$infection_age_dists_diag <- I(infection_age_dists$infection_age_dists_diag)
    #distributions from first sequence
    df$infection_age_dists_seq <- I(infection_age_dists$infection_age_dists_seq)
    #inverse cdf distributions
    df$icdf <- I(icdfs)
    #forwards, real time pdf
    df$pdf <- I(pdf)
    #forwards, real time cdf
    df$cdf <- I(cdf)
  }
  
  return(df)
}

#get the cdf from the pdf for the infection age distribution
find.cdf <- function(pdf, tolerance = 0.01){
  #unnormalized pdf
  cdf_unnorm <- numeric(length = length(pdf$x))
  cdf_unnorm[1] <- 0
  x_diffs <- diff(pdf$x)
  for(i in 2:length(pdf$x)){ 
    cdf_unnorm[i] <- cdf_unnorm[i-1] + mean(pdf$y[i],pdf$y[i-1])*x_diffs[i-1] #linear approximation integral
  }
  if(max(cdf_unnorm) > 1+tolerance || max(cdf_unnorm) < 1-tolerance) stop("CDF error tolerance exceeded")
  cdf_norm <- cdf_unnorm/max(cdf_unnorm)
  return(list(x = pdf$x, y = cdf_norm))
}

#find unique values in numerical cdf
find.unique.cdf.values <- function(cdf.num){
  uniques <- sapply(unique(cdf.num$y), FUN = function(unique_cdf_y, cdf.num){
    which(unique_cdf_y == cdf.num$y)[1]}, 
    cdf.num) #first of each value
  cdf.num.unique <- list(x = cdf.num$x[uniques], y = cdf.num$y[uniques])
}

#function to make an inverse cdf function
find.icdf <- function(cdf){
  cdf.num.unique <- find.unique.cdf.values(cdf)
  icdf <- approxfun(x = cdf.num.unique$y, y = cdf.num.unique$x, 
                    method = "linear", yleft = min(cdf.num.unique$x), yright = max(cdf.num.unique$x))
  return(icdf)
}