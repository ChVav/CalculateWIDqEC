# PMR calculation main work function
# contact: charlottevavourakis@gmail.com
# under development for Shiny App version for TirolPath
# Last update 19 January 2023

calculate_pmr <- function(data,#experimentname,
                          threshold_COL2A1=30,
                          threshold_targets=35,
                          threshold_targets_warning=30,
                          external_curve=FALSE, #these and parameters below are not used
                          curve=NULL,
                          fix_intercept=36.9,
                          fix_slope=-3.4,
                          calib_fixed=FALSE){
  
  # WID-qEC thresholds
  qEC_threshold1 <- 0.03
  #qEC_threshold2 <- 0.63
  
  # initialize list so multiple objects can be returned from function ----
  list_out <- list()
  
  # Check calibration with COL2A1 values ----
  
  # Generate standard curve
  plot <- data %>%
    filter(Target.Name == "COL2A1" & ((Sample.Name %in% c("STD1", "STD2", "STD3","STD4") | 
                                    (Sample.Name %in% c("STD_1", "STD_2", "STD_3", "STD_4")) | 
                                    (Sample.Name %in% c("Std 1", "Std 2", "Std 3", "Std 4")) | 
                                    (Sample.Name %in% c("Std_1", "Std_2", "Std_3", "Std_4")) | 
                                    (Sample.Name %in% c("Std. 1", "Std. 2", "Std. 3", "Std. 4"))))) %>%
    mutate(x = case_when(Sample.Name %in% c("STD1", "STD_1", "Std 1", "Std_1", "Std. 1") ~ 6, # concentration = log10(copy numbers/5uL)
                         Sample.Name %in% c("STD2", "STD_2", "Std 2", "Std_2", "Std. 2") ~ 5,
                         Sample.Name %in% c("STD3", "STD_3", "Std 3", "Std_3", "Std. 3") ~ 4,
                         Sample.Name %in% c("STD4", "STD_4", "Std 4", "Std_4", "Std. 4") ~ 3)) %>%
    mutate(ct = ifelse(CT == "Undetermined", NA, as.numeric(CT)))
  
  fit <- lm(ct ~ x, data = plot)
  anno <- paste0("R2 = ", signif(summary(fit)$r.squared, 3),
                 "\ny = ", round(fit$coefficients[2],2) , "x + ", round(fit$coefficients[1],2))
  y = max(plot$ct) - 0.05*max(plot$ct)
  intercept <- unname(fit$coefficients[1])
  slope <- unname(fit$coefficients[2])
  R2 <- summary(fit)$r.squared
  
  # Plot calibration curve
  curveplot <- plot %>%
    ggplot(aes(x = x,
               y = ct)) +
    geom_smooth(method = "lm",
                formula = y ~ x,
                size = 0.5,
                se = FALSE,
                alpha = 0.3,
                colour = "gray40") +
    geom_point(aes(colour = Sample.Name),
               size = 0.75) +
    theme(legend.title = element_blank(),
          legend.position = "top",
          legend.key = element_blank()) +
    xlab("Concentration \n [ log(copy number/5uL) ]") +
    ylab("Ct value")  +
    annotate("text",
             x = 5,
             y = 28,
             label = anno,
             hjust = 0,
             size = 3,
             colour = ifelse(R2 < 0.99, "red", "black"))
  
  # add plot to list
  list_out[[1]] <- curveplot
  
  rm(plot, curveplot)
  
  # If user choose fixed values overwrite the coefficients from standard curve
  if(calib_fixed==TRUE){
    intercept <- fix_intercept
    slope <- fix_slope
  }
  
  # reformat data in results file based on thresholds and replication ----
  data <- data %>%
    mutate(rep = case_when(grepl("A|C|E|G|I|K|M|O", Well.Position) == TRUE ~ 1,
                           grepl("B|D|F|H|J|L|N|P", Well.Position) == TRUE ~ 2)) %>%
    mutate(ct = as.numeric(ifelse(CT > threshold_targets | Amp.Status=="No Amp", NA, CT))) # set suspicious amplifications to NA for targets
  
  # set suspicious amplifications to NA for COL2A1
  ind <- which(data$Target.Name == "COL2A1" & !data$Sample.Name %in% c("NTC_H2O","NTC", "H2O") & (data$CT > threshold_COL2A1 | data$Amp.Status=="No Amp"))
  if(!is_empty(ind)){
    data[ind,]$ct <- NA
  }
  
  # collect target and sample overviews 
  targets <- unique(data$Target.Name)
  targets2 <- targets[!targets %in% "COL2A1"]
  samples <- unique(data$Sample.Name)
  samplesNoctrl <- setdiff(samples,c("STD1", "STD2", "STD3","STD4","STD_1", "STD_2", "STD_3", "STD_4","Std 1", "Std 2", "Std 3", "Std 4","Std_1", "Std_2", "Std_3", "Std_4","Std. 1", "Std. 2", "Std. 3", "Std. 4","NTC_H2O","NTC", "H2O")) #but include gBlocks C("posCo","PosCo","gBlock", "gBLOCK", "gBlock (+)")
  
  # For each rep evaluate COL2A1 and add a column PASS/FAIL for filtering data
  df <- NULL
  for (s in samples) {
    subset1 <- data %>%
      filter(Sample.Name==s,
             rep==1,
             Target.Name=="COL2A1")
    subset1 <- subset1 %>%
      mutate(COL2A1_check =
               case_when( 
      is.na(ct) ~ "FAIL", TRUE ~"PASS")) %>%
      #is.na(ct) | ct > threshold_COL2A1 ~ "FAIL", TRUE ~ "PASS")) %>%
      dplyr::select(Sample.Name,rep,COL2A1_check)
    df <- rbind(df,subset1)
    
    subset2 <- data %>%
      filter(Sample.Name==s,
             rep==2,
             Target.Name=="COL2A1")
    subset2 <- subset2 %>%
      mutate(COL2A1_check = 
               case_when( 
      is.na(ct) ~ "FAIL", TRUE ~"PASS")) %>%
      #is.na(ct) | ct > threshold_COL2A1 ~ "FAIL", TRUE ~ "PASS")) %>%
      dplyr::select(Sample.Name,rep,COL2A1_check)
    df <- rbind(df,subset2)
  }
  
  data <- left_join(data, df, by=c("Sample.Name","rep"))
  
  # When both reps of a sample are undetermined or above the set threshold for COL2A1,
  # DNA input for this sample was too low, and no PMR could be calculated.
  # If one out of two reps is undetermined or above the set threshold for COL2A1, the sample should be reprocessed.
  # IF one out of two reps is undetermined or above the set threshold for the other targets, the respective amplification is disregarded (no ct value is taken forward) 
  low_input_fail <- NULL # negative control should always endup here
  reprocess_needed <- NULL # one of two COL2A1 reps failed
  
  for (s in samples){
    for(t in targets){
      tmp <- data %>%
        filter(Sample.Name == s & Target.Name == t) %>%
        dplyr::select(Sample.Name, rep, ct, COL2A1_check)
      
      if(any(is.na(tmp$ct)) && any(!is.na(tmp$ct))){
        # set non-replicated amplification to NA
        ind <- which((data$Sample.Name == s & data$Target.Name == t & !is.na(data$ct)))
        data[ind,]$ct <- NA
        
        if(t == "COL2A1"){
          # list sample as reprocessing required
          x <- data.frame(sample = tmp$Sample.Name[is.na(tmp$ct)],
                          rep = tmp$rep[is.na(tmp$ct)],
                          target = t,
                          `COL2A1_check (this replicate)` = tmp$COL2A1_check[is.na(tmp$ct)])
          reprocess_needed <- rbind(reprocess_needed,x)
        }}
      
      if(t == "COL2A1" & all(tmp$COL2A1_check=="FAIL")){
        # if target is COL2A1 and both reps fail, sample should be listed under low_input_fail
        x <- data.frame(sample = tmp$Sample.Name, COL2A1_Ct = tmp$ct)
        low_input_fail <- rbind(low_input_fail, x)
        }
    }
  }
  
  # Flag reactions with high ct (after removing spurious amplification)
  flagged <- data %>% 
    mutate(flag=case_when(ct >= threshold_targets_warning ~ paste0("CT below threshold ", threshold_targets,", but still higher than ",threshold_targets_warning), TRUE ~"")) %>%
    select(Well.Position,Sample.Name,Target.Name,rep,ct,COL2A1_check,flag) %>%
    filter(flag != "")
  
  # Carry on with ct means
  data1 <- data %>%
    pivot_wider(id_cols = Sample.Name,
                names_from = Target.Name,
                values_from = ct,
                values_fn = mean)
  
  data2 <- data %>%
    pivot_wider(id_cols = Sample.Name,
                names_from = Target.Name,
                values_from = ct,
                values_fn = sd)

  # Calculate PMR ----
  
  # Calculate input amounts (using calibration curve COL2A1)
  for (t in targets){
    varname <- paste0("input_",t)
    data1[varname] <- with(data1, 10^((data1[t]-intercept)/slope))
  }
  
  # Normalization: Divide calculated input amount by input amount calculated for COL2A1
  for (t in targets){
    varname <- paste0("ref_",t)
    calcname <- paste0("input_",t)
    data1[varname] <- data1[calcname]/data1["input_COL2A1"]
  }
  
  # Fully methylated DNA (SSS1 or gBlock) - only gBlock implemented
  gblock <- data1 %>%
    filter(Sample.Name %in% c("posCo","PosCo","gBlock", "gBLOCK", "gBlock (+)")) 
  
  # For each ref gene divide normalized input sample by normalized input fully methylated DNA *100
  results <- data1
  
  for (t in targets2){
    varname <- t
    calcname <- paste0("ref_",t)
    div_block <- c(gblock[calcname])
    results[varname] <- (data1[calcname]/div_block)*100
  }
  
  results <- as.data.frame(results)
  results <- results %>%
    select(Sample.Name, all_of(targets2))
  results[is.na(results)] <- 0 # set all NAs to 0
  
  # Calculate concentration input DNA based on COL2A mean Ct
  conc_input <- data1 %>% 
    as.data.frame() %>%
    mutate(input_conc = (data1$COL2A1 - intercept)/slope) %>%
    select(Sample.Name, input_conc) %>% droplevels
  
  # Calculate WID-qEC
  targets_qEC <- c("GYPC1","GYPC2","ZSCAN12")
  if(all(targets_qEC %in% targets)){
    results <- results %>%
    mutate(WIDqEC = round((GYPC1 + GYPC2 + ZSCAN12), digits=3)) # calculate sumPMR for WID-qEC
    results <- results %>% 
      mutate(WIDqEC_interpret = case_when(results$WIDqEC >= qEC_threshold1 ~ "Positiv",
                                     TRUE ~ "Negativ"))
  }
  
  # samples for which both reps COL2A1 failed, PMR should not be 0, but should be NA
  if(!is_empty(low_input_fail)){
    samples_failed <- low_input_fail %>% filter(!sample %in% c("NTC","NTC_H2O")) %>% unique() %>% pull(sample)
    if(!is_empty(samples_failed)){
      for (i in 2:(length(targets2)+3)){
        for (j in 1:length(samples_failed)){
          results[results$Sample.Name==samples_failed[j],i]<- NA
        }
      }
    }
  }
  
  # samples for which one rep of COL2A1 failed, PMR should not be 0, but should be NA
  if(!is_empty(reprocess_needed)){
    samples_failed <- reprocess_needed %>% pull(sample)
    for (i in 2:(length(targets2)+3)){
      for (j in 1:length(samples_failed)){
        results[results$Sample.Name==samples_failed[j],i]<- NA
      }
    }
  }
    
  # collect results mean+SD CT
  data1 <- data1 %>%
    select(Sample.Name, COL2A1, all_of(targets2)) %>%
    as.data.frame()
  
  data2 <- data2 %>%
    select(Sample.Name, COL2A1, all_of(targets2)) %>%
    as.data.frame()
  
  # Additional warning for When COL2A1 passed for both reps but CT stdev > 1.5 CT
  warning = data2 %>% 
    select(Sample.Name,COL2A1) %>%
    mutate(Warning.SD.COL2A1 = case_when(COL2A1 > 1.5 ~ "warning SD > 1.5 CT", TRUE ~ "PASS"))
  
  if(any(warning$Warning.SD.COL2A1== "warning SD > 1.5 CT")){
    warning = warning %>% filter(Warning.SD.COL2A1=="warning SD > 1.5 CT")
  } else{
    warning = NULL
  }
  
  #plot
  ### transform data frames to long format, combine dataframes, drop NAs, select samples only
  Ct_means <- pivot_longer(data1,!Sample.Name, names_to = "Target", values_to = "Ct_mean")
  Ct_sd <- pivot_longer(data2,!Sample.Name, names_to = "Target", values_to = "Ct_sd")
  df <- full_join(Ct_means, Ct_sd) %>%
    drop_na() %>%
    filter(Sample.Name %in% samplesNoctrl)
  
  ### plot output depends on number of 
  # plot processed data all targets and save
  plot <- df %>%
    filter(Target=="COL2A1") %>%
    ggplot(aes(x=Sample.Name, y=Ct_mean)) +
    geom_point()+
    geom_errorbar(aes(ymin=Ct_mean-Ct_sd, ymax=Ct_mean+Ct_sd), width=0.5) +
    #facet_wrap(~target) +
    xlab("") + #rename axis text label
    ylab("Ct (mean)") +
    theme_minimal() +
    theme(axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12),
          axis.text.x=element_text(size=5,angle=60,hjust=1),
          axis.text.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title = element_blank(),
          legend.key.size= unit(6, "mm"))
  
  # make the final summary
  final <- results %>%
    filter(!Sample.Name %in% c("STD_1","STD_2","STD_3","STD_4","posCo","PosCo","NTC_H2O","NTC")) %>%
    select(Sample.Name, WIDqEC, WIDqEC_interpret) %>%
    dplyr::rename(SampleName = Sample.Name)
  final$WIDqEC <- format(final$WIDqEC, nsmall=3, scientific=FALSE)
  
  samples <- samples[!samples %in% c("STD_1","STD_2","STD_3","STD_4","posCo","PosCo","NTC_H2O","NTC")]
  
  QC <- data.frame(SampleName = samples)
  
  QC <- QC %>%
    mutate(QC=case_when(
     SampleName %in% low_input_fail$sample ~ "Insufficient DNA in both reps", #COL2A1 did not amplify in any of the reps 
     SampleName %in% reprocess_needed[,1] ~ "Reprocessing needed insufficient DNA in one rep", # samples for which for only one of two reps COL2A1 failed
     SampleName %in% warning[,1] ~ "STDEV CT COL2A1 exceeds 1.5 cycles",
     TRUE ~ "PASS" #all ok
    ))
  
  final <- full_join(final,QC)
  
  # Save results to list ---
  list_out[[2]] <- results #PMR + WIDqEC + WIDqEC outcome
  list_out[[3]] <- data1 #mean CT
  list_out[[4]] <- data2 #stdev CT
  list_out[[5]] <- conc_input #log(copy number/5uL)
  if(is_empty(low_input_fail)){
    list_out[[6]] <- as.data.frame("you have some amplification in the NTC_H2O")
  } else{
    list_out[[6]] <- low_input_fail #samples for which COL2A1 failed in both reps, should have negative controls
  }
  
  if(!is_empty(reprocess_needed)){ # samples for which one of two reps COL2A1 failed
    list_out[[7]] <- reprocess_needed
  } else {
    list_out[[7]] <- as.data.frame("none")
  }
  
  if(!is_empty(flagged)){ #reactions with high CT, current threshold set at 30
    list_out[[8]] = flagged
  } else {
    list_out[[8]] = as.data.frame("none")
  }
  
  if(!is_empty(warning)){ #marks samples for which STDEV COL2A > 1.5 CT (QC warning)
    list_out[[9]] = warning
  } else {
    list_out[[9]] = as.data.frame("none")
  }
  
  list_out[[10]] <- plot #mean CT COL2A1
  
  list_out[[11]] <- final #information in final csv file
  
  return(list_out)
}
