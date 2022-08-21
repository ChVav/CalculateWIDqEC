# PMR calculation main work function
# contact: charlottevavourakis@gmail.com
# under development for Shiny App version for TirolPath
# Last update 21 August 2022

calculate_pmr <- function(data,#experimentname,
                          threshold_COL2A1=35,
                          external_curve=FALSE,
                          curve=NULL,
                          fix_intercept=36.9,
                          fix_slope=-3.4,
                          calib_fixed=FALSE){
  
  # WID-qEC thresholds
  qEC_threshold1 <- 0.03
  qEC_threshold2 <- 0.63
  
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
  
  # Carry on with CT means
  data <- data %>%
    mutate(rep = case_when(grepl("A|C|E|G|I|K|M|O", Well.Position) == TRUE ~ 1,
                           grepl("B|D|F|H|J|L|N|P", Well.Position) == TRUE ~ 2)) %>%
    mutate(ct = as.numeric(ifelse(CT == "Undetermined", NA, CT)))
  
  # collect target and sample overviews
  targets <- unique(data$Target.Name)
  targets2 <- targets[!targets %in% "COL2A1"]
  samples <- unique(data$Sample.Name)
  
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
      is.na(ct) | ct > threshold_COL2A1 ~ "FAIL", TRUE ~ "PASS")) %>%
      dplyr::select(Sample.Name,rep,COL2A1_check)
    df <- rbind(df,subset1)
    
    subset2 <- data %>%
      filter(Sample.Name==s,
             rep==2,
             Target.Name=="COL2A1")
    subset2 <- subset2 %>%
      mutate(COL2A1_check = 
               case_when( 
      is.na(ct) | ct > threshold_COL2A1 ~ "FAIL", TRUE ~ "PASS")) %>%
      dplyr::select(Sample.Name,rep,COL2A1_check)
    df <- rbind(df,subset2)
  }
  
  data <- left_join(data, df, by=c("Sample.Name","rep"))
  
  # When both reps of a sample are undetermined or above the set threshold for COL2A1,
  # DNA input for this sample was too low, and no PMR could be calculated.
  # If one out of two reps is undetermined or above the set threshold,
  # sample should to be reprocessed (COL2A1) or reprocessing is recommended. 
  # In these cases, one Ct value is taken forward. 
  low_input_fail <- NULL # negative control should always endup here
  reprocess_needed <- NULL
  reprocess <- NULL
  
  for (s in samples){
    for(t in targets){
      tmp <- data %>%
        filter(Sample.Name == s & Target.Name == t) %>%
        dplyr::select(Sample.Name, rep, ct, COL2A1_check)
      
      if(t == "COL2A1" & all(tmp$COL2A1_check=="FAIL")){
        # if target is COL2A1 and both reps fail, sample should be listed under low_input_fail
        x <- data.frame(sample = tmp$Sample.Name, COL2A1_Ct = tmp$ct)
        low_input_fail <- rbind(low_input_fail, x)
      }
      
      if(any(is.na(tmp$ct)) && any(!is.na(tmp$ct))){
        # Carry one Ct value forward
        und <- sum(is.na(tmp$ct))
        ind1 <- na.omit(data$Sample.Name == s & data$Target.Name == t)
        ct <- ifelse(und < 2, na.omit(tmp$ct), NA)
        data[ind1,]$ct <- ct
        
        if(t == "COL2A1"){
          # reprocessing required
          x <- data.frame(sample = tmp$Sample.Name[is.na(tmp$ct)],
                          rep = tmp$rep[is.na(tmp$ct)],
                          target = t,
                          `COL2A1_check (this replicate)` = tmp$COL2A1_check[is.na(tmp$ct)])
          reprocess_needed <- rbind(reprocess_needed,x)
        } else {
          # reprocessing recommended
          x <- data.frame(sample = tmp$Sample.Name[is.na(tmp$ct)],
                          rep = tmp$rep[is.na(tmp$ct)],
                          target = t,
                          `COL2A1_check (this replicate)` = tmp$COL2A1_check[is.na(tmp$ct)])
          reprocess <- rbind(reprocess,x)
        }
      }
    }
  }
  
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
    filter(Sample.Name %in% c("PosCo","gBlock", "gBLOCK", "gBlock (+)")) 
  
  # For each ref gene divide normalized input sample by normalized input fully methylated DNA *100
  results <- data1
  
  for (t in targets2){
    varname <- t
    calcname <- paste0("ref_",t)
    div_block <- c(gblock[calcname])
    results[varname] <- (data1[calcname]/div_block)*100
  }
  
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
    mutate(WIDqEC = (GYPC1 + GYPC2 + ZSCAN12)) # calculate sumPMR for WID-qEC
    results <- results %>% 
      mutate(WIDqEC_test = case_when(results$WIDqEC < qEC_threshold1 ~ "low risk EC/CIN",
                                     results$WIDqEC >= qEC_threshold1 & results$WIDqEC <= qEC_threshold2 ~ "high risk EC/CIN",
                                     results$WIDqEC > qEC_threshold2 ~ "very high risk EC/CIN"))
  }

  # Save results to list ---
  results <- as.data.frame(results)
  
  data1 <- data1 %>%
    select(Sample.Name, COL2A1, all_of(targets2)) %>%
    as.data.frame()
  
  data2 <- data2 %>%
    select(Sample.Name, COL2A1, all_of(targets2)) %>%
    as.data.frame()
  
  list_out[[2]] <- results #PMR + WIDqEC + WIDqEC outcome
  list_out[[3]] <- data1 #mean CT
  list_out[[4]] <- data2 #stdev CT
  list_out[[5]] <- conc_input #log(copy number/5uL)
  list_out[[6]] <- low_input_fail #samples for which COL2A1 failed in both reps, should have negative controls

  if(!is_empty(reprocess_needed)){ # samples for which for only one of two reps COL2A1 failed
    list_out[[7]] <- reprocess_needed
  } else {
    list_out[[7]] <- as.data.frame("none")
  }
  
  if(!is_empty(reprocess)){ # samples for which for only one of two reps target amplified
    list_out[[8]] <- reprocess
  } else {
    list_out[[8]] <- as.data.frame("none")
  }
  
  return(list_out)
}
