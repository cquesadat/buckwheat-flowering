# Author: Flavian Tschurr
# Flowering-phase logistic and abscission-phase linear fitting.
# Original contributions: Flavian Tschurr (2025).
################################################################################
# 1. define maximum
# 1.1 smooth (reduce the effect of outlier to define max) first create 
#         daily values, then smooth
# 1.2 look for max value per plot

# 2. cut period
# 3. fit logistic regression
# 4. get parameters
#     4.1 coefficients
# 4.2 get derivatives (1. and 2.)
# define max speed, and start/stop
# 5. fit a linear model for the abscission phase



# define function
get_parameters_start_to_max <- function(one_genotype, date_range){
  require(ggplot2)
  require(cowplot)
  require(lubridate)
  # Supporting interpolation and derivative functions.
  source(file.path("utils", "FUN_utils.R"))
  
  # output list
  one_genotype_parameters <- list()
  model_outs <- list()
  metrics_list <- list()
  for (plot.ID in unique(one_genotype$Plot)) {
    print(plot.ID)
    one_plot <- subset(one_genotype, Plot == plot.ID)
   
    # step 1: define maximum
   
  # gapfill and smooth
  filled_data <- gapfill_and_smooth_function(one_plot$Total_Area,one_plot$Date)
    filled_data$time_in <- c(1:length(filled_data$value))
    max_value <- which.max(filled_data$value_smoothed)
    
    # step 2: cut period
    
  cut_data <- filled_data[c(1:max_value),]
  data_vect <- cut_data$value_filled
    # step 3: fit logistic model
    
    control_nls <- nls.control(maxiter = 500)
    logistic_model <- try(nls(data_vect ~ SSlogis(time_in,Asym,xmid,scal),
                              cut_data, control= control_nls,trace = TRUE),
                          silent=TRUE)
    
    if(typeof(logistic_model)== "list"){
      fitted_model_values <- fitted(logistic_model)
      cut_data$fitted_model_values <- fitted_model_values
    }else{
      print(plot.ID)
      # we provide more data to try to fit the model (half of the remaining)
      max_value <- (length(filled_data$Date)-max_value)/2
      
      cut_data <- filled_data[c(1:max_value),]
      data_vect <- cut_data$value_filled
      
      logistic_model <- try(nls(data_vect ~ SSlogis(time_in,Asym,xmid,scal),
                                cut_data, control= control_nls,trace = TRUE),
                            silent=TRUE)
      if(typeof(logistic_model)== "list"){
        fitted_model_values <- fitted(logistic_model)
        cut_data$fitted_model_values <- fitted_model_values
      }else{
        message("Logistic model did not converge for plot ", plot.ID)
        next
      }
      
    }
    
    
    
    # step 4: get parameters
    
    # logistic function as expression
    expr <- expression( asym/(1+exp((xmid-x)/ scal)))
    # get coefficients of the logistic model
    coeff <- coefficients(logistic_model)
    
    # put the values into the expression
    asymptote_fun <- do.call(substitute, list(as.list(expr)[[1]], env= 
                                                list(asym=as.numeric(coeff["Asym"]),
                                                     xmid=as.numeric(coeff["xmid"]),
                                                     scal=as.numeric(coeff["scal"]))))
    
    # get second derivative of the asymptote (--> speed: enables us the find the start point)
    asymptote2deriv <- DD(asymptote_fun,"x",2)
    # get data for the 2nd derrivative
    data_2deriv <- eval(asymptote2deriv,list(x=c(1:length(data_vect))))
    # maximal gradient --> max of first derrivative --> mid of the growth phase
    speed_divider_end = 4
    speed_divider_begin = 4
    
    if(coeff[which(names(coeff)=="scal")]>=0){
      max_gradient <- which.max(eval(D(asymptote_fun,"x"),list(x=c(1:length(data_vect)))))
      # get the end of the phase of the slope
      # closest local maximum higher than the max gradient in the data vector
      log_min_speed <- max(which(data_2deriv[max_gradient:length(data_2deriv)]<= 
                                   min(data_2deriv, na.rm = T)/speed_divider_end), 
                           na.rm=TRUE) + max_gradient -1
      end_phase <- log_min_speed
      # get start point from the logistic model second derivative
      # take the first entry which has a higher value than maximal of the second
      # derivative divided by 6 before the max gradient (so before the end, 
      #  where we expect the start point)
      start_phase <- min(which(data_2deriv[1:max_gradient]>= max(data_2deriv,
                                                                 na.rm = T)/speed_divider_begin), na.rm=TRUE)
      
    }else{
      max_gradient <- which.min(eval(D(asymptote_fun,"x"),list(x=c(1:length(data_vect)))))
      # get the end of the phase of the slope
      # closest local maximum higher than the max gradient in the data vector
      log_min_speed <- min(which(data_2deriv[1:max_gradient]<=
                                   min(data_2deriv, na.rm = TRUE)/speed_divider_begin),
                           na.rm = TRUE) #+ max_gradient -1
      
      end_phase <- log_min_speed
      
      # get start point from the logistic model second derivative
      # take the first entry which has a higher value than maximal of the second
      # derivative divided by 6 before the max gradient (so before the end,
      #  where we expect the start point)
      start_phase <- max(which(data_2deriv[max_gradient:length(data_2deriv)]>= 
                                 max(data_2deriv, na.rm = T)/speed_divider_end),
                         na.rm = T) + max_gradient -1
      
    }

    # Flowering onset at the 2% coverage threshold, with derivative fallback.
    threshold_value <- 2
    threshold_index <- NA
    if(any(!is.na(cut_data$value_filled))){
      idxs <- which(cut_data$value_filled >= threshold_value)
      if(length(idxs)>0) threshold_index <- idxs[1]
    }

    # choose start phase: prefer threshold-based if available, otherwise fall back to derivative-based
    if(!is.na(threshold_index)){
      start_phase_final <- threshold_index
    } else if(exists("start_phase") && !is.na(start_phase)){
      start_phase_final <- start_phase
    } else {
      start_phase_final <- NA
    }
    
    
    ## get area under the curve
    #   asym / (1 + exp((xmid - x) / scal))
    # }
    
    logistic_function <- function(x) {
      coeff["Asym"] / (1 + exp((coeff["xmid"] - x) / coeff["scal"]))
    }
    
    #calculate area under the curve
    area_under_curve<-integrate(logistic_function, lower = 0, upper = length(fitted_model_values))
    area_under_curve<-area_under_curve$value
    
    # step 5.1: end phase - cut correct period (from max value to end)
    
    ## end of flowering phase
    cut_after_peak <- filled_data[c(max_value:length(filled_data$value_smoothed)),]
    print(plot.ID)
    
    if(length(cut_after_peak$value)<=1){ # if too short (no lm can be fit, we set all as NA)
      end_phase_intercept_amplitude <-NA # intercept of the model (y value)
      end_phase_slope <- NA # slope of the model
      end_phase_x0 <-NA
      # linear model metrics
      linear_R2 <- NA
      linear_RMSE <- NA
      linear_AIC <- NA
      linear_BIC <- NA
    }else{
      
      # step 5.2: fit lm
      lm_out <- lm(value ~ time_in,cut_after_peak)
      cut_after_peak$fitted_model_values <- predict(lm_out,cut_after_peak)
      end_phase_intercept_amplitude <- coef(lm_out)[1] # intercept of the model (y value)
      end_phase_slope <- coef(lm_out)[2] # slope of the model
      end_phase_x0 <- -coef(lm_out)[1] / coef(lm_out)[2] # where model cuts x axis 
      # linear model metrics
      lm_resid <- residuals(lm_out)
      linear_SSE <- sum(lm_resid^2, na.rm = TRUE)
      linear_SST <- sum((cut_after_peak$value - mean(cut_after_peak$value, na.rm = TRUE))^2, na.rm = TRUE)
      linear_R2 <- ifelse(linear_SST>0, 1 - linear_SSE/linear_SST, NA)
      linear_RMSE <- sqrt(linear_SSE/length(na.omit(lm_resid)))
      linear_AIC <- tryCatch(AIC(lm_out), error=function(e) NA)
      linear_BIC <- tryCatch(BIC(lm_out), error=function(e) NA)
      
      
    }
    
    if(is.na(end_phase_slope)){
      end_phase_date_zero <- NA
      duration_start_flowering_to_end_abscission = NA
      duration_mid_flowering_to_end_abscission = NA
      duration_end_flowering_to_end_abscission = NA
    }else if(end_phase_slope<0){
      end_phase_date_zero <- as.Date(filled_data$Date[1]+ lubridate::days((
        round(end_phase_x0)-1)))
      
      duration_start_flowering_to_end_abscission = if(!is.na(start_phase_final) && !is.na(cut_data$Date[start_phase_final])) abs(as.numeric(
        cut_data$Date[start_phase_final]-end_phase_date_zero)) else NA
      duration_mid_flowering_to_end_abscission = if(!is.na(max_gradient)) abs(as.numeric(
        cut_data$Date[max_gradient]-end_phase_date_zero)) else NA
      duration_end_flowering_to_end_abscission = if(!is.na(end_phase)) abs(as.numeric(
        cut_data$Date[end_phase]-end_phase_date_zero)) else NA
    } else{
      end_phase_date_zero <- NA
      duration_start_flowering_to_end_abscission = NA
      duration_mid_flowering_to_end_abscission = NA
      duration_end_flowering_to_end_abscission = NA
      
    }
    
    ##
    
    output_one_plot <- data.frame(
      "Plot"= one_plot$Plot[1],
      "ID"= one_plot$ID[1],
      "Asym"=coeff[1],
      "xmid"=coeff[2],
      "scal"=coeff[3],
      "area_under_curve" = area_under_curve,
      "start_flowering_date"= if(!is.na(start_phase_final)) cut_data$Date[start_phase_final] else NA,
      "end_flowering_date"= if(!is.na(end_phase)) cut_data$Date[end_phase] else NA,
      "start_flowering_value" = if(!is.na(start_phase_final)) cut_data$value_filled[start_phase_final] else NA,
      "end_flowering_value"= if(!is.na(end_phase)) cut_data$value_filled[end_phase] else NA,
      "max_flowering_rate_date"= if(!is.na(max_gradient)) cut_data$Date[max_gradient] else NA,
      "max_flowering_rate_value" = if(!is.na(max_gradient)) cut_data$value_filled[max_gradient] else NA,
      "duration_flowering_start_end"= if(!is.na(start_phase_final) && !is.na(end_phase)) abs(as.numeric(cut_data$Date[start_phase_final]-cut_data$Date[end_phase])) else NA,
      "duration_flowering_start_mid"=  if(!is.na(start_phase_final) && !is.na(max_gradient)) abs(as.numeric(cut_data$Date[start_phase_final]-cut_data$Date[max_gradient])) else NA,
      "duration_flowering_mid_end"=  if(!is.na(max_gradient) && !is.na(end_phase)) abs(as.numeric(cut_data$Date[max_gradient]-cut_data$Date[end_phase])) else NA,
      "amplitude_flowering_start_end" = if(!is.na(start_phase_final) && !is.na(end_phase)) ( cut_data$value_filled[end_phase]- cut_data$value_filled[start_phase_final]) else NA,
      "amplitude_flowering_start_mid" = if(!is.na(start_phase_final) && !is.na(max_gradient)) ( cut_data$value_filled[max_gradient]- cut_data$value_filled[start_phase_final]) else NA,
      "amplitude_flowering_mid_end" = if(!is.na(max_gradient) && !is.na(end_phase)) ( cut_data$value_filled[end_phase]- cut_data$value_filled[max_gradient]) else NA,
      "abscission_intercept_amplitude"= end_phase_intercept_amplitude,
      "abscission_slope" = end_phase_slope,
      "abscission_x0" = end_phase_x0,
      "abscission_projected_end_date" = end_phase_date_zero,
      "duration_start_flowering_to_end_abscission" = duration_start_flowering_to_end_abscission,
      "duration_mid_flowering_to_end_abscission" = duration_mid_flowering_to_end_abscission,
      "duration_end_flowering_to_end_abscission" = duration_end_flowering_to_end_abscission
    )
    
    # --- Goodness-of-fit metrics for logistic model ---
    # compute residual-based metrics
    logistic_R2 <- NA; logistic_RMSE <- NA; logistic_AIC <- NA; logistic_BIC <- NA
    if(typeof(logistic_model) == "list"){
      log_resid <- residuals(logistic_model)
      SSE <- sum(log_resid^2, na.rm = TRUE)
      SST <- sum((data_vect - mean(data_vect, na.rm = TRUE))^2, na.rm = TRUE)
      logistic_R2 <- ifelse(SST>0, 1 - SSE/SST, NA)
      logistic_RMSE <- sqrt(SSE/length(na.omit(log_resid)))
      logistic_AIC <- tryCatch(AIC(logistic_model), error=function(e) NA)
      logistic_BIC <- tryCatch(BIC(logistic_model), error=function(e) NA)
    }

    # metrics row
    metrics_row <- data.frame(
      Plot = one_plot$Plot[1],
      ID = one_plot$ID[1],
      logistic_R2 = logistic_R2,
      logistic_RMSE = logistic_RMSE,
      logistic_AIC = logistic_AIC,
      logistic_BIC = logistic_BIC,
      linear_R2 = if(exists("linear_R2")) linear_R2 else NA,
      linear_RMSE = if(exists("linear_RMSE")) linear_RMSE else NA,
      linear_AIC = if(exists("linear_AIC")) linear_AIC else NA,
      linear_BIC = if(exists("linear_BIC")) linear_BIC else NA,
      threshold_used = ifelse(!is.na(threshold_index), TRUE, FALSE),
      threshold_value = threshold_value
    )
    metrics_list[[plot.ID]] <- metrics_row
    
    one_genotype_parameters[[plot.ID]] <- output_one_plot
    filled_data$Plot <- plot.ID
    filled_data$fitted_model_values <- c(cut_data$fitted_model_values,
                                         cut_after_peak$fitted_model_values[-1])[1:length(filled_data$Plot)]
    
    model_outs[[plot.ID]] <- filled_data
    
  }
  
  
  # prepare output
  output_list <- list()
  
  # create overview figure
  model_out_df <- do.call("rbind",model_outs)
  one_genotype_parameters_df <- do.call("rbind",one_genotype_parameters)
  
  one_genotype_parameters_df$start_flowering_date_jittered <- 
    one_genotype_parameters_df$start_flowering_date + 
    runif(nrow(one_genotype_parameters_df), min=-0.3, max=0.2)
  one_genotype_parameters_df$end_flowering_date_jittered <- 
    one_genotype_parameters_df$end_flowering_date + 
    runif(nrow(one_genotype_parameters_df), min=-0.3, max=0.2)
  
  one_genotype_parameters_df$end_abscission_date_jittered <- 
    one_genotype_parameters_df$abscission_projected_end_date + 
    runif(nrow(one_genotype_parameters_df), min=-0.3, max=0.2)
  
  output_list[["figure"]]<- ggplot(model_out_df, aes(x=Date, color=Plot ))+
    geom_line(aes(y=fitted_model_values), size = 1.5, show.legend = F)+
    geom_vline(one_genotype_parameters_df, mapping= aes(xintercept = 
                                                          start_flowering_date_jittered, color=Plot),
               size = 1.2, show.legend = F)+
    geom_vline(one_genotype_parameters_df, mapping= aes(xintercept= 
                                                          end_flowering_date_jittered, color=Plot), linetype=2 ,
               size = 1.2, show.legend = F)+
    geom_vline(one_genotype_parameters_df, mapping= aes(xintercept= 
                                                          end_abscission_date_jittered, color=Plot), linetype=3 ,
               size = 1.2, show.legend = F)+
    geom_point(aes(y=value, fill = Plot), shape = 21, color = "black", size = 3, show.legend = F)+
    ylab("% Flowering")+ xlab("Date")+
    ggtitle(paste0("Genotype: ",one_genotype$ID[1]))+
    theme_cowplot() +
    scale_x_date(limits = as.Date(date_range)) +
    theme(legend.position = "none") 
  
  
  output_list[["parameters"]] <- do.call("rbind",one_genotype_parameters)
  # attach metrics if present
  if(length(metrics_list)>0){
    output_list[["metrics"]] <- do.call("rbind", metrics_list)
  } else {
    output_list[["metrics"]] <- NULL
  }
  
  return(output_list)
  
}

