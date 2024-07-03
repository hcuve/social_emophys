# Load necessary libraries
library(dplyr)
library(ggplot2)
library(zoo)  # For na.approx function

# Helper function to interpolate missing data
interpolate <- function(data, method = "linear") {
  na.approx(data, rule = 2)
}

# Helper function to calculate speed
calculate_speed <- function(x, y, samplerate) {
  dx <- diff(x, lag = 1)
  dy <- diff(y, lag = 1)
  speed <- sqrt(dx^2 + dy^2) * samplerate
  c(speed, NA) # Append NA to match length of original vectors
}

# Fixation and Saccade Classification Function
fixationANDsaccade <- function(speed, thres_vel, thres_dur, Hz){
  thres_dur <- thres_dur * (Hz / 1000)
  fixsac <- ifelse(speed > thres_vel, 's', 'f')
  rle <- rle(fixsac)
  fixsac[cumsum(rle$lengths)[which(rle$lengths < 10 * (Hz / 1000) & rle$values == 's')]] <- 'f'
  
  while(length(which(rle$lengths < 10 * (Hz / 1000) & rle$values == 's')) != 0){
    rle <- rle(fixsac)
    fixsac[cumsum(rle$lengths)[which(rle$lengths < 10 * (Hz / 1000) & rle$values == 's')]] <- 'f'
  }
  classify <- numeric()
  for(i in 1:length(rle$values)){
    if(is.na(rle$values[i])){
      classify <- c(classify, rep(NA, rle$lengths[i]))
    } else{
      if(rle$values[i] == 'f' & rle$lengths[i] >= thres_dur){
        classify <- c(classify, rep('f', rle$lengths[i]))
      }
      if(rle$values[i] == 'f' & rle$lengths[i] < thres_dur){
        classify <- c(classify, rep('u', rle$lengths[i]))
      }
      if(rle$values[i] == 's'){
        classify <- c(classify, rep('s', rle$lengths[i]))
      }
    }
  }
  return(classify)
}

# Main gazepath function
gazepath <- function(data, x1, y1, x2 = NULL, y2 = NULL, d1, d2 = NULL, trial, height_px, height_mm, width_px, width_mm, extra_var = NULL, res_x = 1280, res_y = 1024, samplerate = 500, method = 'Mould', posthoc = FALSE, thres_vel = 35, thres_dur = 100, min_dist = 250, in_thres = 150){
  if(!is.data.frame(data)) {
    stop('please insert a data frame and define the column numbers of the variables')
  }
  
  if(is.numeric(trial)){
    trial <- colnames(data)[trial]
  }
  
  data[[trial]] <- as.character(data[[trial]])
  
  unique_trials <- unique(data[[trial]])
  rle_lengths <- rle(as.numeric(factor(data[[trial]])))$lengths
  
  if(length(unique_trials) != length(rle_lengths)){
    TRIAL_NEW <- rep.int(1:length(rle_lengths), rle_lengths)
    data <- data.frame(data, TRIAL_NEW)
    if(is.numeric(trial)){
      names(data)[trial] <- 'TRIAL_OLD'
      trial <- 'trial'
    } else {
      names(data)[which(names(data) == trial)] <- 'TRIAL_OLD'
    }
    names(data)[ncol(data)] <- trial
    warning('The trial index in the data frame was not unique, therefore trials are renamed to be unique and the old trial index is stored in the data as TRIAL_OLD')
  }
  
  extra <- list()
  
  if(!is.null(extra_var)){
    if(sum(extra_var %in% names(data)) == length(extra_var)){
      for(i in unique(data[[trial]])){
        extra[[i]] <- sapply(extra_var, function(var) {
          if(is.factor(data[[var]])) {
            return(as.character(data[data[[trial]] == i, which(names(data) == var)][1]))
          } else {
            return(data[data[[trial]] == i, which(names(data) == var)][1])
          }
        })
      }
    } else {
      extra_var <- NULL
      print('Please make sure the variables to pass through have the correct names')
    }
  }
  
  data[[d1]] <- ifelse(data[[d1]] < min_dist, NA, data[[d1]])
  if(is.null(d2)){
    D <- by(data[[d1]], data[[trial]], data.frame)
  } else {
    data[[d2]] <- ifelse(data[[d2]] < min_dist, NA, data[[d2]])
    D <- by((data[[d1]] + data[[d2]]) / 2, data[[trial]], data.frame)
  }
  
  if(!is.null(x2) & !is.null(y2)){
    data[[x1]] <- ifelse(is.na(data[[x1]]), data[[x2]], data[[x1]])
    data[[y1]] <- ifelse(is.na(data[[y1]]), data[[y2]], data[[y1]])
  }
  
  # Interpolate missing data
  data[[x1]] <- interpolate(data[[x1]])
  data[[y1]] <- interpolate(data[[y1]])
  
  # Calculate speed
  data$speed <- calculate_speed(data[[x1]], data[[y1]], samplerate)
  
  # Identify fixations and saccades
  data$class <- fixationANDsaccade(data$speed, thres_vel, thres_dur, samplerate)
  
  # Return the processed data
  return(list(data = data, extra = extra))
}

# Example usage:
# data <- read.csv('path_to_your_data.csv')
# result <- gazepath(data, 
#                    x1 = "x_column_name", 
#                    y1 = "y_column_name", 
#                    d1 = "distance_column_name", 
#                    trial = "trial_column_name", 
#                    height_px = 1080, 
#                    height_mm = 295, 
#                    width_px = 1920, 
#                    width_mm = 525, 
#                    samplerate = 150, 
#                    method = 'gazepath')

# Access the processed data
# processed_data <- result$data
# extra_info <- result$extra