
# this kinda works but there are steps missing
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
  c(speed, NA)  # Append NA to match length of original vectors
}

# Fixation and Saccade Classification Function
fixationANDsaccade <- function(speed, thres_vel, thres_dur, Hz) {
  thres_dur <- thres_dur * (Hz / 1000)
  fixsac <- ifelse(speed > thres_vel, 's', 'f')
  rle <- rle(fixsac)
  fixsac[cumsum(rle$lengths)[which(rle$lengths < 10 * (Hz / 1000) & rle$values == 's')]] <- 'f'
  
  while (length(which(rle$lengths < 10 * (Hz / 1000) & rle$values == 's')) != 0) {
    rle <- rle(fixsac)
    fixsac[cumsum(rle$lengths)[which(rle$lengths < 10 * (Hz / 1000) & rle$values == 's')]] <- 'f'
  }
  classify <- numeric()
  for (i in 1:length(rle$values)) {
    if (is.na(rle$values[i])) {
      classify <- c(classify, rep(NA, rle$lengths[i]))
    } else {
      if (rle$values[i] == 'f' & rle$lengths[i] >= thres_dur) {
        classify <- c(classify, rep('f', rle$lengths[i]))
      }
      if (rle$values[i] == 'f' & rle$lengths[i] < thres_dur) {
        classify <- c(classify, rep('u', rle$lengths[i]))
      }
      if (rle$values[i] == 's') {
        classify <- c(classify, rep('s', rle$lengths[i]))
      }
    }
  }
  return(classify)
}

# Simplify Function
simplify <- function(classification, x, y, Hz, D, width_px, width_mm, extra, extra_var){
  class <- rle(classification)
  simple <- data.frame(class$values, class$lengths, 
                       c(1, cumsum(class$lengths) + 1)[-(length(class$values) + 1)], 
                       cumsum(class$lengths))
  if(length(which(class$values == 'f')) > 0){
    x_start <- y_start <- x_end <- y_end <- mean_x <- mean_y <- POGvar <- RMS <- numeric()
    for(i in 1:nrow(simple)){
      x_start <- c(x_start, x[simple[i, 3]])
      y_start <- c(y_start, y[simple[i, 3]])
      x_end <- c(x_end, x[simple[i, 4]])
      y_end <- c(y_end, y[simple[i, 4]])
      mean_x <- c(mean_x, mean(x[simple[i,3] : simple[i,4]]))
      mean_y <- c(mean_y, mean(y[simple[i,3] : simple[i,4]]))
      m <- as.matrix(dist(cbind(c(mean_x[length(mean_x)], x[simple[i,3] : simple[i,4]]), c(mean_y[length(mean_y)], y[simple[i,3] : simple[i,4]]))))
      RMS <- c(RMS, sqrt(mean((atan((diag(m[-1,-c(1,2)]) / 2) / mean(D, na.rm = TRUE)) * (180 / pi) * (width_mm / width_px) * 2)^2)))
      POGvar <- c(POGvar, mean(m[-1,1]))
    }
    ## Calculate saccade amplitude and transform POGvar from pixels to degrees of visual angle and to sd
    ss <- which(class$values == 's')
    POGvar[ss] <- sqrt((x_start[ss] - x_end[ss]) ^ 2 + (y_start[ss] - y_end[ss]) ^ 2)
    POGsdSacAmp <- atan((POGvar / 2) / mean(D, na.rm = TRUE)) * (180 / pi) * (width_mm / width_px) * 2
    POGsdSacAmp[!ss] <- sqrt(POGsdSacAmp)
    RMS[ss] <- NA
    
    simple <- data.frame(class$values, class$lengths * (1000/Hz), 
                         c(1, cumsum(class$lengths * (1000/Hz)) + 1)[-(length(class$values) + 1)], 
                         cumsum(class$lengths * (1000/Hz)), x_start, y_start, x_end, y_end, mean_x, mean_y, POGsdSacAmp, RMS)
    names(simple)[1:4] <- c('Value', 'Dur', 'Start', 'End')
    if(!is.null(extra_var)){
      for(i in seq_along(extra_var)){
        simple <- data.frame(simple, extra[i])
        names(simple)[ncol(simple)] <- extra_var[i]
      }
    }
  }
  # Remove NA values
  if(length(which(is.na(simple[,1]))) != 0){
    simple <- simple[-which(is.na(simple[,1])),]
  }
  return(simple)
}

# Main gazepath function
gazepath <- function(data, x1, y1, x2 = NULL, y2 = NULL, d1, d2 = NULL, trial, height_px, height_mm, width_px, width_mm, extra_var = NULL, res_x = 1280, res_y = 1024, samplerate = 500, method = 'Mould', posthoc = FALSE, thres_vel = 35, thres_dur = 100, min_dist = 250, in_thres = 150) {
  if (!is.data.frame(data)) {
    stop('please insert a data frame and define the column numbers of the variables')
  }
  
  if (is.numeric(trial)) {
    trial <- colnames(data)[trial]
  }
  
  data[[trial]] <- as.character(data[[trial]])
  
  unique_trials <- unique(data[[trial]])
  rle_lengths <- rle(as.numeric(factor(data[[trial]])))$lengths
  
  if (length(unique_trials) != length(rle_lengths)) {
    TRIAL_NEW <- rep.int(1:length(rle_lengths), rle_lengths)
    data <- data.frame(data, TRIAL_NEW)
    if (is.numeric(trial)) {
      names(data)[trial] <- 'TRIAL_OLD'
      trial <- 'trial'
    } else {
      names(data)[which(names(data) == trial)] <- 'TRIAL_OLD'
    }
    names(data)[ncol(data)] <- trial
    warning('The trial index in the data frame was not unique, therefore trials are renamed to be unique and the old trial index is stored in the data as TRIAL_OLD')
  }
  
  extra <- list()
  
  if (!is.null(extra_var)) {
    if (sum(extra_var %in% names(data)) == length(extra_var)) {
      for (i in unique(data[[trial]])) {
        extra[[i]] <- sapply(extra_var, function(var) {
          if (is.factor(data[[var]])) {
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
  if (is.null(d2)) {
    D <- by(data[[d1]], data[[trial]], data.frame)
  } else {
    data[[d2]] <- ifelse(data[[d2]] < min_dist, NA, data[[d2]])
    D <- by((data[[d1]] + data[[d2]]) / 2, data[[trial]], data.frame)
  }
  
  if (!is.null(x2) & !is.null(y2)) {
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
  
  # Simplify the results for each trial
  simplified_results <- by(data, data[[trial]], function(trial_data) {
    simplify(
      trial_data$class, 
      trial_data[[x1]], 
      trial_data[[y1]], 
      samplerate, 
      trial_data[[d1]], 
      width_px, 
      width_mm, 
      extra[[unique(trial_data[[trial]])]], 
      extra_var
    )
  })
  
  # Create list to store results
  results <- list(
    Classifications = by(data$class, data[[trial]], identity),
    x_coordinates = by(data[[x1]], data[[trial]], identity),
    y_coordinates = by(data[[y1]], data[[trial]], identity),
    Method = method,
    Robustness = by(data$speed, data[[trial]], function(s) mean(!is.na(s) * 1000 / samplerate)),
    Precision = by(data$speed, data[[trial]], function(s) mean(abs(s - mean(s, na.rm = TRUE)), na.rm = TRUE)),
    Velocity_thresholds = thres_vel,
    Duration_threshold = thres_dur,
    Speed = by(data$speed, data[[trial]], identity),
    Samplerate = samplerate,
    Head_target_distance = D,
    Height_px = height_px,
    Height_mm = height_mm,
    Width_px = width_px,
    Width_mm = width_mm,
    Fixations_and_saccades = simplified_results,
    Extra = extra
  )
  
  return(results)
}



# Plot function
plot.gazepath <- function(x, ..., trial_index = 1){
  i <- trial_index
  object <- x
  if(length(dim(object$Fixations_and_saccades[[i]])) == 2 && dim(object$Fixations_and_saccades[[i]])[2] == 4){
    warning('There is not enough data to identify fixations and saccades in this trial')
  } else {
    layout(matrix(c(1, 1:3), 2, 2))
    
    # Plot X and Y coordinates
    plot(object$x_coordinates[[i]], object$y_coordinates[[i]], xlab = "X", ylab = 'Y', type = 'l', las = 1, ylim = c(max(object$y_coordinates[[i]], na.rm = TRUE), min(object$y_coordinates[[i]], na.rm = TRUE)))
    sim <- object$Fixations_and_saccades[[i]]
    points(sim$Fixations[,2:3], pch = 19, cex = 1, col = 'blue')
    points(sim$Saccades[,2:3], pch = 1, cex = 1, col = 'red')
    
    # Plot position over time
    plot(object$x_coordinates[[i]], ylim = c(-50, max(object$y_coordinates[[i]], na.rm = TRUE) + 100), type = 'l', xlab = 'Time (msec)', ylab = 'Position', las = 1, xaxt = 'n', col = 5)
    lines(object$y_coordinates[[i]], col = 4)
    axis(1, at = seq(0, length(object$x_coordinates[[i]]), length.out = 6), labels = round(seq(0, length(object$x_coordinates[[i]]) * (1000 / object$Samplerate), length.out = 6)))
    fix <- object$Fixations_and_saccades[[i]]$Fixations[,3:4] / (1000 / object$Samplerate)
    rect(fix[,1], -50, fix[,2], 0, col = 'blue', border = NA)
    legend('topleft', c('X-coordinates', 'Y-coordinates', 'Fixations', 'Saccades'), col = c(5, 4, 'blue', 'red'), lwd = 2, bty = 'n', horiz = TRUE)
    
    # Plot speed over time
    if(object$Method != 'Tobii' & object$Method != 'Eyelink'){
      plot(object$Speed[[i]], type = 'l', xlab = 'Time (msec)', ylab = 'Speed (deg/s)', las = 1, xaxt = 'n')
      axis(1, at = seq(0, length(object$Speed[[i]]), length.out = 6), labels = round(seq(0, length(object$Speed[[i]]) * (1000 / object$Samplerate), length.out = 6)))
      if(object$Method == 'Mould.all' | object$Method == 'Mould.allDur'){
        segments(0, object$Velocity_thresholds, length(object$Speed[[i]]), object$Velocity_thresholds, col = 2, lwd = 2)
      } else {
        segments(0, object$Velocity_thresholds[[i]], length(object$Speed[[i]]), object$Velocity_thresholds[[i]], col = 2, lwd = 2)
      }
    }
  }
  layout(1)
}



# Plot function

plot.gazepath <- function(x, ..., trial_index = 1){
  i <- trial_index
  object <- x
  
  if(length(dim(object$Fixations_and_saccades[[i]])) == 2 && dim(object$Fixations_and_saccades[[i]])[2] == 4){
    warning('There is not enough data to identify fixations and saccades in this trial')
  } else {
    layout(matrix(c(1, 1:3), 2, 2))
    
    # Plot X and Y coordinates
    # Plot X and Y coordinates
    plot(object$x_coordinates[[i]], object$y_coordinates[[i]], xlab = "X", ylab = 'Y', type = 'l', las = 1, ylim = c(max(object$y_coordinates[[i]], na.rm = TRUE), min(object$y_coordinates[[i]], na.rm = TRUE)))
   
    # Extract the fixation and saccade data for the current trial
    sim <- object$Fixations_and_saccades[[i]]
    
    # Add fixation points (blue) and saccade points (red) to the plot
    points(sim$Fixations[,2:3], pch = 19, cex = 1, col = 'blue')
    points(sim$Saccades[,2:3], pch = 1, cex = 1, col = 'red')
    
    
    # Plot position over time without interpolation
    plot(object$x_coordinates[[i]], 
         ylim = c(-50, max(object$y_coordinates[[i]], na.rm = TRUE) + 100), 
         type = 'n', xlab = 'Time (msec)', ylab = 'Position', las = 1, xaxt = 'n', col = 5)
    lines(1:length(object$x_coordinates[[i]]), object$x_coordinates[[i]], col = 5)
    lines(1:length(object$y_coordinates[[i]]), object$y_coordinates[[i]], col = 4)
    axis(1, at = seq(0, length(object$x_coordinates[[i]]), length.out = 6), 
         labels = round(seq(0, length(object$x_coordinates[[i]]) * (1000 / object$Samplerate), 
                            length.out = 6)))
    
    if(!is.null(sim$Fixations)){
      fix <- sim$Fixations[, c("Start", "End")] / (1000 / object$Samplerate)
      rect(fix[,1], -50, fix[,2], 0, col = 'blue', border = NA)
    }
    legend('topleft', c('X-coordinates', 'Y-coordinates', 'Fixations', 'Saccades'), 
           col = c(5, 4, 'blue', 'red'), lwd = 2, bty = 'n', horiz = TRUE)
    
    # Plot speed over time without interpolation
    if(object$Method != 'Tobii' & object$Method != 'Eyelink'){
      plot(1:length(object$Speed[[i]]), object$Speed[[i]], type = 'n', xlab = 'Time (msec)', ylab = 'Speed (deg/s)', las = 1, xaxt = 'n')
      lines(1:length(object$Speed[[i]]), object$Speed[[i]])
      axis(1, at = seq(0, length(object$Speed[[i]]), length.out = 6), labels = round(seq(0, length(object$Speed[[i]]) * (1000 / object$Samplerate), length.out = 6)))
      if(object$Method == 'Mould.all' | object$Method == 'Mould.allDur'){
        segments(0, object$Velocity_thresholds, length(object$Speed[[i]]), object$Velocity_thresholds, col = 2, lwd = 2)
      } else {
        segments(0, object$Velocity_thresholds[[i]], length(object$Speed[[i]]), object$Velocity_thresholds[[i]], col = 2, lwd = 2)
      }
    }
  }
  layout(1)
}




# Function to plot segments without interpolation over NA values
plot_segments <- function(x, y, col = "black") {
  non_na_indices <- which(!is.na(y))
  if (length(non_na_indices) > 1) {
    segments <- split(non_na_indices, cumsum(c(1, diff(non_na_indices) != 1)))
    for (segment in segments) {
      lines(x[segment], y[segment], col = col)
    }
  }
}




plot.gazepath <-
  function(x, ..., trial_index = 1){
    i <- trial_index
    object <- x
    if(dim(object[[16]][[i]])[2] == 4){
      warning('There is not enough data to identify fixations and saccades in this trial')
    } else {
      layout(matrix(c(1, 1:3), 2, 2))
      plot(object[[2]][[i]], object[[3]][[i]], xlab = "X", ylab = 'Y', type = 'l', las = 1, ylim = c(max(object[[3]][[i]], na.rm = T), min(object[[3]][[i]], na.rm = T)))
      sim <- object[[16]][[i]]
      points(sim[sim[,1] == 'f', 9:10], pch = letters, cex = 3, col = 4)
      
      plot(object[[2]][[i]], ylim = c(-50, max(object[[14]][[i]], object[[12]][[i]]) + 100), type = 'l', xlab = 'Time (msec)', ylab = 'position', las = 1, xaxt = 'n', col = 5)
      lines(object[[3]][[i]], col = 4)
      axis(1, at = seq(0, length(object[[2]][[i]]), length.out = 6), labels = round(seq(0, length(object[[2]][[i]]) * (1000 / object[[10]]), length.out = 6)))
      fix <- sim[sim[,1] == 'f',3:4] / (1000 / object[[10]])
      rect(fix[,1], -50, fix[,2], 0, col = 3)
      legend('topleft', c('X-coordinates', 'Y-coordinates', 'Fixations'), col = 5:3, lwd = 2, bty = 'n', horiz = TRUE)
      
      if(object[[4]] != 'Tobii' & object[[4]] != 'Eyelink'){
        plot(object[[9]][[i]], typ = 'l', xlab = 'Time (msec)', ylab = 'Speed (deg/s)', las = 1, xaxt = 'n')
        axis(1, at = seq(0, length(object[[9]][[i]]), length.out = 6), labels = round(seq(0, length(object[[9]][[i]]) * (1000 / object[[10]]), length.out = 6)))
        if(object[[4]] == 'Mould.all' | object[[4]] == 'Mould.allDur'){
          segments(0, object[[7]], length(object[[9]][[i]]), object[[7]], col = 2, lwd= 2)
        } else {
          segments(0, object[[7]][[i]], length(object[[9]][[i]]), object[[7]][[i]], col = 2, lwd= 2)
        }
      }
    }
    layout(1)
  }


# Load your data
dta_142<- dta_4gazepath_stim%>%
  subset(ssid == 122)

# Ensure the trial column is a character or factor
dta_140$trial_no <- as.character(dta_140$trial_no)

# Run the gazepath function
ppt142 <- gazepath(dta_142, x1 = 17, y1 = 19, x2 = 18, y2 = 20,
                   d1 = 23, d2 = 23,
                   trial = 28, # Ensure this is a character or factor
                   height_px = 1080, height_mm = 295,
                   width_px = 1920, width_mm = 525, 
                   # extra_var = c("text","participant", "trial_no_gpid_stim","trial_no"),
                   method = 'gazepath', 
                   samplerate = 150)

# View unique trials
unique(dta_140$trial_no)

# Access extra variables
ppt140$Extra



# Generate a summary of the results
summary_result <- summary.gazepath(ppt142)

# Plot the results for a specific trial
plot.gazepath(ppt142, trial_index = 1)

