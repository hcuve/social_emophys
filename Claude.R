# Helper Functions

lomax <- function(x){
  which(diff(c(TRUE, diff(x) >= 0, FALSE)) < 0)
}

Mould_vel <- function(speed, Hz, plot = F){
  lmax <- speed[lomax(speed)]
  if(length(lmax) < 10){
    return(NA)
    warning('There are not enough data points to estimate a velocity threshold')
  } else {
    thresholds <- seq(min(lmax), max(lmax), length.out = Hz)
    set <- sapply(thresholds, function(x) {length(which(lmax > x))})
    uni <- seq(length(lmax), 0, length.out = Hz)
    gap <- uni - set
    
    if(Hz < 250) {h <- .1} else {h <- .05}
    
    gap <- predict(loess(gap ~ log(thresholds), span = h))
    while (length(lomax(gap)) > 1 & h < 1) {
      h <- h + .01
      gap <- predict(loess((uni - set) ~ log(thresholds), span = h, surface = 'direct', cell = 1))
    }
    
    if(plot == T) plotMould(uni, set, gap, thresholds, lmax, Hz)
    if(h != 1) return(thresholds[which.max(gap)]) else return(NA)
  }
}

comhull <- function(d, classification, dat_x, dat_y, in_thres, Hz = Hz, M, D, res_x = res_x, width_mm = width_mm){
  d <- d[d$dur > 1,]
  fix <- tail(which(d$index == 'fixation'), 1)
  count <- length(which(d$index == 'fixation')) - 1
  while(count >= 1){
    fix2 <- which(d$index == 'fixation')[count]
    cvhull <- chull(cbind(dat_x, dat_y)[d[fix,3] : d[fix,4],])
    POLY_FIX <- cbind(dat_x[d[fix,3] : d[fix,4]][cvhull], dat_y[d[fix,3] : d[fix,4]][cvhull])
    PNT <- sum(pnt.in.poly(cbind(dat_x, dat_y)[d[fix2,3] : d[fix2,4],], POLY_FIX)[,3])
    if(PNT != 0){
      dis <- dist(rbind(t(apply(cbind(dat_x, dat_y)[d[fix,3] : d[fix,4],], 2, mean)),
                        t(apply(cbind(dat_x, dat_y)[d[fix2,3] : d[fix2,4],], 2, mean))))
      thres_d <- atan((width_mm/2)/D) * (180/pi) * 2 *(dis/res_x)
      if(thres_d < M & (d[fix,3] - d[fix2,4]) < in_thres * (Hz / 1000)){
        classification[d[fix2,3] : d[fix,4]] <- 'fixation'
        CL <- rle(classification)
        index <- rep.int(1:length(CL$value), CL$lengths)
        POG <- sapply(unique(index[!is.na(index)]), function(i) mean(dist(cbind(dat_x[index == i], dat_y[index == i])), na.rm = T))
        POG[is.na(POG)] <- 0
        mean_x <- as.vector(by(dat_x, index, function(i) mean(i, na.rm = T)))
        mean_y <- as.vector(by(dat_y, index, function(i) mean(i, na.rm = T)))
        
        dat_x[d[fix2,3] : d[fix,4]] <- na.approx(dat_x[d[fix2,3] : d[fix,4]])
        dat_y[d[fix2,3] : d[fix,4]] <- na.approx(dat_y[d[fix2,3] : d[fix,4]])
        
        d <- data.frame(CL$value, CL$length, c(1, cumsum(CL$length)[-length(CL$length)] + 1), cumsum(CL$length), POG, mean_x, mean_y)
        names(d)[1:4] <- c('index', 'dur', 'start', 'end')
        d <- d[d$dur > 1,]
      } 
    }
    fix <- fix2
    count <- count - 1
  }
  return(list(classification, dat_x, dat_y))
}

MFW <- function(x){
  if(sum(is.na(x)) / length(x) > .5){
    return(rep(NA, length(x)))
  } else {
    return(rep(median(x, na.rm = T), length(x)))
  }
}

Boundary <- function(X, min, max){
  X <- ifelse(X < min | X > max,  NA, X)
  return(X)
}

Speed_Deg <- function(X, Y, distance, height_mm, width_mm, height_px, width_px, Hz){
  hor <- atan((width_mm / 2) / distance) * (180 / pi) * 2 / width_px * X
  ver <- atan((height_mm / 2) / distance) * (180 / pi) * 2 / height_px * Y
  speed <- sqrt(diff(hor,2) ^ 2 + diff(ver,2) ^ 2) * (Hz/2)
  return(c(.001, speed, 0.001))
}

Speed <- function(X, Y, distance, height_mm, width_mm, height_px, width_px, res_x = 1280, res_y = 1024, Hz){
  d2p_px <- sqrt(abs(X - (res_x / 2))**2 + abs(Y - (res_y / 2))**2 + (distance * width_px / width_mm)**2)
  dbp_px <- c(.001, sqrt(diff(X, 2)**2 + diff(Y, 2)**2), .001)
  speed <- atan((dbp_px / 2) / d2p_px) * (180 / pi) * 2 * (Hz / 2)
  return(speed)
}

robust <- function(x, Hz) {
  mis <- rle(ifelse(is.na(x), 0, 1))
  return(mean(mis$length[mis$values == 1]) / Hz)
}

precision <- function(X, Hz){
  end <- length(X) / (Hz / 10)
  window <- rep(1:end, each = (Hz / 10))
  Smooth <- as.vector(unlist(by(X[1:length(window)], window, MFW)))
  
  return(mean(abs(X[1:length(window)] - Smooth), na.rm = T))
}

simplify <- function(classification, x, y, Hz, D, width_px, width_mm, extra, extra_var){
  class <- rle(classification)
  simple <- data.frame(class$values, class$lengths, 
                       c(1, cumsum(class$lengths) + 1)[-(length(class$values) + 1)], 
                       cumsum(class$lengths))
  if(length(which(class$values == 'f')) > 0){
    x_start <- y_start <- x_end <- y_end <- mean_x <- mean_y <- POGvar <- RMS <- numeric()
    for(i in 1:dim(simple)[1]){
      x_start <- c(x_start, x[simple[i, 3]])
      y_start <- c(y_start, y[simple[i, 3]])
      x_end <- c(x_end, x[simple[i, 4]])
      y_end <- c(y_end, y[simple[i, 4]])
      mean_x <- c(mean_x, mean(x[simple[i,3] : simple[i,4]]))
      mean_y <- c(mean_y, mean(y[simple[i,3] : simple[i,4]]))
      m <- as.matrix(dist(cbind(c(mean_x[length(mean_x)], x[simple[i,3] : simple[i,4]]), c(mean_y[length(mean_y)], y[simple[i,3] : simple[i,4]]))))
      RMS <- c(RMS, sqrt(mean((atan((diag(m[-1,-c(1,2)]) / 2) / mean(D, na.rm = T)) * (180 / pi) * (width_mm / width_px) * 2)**2)))
      POGvar <- c(POGvar, mean(m[-1,1]))
    }
    ss <- which(class$values == 's')
    POGvar[ss] <- sqrt((x_start[ss] - x_end[ss]) ^ 2 + (y_start[ss] - y_end[ss]) ^ 2)
    POGsdSacAmp <- atan((POGvar / 2) / mean(D, na.rm = T)) * (180 / pi) * (width_mm / width_px) * 2
    POGsdSacAmp[!ss] <- sqrt(POGsdSacAmp)
    RMS[ss] <- NA
    
    simple <- data.frame(class$values, class$lengths * (1000/Hz), 
                         c(1, cumsum(class$lengths * (1000/Hz)) + 1)[-(length(class$values) + 1)], 
                         cumsum(class$lengths * (1000/Hz)), x_start, y_start, x_end, y_end, mean_x, mean_y, POGsdSacAmp, RMS)
    names(simple)[1:4] <- c('Value', 'Dur', 'Start', 'End')
    if(!is.null(extra_var)){
      for(i in 1:length(extra_var)){
        simple <- data.frame(simple, extra[i])
        names(simple)[dim(simple)[2]] <- extra_var[i]
      }
    }
  }
  if(length(which(is.na(simple[,1]))) != 0){
    simple <- simple[-which(is.na(simple[,1])),]
  }
  return(simple)
}

Interpolate <- function(X, Y, D, height_mm, width_mm, height_px, width_px, res_x = res_x, res_y = res_y, Hz = Hz, in_thres = in_thres, thres_dur = thres_dur){
  
  s <- Speed(X, Y, D, height_mm, width_mm, height_px, width_px, res_x = res_x, res_y = res_y, Hz)
  s <- ifelse(s > 1000, NA, s)
  if(length(lomax(s)) < 10){
    return(list('No Return', 'No Return','No Return','No Return','No Return','No Return','No Return','No Return'))
  } else {
    M <- Mould_vel(s, Hz)
    
    classification <- ifelse(s > M, 'saccade', 'fixation')
    classification[is.na(classification)] <- 'missing'
    
    CL <- rle(classification)
    d <- data.frame(CL$value, CL$length, c(1, cumsum(CL$length)[-length(CL$length)] + 1), cumsum(CL$length))
    names(d) <- c('index', 'dur', 'start', 'end')
    
    dat_x <- X
    dat_y <- Y
    dat_d <- D
    
    for(i in which(d$index == 'missing')){
      if(i > 1 & i < dim(d)[1] & d[i, 2] < (in_thres * (Hz / 1000))){
        if(d[i + 1, 1] == 'fixation' & d[i - 1, 1] == 'fixation'){
          ii_s <- d[i - 1, 4]
          ii_e <- d[i + 1, 3]
          speed <- Speed(c(dat_x[ii_s], dat_x[ii_s], dat_x[ii_e]), c(dat_y[ii_s], dat_y[ii_s], dat_y[ii_e]), c(dat_d[ii_s], dat_d[ii_s], dat_d[ii_e]), height_mm, width_mm, height_px, width_px, res_x = res_x, res_y = res_y, Hz)
          if(speed[2] < M){
            dat_x[d[i, 3] : d[i, 4]] <- dat_x[ii_s]
            dat_y[d[i, 3] : d[i, 4]] <- dat_y[ii_s]
            dat_d[d[i, 3] : d[i, 4]] <- dat_d[ii_s]
          }
        }
      }
    }
    
    s <- Speed(dat_x, dat_y, dat_d, height_mm, width_mm, height_px, width_px, res_x = res_x, res_y = res_y, Hz)
    s <- ifelse(s > 1000, NA, s)
    
    classification <- ifelse(s > M, 'saccade', 'fixation')
    classification[is.na(classification)] <- 'missing'
    CL <- rle(classification)
    
    index <- rep.int(1:length(CL$value), CL$lengths)
    POG <- sapply(unique(index[!is.na(index)]), function(i) mean(dist(cbind(dat_x[index == i], dat_y[index == i])), na.rm = T))
    POG[is.na(POG)] <- 0
    mean_x <- as.vector(by(dat_x, index, function(i) mean(i, na.rm = T)))
    mean_y <- as.vector(by(dat_y, index, function(i) mean(i, na.rm = T)))
    
    d <- data.frame(CL$value, CL$length, c(1, cumsum(CL$length)[-length(CL$length)] + 1), cumsum(CL$length), POG, mean_x, mean_y)
    names(d)[1:4] <- c('index', 'dur', 'start', 'end')
    
    dimd_new <- dim(d)[1] + 1
    while(dimd_new != dim(d)[1]){
      dimd_new <- dim(d)[1]
      classif <- comhull(d, classification, dat_x, dat_y, in_thres, Hz, M = sqrt(M)/10, mean(dat_d, na.rm = T), res_x = res_x, width_mm = width_mm)
      
      CL <- rle(classif[[1]])
      classification <- classif[[1]]
      index <- rep.int(1:length(CL$value), CL$lengths)
      dat_x <- classif[[2]]
      dat_y <- classif[[3]]
      POG <- sapply(unique(index[!is.na(index)]), function(i) mean(dist(cbind(dat_x[index == i], dat_y[index == i])), na.rm = T))
      POG[is.na(POG)] <- 0
      mean_x <- as.vector(by(dat_x, index, function(i) mean(i, na.rm = T)))
      mean_y <- as.vector(by(dat_y, index, function(i) mean(i, na.rm = T)))
      
      d <- data.frame(CL$value, CL$length, c(1, cumsum(CL$length)[-length(CL$length)] + 1), cumsum(CL$length), POG, mean_x, mean_y)
      names(d)[1:4] <- c('index', 'dur', 'start', 'end')
    }
    
    for(i in which(CL$value == 'fixation' & CL$length < (Hz / 1000 * thres_dur))){
      classification[((cumsum(CL$length) - CL$length) + 1)[i] : cumsum(CL$length)[i]] <- 'saccade' 
    }
    
    CL <- rle(classification)
    index <- rep.int(1:length(CL$value), CL$lengths)
    POG <- sapply(unique(index[!is.na(index)]), function(i) mean(dist(cbind(dat_x[index == i], dat_y[index == i])), na.rm = T))
    POG[is.na(POG)] <- 0
    mean_x <- as.vector(by(dat_x, index, function(i) mean(i, na.rm = T)))
    mean_y <- as.vector(by(dat_y, index, function(i) mean(i, na.rm = T)))
    
    d <- data.frame(CL$value, CL$length, c(1, cumsum(CL$length)[-length(CL$length)] + 1), cumsum(CL$length), POG, mean_x, mean_y)
    names(d)[1:4] <- c('index', 'dur', 'start', 'end')
    
    classification <- comhull(d, classification, dat_x, dat_y, in_thres, Hz,  M = sqrt(M)/10, D = mean(dat_d, na.rm = T), res_x = res_x, width_mm = width_mm)
    
    clas <- classification[[1]]
    CL <- rle(clas)
    dat_x <- classification[[2]]
    dat_y <- classification[[3]]
    
    index <- rep.int(1:length(CL$value), CL$lengths)
    mean_x <- as.vector(by(dat_x, index, function(i) mean(i, na.rm = T)))
    mean_y <- as.vector(by(dat_y, index, function(i) mean(i, na.rm = T)))
    
    index <- CL$value
    end <- cumsum(CL$length) * (1000 / Hz)
    dur <- CL$length * (1000 / Hz)
    start <- (end - dur) + 1
    
    d <- data.frame(index, dur, start, end, mean_x, mean_y)
    d <- data.frame(d, order=1:dim(d)[1])
    
    return(list(dat_x, dat_y, dat_d, d, M, s, clas, 'Return'))
  }
}

# Main gazepath function
gazepath <- function(data, x1, y1, x2 = NULL, y2 = NULL, d1, d2 = NULL, trial, height_px, height_mm, width_px, width_mm, extra_var = NULL, res_x = 1280, res_y = 1024, samplerate = 500, method = 'Mould', posthoc = FALSE, thres_vel = 35, thres_dur = 100, min_dist = 250, in_thres = 150){
  if(!is.data.frame(data)) {
    stop('please insert a data frame and define the column numbers of the variables')
  }
  
  # Ensure trial is a vector, not a data frame
  if(is.data.frame(data[,trial])) {
    stop('The trial column appears to be a data frame. Please ensure it is a vector.')
  }
  
  # Convert trial column to factor if it's not already
  if(!is.factor(data[,trial])) {
    data[,trial] <- as.factor(data[,trial])
  }
  
  # Check for unique trials
  trial_levels <- levels(data[,trial])
  trial_counts <- table(data[,trial])
  
  if(length(trial_levels) != length(trial_counts)){
    TRIAL_NEW <- as.numeric(data[,trial])
    data$TRIAL_NEW <- TRIAL_NEW
    if(is.numeric(trial)){
      names(data)[trial] <- 'TRIAL_OLD'
      trial <- which(names(data) == 'TRIAL_NEW')
    } else {
      data$TRIAL_OLD <- data[,trial]
      trial <- which(names(data) == 'TRIAL_NEW')
    }
    warning('The trial index in the data frame was not unique, therefore trials are renamed to be unique and the old trial index is stored in the data as TRIAL_OLD')
  }
  
  extra <- list()
  
  if(!is.null(extra_var)){
    if(sum(extra_var %in% names(data)) == length(extra_var)){
      for(i in unique(data[,trial])){
        extra[[i]] <- sapply(1:length(extra_var), function(j) head(as.character(data[data[,trial] == i, which(names(data) == extra_var[j])]), 1))
      }
    } else {
      extra_var <- NULL
      print('Please make sure the variables to pass through have the correct names')
    }
  }
  data[,d1] <- ifelse(data[,d1] < min_dist, NA, data[,d1])
  if(is.null(d2)){
    D <- by(data[,d1], data[,trial], data.frame)
  } else {
    data[,d2] <- ifelse(data[,d2] < min_dist, NA, data[,d2])
    D <- by((data[,d1] + data[,d2]) / 2, data[,trial], data.frame)
  }
  
  if(!is.null(x2) & !is.null(y2)){
    data[,x1] <- ifelse(is.na(data[,x1]), data[,x2], data[,x1])
    data[,y1] <- ifelse(is.na(data[,y1]), data[,y2], data[,y1])
    data[,x2] <- ifelse(is.na(data[,x2]), data[,x1], data[,x2])
    data[,y2] <- ifelse(is.na(data[,y2]), data[,y1], data[,y2])
    X <- by((data[,x1] + data[,x2]) / 2, data[,trial], data.frame)
    Y <- by((data[,y1] + data[,y2]) / 2, data[,trial], data.frame)
  } else {
    X <- by(data[,x1], data[,trial], data.frame)
    Y <- by(data[,y1], data[,trial], data.frame)
  }
  
  Rob_x <- sapply(1:length(X), function(i) robust(X[[i]], samplerate))
  Rob_y <- sapply(1:length(X), function(i) robust(Y[[i]], samplerate))
  Robustness <- (Rob_x + Rob_y) / 2
  Pre_x <- sapply(1:length(X), function(i) precision(X[[i]], samplerate))
  Pre_y <- sapply(1:length(X), function(i) precision(Y[[i]], samplerate))
  Precision <- (Pre_x + Pre_y) / 2
  
  if(length(height_px) == 1) height_px <- rep(height_px, length(unique(data[,trial])))
  if(length(height_mm) == 1) height_mm <- rep(height_mm, length(unique(data[,trial])))
  if(length(width_px) == 1) width_px <- rep(width_px, length(unique(data[,trial])))
  if(length(width_mm) == 1) width_mm <- rep(width_mm, length(unique(data[,trial])))
  
  final <- 'Please insert a correct method'
  s <- NA
  
  if(method == 'gazepath'){
    fix <- thres_vel <- numeric()
    final <- s <- list()
    for(i in 1:length(unique(data[,trial]))){
      X[[i]] <- Boundary(X[[i]], (res_x - width_px[i]) / 2, res_x - (res_x - width_px[i]) / 2)
      Y[[i]] <- Boundary(Y[[i]], (res_y - height_px[i]) / 2, res_y - (res_y - height_px[i]) / 2)
      if(length(which(!is.na(X[[i]]))) > samplerate & length(which(!is.na(Y[[i]]))) > samplerate & length(which(!is.na(D[[i]]))) > samplerate){
        interpol <- Interpolate(X[[i]], Y[[i]], D[[i]], height_mm[i], width_mm[i], height_px[i], width_px[i], res_x = res_x, res_y = res_y, Hz = samplerate, in_thres = in_thres, thres_dur = thres_dur)
        if(interpol[[8]] == 'Return'){
          final[[i]] <- ifelse(interpol[[7]] == 'missing', NA, ifelse(interpol[[7]] == 'fixation', 'f', 's'))
          thres_vel[i] <- interpol[[5]]
          s[[i]] <- interpol[[6]]
          X[[i]] <- interpol[[1]]
          Y[[i]] <- interpol[[2]]
        } else {
          s[[i]] <- NA; thres_vel[i] <- NA; final[[i]] <- NA
        }
      } else {
        s[[i]] <- NA; thres_vel[i] <- NA; final[[i]] <- NA
      }
    }
  }
  
  if(posthoc == TRUE){
    for(i in 1:length(X)) {
      PH <- posthocCheck(final[[i]], X[[i]], Y[[i]])
      final[[i]] <- PH[[1]]
      X[[i]] <- PH[[2]]
      Y[[i]] <- PH[[3]]
    }
  }
  sim <- list()
  for(i in 1:length(X)){
    sim[[i]] <- simplify(final[[i]], X[[i]], Y[[i]], samplerate, D[[i]], width_px[i], width_mm[i], extra[[i]], extra_var)
  }
  
  output <- list(final, X, Y, method, Robustness, Precision, thres_vel, thres_dur, s, samplerate, D, height_px, height_mm, width_px, width_mm, sim)
  names(output) <- c('classifications', 'x-coor', 'y-coor', 'method', 'robustness',
                     'precision', 'vel_thres', 'dur_thres', 'speed', 'samplerate',
                     'distance', 'height_px', 'height_mm', 'width_px', 'width_mm',
                     'fixations')
  class(output) <- 'gazepath'
  return(output)
}



# Load your data
dta_140<- dta_4gazepath_stim%>%
  subset(ssid == 140)

# Ensure the trial column is a character or factor
dta_140$trial_no <- as.character(dta_140$trial_no)

dta_140$trial_no<-as.factor(as.numeric(as.character(dta_140$trial_no)))

dta_140$trial_no<- as.vector(dta_140$trial_no)

install.packages("sp")
install.packages("sdm")
length(unique(dta_140[, 28])) == length(levels(dta_140[, 28]))
colnames(dta_140)
# Run the gazepath function
ppt140 <- gazepath(dta_140, x1 = 17, y1 = 19, x2 = 18, y2 = 20,
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
summary_result <- summary.gazepath(ppt140)

# Plot the results for a specific trial
plot.gazepath(ppt142)