##This is a tiny function for matrix which is for drawing heatmap...

#This returns top values for each peak
#ARGV1: matrix with each peak as row, each window as column.
return_top <- function (m) {
  r <- rep(-10000, nrow(m)) ##vector of minimal values with the length of peak number
  for (i in 1:nrow(m)) {
    for (j in 1:ncol(m)) {
       if (r[i] < m[i,j]) {r[i] <- m[i,j]}
    }
  }

  return (r)
}

#This returns the window of top values for each peak
#ARGV1: matrix with each peak as row, each window as column.
return_top_win <- function (m) {
  win <- vector("numeric", length=nrow(m))
  for (i in 1:nrow(m)) {
    r <- -10000 #minimal value
    for (j in 1:ncol(m)) {
       if (r < m[i,j]) {
         r <- m[i,j]
         win[i] <- j
       }
    }
  }

  return (win)
}

#This returns the matrix (ncol=2, nrow=nrow(m)) composed of x position and y as rows. For each row of input matrix, this scans the col and put the right most col number which satisfy the threshold.
#ARGV1: matrix with each peak as row, each window as column.
return_thresh_mat <- function (m, thre=1, continuous=F) {
  mat <- matrix(0, nrow=nrow(m), ncol=2)
  for (i in 1:nrow(m)) {
    s <- 1
    c <- 1 #left most pos
    if (continuous) {
      for (j in 1:ncol(m)) {
        if (thre <= m[i,j] & s >= continuous) { ##if s is eq or more than continuous
          c <- j
          s <- s + 1
        } else if (thre <= m[i,j]) {
          s <- s + 1
        } else {s <- 1} ##reset s
      }
    } else {
      for (j in 1:ncol(m)) {
        if (thre <= m[i,j]) {
          c <- j
        }
      }
    }

    mat[i,1] <- i
    mat[i,2] <- c
  }

  return (mat)
}

#This returns the window of bottom values for each peak
#ARGV1: matrix with each peak as row, each window as column.
return_bottom_win <- function (m) {
  win <- vector("numeric", length=nrow(m))
  for (i in 1:nrow(m)) {
    r <- 10000 #minimal value
    for (j in 1:ncol(m)) {
       if (r > m[i,j]) {
         r <- m[i,j]
         win[i] <- j
       }
    }
  }

  return (win)
}

#This returns the log2 ratio of specific windows for each peak
#ARGV1: matrix with each peak as row, each window as column.
return_ratio_win <- function (m, s1, e1, s2, e2) {
  ratio <- vector("numeric", length=nrow(m))
  for (i in 1:nrow(m)) {
    ratio[i] <- log2 (mean(as.numeric(m[i,s2:e2])) / mean(as.numeric(m[i,s1:e1])))
  }

  return (ratio)
}

#This returns the difference {(s2,e2) - (s1,e1)} of specific windows for each peak
#ARGV1: matrix with each peak as row, each window as column.
return_diff_win <- function (m, s1, e1, s2, e2) {
  diff <- vector("numeric", length=nrow(m))
  for (i in 1:nrow(m)) {
    diff[i] <- mean(as.numeric(m[i,s2:e2])) - mean(as.numeric(m[i,s1:e1]))
  }

  return (diff)
}

#This returns the normalized matrix with value range from 0 to 1 for each summit.
#ARGV1: matrix with each peak as row, each window as column.
return_mat_norm <- function (m) {
  mat <- m
  for (i in 1:nrow(m)) {
    val_min <- 1e5
    val_max <- 1e-5
    for (j in 1:ncol(m)) {
       if (is.na(m[i,j])) {next}
       if (val_min > m[i,j]) {val_min <- m[i,j]}
       if (val_max < m[i,j]) {val_max <- m[i,j]}
    }
    x <- val_max - val_min ##difference of max - min
    for (j in 1:ncol(m)) { #normalizing
      mat[i,j] <- (m[i,j] - val_min)/ x
    }
  }

  return (mat)
}


#This returns the normalized matrix with value range from 0 to 1 for all summits.
#ARGV1: matrix with each peak as row, each window as column.
#ARGV2: 0 < thresh < 1. Values above this threshold is considered as 1.
return_mat_norm_all <- function (m, thresh=1) {
  mat <- m
  val_min <- 0
  vals <- c()
  for (i in 1:nrow(m)) {
    vals <- c(vals,as.numeric(m[i,]))
  }
  val_max <- as.numeric(quantile(vals,probs=thresh)) ##val_max is quantile of vals
  x <- val_max - val_min ##difference of max - min

  for (i in 1:nrow(m)) {
    for (j in 1:ncol(m)) { #normalizing
      if (m[i,j] >= val_max) {mat[i,j] <- 1}
      else { mat[i,j] <- (m[i,j] - val_min)/ x }
    }
  }

  return (mat)
}


#This returns the matrix where NA was replaced to zero.
#ARGV1: matrix with each peak as row, each window as column.
return_mat_NA2zero <- function (m) {
  mat <- m
  for (i in 1:nrow(m)) {
    for (j in 1:ncol(m)) {
       if (is.na(m[i,j])) {
         mat[i,j] <- 0
       }
    }
  }

  return (mat)
}





