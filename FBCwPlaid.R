FBCwPlaid <- function(data_input, weight_input, iter.startup, 
                      max.layers, iter.layer, iter.bin) {
  initial_dec <- function(data, side, iter.startup) {
    if (side == 1) {
      row_test <- try(kmeans(data, 2, iter.startup), TRUE)
      row_index_all <- row_test$cluster
      row_index_one <- which(row_index_all == 1)
      row_index_two <- which(row_index_all == 2)
      row_one_len <- length(row_index_one)
      row_two_len <- length(row_index_two)
      if (row_two_len <= row_one_len) {
        initial_row_num <- row_two_len
        initial_row <- row_index_two
      } else {
        initial_row_num <- row_one_len
        initial_row <- row_index_one
      }
      row_plus <- runif(initial_row_num, 0, 0.1)
      row_reduce <- runif((nrow(data) - initial_row_num), -0.1, 0)
      row_dec <- matrix(nrow = nrow(data), ncol = ncol(data))
      post_plus <- 1
      post_reduce <- 1
      for (var in 1:nrow(data)) {
        if (var %in% initial_row) {
          row_dec[var, ] <- 0.5 + row_plus[post_plus]
          post_plus <- post_plus + 1
        } else {
          row_dec[var, ] <- 0.5 + row_reduce[post_reduce]
          post_reduce <- post_reduce + 1
        }
      }
      return(row_dec)
    }
    
    if (side == 2) {
      col_test <- try(kmeans(t(data), 2, iter.startup), TRUE)
      col_index_all <- col_test$cluster
      col_index_one <- which(col_index_all == 1)
      col_index_two <- which(col_index_all == 2)
      col_one_len <- length(col_index_one)
      col_two_len <- length(col_index_two)
      if (col_two_len <= col_one_len) {
        initial_col_num <- col_two_len
        initial_col <- col_index_two
      } else {
        initial_col_num <- col_one_len
        initial_col <- col_index_one
      }
      col_plus <- runif(initial_col_num, 0, 0.1)
      col_reduce <- runif((ncol(data) - initial_col_num), -0.1, 0)
      col_dec <- matrix(nrow = nrow(data), ncol = ncol(data))
      post_plus <- 1
      post_reduce <- 1
      for (var in 1:ncol(data)) {
        if (var %in% initial_col) {
          col_dec[, var] <- 0.5 + col_plus[post_plus]
          post_plus <- post_plus + 1
        } else {
          col_dec[, var] <- 0.5 + col_reduce[post_reduce]
          post_reduce <- post_reduce + 1
        }
      }
      return(col_dec)
    }
  }
  
  
  updata_theta <- function(data, row_fun, col_fun) {
    row_ss <- sum(row_fun[, 1] ^ 2)
    col_ss <- sum(col_fun[1, ] ^ 2)
    if ((row_ss * col_ss) == 0) {
      mu <- 0
      alpha <- vector()
      beta <- vector()
      alpha[1:nrow(data)] <- 0
      beta[1:ncol(data)] <- 0
      mu_matrix <- matrix(mu, nrow = nrow(data), ncol = ncol(data))
    } else {
      mu <- sum(row_fun * col_fun * data) / (row_ss * col_ss)
      mu_matrix <- matrix(mu, nrow = nrow(data), ncol = ncol(data))
      data_mu <- data - (mu * row_fun * col_fun)
      alpha <- vector()
      for (var in 1:nrow(data)) {
        if (row_fun[var, 1] < 0.001) {
          alpha[var] <- 0
        } else {
          alpha[var] <- sum(data_mu[var, ] * col_fun[var, ]) / (row_fun[var, 1] * col_ss)
        }
      }
      beta <- vector()
      for (var in 1:ncol(data)) {
        if (col_fun[1, var] < 0.001) {
          beta[var] <- 0
        } else {
          beta[var] <- sum(data_mu[, var] * row_fun[, var]) / (col_fun[1, var] * row_ss)
        }
      }
    }
    alpha_matrix <- matrix(rep(alpha, ncol(a)), ncol = ncol(a))
    beta_matrix <- t(matrix(rep(beta, nrow(a)), nrow = ncol(a)))
    theta <- mu_matrix + alpha_matrix + beta_matrix
    out <- list(mu_out = mu, alpha_out = alpha_matrix, beta_out = beta_matrix, theta_out = theta)
    return(out)
  }
  
  
  updata_dec <- function(data, theta_fun, row_fun, col_fun, layer_cycle, 
                         iter.layer, iter.bin, side) {
    if (side == 1) {
      row_dec_test <- vector()
      for (var in 1:nrow(data)) {
        row_dec_test[var] <- sum(theta_fun[var, ] * col_fun[var, ] * data[var, ]) /
          sum((theta_fun[var, ] ^ 2) * (col_fun[var, ] ^ 2))
        if (layer_cycle < iter.layer - iter.bin) {
          if (row_dec_test[var] > 0.5) {
            row_dec_test[var] <- 0.5 + (layer_cycle / (2  * (iter.layer - iter.bin)))
          }
          if (row_dec_test[var] < 0.5) {
            row_dec_test[var] <- 0.5 - (layer_cycle / (2  * (iter.layer - iter.bin)))
          }
        }
        if (layer_cycle >= iter.layer - iter.bin) {
          if (row_dec_test[var] > 0.5) {
            row_dec_test[var] <- 1
          }
          if (row_dec_test[var] < 0.5) {
            row_dec_test[var] <- 0
          }
        }
      }
      return(row_dec_test)
    }
    if (side == 2) {
      col_dec_test <- vector()
      for (var in 1:ncol(data)) {
        col_dec_test[var] <- sum(theta_fun[, var] * row_fun[, var] * data[, var]) /
          sum((theta_fun[, var] ^ 2) * (row_fun[, var] ^ 2))
        if (layer_cycle < iter.layer - iter.bin) {
          if (col_dec_test[var] > 0.5) {
            col_dec_test[var] <- 0.5 + (layer_cycle / (2  * (iter.layer - iter.bin)))
          }
          if (col_dec_test[var] < 0.5) {
            col_dec_test[var] <- 0.5 - (layer_cycle / (2  * (iter.layer - iter.bin)))
          }
        }
        if (layer_cycle >= iter.layer - iter.bin) {
          if (col_dec_test[var] > 0.5) {
            col_dec_test[var] <- 1
          }
          if (col_dec_test[var] < 0.5) {
            col_dec_test[var] <- 0
          }
        }
      }
      return(col_dec_test)
    }
  }
  
  
  SDwC_cal <- function(data) {
    SDwC_temp <- sum((data - mean(data)) ^ 2) / (nrow(data) * ncol(data))
    return(SDwC_temp)
  }
  
  
  drop_release_row <- function(data, row_fun, col_fun, iter.startup) {
    row_second <- row_fun
    col_second <- col_fun
    row_index <- which(row_fun[, 1] == 1)
    col_index <- which(col_fun[1, ] == 1)
    data_min <- matrix(nrow = length(row_index), ncol = length(col_index))
    data_min <- data[row_index, col_index]
    SDwC_original <- SDwC_cal(data_min)
    mean_original <- mean(data_min)
    if (nrow(data_min) > 2) {
      row_test <- try(kmeans(data_min, 2, iter.startup), TRUE)
      row_index_all <- row_test$cluster
      row_index_one <- which(row_index_all == 1)
      row_index_two <- which(row_index_all == 2)
      row_one <- as.numeric(row_index[-c(row_index_one)])
      row_two <- as.numeric(row_index[-c(row_index_two)])
      data_min_one <- as.matrix(data[row_one, col_index])
      data_min_two <- as.matrix(data[row_two, col_index])
      SDwC_one <- SDwC_cal(data_min_one)
      SDwC_two <- SDwC_cal(data_min_two)
      mean_one <- mean(data_min_one)
      mean_two <- mean(data_min_two)
      row_final <- vector()
      
      if (mean_one >= mean_original && SDwC_one <= SDwC_original) {
        row_final <- row_one
        SDwC_min <- SDwC_one
      }
      if (mean_two >= mean_original && SDwC_two <= SDwC_original) {
        row_final <- row_two
        SDwC_min <- SDwC_two
      }
      
      if (length(row_final) != 0) {
        for (var in 1:nrow(data)) {
          if (var %in% row_final) {
            row_second[var, ] <- 1
          } else {
            row_second[var, ] <- 0
          }
        }
      } 
    }
    out <- list(row = row_second, col = col_second)
    return(out)
  }
  
  
drop_release_col <- function(data, row_fun, col_fun, iter.startup) {
    row_second <- row_fun
    col_second <- col_fun
    row_index <- which(row_fun[, 1] == 1)
    col_index <- which(col_fun[1, ] == 1)
    data_min <- matrix(nrow = length(row_index), ncol = length(col_index))
    data_min <- data[row_index, col_index]
    SDwC_original <- SDwC_cal(data_min)
    mean_original <- mean(data_min)
    if (nrow(t(data_min)) > 2) {
      col_test <- try(kmeans(t(data_min), 2, iter.startup), TRUE)
      col_index_all <- col_test$cluster
      col_index_one <- which(col_index_all == 1)
      col_index_two <- which(col_index_all == 2)
      col_one_len <- length(col_index_one)
      col_two_len <- length(col_index_two)
      
      col_one <- as.numeric(col_index[-c(col_index_one)])
      col_two <- as.numeric(col_index[-c(col_index_two)])
      data_min_one <- as.matrix(data[row_index, col_one])
      data_min_two <- as.matrix(data[row_index, col_two])
      SDwC_one <- SDwC_cal(data_min_one)
      SDwC_two <- SDwC_cal(data_min_two)
      mean_one <- mean(data_min_one)
      mean_two <- mean(data_min_two)
      col_final <- vector()
      
      if (length(col_one) != 1 && length(col_two) != 1) {
        if (mean_one > mean_original && SDwC_one < SDwC_original) {
          col_final <- col_one
          SDwC_min <- SDwC_one
        }
        if (mean_two > mean_original && SDwC_two < SDwC_original) {
          col_final <- col_two
          SDwC_min <- SDwC_two
        }
      }
      
      if (length(col_final) != 0) {
        for (var in 1:ncol(data)) {
          if (var %in% col_final) {
            col_second[, var] <- 1
          } else {
            col_second[, var] <- 0
          }
        }
      }
    }
    col_out <- list(row = row_second, col = col_second)
    return(col_out)
  }
  
  backfitting <- function(data, row_fun_dec, col_fun_dec, layer_num, back_num, mu_bk,
                          row_eff, col_eff) {
    if (layer_num != 1) {
      data_back <- data - mean(data)
      for (var_bk_1 in 1:back_num) {
        for (var_bk_2 in 1:layer_num) {
          data_back_temp <- data_back
          cyc_num <- seq(1, layer_num, 1)
          cyc_drop_index <- which(cyc_num == var_bk_2)
          cyc_num <- cyc_num[-cyc_drop_index]
          for (var_bk_3 in cyc_num) {
            a_row <- as.matrix(row_fun_dec[, var_bk_3])
            a_col <- t(as.matrix(col_fun_dec[var_bk_3, ]))
            region <- a_row %*% a_col
            a_row <- as.matrix(row_eff[, var_bk_3])
            a_col <- t(as.matrix(col_eff[var_bk_3, ]))
            value <- a_row %*% a_col
            bi_value <- (value + mu_bk[var_bk_3]) * region
            data_back_temp <- data_back_temp - bi_value
          }
          row_back <- row_fun_dec[, var_bk_2]
          col_back <- col_fun_dec[var_bk_2, ]
          row_back <- matrix(rep(row_back, ncol(data)), ncol = ncol(data))
          col_back <- t(matrix(rep(col_back), nrow(data), nrow = ncol(data)))
          out_temp <- updata_theta(data_back_temp, row_fun = row_back, col_fun = col_back)
          mu <- out_temp$mu_out
          alpha_matrix <- out_temp$alpha_out
          beta_matrix <- out_temp$beta_out
          theta <- out_temp$theta_out
          mu_bk[var_bk_2] <- mu
          row_eff[, var_bk_2] <- alpha_matrix[, 1]
          col_eff[var_bk_2, ] <- beta_matrix[1, ]
        }
      }
    }
    bk_out <- list(mu = mu_bk, row = row_eff, col = col_eff)
    return(bk_out)
  }
  
  bc_filter <- function(row, col) {
    row_temp <- row
    col_temp <- col
    filter_row <- vector()
    post <- 1
    for (var in 1:ncol(row_temp)) {
      temp <- row_temp[, var]
      filter_temp <- which(temp == 1)
      if (length(filter_temp) == 0) {
        filter_row[post] <- var
        post <- post + 1
      }
    }
    filter_col <- vector()
    post <- 1
    for (var in 1:nrow(col_temp)) {
      temp <- col_temp[var, ]
      filter_temp <- which(temp == 1)
      if (length(filter_temp) == 0) {
        filter_col[post] <- var
        post <- post + 1
      }
    }
    filter <- c(filter_row, filter_col)
    filter <- unique(filter)
    if (length(filter) != 0) {
      row_temp <- row_temp[, -filter]
      col_temp <- col_temp[-filter, ]
    }
    out <- list(row_dec = row_temp, col_dec = col_temp)
    return(out)
  }
  
  iter.epoch <- 2 * max.layers
  
  a <- data_input
  a <- as.matrix(a)
  a <- apply(a, 2, as.numeric)
  mu_0 <- mean(a)
  data_temp <- a - mu_0
  Q <- sum(data_temp ^ 2)
  
  weight_input <- as.matrix(weight_input)
  weight_input <- apply(weight_input, 2, as.numeric)
  weight <- weight_input
  weight_temp <- weight_input
  

  mu_rec <- vector()
  CPS_rec <- vector()
  alpha_rec <- matrix(nrow = nrow(a), ncol = iter.epoch)
  beta_rec <- matrix(nrow = iter.epoch, ncol = ncol(a))
  theta_rec <- matrix(nrow = iter.epoch, ncol = ncol(a))
  row_dec_rec <- matrix(nrow = nrow(a), ncol = iter.epoch)
  col_dec_rec <- matrix(nrow = iter.epoch, ncol = ncol(a))
  
  
  

  layer_find <- 0
  for (bi_num in 1:iter.epoch) {
    if (layer_find < max.layers) {
      layer_log <- FALSE
      print(layer_find)
      center_min <- min(data_temp)
      center_diff <- max(data_temp) - min(data_temp)
      data_norm <- (data_temp - center_min) / center_diff
      
      center_min <- min(weight_temp)
      center_diff <- max(weight_temp) - min(weight_temp)
      weight_norm <- (weight_temp - center_min) / center_diff
      
      data_weight <- data_norm * (weight_norm ^ 3)
      center_mean <- mean(data_weight)
      center_diff <- max(data_weight) - min(data_weight)
      data_weight_stan <- (data_weight - center_mean) / center_diff
      
  
      row_dec <- initial_dec(data_weight_stan, 1, iter.startup)
      col_dec <- initial_dec(data_weight_stan, 2, iter.startup)
      
      for (layer_cycle in 1:iter.layer) {
        out_temp <- updata_theta(data_weight_stan, row_fun = row_dec, col_fun = col_dec)
        mu <- out_temp$mu_out
        alpha_matrix <- out_temp$alpha_out
        beta_matrix <- out_temp$beta_out
        theta <- out_temp$theta_out
        row_dec_test <- updata_dec(data_weight_stan, theta_fun = theta, 
                                   col_fun = col_dec, layer_cycle = layer_cycle,
                                   iter.layer = iter.layer, iter.bin = iter.bin, side = 1)
        col_dec_test <- updata_dec(data_weight_stan, theta_fun = theta,
                                   row_fun = row_dec, layer_cycle = layer_cycle,
                                   iter.layer = iter.layer, iter.bin = iter.bin, side = 2)
        row_dec <- matrix(rep(row_dec_test, ncol(a)), ncol = ncol(a))
        col_dec <- t(matrix(rep(col_dec_test), nrow(a), nrow = ncol(a)))
      }
      
      drop_row_dec <- drop_release_row(data = data_weight_stan,
                                       row_fun = row_dec, col_fun = col_dec,
                                       iter.startup)
      if (length(drop_row_dec) != 0) {
        row_dec <- drop_row_dec$row
        col_dec <- drop_row_dec$col
      }
      drop_col_dec <- drop_release_col(data = data_weight_stan,
                                       row_fun = row_dec, col_fun = col_dec,
                                       iter.startup)
      if (length(drop_col_dec) != 0) {
        row_dec <- drop_col_dec$row
        col_dec <- drop_col_dec$col
      }
      out_temp <- updata_theta(data_temp, row_fun = row_dec, col_fun = col_dec)
      mu <- out_temp$mu_out
      alpha_matrix <- out_temp$alpha_out
      beta_matrix <- out_temp$beta_out
      theta <- out_temp$theta_out
      CPS <- sum(theta * row_dec * col_dec)
      print(CPS)
      
      if (CPS > 0) {
        CPS_rec[bi_num] <- CPS
        row_dec_rec[, bi_num] <- row_dec[, 1]
        col_dec_rec[bi_num, ] <- col_dec[1, ]
        mu_rec[bi_num] <- mu
        alpha_rec[, bi_num] <- alpha_matrix[, 1]
        beta_rec[bi_num, ] <- beta_matrix[1, ]
        
        aaa_test <- backfitting(data = a, row_fun_dec = row_dec_rec, col_fun_dec = col_dec_rec,
                                layer_num = bi_num, back_num = 1, mu_bk = mu_rec,
                                row_eff = alpha_rec, col_eff = beta_rec)
        mu_rec <<- aaa_test$mu
        alpha_rec <<- aaa_test$row
        beta_rec <<- aaa_test$col
        
        for (var_cancha in 1:bi_num) {
          mu <- mu_rec[var_cancha]
          alpha_matrix <- matrix(rep(alpha_rec[, var_cancha], ncol(a)), ncol = ncol(a))
          beta_matrix <- t(matrix(rep(beta_rec[var_cancha, ], nrow(a)), nrow = ncol(a)))
          theta <- mu + alpha_matrix + beta_matrix
          row_dec_temp <- row_dec_rec[, var_cancha]
          col_dec_temp <- col_dec_rec[var_cancha, ]
          row_dec_temp <- matrix(rep(row_dec_temp, ncol(a)), ncol = ncol(a))
          col_dec_temp <- t(matrix(rep(col_dec_temp), nrow(a), nrow = ncol(a)))
          if (var_cancha == 1) {
            data_temp <- data_temp - theta * row_dec_temp * col_dec_temp
            row_index <- which(row_dec_temp[, 1] == 1)
            col_index <- which(col_dec_temp[1, ] == 1)
            weight_min <- weight_temp[row_index, col_index]
            weight_temp[row_index, col_index] <- weight_temp[row_index, col_index] - mean(weight_min)
          } else {
            data_temp <- data_temp - theta * row_dec_temp * col_dec_temp
            row_index <- which(row_dec_temp[, 1] == 1)
            col_index <- which(col_dec_temp[1, ] == 1)
            weight_min <- weight_temp[row_index, col_index]
            weight_temp[row_index, col_index] <- weight_temp[row_index, col_index] - mean(weight_min)
          }
        }
      } else {
        break()
      }
      if (layer_log == FALSE) {
        layer_find <- layer_find + 1
      }
      
    }
  }
  
  result <- bc_filter(row_dec_rec, col_dec_rec)
  row_dec_rec <- result$row_dec
  col_dec_rec <- result$col_dec
  
  setClass("Biclust", slots = list(RowxNumber = "matrix", NumberxCol = "matrix",
                                   Number = "numeric"))
  bicluster_num <- length(CPS_rec)
  if (bicluster_num != 0) {
    if (bicluster_num == 1) {
      row_dec_rec_log <- (row_dec_rec == 1)
      col_dec_rec_log <- (col_dec_rec == 1)
    }
    if (bicluster_num > 1) {
      row_dec_rec_log <- apply(row_dec_rec, 2, as.logical)
      col_dec_rec_log <- apply(col_dec_rec, 2, as.logical)
    }
    row_dec_rec_log <- as.matrix(as.data.frame(row_dec_rec_log))
    col_dec_rec_log <- as.matrix(as.data.frame(col_dec_rec_log))
    bc <- new("Biclust", RowxNumber = row_dec_rec_log, NumberxCol = col_dec_rec_log,
               Number = bicluster_num)
    cat("The number of patterns found:", bicluster_num, "\n")
    alpha_rec <<- alpha_rec[, 1:bicluster_num]
    beta_rec <<- beta_rec[1:bicluster_num, ]
    CPS_rec <<- CPS_rec
  }
  
  if (bicluster_num == 0) {
    alpha_rec <<- 0
    beta_rec <<- 0
    mu_rec <<- 0
    CPS_rec <<- 0
    print("No pattern was found!!")
  }
  return(bc)
}
