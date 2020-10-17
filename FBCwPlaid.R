FBCwPlaid <- function(FPKM.IP, FPKM.input, Methylation.level, Expression.level, 
                      optimization, GENES.CORRES.SITES, GENE.ID.TYPES, exponent.num, 
                      kmeans.startup, max.layers, iter.layer, iter.bin, backfitting.num, 
                      verbose) {
  
  # Detect whether the data is missing. 
  if ((missing(FPKM.IP) || missing(FPKM.input))) {
    if ((missing(Methylation.level) || missing(Expression.level))) {
      stop("Missing input data!")
    } else {
      if (nrow(Methylation.level) != nrow(Expression.level) || 
          ncol(Methylation.level) != ncol(Expression.level)) {
        stop("Please check whether the dimensions of the two sets of input data are the same!")
      }
    }
  }
  
  # Methylation level matrix and expression level matrix are generated 
  # by using IP and input samples. 
  if ((missing(Methylation.level) || missing(Expression.level))) {
    if (nrow(FPKM.IP) != nrow(FPKM.input) || ncol(FPKM.IP) != ncol(FPKM.input)) {
      stop("Please check whether the dimensions of the two sets of input data are the same!")
    }
    FPKM.IP = apply(FPKM.IP, 2, as.numeric)
    FPKM.input = apply(FPKM.input, 2, as.numeric)
    FPKM.IP <- as.matrix(FPKM.IP) + 0.0001
    FPKM.input <- as.matrix(FPKM.input) + 0.0001
    data.sum <- FPKM.IP + FPKM.input
    Methylation.level <- FPKM.IP / data.sum
    Expression.level <- log2(data.sum + 1 - 0.0002)
  }
  
  # Set default parameters for optimization and verbose. 
  if (missing(optimization)) {
    optimization <- FALSE
  }
  if (missing(verbose)) {
    verbose <- FALSE
  }
  
  # Whether to open the enrichment constraint module. 
  if (optimization == FALSE) {
    if (missing(exponent.num)) {
      power.num.min <- 1
      power.num.max <- 1
    } else {
      power.num.min <- exponent.num
      power.num.max <- exponent.num
    }
    
  } else {
    if (missing(GENES.CORRES.SITES) || missing(GENE.ID.TYPES)) {
      stop("Lack of information on genes corresponding to the sites!")
    }
    power.num.min <- 1
    power.num.max <- 5
    # Import dependent packages. 
    library(clusterProfiler)
    library(org.Hs.eg.db)
  }
  
  
  # Initialization of row and column indicators. 
  initial_dec <- function(data, side, kmeans.startup) {
    # Initialization of row indicators. 
    if (side == 1) {
      row_test <- try(kmeans(data, 2, kmeans.startup), TRUE)
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
    # Initialization of column indicators. 
    if (side == 2) {
      col_test <- try(kmeans(t(data), 2, kmeans.startup), TRUE)
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
  
  
  # Update theta. 
  updata_theta <- function(data, row_fun, col_fun) {
    # Start iterative update. 
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
      # Updata alpha and beta. 
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
    # Generate theta, alpha and beta matrices. 
    alpha_matrix <- matrix(rep(alpha, ncol(a)), ncol = ncol(a))
    beta_matrix <- t(matrix(rep(beta, nrow(a)), nrow = ncol(a)))
    theta <- mu_matrix + alpha_matrix + beta_matrix
    out <- list(mu_out = mu, alpha_out = alpha_matrix, beta_out = beta_matrix, theta_out = theta)
    return(out)
  }
  
  
  # Update row and column indicators. 
  updata_dec <- function(data, theta_fun, row_fun, col_fun, layer_cycle, 
                         iter.layer, iter.bin, side) {
    # Update row indicators (rou).
    if (side == 1) {
      row_dec_test <- vector()
      for (var in 1:nrow(data)) {
        if (sum((theta_fun[var, ] ^ 2) * (col_fun[var, ] ^ 2)) == 0) {
          row_dec_test[var] <- 0
        } else {
          row_dec_test[var] <- sum(theta_fun[var, ] * col_fun[var, ] * data[var, ]) /
            sum((theta_fun[var, ] ^ 2) * (col_fun[var, ] ^ 2))
        }
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
    # Update column indicators (kapa). 
    if (side == 2) {
      col_dec_test <- vector()
      for (var in 1:ncol(data)) {
        if (sum((theta_fun[, var] ^ 2) * (row_fun[, var] ^ 2)) == 0) {
          col_dec_test[var] <- 0
        } else {
          col_dec_test[var] <- sum(theta_fun[, var] * row_fun[, var] * data[, var]) /
            sum((theta_fun[, var] ^ 2) * (row_fun[, var] ^ 2))
        }
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
  
  
  # Calculate intra-cluster variance. 
  SDwC_cal <- function(data) {
    SDwC_temp <- sum((data - mean(data)) ^ 2) / (nrow(data) * ncol(data))
    return(SDwC_temp)
  }
  
  
  # Secondary contraction of rows using k-means. 
  drop_release_row <- function(data, row_fun, col_fun, kmeans.startup) {
    row_second <- row_fun
    col_second <- col_fun
    row_index <- which(row_fun[, 1] == 1)
    col_index <- which(col_fun[1, ] == 1)
    data_min <- matrix(nrow = length(row_index), ncol = length(col_index))
    data_min <- as.matrix(data[row_index, col_index])
    SDwC_original <- SDwC_cal(data_min)
    mean_original <- mean(data_min)
    if (nrow(data_min) > 2) {
      row_test <- try(kmeans(data_min, 2, kmeans.startup), TRUE)
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
      
      length.small <- min(length(row_one), length(row_two))
      if (length.small / length(row_index) >= 0.15) {
        if (mean_one >= mean_original && SDwC_one <= 1.5 * SDwC_original && 
            mean_one >= 1.1 * mean_two) {
          row_final <- row_one
          SDwC_min <- SDwC_one
        }
        if (mean_two >= mean_original && SDwC_two <= 1.5 * SDwC_original && 
            mean_two >= 1.1 * mean_one) {
          row_final <- row_two
          SDwC_min <- SDwC_two
        }
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
  
  
  # Secondary contraction of columns using k-means.
  drop_release_col <- function(data, row_fun, col_fun, kmeans.startup) {
    row_second <- row_fun
    col_second <- col_fun
    row_index <- which(row_fun[, 1] == 1)
    col_index <- which(col_fun[1, ] == 1)
    data_min <- matrix(nrow = length(row_index), ncol = length(col_index))
    data_min <- as.matrix(data[row_index, col_index])
    SDwC_original <- SDwC_cal(data_min)
    mean_original <- mean(data_min)
    if (nrow(t(data_min)) > 2) {
      col_test <- try(kmeans(t(data_min), 2, kmeans.startup), TRUE)
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
        if ((mean_one > mean_original && SDwC_one < 1.3 * SDwC_original) || mean_two < 0) {
          col_final <- col_one
          SDwC_min <- SDwC_one
        }
        if ((mean_two > mean_original && SDwC_two < 1.3 * SDwC_original) || mean_one < 0) {
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
  
  
  # Backfitting. 
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
  
  
  # Filter the final result. 
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
  
  
  # Main part. 
  iter.epoch <- 3 * max.layers
  WE_score_rec <- vector()
  # If optimization is TRUE, start enrichment constraint. 
  for (var.optimization in power.num.min:power.num.max) {
    if (optimization) {
      cat("The value of the current exponential power:", var.optimization, "\n")
    }
    
    # Data recording and preprocessing. 
    a <- Methylation.level
    a <- as.matrix(a)
    a <- apply(a, 2, as.numeric)
    mu_0 <- mean(a)
    
    if (!optimization) {
      mu.all <<- mu_0
    } else {
      mu.all <- mu_0
    }
    
    data_temp <- a - mu_0
    Q <- sum(data_temp ^ 2)
    weight <- Expression.level
    weight_temp <- Expression.level
    
    # Create a location for the storage of parameters. 
    mu.rec <- vector()
    CPS.rec <- vector()
    alpha.rec <- matrix(nrow = nrow(a), ncol = iter.epoch)
    beta.rec <- matrix(nrow = iter.epoch, ncol = ncol(a))
    theta_rec <- matrix(nrow = iter.epoch, ncol = ncol(a))
    row_dec_rec <- matrix(nrow = nrow(a), ncol = iter.epoch)
    col_dec_rec <- matrix(nrow = iter.epoch, ncol = ncol(a))
    
    # Parameters optimization process. 
    layer_find <- 0
    for (bi_num in 1:iter.epoch) {
      if (layer_find < max.layers) {
        layer_log <- FALSE
        if (!optimization && verbose) {
          cat("Number of patterns found:", layer_find, "\t\t")
        }
        center_min <- min(data_temp)
        center_diff <- max(data_temp) - min(data_temp)
        data_norm <- (data_temp - center_min) / center_diff
        
        center_min <- min(weight_temp)
        center_diff <- max(weight_temp) - min(weight_temp)
        weight_norm <- (weight_temp - center_min) / center_diff
        
        data_weight <- data_norm * (weight_norm ^ var.optimization)
        
        center_mean <- mean(data_weight)
        center_diff <- max(data_weight) - min(data_weight)
        data_weight_stan <- (data_weight - center_mean) / center_diff
        
        if (bi_num != 1) {
          data_weight_stan_mean <- mean(data_weight_stan)
          data_weight_stan_sd <- sqrt(SDwC_cal(data_weight_stan))
          row.dec.index <- vector()
          col.dec.index <- vector()
          for (random.selected in 1:(bi_num - 1)) {
            row_dec <- row_dec_rec[, random.selected]
            col_dec <- col_dec_rec[random.selected, ]
            row.dec.index <- which(row_dec == 1)
            col.dec.index <- which(col_dec == 1)
            if (length(row.dec.index) > 0 && length(col.dec.index) > 0) {
              selected.num <- length(row.dec.index) * length(col.dec.index)
              data_weight_stan[row.dec.index, col.dec.index] <- rnorm(selected.num, 
                                                                      data_weight_stan_mean,
                                                                      data_weight_stan_sd)
            }
          }
        }
        # drawHeatmap(data_weight_stan)
        
        # Initialize the row and column indicators. 
        row_dec <- initial_dec(data_weight_stan, 1, kmeans.startup)
        col_dec <- initial_dec(data_weight_stan, 2, kmeans.startup)
        # Start iterative update. 
        for (layer_cycle in 1:iter.layer) {
          # Update theta, alpha and beta. 
          out_temp <- updata_theta(data_weight_stan, row_fun = row_dec, col_fun = col_dec)
          mu <- out_temp$mu_out
          alpha_matrix <- out_temp$alpha_out
          beta_matrix <- out_temp$beta_out
          theta <- out_temp$theta_out
          # Update row and column indicators (rou and kapa). 
          row_dec_test <- updata_dec(data_weight_stan, theta_fun = theta, 
                                     col_fun = col_dec, layer_cycle = layer_cycle,
                                     iter.layer = iter.layer, iter.bin = iter.bin, side = 1)
          col_dec_test <- updata_dec(data_weight_stan, theta_fun = theta,
                                     row_fun = row_dec, layer_cycle = layer_cycle,
                                     iter.layer = iter.layer, iter.bin = iter.bin, side = 2)
          row_dec <- matrix(rep(row_dec_test, ncol(a)), ncol = ncol(a))
          col_dec <- t(matrix(rep(col_dec_test), nrow(a), nrow = ncol(a)))
        }
        # drawHeatmap(data_weight_stan)
        # drawHeatmap(row_dec * col_dec)
        # Use k-means to shrink rows and columns. 
        not.drop.row <- FALSE
        drop_row_dec <- drop_release_row(data = data_weight_stan,
                                         row_fun = row_dec, col_fun = col_dec,
                                         kmeans.startup)
        if (length(drop_row_dec) != 0) {
          row_dec_drop <- drop_row_dec$row
        }
        if (sum(row_dec[, 1]) == sum(row_dec_drop[, 1])) {
          not.drop.row <- TRUE
        } else {
          row_dec <- row_dec_drop
        }
        drop_col_dec <- drop_release_col(data = data_weight_stan, row_fun = row_dec,
                                         col_fun = col_dec, kmeans.startup)
        if (length(drop_col_dec) != 0) {
          col_dec <- drop_col_dec$col
        }
        if (not.drop.row) {
          drop_row_dec <- drop_release_row(data = data_weight_stan,
                                           row_fun = row_dec, col_fun = col_dec,
                                           kmeans.startup)
          if (length(drop_row_dec) != 0) {
            row_dec <- drop_row_dec$row
          }
        }
        
        # Update one more round of theta, alpha and beta. 
        out_temp <- updata_theta(data_temp, row_fun = row_dec, col_fun = col_dec)
        mu <- out_temp$mu_out
        alpha_matrix <- out_temp$alpha_out
        beta_matrix <- out_temp$beta_out
        theta <- out_temp$theta_out
        CPS <- sum(theta * row_dec * col_dec)
        if (!optimization && verbose) {
          cat("CPS:", CPS, "\n")
        }
        
        if (CPS > 0) {
          CPS.rec[bi_num] <- CPS
          row_dec_rec[, bi_num] <- row_dec[, 1]
          col_dec_rec[bi_num, ] <- col_dec[1, ]
          mu.rec[bi_num] <- mu
          alpha.rec[, bi_num] <- alpha_matrix[, 1]
          beta.rec[bi_num, ] <- beta_matrix[1, ]
          
          back.test <- backfitting(data = a, row_fun_dec = row_dec_rec, 
                                   col_fun_dec = col_dec_rec, layer_num = bi_num, 
                                   back_num = backfitting.num, mu_bk = mu.rec,
                                   row_eff = alpha.rec, col_eff = beta.rec)
          if (!optimization) {
            mu.rec <<- back.test$mu
            alpha.rec <<- back.test$row
            beta.rec <<- back.test$col
          } else {
            mu.rec <- back.test$mu
            alpha.rec <- back.test$row
            beta.rec <- back.test$col
          }
          # After the backfitting is complete, update data_temp. 
          for (var_cancha in 1:bi_num) {
            mu <- mu.rec[var_cancha]
            alpha_matrix <- matrix(rep(alpha.rec[, var_cancha], ncol(a)), ncol = ncol(a))
            beta_matrix <- t(matrix(rep(beta.rec[var_cancha, ], nrow(a)), nrow = ncol(a)))
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
          CPS <- 0
          mu <- 0
          mu.rec[bi_num] <- mu
          alpha.rec[, bi_num] <- 0
          beta.rec[bi_num, ] <- 0
          row_dec_rec[, bi_num] <- 0
          col_dec_rec[bi_num, ] <- 0
          CPS.rec[bi_num] <- CPS
          data_temp <- data_temp
          weight_temp <- weight_temp
        }
        if (layer_log == FALSE) {
          layer_find <- layer_find + 1
        }
      }
    }
    
    result <- bc_filter(row_dec_rec, col_dec_rec)
    row_dec_rec <- result$row_dec
    col_dec_rec <- result$col_dec
    
    setClass("biclust", slots = list(RowxNumber = "matrix", NumberxCol = "matrix",
                                     Number = "numeric"))
    bicluster_num <- sum(CPS.rec != 0)
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
      pattern.result <- new("biclust", RowxNumber = row_dec_rec_log, NumberxCol = col_dec_rec_log,
                Number = bicluster_num)
      if (!optimization) {
        if (verbose) {
          cat("The number of patterns found:", bicluster_num, "\n")
        }
        alpha.rec <<- alpha.rec[, 1:bicluster_num]
        beta.rec <<- beta.rec[1:bicluster_num, ]
        CPS.rec <<- CPS.rec
      } else {
        alpha.rec <- alpha.rec[, 1:bicluster_num]
        beta.rec <- beta.rec[1:bicluster_num, ]
        CPS.rec <- CPS.rec
      }
      
      
    }
    if (bicluster_num == 0) {
      if (!optimization) {
        alpha.rec <<- 0
        beta.rec <<- 0
        mu.rec <<- 0
        CPS.rec <<- 0
        if (verbose) {
          cat("No pattern was found!!", "\n")
        }
      } else {
        alpha.rec <- 0
        beta.rec <- 0
        mu.rec <- 0
        CPS.rec <- 0
      }
    }
    
    if (optimization) {
      # Extract the corresponding genes in each pattern. 
      RowxNumber <- pattern.result@RowxNumber
      pattern.num <- as.numeric(pattern.result@Number)
      we.power <- vector()
      gene.obtained <- data.frame()
      for (pattern.temp in 1:pattern.num) {
        sites.index <- which(RowxNumber[, pattern.temp] == TRUE)
        gene.obtained[1:length(sites.index), pattern.temp] <- GENES.CORRES.SITES[sites.index]
      }
      
      cat("Number of patterns found:", ncol(gene.obtained), "\n")
      for (var.enrich in 1:ncol(gene.obtained)) {
        cat("Current pattern of GO analysis:", var.enrich, "\n")
        gene.temp <- gene.obtained[, var.enrich]
        gene.temp <- na.omit(gene.temp) 
        gene.temp <- as.character(as.vector(gene.temp))
        # Gene ontology analysis. 
        go <- enrichGO(gene.temp, "org.Hs.eg.db", ont = "ALL", 
                       pAdjustMethod = "BH", pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.2, keyType = GENE.ID.TYPES, readable = TRUE)
        if (nrow(go@result) == 0) {
          we.power[var.enrich] <- 0
        } else {
          go.result <- matrix(nrow = nrow(go@result), ncol = 2)
          go.result[, 1] <- go@result[, 10]
          go.result[, 2] <- as.numeric(go@result[, 6])
          temp.count <- go.result[1, 1]
          temp.ratio <- go@result[1, 4]
          temp.ratio <- eval(parse(text = temp.ratio))
          not.included.num <- length(gene.temp) - (temp.count / temp.ratio)
          go.result[, 1] <- go.result[, 1] / length(gene.temp)
          go.result[, 2] <- -log(go.result[, 2], 10)
          denominator <- (apply(go.result, 2, sum))[1]
          denominator <- denominator + (not.included.num / length (gene.temp))
          molecule <- sum(go.result[, 1] * go.result[, 2])
          we.power[var.enrich] <- molecule / denominator
        }
      }
      WE_score_rec[var.optimization] <- mean(we.power)
    }
  }
  if (optimization) {
    exponent <<- which(WE_score_rec == max(WE_score_rec))
    cat("The optimal choice of exponential power is:", exponent, "\n")
    return(exponent)
  } else {
    return(pattern.result)
  }
}
