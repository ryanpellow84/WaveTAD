#### Files ####

input_args <- commandArgs(TRUE)
left_coverage_file <- input_args[1]
right_coverage_file <- input_args[2]
results_file <- input_args[3]
topdom_files <- c(input_args[4],
                  input_args[5],
                  input_args[6],
                  input_args[7])
loops_files <- c(input_args[8],
                 input_args[9],
                 input_args[10],
                 input_args[11])
chr <- input_args[12]
genome_size <- input_args[13]

#### Libraries####
library(ggplot2)
library(scales)
library(wavelets)
library(matrixStats)
library(reshape2)
library(IRanges)
library(plyr)
library(dplyr)
library(stringi)
library(iotools)
options(scipen = 10)

#### WUBBA Pairing Function ####
wubba_pairing <- function(left_cov,
                          right_cov,
                          chromosome = "X",
                          min_pair_distance = 40000,
                          max_pair_distance = 2000000,
                          wubba_alpha = 0.05, 
                          wubba_from_scale = 6,
                          wubba_to_scale = -2,
                          wubba_filt = "c6",
                          wubba_p_method = "holm",
                          resolution_size = 5000,
                          topdom_input_files,
                          loops_input_files,
                          bin_vec = c(0, 7500, 17500, 100000),
                          resolution_vec = c(5000, 10000, 25000),
                          error_val = 4000){
  #### Functions ####
  # comp_wt Function 
  comp_wt <- function(input_data, 
                      alpha = 0.05, 
                      from_scale = 5,
                      to_scale = "max",
                      prox = FALSE,
                      filt = "c6",
                      p_method = "holm"){
    
    #### normalizes data ####
    norm_data <- function(d){
      n_d <- (d - mean(d)) / sd(d)
      return(n_d)
    }
    #### returns p-values ####
    probs <- function(w1,w2){
      sample_mean <- mean(w1)
      sd_of_previous_scale <- sd(w2)
      z_scores <- (w1 - sample_mean)/(sd_of_previous_scale)
      # p_vals <- 2*pnorm(-abs(z_scores))
      p_vals <- 1 - pnorm(z_scores)
      return(p_vals)
      
    }
    #### complete p-value function ####
    pcom <- function(w, locs){
      p <- data.frame(stringsAsFactors = FALSE)
      for (i in 1:length(w)){
        if (i == 1){
          j <- 1
        } else{
          j <- i-1
        }
        get_pvals <- probs(w[[i]], w[[j]])
        if (p_method == "bonferroni"){
          correction_value <- 2^(i-1)
          adj_pvals <- get_pvals * (length(get_pvals)/correction_value)
        }
        
        if (p_method == "holm"){
          correction_value <- ceiling(length(get_pvals)/(2^(i-1)))
          correction_vec <- rep(1:correction_value, each = 2^(i-1))
          correction_vec <- correction_vec[1:length(get_pvals)]
          correction_vec <- (correction_value + 1) - correction_vec
          holm_df <- data.frame(Index = 1:length(get_pvals), 
                                Pvals = get_pvals)
          holm_df <- holm_df[order(holm_df$Pvals),]
          holm_df <- data.frame(Index = holm_df$Index, 
                                Pvals = holm_df$Pvals,
                                Adj = holm_df$Pvals * correction_vec)
          holm_df <- holm_df[order(holm_df$Index),]
          adj_pvals <- holm_df$Adj
        }
        
        print(paste("Hits at scale ", 
                    i, 
                    ": ", 
                    sum(adj_pvals <= .05),
                    sep = ""))
        temp_df <- data.frame(s = rep(2^(length(w)-i+2),length(locs)),
                              l = locs,
                              p = adj_pvals)
        temp_df = temp_df[temp_df$p <= 0.05,]
        p <- rbind(p, temp_df)
      }
      return(p)
    }
    
    #### creates list of data sets for each scale for prox function ####
    prox_list <- function(sig_loc){
      prox_tab_1 <- sig_loc
      prox_tab_2 <- data.frame(scale = double(), oligo = double(), pval = double())
      prox_tab_2 <- prox_tab_1[prox_tab_1$pval <= alpha,]
      prox_scales_1 <- c()
      for (i in 1:length(w)){
        prox_scales_1[i] <- 2^(i+1)
      }
      prox_scales_1 <- rev(prox_scales_1)
      n <- length(prox_scales_1)
      scales_picked <- from_scale:to_scale
      prox_tab_final <- 
        data.frame(scale = double(), oligo = double(), pval = double())
      for (i in scales_picked){
        prox_tab <- prox_tab_2[prox_tab_2$scale == prox_scales_1[i],]
        prox_tab_final <- rbind(prox_tab_final, prox_tab)
      }
      prox_scales_2 <- as.numeric(levels(factor(prox_tab_final$scale)))
      prox_list <- c()
      if (length(prox_scales_2) != 0){
        for (i in 1:length(prox_scales_2)){
          prox_list[[i]] <- 
            (prox_tab_final[prox_tab_final$scale == prox_scales_2[i],2:3])
        }
        names(prox_list) <- levels(factor(prox_tab_final$scale)) 
      }
      return(prox_list)
    }
    
    #### determines proximate location of significant hits ####
    prox_func <- function(l){
      sig_list <- c()
      if (!is.null(l)){
        for (i in 1:length(l)){
          x1 <- l[[i]]
          sig_table <- c()
          j <- 1
          while (j <= length(x1$oligo)){
            start_val <- x1$oligo[j]
            end_val <- NA
            while(is.na(end_val) & (j+1 <= length(x1$oligo))){
              if (x1$oligo[j]+1 != x1$oligo[j+1]){
                end_val <- x1$oligo[j]
              } else{
                end_val <- NA
                j <- j + 1
              }
            }
            if (j == length(x1$oligo)){
              end_val <- x1$oligo[j]
            }
            sig_table <- rbind(sig_table, cbind(start_val, end_val, x1$pval[j]))
            colnames(sig_table) <- c("start", "end", "pval")
            j <- j + 1
          }
          sig_list[[i]] <- sig_table
        }
        for (k in 1:length(sig_list)){
          avg <- round((sig_list[[k]][,1] + sig_list[[k]][,2])/2)
          sig_list[[k]] <- cbind(sig_list[[k]], avg = avg)
        }
        names(sig_list) <- names(l)
      }
      return(sig_list)
    }
    
    #### complete function ####
    input_data$Value <- norm_data(input_data$Value)
    samp_dwt <- align(modwt(input_data$Value, filter = filt))
    samp_dwt <- align(modwt(samp_dwt@V[[5]], filter = filt)) # Denoised data
    if(to_scale == "max"){
      to_scale  <- length(samp_dwt@W)
    } else if(to_scale < 0){
      to_scale  <- length(samp_dwt@W) + to_scale
    }
    w <- samp_dwt@W[]
    x <- input_data$Location
    df <- pcom(w, x)
    sig_loc <- df
    colnames(sig_loc) <- c("scale", "oligo", "pval")
    if (prox){
      print("Start prox_list function")
      print(paste("System Time:", Sys.time()))
      sig_loc <- prox_list(sig_loc)
      print("Start prox function")
      print(paste("System Time:", Sys.time()))
      sig_loc <- prox_func(sig_loc)
    } 
    return(sig_loc)
  }
  
  # WUBBA Loop Confirmation Function
  loop_confirm <- function(pairs_input,
                           loops_input_files,
                           chromosome,
                           resolution_vec){
    loop_domains <- data.frame(stringsAsFactors = FALSE)
    for (i in 1:length(loops_input_files)){
      for (j in 1:length(resolution_vec)) {
        if (j >= i){
          loop_domains_temp <- read.delim(loops_input_files[i],
                                          header = FALSE,
                                          col.names = c("Chromosome1", "Start1", "End1", "Chromosome2", "Start2", "End2", "AncPval"),
                                          stringsAsFactors = FALSE)
          loop_domains_temp <- loop_domains_temp[loop_domains_temp$Chromosome1 == loop_domains_temp$Chromosome2,]
          loop_domains_temp <- loop_domains_temp[loop_domains_temp$Chromosome1 == chromosome,]
          loop_domains_temp$Start1 <- round_any(loop_domains_temp$Start1, resolution_vec[j], f = floor)
          loop_domains_temp$End1 <- round_any(loop_domains_temp$End1, resolution_vec[j], f = ceiling)
          loop_domains_temp$Start2 <- round_any(loop_domains_temp$Start2, resolution_vec[j], f = floor)
          loop_domains_temp$End2 <- round_any(loop_domains_temp$End2, resolution_vec[j], f = ceiling)
          possible_boundaries <- data.frame(Chromosome1 = rep(chromosome, nrow(pairs_input)),
                                            Start1 = round_any(pairs_input$Start, resolution_vec[j], f = floor),
                                            End1 = round_any(pairs_input$Start, resolution_vec[j], f = ceiling),
                                            Chromosome2 = rep(chromosome, nrow(pairs_input)),
                                            Start2 = round_any(pairs_input$End, resolution_vec[j], f = floor),
                                            End2 = round_any(pairs_input$End, resolution_vec[j], f = ceiling),
                                            StartSave = pairs_input$Start,
                                            EndSave = pairs_input$End,
                                            Pval1 = pairs_input$Pval1,
                                            Pval2 = pairs_input$Pval2,
                                            Pval = pairs_input$Pval,
                                            stringsAsFactors = FALSE)
          boundary_hits <- full_join(possible_boundaries, 
                                     loop_domains_temp, 
                                     by = c("Chromosome1", "Start1", "End1", "Chromosome2", "Start2", "End2"))
          boundary_hits <- boundary_hits[!is.na(boundary_hits$AncPval),]
          boundary_hits <- boundary_hits[!is.na(boundary_hits$Pval),]
          loop_domains <- rbind(loop_domains, boundary_hits)
        }
      }
    }
    
    
    boundary_keeps <- data.frame(Start = loop_domains$StartSave,
                                 End = loop_domains$EndSave,
                                 Pval1 = loop_domains$Pval1,
                                 Pval2 = loop_domains$Pval2,
                                 Pval = loop_domains$Pval,
                                 AncPval = loop_domains$AncPval,
                                 stringsAsFactors = FALSE)
    boundary_keeps <- arrange(boundary_keeps, Start, End, AncPval)
    
    # Get unique domains with lowest p-value
    rownames(boundary_keeps) <- NULL
    unique_starts <- unique(boundary_keeps[,1:2])
    unique_index <- as.numeric(rownames(unique_starts))
    boundary_keeps <- boundary_keeps[unique_index,]
    
    return(boundary_keeps)
  }
  # WUBBA TopDom Pre-Filter
  topdom_collab_pre <- function(anchor_pairs_input,
                                chromosome,
                                topdom_input_files,
                                bin_vec = c(0, 7500, 17500, 100000),
                                resolution_vec = c(5000, 10000, 25000),
                                error_val = 4000){
    #### Functions ####
    # Takes two vectors and returns the value in the reference vec closest to the nth subject
    get_closest <- function(subject_vec, reference_vec){
      abs_diffs <- abs(outer(subject_vec, reference_vec, "-"))
      min_index <- apply(abs_diffs, 1, which.min)
      closest_pair <- reference_vec[min_index]
    }
    
    #### Prep Anchor Pairs ####
    # Add Chromosome Feature
    anchor_pairs_input$Chromosome <- rep(chromosome, nrow(anchor_pairs_input))
    # Add Size Feature
    anchor_pairs_input$Size <- anchor_pairs_input$End - anchor_pairs_input$Start
    # Add Resolution Feature
    anchor_pairs_input$Resolution <- round(anchor_pairs_input$Size / 40)
    # Group TAD domains by resolution
    anchor_pairs_input$Resolution <-  cut(anchor_pairs_input$Resolution, bin_vec)
    # Rename Groups
    levels(anchor_pairs_input$Resolution) <- resolution_vec 
    # Duplicate Start Feature (Important Later)
    anchor_pairs_input$StartSave <- anchor_pairs_input$Start
    # Duplicate End Feature (Important Later)
    anchor_pairs_input$EndSave <- anchor_pairs_input$End
    #### Prep TopDom Files ####
    # Load TopDom files
    topdom_domains <- data.frame(stringsAsFactors = FALSE)
    topdom_transfer_big_domains <- data.frame(stringsAsFactors = FALSE)
    topdom_transfer_small_domains <- data.frame(stringsAsFactors = FALSE)
    for (i in 1:length(topdom_files)){
      topdom_domains_temp <- read.delim(topdom_input_files[i],
                                        header = FALSE,
                                        stringsAsFactors = FALSE)
      topdom_domains_temp <- topdom_domains_temp[,1:3]
      names(topdom_domains_temp) <- c("Chromosome", "Start", "End")
      topdom_domains_temp <- topdom_domains_temp[topdom_domains_temp$Chromosome == chromosome,]
      # Add Size and Resolution Features
      topdom_domains_temp$Size <- topdom_domains_temp$End - topdom_domains_temp$Start
      topdom_domains_temp$Resolution <- topdom_domains_temp$Size / 40
      # Group TAD domains by resolution
      topdom_domains_temp$Resolution <-  cut(topdom_domains_temp$Resolution, bin_vec)
      # Rename Groups
      levels(topdom_domains_temp$Resolution) <- resolution_vec 
      topdom_domains_temp$Resolution <- as.numeric(levels(topdom_domains_temp$Resolution))[topdom_domains_temp$Resolution]
      # Add to transfer df
      topdom_transfer_big_domains <- rbind(topdom_transfer_big_domains, 
                                           topdom_domains_temp[topdom_domains_temp$Resolution > resolution_vec[i],])
      topdom_transfer_small_domains <- rbind(topdom_transfer_small_domains, 
                                             topdom_domains_temp[topdom_domains_temp$Resolution < resolution_vec[i],])
      topdom_domains_temp <- rbind(topdom_domains_temp, 
                                   topdom_transfer_big_domains[topdom_transfer_big_domains$Resolution == resolution_vec[i],])
      topdom_domains_temp <- rbind(topdom_domains_temp, 
                                   topdom_transfer_small_domains[topdom_transfer_small_domains$Resolution == resolution_vec[i],])
      # Convert to boundaries
      topdom_domains_temp <- data.frame(Chromosome = c(topdom_domains_temp$Chromosome, topdom_domains_temp$Chromosome),
                                        Boundary = c(topdom_domains_temp$Start, topdom_domains_temp$End))
      # Remove Duplicates
      topdom_domains_temp <- topdom_domains_temp[!duplicated(topdom_domains_temp),]
      # Add Resolution Feature
      topdom_domains_temp$Resolution <- rep(resolution_vec[i], nrow(topdom_domains_temp))
      topdom_domains <- rbind(topdom_domains, topdom_domains_temp)
    }
    # Transfer small domains
    topdom_transfer_small_boundaries <- data.frame(Chromosome = c(topdom_transfer_small_domains$Chromosome, 
                                                                  topdom_transfer_small_domains$Chromosome),
                                                   Boundary = c(topdom_transfer_small_domains$Start, 
                                                                topdom_transfer_small_domains$End),
                                                   Resolution = c(topdom_transfer_small_domains$Resolution,
                                                                  topdom_transfer_small_domains$Resolution))
    topdom_domains <- rbind(topdom_domains, topdom_transfer_small_boundaries)
    
    #### Pair WUBBA and TopDom Boundaries ####
    
    # Within each group pair WUBBA called TAD with best TopDom counterpart
    wubba_topdom_domains <- data.frame(stringsAsFactors = FALSE)
    error_vec <- resolution_vec * 2
    # error_vec <- rep(resolution_vec[1] * 4, 4)
    # error_vec <- rep(error_val, 4)
    for (i in 1:length(resolution_vec)){
      print(paste("      Wubba and TopDom Pairing Resolution: ", resolution_vec[i]))
      print(paste("      System Time:", Sys.time()))
      # Get TopDom group
      topdom_domains_sub <- topdom_domains[topdom_domains$Resolution == resolution_vec[i] &
                                             topdom_domains$Chromosome == chromosome,]
      # Get Wubba group
      anchor_pairs_input_sub <- anchor_pairs_input[anchor_pairs_input$Resolution == resolution_vec[i] & 
                                                     anchor_pairs_input$Chromosome == chromosome,]
      # Get closest TopDom boundary to each WUBBA boundary
      start_vec <- c()
      end_vec <- c()
      count <- 0
      est_count <- ceiling(length(anchor_pairs_input_sub$Start) / 30000)
      count_left <- 1
      count_right <- 30000
      if(count_right > length(anchor_pairs_input_sub$Start)){
        count_right <- length(anchor_pairs_input_sub$Start)
      }
      print(paste(" Wubba and TopDom Pairing Pass Estimated Count:", est_count))
      while(length(start_vec) < length(anchor_pairs_input_sub$Start)){
        count <- count + 1
        #print(paste("Wubba and TopDom Pairing Pass:", count, "out of", est_count))
        start_vec <- c(start_vec, get_closest(subject_vec = anchor_pairs_input_sub$Start[count_left:count_right], 
                                              reference_vec = topdom_domains_sub$Boundary))
        end_vec <- c(end_vec, get_closest(subject_vec = anchor_pairs_input_sub$End[count_left:count_right], 
                                          reference_vec = topdom_domains_sub$Boundary))
        count_left <- count_left + 30000
        count_right <- count_right + 30000
        if(count_right > length(anchor_pairs_input_sub$Start)){
          count_right <- length(anchor_pairs_input_sub$Start)
        }
      }
      anchor_pairs_input_sub$Start <- start_vec
      anchor_pairs_input_sub$End <- end_vec
      # Calculate differences
      anchor_pairs_input_sub$StartDiff <- abs(anchor_pairs_input_sub$Start - anchor_pairs_input_sub$StartSave)
      anchor_pairs_input_sub$EndDiff <- abs(anchor_pairs_input_sub$End - anchor_pairs_input_sub$EndSave)
      # Remove Far Away hits
      
      
      
      anchor_pairs_input_sub <- anchor_pairs_input_sub[anchor_pairs_input_sub$StartDiff <= error_vec[i] &
                                                         anchor_pairs_input_sub$EndDiff <= error_vec[i],]
      
      
      # Add to hits
      wubba_topdom_domains <- rbind(wubba_topdom_domains, anchor_pairs_input_sub)
    }
    # Remove bad domains where Start is equal to End
    wubba_topdom_domains <- wubba_topdom_domains[wubba_topdom_domains$Start < wubba_topdom_domains$End,]
    
    
    # Sort based on Error, Start, and End
    wubba_topdom_domains$Error <- wubba_topdom_domains$StartDiff + wubba_topdom_domains$EndDiff
    wubba_topdom_domains <- arrange(wubba_topdom_domains, Start, End, Error)
    
    # Sort based on Pval, Start and End
    # wubba_topdom_domains <- arrange(wubba_topdom_domains, Start, End, Pval)
    
    
    # Unfactor Resolution feature and add difference features
    wubba_topdom_domains$Resolution <- as.numeric(levels(wubba_topdom_domains$Resolution))[wubba_topdom_domains$Resolution]
    # Get unique domains with lowest error
    rownames(wubba_topdom_domains) <- NULL
    unique_starts <- unique(wubba_topdom_domains[,1:2])
    unique_index <- as.numeric(rownames(unique_starts))
    wubba_topdom_domains <- wubba_topdom_domains[unique_index,]
    # Return Hits
    return(wubba_topdom_domains)
  }
  #### Finds Significant Boundaries, Pairs and Filters ####
  # Performs wavelet analysis on left contact coverage and coverts output to proper dataframe 
  left_wubba <- comp_wt(input_data = left_cov, 
                        prox = TRUE,
                        alpha = wubba_alpha, 
                        from_scale = wubba_from_scale,
                        to_scale = wubba_to_scale,
                        filt = wubba_filt,
                        p_method = wubba_p_method)
  left_df <- data.frame(stringsAsFactors = FALSE)
  for (i in 1:length(left_wubba)){
    left_df_sub <- data.frame(Start = left_wubba[[i]][,1],
                              End = left_wubba[[i]][,2],
                              Pval = left_wubba[[i]][,3],
                              Avg = left_wubba[[i]][,4],
                              Scale = rep(names(left_wubba)[i], nrow(left_wubba[[i]])),
                              Side = rep("left", nrow(left_wubba[[i]])),
                              stringsAsFactors = FALSE)
    left_df <- rbind(left_df, left_df_sub)
  }
  # Performs wavelet analysis on right contact coverage and coverts output to proper dataframe
  right_wubba <- comp_wt(input_data = right_cov, 
                         prox = TRUE,
                         alpha = wubba_alpha, 
                         from_scale = wubba_from_scale,
                         to_scale = wubba_to_scale,
                         filt = wubba_filt,
                         p_method = wubba_p_method)
  right_df <- data.frame(stringsAsFactors = FALSE)
  for (i in 1:length(right_wubba)){
    right_df_sub <- data.frame(Start = right_wubba[[i]][,1],
                               End = right_wubba[[i]][,2],
                               Pval = right_wubba[[i]][,3],
                               Avg = right_wubba[[i]][,4],
                               Scale = rep(names(right_wubba)[i], nrow(right_wubba[[i]])),
                               Side = rep("right", nrow(right_wubba[[i]])),
                               stringsAsFactors = FALSE)
    right_df <- rbind(right_df, right_df_sub)
  }
  
  # Creates all possible pairs of significant breakpoints between min and max distance
  print(paste("Right Sigs:", nrow(right_df)))
  print(paste("Left Sigs:", nrow(left_df)))
  print(paste("System Time:", Sys.time()))
  pairs_df <- data.frame(stringsAsFactors = FALSE)
  anchor_pairs_df <- data.frame(stringsAsFactors = FALSE)
  primer <- TRUE
  start_index <- 1
  index_bin <- floor(1000000000 / nrow(right_df))
  end_index <- index_bin
  count <- 0
  est_count <- ceiling(nrow(left_df) / index_bin)
  while(primer){
    if(end_index > nrow(left_df)){
      end_index <- nrow(left_df)
    }
    count <- count + 1
    print(paste("Loop:", count, "out of", est_count))
    print(paste("      System Time:", Sys.time()))
    ind <- start_index:end_index
    pairs_expand <- expand.grid(paste(right_df$Avg, right_df$Pval, sep = "_"), 
                                paste(left_df$Avg[ind], left_df$Pval[ind], sep = "_"),
                                stringsAsFactors = FALSE)
    print(paste("      System Time Pairs Expanded Completed:", Sys.time()))
    
    var1_split <- mstrsplit(pairs_expand$Var1, "_")
    var2_split <- mstrsplit(pairs_expand$Var2, "_")
    print(paste("      System Time String Split Completed:", Sys.time()))
    pairs_df_sub <- data.frame(Start = as.numeric(var1_split[,1]),
                               End = as.numeric(var2_split[,1]),
                               Pval1 = as.numeric(var1_split[,2]),
                               Pval2 = as.numeric(var2_split[,2]),
                               stringsAsFactors = FALSE)
    
    print(paste("      System Time First Pairs DF Completed:", Sys.time()))
    pairs_df_sub <- pairs_df_sub[(pairs_df_sub$End - pairs_df_sub$Start) > 
                                   min_pair_distance & 
                                   (pairs_df_sub$End - pairs_df_sub$Start) < 
                                   max_pair_distance,]
    print(paste("      System Time Filtered Pairs DF:", Sys.time()))
    pairs_df_sub$Pval <- pairs_df_sub$Pval1 * pairs_df_sub$Pval2 
    if(end_index == nrow(left_df)){
      primer <- FALSE
    }
    start_index <- start_index + index_bin
    end_index <- end_index + index_bin
    pairs_df <- pairs_df_sub
    pairs_df <- pairs_df[!duplicated(pairs_df[,c(1,2)]),]
    print(paste("      Pairs Dataframe Size:", nrow(pairs_df)))
    print(paste("      System Time:", Sys.time()))
    anchor_pairs_df_sub <- loop_confirm(pairs_input = pairs_df,
                                        chromosome = chr,
                                        loops_input_files = loops_input_files,
                                        resolution_vec = resolution_vec)
    print(paste("      Anchor Pairs Dataframe Size:", nrow(anchor_pairs_df_sub)))
    print(paste("      System Time:", Sys.time()))
    
    anchor_topdom_df <- topdom_collab_pre(anchor_pairs_input = anchor_pairs_df_sub,
                                          chromosome = chr,
                                          topdom_input_files = topdom_input_files,
                                          bin_vec = bin_vec,
                                          resolution_vec = resolution_vec,
                                          error_val = error_val)
    anchor_topdom_df <- anchor_topdom_df[,c(10,11,3:6)]
    colnames(anchor_topdom_df) <- c("Start", "End", "Pval1", "Pval2", "Pval", "AncPval")
    anchor_pairs_df <- rbind(anchor_pairs_df, anchor_topdom_df)
    print(paste("      Anchor TopDom Pairs Dataframe Size:", nrow(anchor_pairs_df)))
    print(paste("      System Time:", Sys.time()))
  }
  return(anchor_pairs_df)
}

#### TopDom Collab ####
topdom_collab <- function(anchor_pairs_input,
                          chromosome,
                          topdom_input_files,
                          bin_vec = c(0, 7500, 17500, 75000),
                          resolution_vec = c(5000, 10000, 25000),
                          error_val = 4000){
  #### Functions ####
  # Takes two vectors and returns the value in the reference vec closest to the nth subject
  get_closest <- function(subject_vec, reference_vec){
    abs_diffs <- abs(outer(subject_vec, reference_vec, "-"))
    min_index <- apply(abs_diffs, 1, which.min)
    closest_pair <- reference_vec[min_index]
  }
  
  #### Prep Anchor Pairs ####
  # Add Chromosome Feature
  anchor_pairs_input$Chromosome <- rep(chromosome, nrow(anchor_pairs_input))
  # Add Size Feature
  anchor_pairs_input$Size <- anchor_pairs_input$End - anchor_pairs_input$Start
  # Add Resolution Feature
  anchor_pairs_input$Resolution <- round(anchor_pairs_input$Size / 40)
  # Group TAD domains by resolution
  anchor_pairs_input$Resolution <-  cut(anchor_pairs_input$Resolution, bin_vec)
  # Rename Groups
  levels(anchor_pairs_input$Resolution) <- resolution_vec 
  # Duplicate Start Feature (Important Later)
  anchor_pairs_input$StartSave <- anchor_pairs_input$Start
  # Duplicate End Feature (Important Later)
  anchor_pairs_input$EndSave <- anchor_pairs_input$End
  #### Prep TopDom Files ####
  # Load TopDom files
  topdom_domains <- data.frame(stringsAsFactors = FALSE)
  topdom_transfer_big_domains <- data.frame(stringsAsFactors = FALSE)
  topdom_transfer_small_domains <- data.frame(stringsAsFactors = FALSE)
  for (i in 1:length(topdom_files)){
    topdom_domains_temp <- read.delim(topdom_input_files[i],
                                      header = FALSE,
                                      stringsAsFactors = FALSE)
    topdom_domains_temp <- topdom_domains_temp[,1:3]
    names(topdom_domains_temp) <- c("Chromosome", "Start", "End")
    topdom_domains_temp <- topdom_domains_temp[topdom_domains_temp$Chromosome == chromosome,]
    # Add Size and Resolution Features
    topdom_domains_temp$Size <- topdom_domains_temp$End - topdom_domains_temp$Start
    topdom_domains_temp$Resolution <- topdom_domains_temp$Size / 40
    # Group TAD domains by resolution
    topdom_domains_temp$Resolution <-  cut(topdom_domains_temp$Resolution, bin_vec)
    # Rename Groups
    levels(topdom_domains_temp$Resolution) <- resolution_vec 
    topdom_domains_temp$Resolution <- as.numeric(levels(topdom_domains_temp$Resolution))[topdom_domains_temp$Resolution]
    # Add to transfer df
    topdom_transfer_big_domains <- rbind(topdom_transfer_big_domains, 
                                         topdom_domains_temp[topdom_domains_temp$Resolution > resolution_vec[i],])
    topdom_transfer_small_domains <- rbind(topdom_transfer_small_domains, 
                                           topdom_domains_temp[topdom_domains_temp$Resolution < resolution_vec[i],])
    topdom_domains_temp <- rbind(topdom_domains_temp, 
                                 topdom_transfer_big_domains[topdom_transfer_big_domains$Resolution == resolution_vec[i],])
    topdom_domains_temp <- rbind(topdom_domains_temp, 
                                 topdom_transfer_small_domains[topdom_transfer_small_domains$Resolution == resolution_vec[i],])
    # Convert to boundaries
    topdom_domains_temp <- data.frame(Chromosome = c(topdom_domains_temp$Chromosome, topdom_domains_temp$Chromosome),
                                      Boundary = c(topdom_domains_temp$Start, topdom_domains_temp$End))
    # Remove Duplicates
    topdom_domains_temp <- topdom_domains_temp[!duplicated(topdom_domains_temp),]
    # Add Resolution Feature
    topdom_domains_temp$Resolution <- rep(resolution_vec[i], nrow(topdom_domains_temp))
    topdom_domains <- rbind(topdom_domains, topdom_domains_temp)
  }
  # Transfer small domains
  topdom_transfer_small_boundaries <- data.frame(Chromosome = c(topdom_transfer_small_domains$Chromosome, 
                                                                topdom_transfer_small_domains$Chromosome),
                                                 Boundary = c(topdom_transfer_small_domains$Start, 
                                                              topdom_transfer_small_domains$End),
                                                 Resolution = c(topdom_transfer_small_domains$Resolution,
                                                                topdom_transfer_small_domains$Resolution))
  topdom_domains <- rbind(topdom_domains, topdom_transfer_small_boundaries)
  
  # Remove large domains
  topdom_domains <- topdom_domains[!is.na(topdom_domains$Boundary),]
  
  #### Pair WUBBA and TopDom Boundaries ####
  
  # Within each group pair WUBBA called TAD with best TopDom counterpart
  wubba_topdom_domains <- data.frame(stringsAsFactors = FALSE)
  error_vec <- resolution_vec * 2
  # error_vec <- rep(resolution_vec[1] * 4, 4)
  # error_vec <- rep(error_val, 4)
  for (i in 1:length(resolution_vec)){
    print(paste("Wubba and TopDom Pairing Resolution: ", resolution_vec[i]))
    # Get TopDom group
    topdom_domains_sub <- topdom_domains[topdom_domains$Resolution == resolution_vec[i] &
                                           topdom_domains$Chromosome == chromosome,]
    # Get Wubba group
    anchor_pairs_input_sub <- anchor_pairs_input[anchor_pairs_input$Resolution == resolution_vec[i] & 
                                                   anchor_pairs_input$Chromosome == chromosome,]
    # Get closest TopDom boundary to each WUBBA boundary
    start_vec <- c()
    end_vec <- c()
    count <- 0
    est_count <- ceiling(length(anchor_pairs_input_sub$Start) / 40000)
    count_left <- 1
    count_right <- 40000
    if(count_right > length(anchor_pairs_input_sub$Start)){
      count_right <- length(anchor_pairs_input_sub$Start)
    }
    print(paste(" Wubba and TopDom Pairing Pass Estimated Count:", est_count))
    while(length(start_vec) < length(anchor_pairs_input_sub$Start)){
      count <- count + 1
      start_vec <- c(start_vec, get_closest(anchor_pairs_input_sub$Start[count_left:count_right], topdom_domains_sub$Boundary))
      end_vec <- c(end_vec, get_closest(anchor_pairs_input_sub$End[count_left:count_right], topdom_domains_sub$Boundary))
      count_left <- count_left + 40000
      count_right <- count_right + 40000
      if(count_right > length(anchor_pairs_input_sub$Start)){
        count_right <- length(anchor_pairs_input_sub$Start)
      }
    }
    anchor_pairs_input_sub$Start <- start_vec
    anchor_pairs_input_sub$End <- end_vec
    # Calculate differences
    anchor_pairs_input_sub$StartDiff <- abs(anchor_pairs_input_sub$Start - anchor_pairs_input_sub$StartSave)
    anchor_pairs_input_sub$EndDiff <- abs(anchor_pairs_input_sub$End - anchor_pairs_input_sub$EndSave)
    # Remove Far Away hits
    
    
    
    anchor_pairs_input_sub <- anchor_pairs_input_sub[anchor_pairs_input_sub$StartDiff <= error_vec[i] &
                                                       anchor_pairs_input_sub$EndDiff <= error_vec[i],]
    
    
    # Add to hits
    wubba_topdom_domains <- rbind(wubba_topdom_domains, anchor_pairs_input_sub)
  }
  # Remove bad domains where Start is equal to End
  wubba_topdom_domains <- wubba_topdom_domains[wubba_topdom_domains$Start < wubba_topdom_domains$End,]
  
  
  # Sort based on Error, Start, and End
  wubba_topdom_domains$Error <- wubba_topdom_domains$StartDiff + wubba_topdom_domains$EndDiff
  wubba_topdom_domains <- arrange(wubba_topdom_domains, Start, End, Error)
  
  # Sort based on Pval, Start and End
  # wubba_topdom_domains <- arrange(wubba_topdom_domains, Start, End, Pval)
  
  
  # Unfactor Resolution feature and add difference features
  wubba_topdom_domains$Resolution <- as.numeric(levels(wubba_topdom_domains$Resolution))[wubba_topdom_domains$Resolution]
  # Get unique domains with lowest p-value
  rownames(wubba_topdom_domains) <- NULL
  unique_starts <- unique(wubba_topdom_domains[,1:2])
  unique_index <- as.numeric(rownames(unique_starts))
  wubba_topdom_domains <- wubba_topdom_domains[unique_index,]
  # Return Hits
  return(wubba_topdom_domains)
}

#### WUBBA Filtering Function ####
wubba_filtering <- function(domains_df,
                            max_passes = 5,
                            min_side_distance = 0,
                            metric = "Pval"){
  
  primer <- TRUE
  best_domains_df <- data.frame(stringsAsFactors = FALSE)
  count <- 1
  while(primer){
    print(paste("pass: ", count, sep = ""))
    too_close_pairs <- data.frame(stringsAsFactors = FALSE)
    for (i in 1:nrow(best_domains_df)){
      check1 <- abs(domains_df$Start - best_domains_df$Start[i]) <= min_side_distance &
        abs(domains_df$End - best_domains_df$End[i]) <= min_side_distance
      too_close_pairs <- rbind(too_close_pairs, domains_df[check1,])
    }
    too_close_pairs <- too_close_pairs[!duplicated(too_close_pairs),]
    
    best_pass_domains_df <- data.frame(stringsAsFactors = FALSE)
    pass_domains_df <- domains_df[!(domains_df$Pval %in% too_close_pairs$Pval),]
    while(nrow(pass_domains_df) != 0){
      pass_domains_df <- pass_domains_df[order(pass_domains_df[[metric]]),]
      best_pass_domains_df <- rbind(pass_domains_df[1,], best_pass_domains_df)
      check1 <- (pass_domains_df$Start[1] >= pass_domains_df$Start) & 
        (pass_domains_df$Start[1] < pass_domains_df$End)
      check2 <- (pass_domains_df$End[1] > pass_domains_df$Start) & 
        (pass_domains_df$End[1] <= pass_domains_df$End)
      check3 <- (pass_domains_df$Start[1] <= pass_domains_df$Start) & 
        (pass_domains_df$End[1] >= pass_domains_df$End)
      pass_domains_df <- pass_domains_df[!(check1 | check2 | check3),]
    }
    if (nrow(best_pass_domains_df) == 0){
      primer <- FALSE
    }
    Pass <- rep(count, nrow(best_pass_domains_df))
    count <- count + 1
    best_pass_domains_df <- cbind(best_pass_domains_df, Pass)
    best_domains_df <- rbind(best_domains_df, best_pass_domains_df)
    if (count > max_passes){
      primer <- FALSE
    }
  }
  return(best_domains_df)
  
  
}

#### Application ####
# Adjust for genome size
if(genome_size == "big"){
  min_pair_distance <- 40000
  max_pair_distance <- 5000000
  resolution_size <- 5000
  bin_vec <- c(0, 1625, 3000, 25000, (max_pair_distance/40))
  resolution_vec <- c(5000, 10000, 25000, 50000)
  error_val <- 8000
} else {
  min_pair_distance <- 40000
  max_pair_distance <- 5000000
  resolution_size <- 1000
  # bin_vec <- c(0, 400, 1625, 3000, 25000)
  bin_vec <- c(0, 400, 1625, 3000, (max_pair_distance/40))
  resolution_vec <- c(1000, 5000, 10000, 25000)
  error_val <- 4000
}


# Load Coverage and Contact Files
min_cov <- 1
left_cov <- read.delim(left_coverage_file,
                       stringsAsFactors = FALSE,
                       header = FALSE)
names(left_cov) <- c("Chromosome", "Location", "Value")
left_cov_sub <- left_cov[left_cov[,1] == chr,]
left_cov_sub <- left_cov_sub[left_cov_sub[,3] >= min_cov,]
left_cov_sub[,3] <- log(left_cov_sub[,3])
rm(left_cov)


right_cov <- read.delim(right_coverage_file,
                        stringsAsFactors = FALSE,
                        header = FALSE)
names(right_cov) <- c("Chromosome", "Location", "Value")
right_cov_sub <- right_cov[right_cov[,1] == chr,]
right_cov_sub <- right_cov_sub[right_cov_sub[,3] >= min_cov,]
right_cov_sub[,3] <- log(right_cov_sub[,3])
rm(right_cov)

print(paste("Chromosome: ", 
            chr,
            sep = ""))

# Begin Pairing
print("Begin Pairing")
print(paste("System Time:", Sys.time()))
pairs_df <- wubba_pairing(left_cov = left_cov_sub,
                          right_cov = right_cov_sub,
                          chromosome = chr,
                          min_pair_distance = min_pair_distance,
                          max_pair_distance = max_pair_distance,
                          wubba_alpha = 0.05, 
                          wubba_from_scale = 5,
                          wubba_to_scale = -2,
                          wubba_filt = "c6",
                          wubba_p_method = "holm",
                          resolution_size = resolution_size,
                          topdom_input_files = topdom_files,
                          loops_input_files = loops_files,
                          bin_vec = bin_vec,
                          resolution_vec = resolution_vec,
                          error_val = error_val)

# Begin TopDom Collab
print("Begin TopDom Collab")
print(paste("System Time:", Sys.time()))
collab_pairs <- topdom_collab(anchor_pairs_input = pairs_df,
                              chromosome = chr,
                              topdom_input_files = topdom_files,
                              bin_vec = bin_vec,
                              resolution_vec = resolution_vec,
                              error_val = error_val)

# Begin Filtering
print("Begin Filtering")
print(paste("System Time:", Sys.time()))
results <- wubba_filtering(domains_df = collab_pairs,
                           max_passes = 10000,
                           min_side_distance = 40000,
                           metric = "Pval")

print(paste("Results Data Size: ", 
            nrow(results),
            sep = ""))



#### Create Final File ####
options(scipen=10)

results$FinalPval <- results$Pval * results$AncPval
results <- results[,c(7,10,11,3,4,6,16)]
colnames(results) <- c("Chromosome", "Start", "End", "Boundary_5_Prime_Pval", "Boundary_3_Prime_Pval", "Loop_Pval", "Pval")

write.table(results,
            results_file,
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")




