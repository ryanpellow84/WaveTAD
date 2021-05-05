#### Files ####
input_args <- commandArgs(TRUE)

#### Libraries ####
library(TopDom)
library(reshape2)
library(zoo)
#### Correct Missing Bins ####
correct_missing_bins <- function(bad_contact_matrix, 
                                 bin_size,
                                 chr_length){
  full_start <- seq(0, chr_length, by = bin_size)
  full_start <- full_start[1:(length(full_start)-1)]
  full_end <- seq(bin_size, chr_length, by = bin_size)
  bad_bins_start1 <- as.numeric(levels(factor(bad_contact_matrix$start1)))
  bad_bins_end1 <- as.numeric(levels(factor(bad_contact_matrix$end1)))
  bad_bins_start2 <- colnames(bad_contact_matrix)
  bad_bins_start2 <- bad_bins_start2[-c(1:3)]
  bad_bins_start2 <- unlist(strsplit(bad_bins_start2, "_"))
  bad_bins_start2_seq <- seq(2, length(bad_bins_start2), by = 3)
  bad_bins_start2 <- as.numeric(bad_bins_start2[bad_bins_start2_seq])
  bad_bins_end2 <- colnames(bad_contact_matrix)
  bad_bins_end2 <- bad_bins_end2[-c(1:3)]
  bad_bins_end2 <- unlist(strsplit(bad_bins_end2, "_"))
  bad_bins_end2_seq <- seq(3, length(bad_bins_end2), by = 3)
  bad_bins_end2 <- as.numeric(bad_bins_end2[bad_bins_end2_seq])
  bad_start1 <- full_start[!(full_start %in% bad_bins_start1)]
  bad_end1 <- full_end[!(full_end %in% bad_bins_end1)]
  bad_start2 <- full_start[!(full_start %in% bad_bins_start2)]
  bad_end2 <- full_end[!(full_end %in% bad_bins_end2)]
  
  
  temp_matrix <- bad_contact_matrix[,-c(1:3)]
  colnames(temp_matrix) <- bad_bins_start2
  rownames(temp_matrix) <- bad_bins_start1
  na_matrix_rows <- matrix(NA, nrow = length(bad_start1), ncol = ncol(temp_matrix))
  rownames(na_matrix_rows) <- bad_start1
  colnames(na_matrix_rows) <- bad_bins_start2
  temp_matrix <- rbind(temp_matrix, na_matrix_rows)
  na_matrix_cols <- matrix(NA, nrow = nrow(temp_matrix), ncol = length(bad_start2))
  rownames(na_matrix_cols) <- rownames(temp_matrix)
  colnames(na_matrix_cols) <- bad_start2
  temp_matrix <- cbind(temp_matrix, na_matrix_cols)
  temp_matrix <- temp_matrix[order(as.numeric(rownames(temp_matrix))),]
  temp_matrix <- temp_matrix[,order(as.numeric(colnames(temp_matrix)))]
  temp_matrix[is.na(temp_matrix)] <- -100
  temp_matrix_long <- melt(temp_matrix)
  temp_matrix_long$value[temp_matrix_long$value == -100] <- NA
  if(is.na(temp_matrix_long[1,2])){
    temp_matrix_long[1,2] <- 0
  }
  if(is.na(temp_matrix_long[nrow(temp_matrix_long),2])){
    temp_matrix_long[nrow(temp_matrix_long),2] <- 0
  }
  temp_matrix_long$value <- na.approx(temp_matrix_long$value)
  temp_matrix_no_na <- matrix(temp_matrix_long$value, 
                              nrow = nrow(temp_matrix), 
                              ncol = ncol(temp_matrix))
  temp_upper_tri_w_diag <- upper.tri(temp_matrix_no_na, diag = TRUE)*temp_matrix_no_na
  temp_upper_tri <- upper.tri(temp_matrix_no_na, diag = FALSE)*temp_matrix_no_na
  temp_lower_tri <- t(temp_upper_tri)
  temp_matrix_sym <- temp_upper_tri_w_diag + temp_lower_tri
  return(temp_matrix_sym)
}
#### Application ####
contacts <- read.table(input_args[1], header = TRUE)
chr_vec <- read.table(input_args[5],
                      header = FALSE,
                      stringsAsFactors = FALSE)[,1]
bin_size <- as.numeric(input_args[4])
topdom_domains <- data.frame(stringsAsFactors = FALSE)
for (i in 1:length(chr_vec)){
  print(paste("Chromosome: ",
              chr_vec[i],
              sep = ""))
  contacts_sub <- contacts[contacts$chrom1 == chr_vec[i] & 
                             contacts$chrom2 == chr_vec[i],]
  max_chr_length <- max(c(contacts_sub$end1, contacts_sub$end2))
  ## new section to ensure proper dimensions ##
  contacts_sub <- contacts_sub[contacts_sub$start1 %in% contacts_sub$start2 & 
                                 contacts_sub$start2 %in% contacts_sub$start1,]
  contacts_sub <- contacts_sub[contacts_sub$end1 %in% contacts_sub$end2 & 
                                 contacts_sub$end2 %in% contacts_sub$end1,]
  ## ##
  print("Contacts Loaded")
  contacts_sub <- dcast(contacts_sub, 
                        chrom1+start1+end1 ~ chrom2+start2+end2, 
                        value.var = "count")
  contacts_sub <- correct_missing_bins(bad_contact_matrix = contacts_sub,
                                       bin_size = bin_size,
                                       chr_length = max_chr_length)
  print(paste("Dimension Check 1: ",
              dim(contacts_sub),
              sep = ""))
  
  write.table(contacts_sub, 
              input_args[2],
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE, 
              sep = "\t")
  contacts_sub <- readHiC(input_args[2],
                          chr = chr_vec[i],
                          binSize = bin_size)
  topdom <- TopDom(data = contacts_sub,
                   window.size = 5)
  
  topdom_domains_temp <- topdom$bed[topdom$bed$name == "domain",]
  topdom_domains <- rbind(topdom_domains, topdom_domains_temp)
  print(paste("Domain Check: ",
              nrow(topdom_domains),
              sep = ""))
}

topdom_domains_bed <- data.frame(Chr = topdom_domains$chrom,
                                 Start = as.numeric(topdom_domains$chromStart),
                                 End = as.numeric(topdom_domains$chromEnd),
                                 ID = 1:nrow(topdom_domains),
                                 Pval = rep(.5, nrow(topdom_domains)),
                                 Blank = rep(".", nrow(topdom_domains)),
                                 Start1 = as.numeric(topdom_domains$chromStart),
                                 End1 = as.numeric(topdom_domains$chromEnd),
                                 Coordinates = rep(c("31,120,180", "51,160,44"), 
                                                   length.out = nrow(topdom_domains)),
                                 stringsAsFactors = FALSE)
topdom_domains_bed$Start <- format(topdom_domains_bed$Start, scientific = FALSE)
topdom_domains_bed$End <- format(topdom_domains_bed$End, scientific = FALSE)
topdom_domains_bed$Start1 <- format(topdom_domains_bed$Start1, scientific = FALSE)
topdom_domains_bed$End1 <- format(topdom_domains_bed$End1, scientific = FALSE)

write.table(topdom_domains_bed, 
            input_args[3],
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE, 
            sep = "\t")
