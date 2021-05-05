#### Files ####
input_args <- commandArgs(TRUE)
best_files <- input_args[1]
result_file <- input_args[2]

#### Work ####
options(scipen=10)

best_files_vec <- unlist(strsplit(best_files, ".bed"))
for (i in 1:length(best_files_vec)){
  best_files_vec[i] <- paste(best_files_vec[i], ".bed", sep = "")
}

best_df <- data.frame(stringsAsFactors = FALSE)
for (i in 1:length(best_files_vec)) {
  if (file.exists(best_files_vec[i])) {
    best_df_temp <- read.delim(best_files_vec[i],
                               header = TRUE)
    best_df <- rbind(best_df, best_df_temp)
  }
}
write.table(best_df,
            result_file,
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")