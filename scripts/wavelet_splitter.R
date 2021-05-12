options(scipen = 10)
input_args <- commandArgs(TRUE)
input_bed1 <- input_args[1]
input_bed2 <- input_args[2]
output_bed1 <- input_args[3]
output_bed2 <- input_args[4]
read_length <- as.numeric(input_args[5])
chr_sizes <- input_args[6]

chr_df <- read.delim(chr_sizes,
                     stringsAsFactors = FALSE,
                     header = FALSE,
                     col.names = c("Chromosome", "Size"))

chr_vec <- chr_df$Chromosome

# First file edits
st1 <- read.delim(input_bed1, 
                  stringsAsFactors = FALSE, 
                  header = FALSE,
                  col.names = c("ID", "Flag", "Chromosome", "Position"))

st1_int_1 <- st1[,-2]
rm(st1)
st1_int_2 <- st1_int_1[st1_int_1$Chromosome %in% chr_vec,]
rm(st1_int_1)

# Second file edits
st2 <- read.delim(input_bed2, stringsAsFactors = FALSE, header = FALSE,
                  col.names = c("ID", "Flag", "Chromosome", "Position"))

st2_int_1 <- st2[,-2]
rm(st2)
st2_int_2 <- st2_int_1[st2_int_1$Chromosome %in% chr_vec,]
rm(st2_int_1)


# Ensure proper pairing
pst1_int_1 <- st1_int_2[st1_int_2$ID %in% st2_int_2$ID,]
pst2_int_1 <- st2_int_2[st2_int_2$ID %in% st1_int_2$ID,]

rm(st1_int_2)
rm(st2_int_2)

# Ensure seperating distance
pst1_int_2 <- pst1_int_1[(pst1_int_1$Chromosome == pst2_int_1$Chromosome) &
                           (abs(as.numeric(pst1_int_1$Position) - 
                                  as.numeric(pst2_int_1$Position)) > 500),]
pst2_int_2 <- pst2_int_1[(pst2_int_1$Chromosome == pst1_int_1$Chromosome) & 
                           (abs(as.numeric(pst1_int_1$Position) - 
                                  as.numeric(pst2_int_1$Position)) > 500),]

rm(pst1_int_1)
rm(pst2_int_1)

pst1_int_3 <- pst1_int_2[(pst1_int_2$Chromosome == pst2_int_2$Chromosome) &
                           (abs(as.numeric(pst1_int_2$Position) - 
                                  as.numeric(pst2_int_2$Position)) <= 5000000),]
pst2_int_3 <- pst2_int_2[(pst2_int_2$Chromosome == pst1_int_2$Chromosome) & 
                           (abs(as.numeric(pst1_int_2$Position) - 
                                  as.numeric(pst2_int_2$Position)) <= 5000000),]
rm(pst1_int_2)
rm(pst2_int_2)

right_contacts_sam <- rbind(pst1_int_3[pst1_int_3$Position < pst2_int_3$Position,],
                            pst2_int_3[!(pst1_int_3$Position < pst2_int_3$Position),])
left_contacts_sam<- rbind(pst2_int_3[!(pst1_int_3$Position > pst2_int_3$Position),],
                          pst1_int_3[pst1_int_3$Position > pst2_int_3$Position,])

rm(pst1_int_3)
rm(pst2_int_3)

right_contacts_sam <- right_contacts_sam[,-1]
left_contacts_sam <- left_contacts_sam[,-1]

right_contacts_sam$End <- right_contacts_sam$Position + read_length - 1
left_contacts_sam$End <- left_contacts_sam$Position + read_length - 1

right_contacts_sam <- right_contacts_sam[order(right_contacts_sam$Position),]
right_contacts_sam <- right_contacts_sam[order(right_contacts_sam$Chromosome),]
left_contacts_sam <- left_contacts_sam[order(left_contacts_sam$Position),]
left_contacts_sam <- left_contacts_sam[order(left_contacts_sam$Chromosome),]

write.table(right_contacts_sam, output_bed1, col.names = FALSE,
            row.names = FALSE, quote = FALSE, na = "", sep = "\t")
write.table(left_contacts_sam, output_bed2, col.names = FALSE,
            row.names = FALSE, quote = FALSE, na = "", sep = "\t")

