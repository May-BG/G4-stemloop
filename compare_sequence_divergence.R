library(tidyverse)
library(stringi)
library(data.table)
rc <- function(nucSeq)
  return(stri_reverse(chartr("acgtACGT", "tgcaTGCA", nucSeq)))

args <- commandArgs(trailingOnly = TRUE)
input <- paste0("hs1.chr", args[1], ".2023v2.processed.analysed.quadron.dat")
human_output <- paste0("shared_g4/hsa", args[1], ".human.chimp.tri.allG4.noanc.csv")
orang_output <- paste0("shared_g4/hsa", args[1], ".human.orang.tri.allG4.noanc.csv")
human_output_all <- paste0("shared_g4/hsa", args[1], ".human.chimp.snp.allG4.noanc.csv")
orang_output_all <- paste0("shared_g4/hsa", args[1], ".human.orang.snp.allG4.noanc.csv")
g4_table_outfile <- paste0("shared_g4/hsa", args[1], ".g4table.allG4.txt")



lines <- readLines(input)



data_list <- list()


current_number1 <- NA
current_number2 <- NA


process_line <- function(line) {
  if (grepl("^##", line)) {
    return(NULL)
  } else if (grepl("^#", line)) {
    current_number1 <<- as.numeric(gsub("#", "", line))
    current_number2 <<- NA
    return(NULL)
  } else if (grepl("^[0-9]+$", line)) {
    current_number2 <<- as.numeric(line)
    return(NULL)
  } else {
    id <- as.numeric(sub("@.*", "", strsplit(line, " ")[[1]][1]))
    if (id %in% c(1, 3, 5, 11)) {
      parts <- strsplit(gsub("\\s+", " ", trimws(line)), " ")[[1]]
      if (length(parts) >= 8) {
        list(Number1 = current_number1, Number2 = current_number2,
             ID = id, Start = as.numeric(parts[2]),
             Query_Length = ifelse(parts[3] == ".", NA, as.numeric(parts[3])),
             Type = parts[4], Score = ifelse(parts[5] == ".", NA, as.numeric(parts[5])),
             Strand = parts[6], Alignment_Length = as.numeric(parts[7]),
             Sequence = parts[8])
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
  }
}

data_list <- lapply(lines, process_line)
data_list <- Filter(Negate(is.null), data_list)

# Convert the list of data rows to a data frame
data <- do.call(rbind, data_list)
colnames(data) <- c("Number1", "Number2", "SpeciesID_ChrID", "Start", "Query_Length",
                    "Type", "Score", "Strand", "Alignment_Length", "Sequence")
filtered = as.data.table(data)
head(filtered)

filtered$Number1 <- sapply(filtered$Number1, `[`, 1)  # Assuming you want the first element
filtered$Number2 <- sapply(filtered$Number2, `[`, 1)  # Adjust the index if needed
filtered$SpeciesID_ChrID <- sapply(filtered$SpeciesID_ChrID, `[`, 1)  # Assuming you want the first element
filtered$Start <- sapply(filtered$Start, `[`, 1)  # Adjust the index if needed
filtered$Type <- sapply(filtered$Type, `[`, 1)  # Assuming you want the first element
filtered$Sequence <- sapply(filtered$Sequence, `[`, 1)  # Adjust the index if needed
filtered$Score <- sapply(filtered$Score, `[`, 1)  # Assuming you want the first element
filtered$Query_Length <- sapply(filtered$Query_Length, `[`, 1)  # Adjust the index if needed
# Convert list columns to their first element (if appropriate)
filtered$Strand <- sapply(filtered$Strand, `[`, 1)
filtered$Alignment_Length <- sapply(filtered$Alignment_Length, `[`, 1)
head(filtered)

DT.c2 <- dcast(
  filtered, 
  Number1 + Number2 ~ SpeciesID_ChrID, 
  value.var = c("Start", "Type", "Score", "Query_Length", "Sequence","Strand"),
  fun.aggregate = function(x) if (is.numeric(x)) mean(x, na.rm = TRUE) else first(na.omit(x))
)
#DT.c2[Type_1=="fullMAF" & Type_3=="fullMAF" &Type_5=="fullMAF",]
g4_table <- DT.c2[Type_1=="fullMAF" & Type_3=="fullMAF" &Type_5=="fullMAF",.(Strand_1,Start_1, Start_3, Query_Length_1, Score_1,Sequence_1, Sequence_3, Sequence_5, Sequence_11)]
#colnames(g4_table) <- c("human_start","chimp_start","human_length","human_score","human_seq","chimp_seq","orang_seq","anc_seq")
colnames(g4_table) <- c("Strand_1","human_start","chimp_start","human_length","human_score","human_seq","chimp_seq","orang_seq","anc_seq")
g4_table$Strand_1 <- as.character(g4_table$Strand_1)
g4_table <- na.omit(g4_table)
g4_table$human_start <- as.numeric(as.character(g4_table$human_start))
g4_table$chimp_start <- as.numeric(as.character(g4_table$chimp_start))
g4_table$human_score <- as.numeric(as.character(g4_table$human_score))
g4_table$human_seq <- as.character(g4_table$human_seq)
g4_table$chimp_seq <- as.character(g4_table$chimp_seq)
g4_table$orang_seq <- as.character(g4_table$orang_seq)
g4_table$anc_seq <- as.character(g4_table$anc_seq)
#write.table(g4_table, "g4_table_outfile.txt", sep="\t", quote=FALSE, row.names=FALSE)

write.table(g4_table, g4_table_outfile, sep="\t", quote=FALSE, row.names=FALSE)
df <- as.data.frame(g4_table)
human_chimp_diffs <- data.frame(position = integer(), context_human = character(), context_chimp = character(), stringsAsFactors = FALSE)
human_orang_diffs <- data.frame(position = integer(), context_human = character(), context_orang = character(), stringsAsFactors = FALSE)
for (i in 1:nrow(df)) {
human <- tolower(strsplit(df$human_seq[i], "")[[1]])
 print(human)
if (human[1]=="c") {
    new_df <- df
    new_df$human_seq[i] <- rc(df$human_seq[i])
    new_df$chimp_seq[i] <- rc(df$chimp_seq[i])
    new_df$orang_seq[i] <- rc(df$orang_seq[i])
    new_df$anc_seq[i] <- rc(df$anc_seq[i])
    human <- tolower(strsplit(new_df$human_seq[i], "")[[1]])
  chimp <- tolower(strsplit(new_df$chimp_seq[i], "")[[1]])
  orang <- tolower(strsplit(new_df$orang_seq[i], "")[[1]])
  anc <- tolower(strsplit(new_df$anc_seq[i], "")[[1]])
    if (df$Strand_1[i]=="+"){
     for (j in 1:length(human)) {  
  if (j > 1 && j < length(human)) {
      print(c(j, human[j],orang[j]))# Ensure j-1 and j+1 are valid indices
    if (human[j] != chimp[j]) {
      # Extract three-nucleotide contexts
      context_human <- paste(human[(j-1):(j+1)], collapse = "")
      context_chimp <- paste(chimp[(j-1):(j+1)], collapse = "")
        gaps_human <- sum(human[1:j]=="-")
      #print(c(j, df$human_start[i]+j,context_human,context_anc))
        # Add to data frame
      human_chimp_diffs <- rbind(human_chimp_diffs, data.frame(position = df$human_start[i]+j-gaps_human, context_human = context_human, context_chimp = context_chimp))
    }
    if (!is.na(human[j]) && !is.na(orang[j]) && human[j] != orang[j]) {
      # Extract three-nucleotide contexts
      context_human <- paste(human[(j-1):(j+1)], collapse = "")
      context_orang <- paste(orang[(j-1):(j+1)], collapse = "")
        gaps_human <- sum(human[1:j]=="-")
      # Add to data frame
      human_orang_diffs <- rbind(human_orang_diffs, data.frame(position = df$human_start[i]+j-gaps_human, context_human = context_human, context_orang = context_orang))
    }
  }
    }   
    }
    if (df$Strand_1[i]=="-"){
    for (j in 1:length(human)) {  
  if (j > 1 && j < length(human)) {
      k = length(human)+1-j
      print(c(j, human[j],orang[j]))# Ensure j-1 and j+1 are valid indices
    if (human[j] != chimp[j]) {
      # Extract three-nucleotide contexts
      context_human <- paste(human[(j-1):(j+1)], collapse = "")
      context_chimp <- paste(chimp[(j-1):(j+1)], collapse = "")
        gaps_human <- sum(human[1:j]=="-")
      #print(c(j, df$human_start[i]+j,context_human,context_anc))
        # Add to data frame
      human_chimp_diffs <- rbind(human_chimp_diffs, data.frame(position = df$human_start[i]+k-gaps_human, context_human = context_human, context_chimp = context_chimp))
    }
    if (!is.na(human[j]) && !is.na(orang[j]) && human[j] != orang[j]) {
      # Extract three-nucleotide contexts
      context_human <- paste(human[(j-1):(j+1)], collapse = "")
      context_orang <- paste(orang[(j-1):(j+1)], collapse = "")
        gaps_human <- sum(human[1:j]=="-")
      # Add to data frame
      human_orang_diffs <- rbind(human_orang_diffs, data.frame(position = df$human_start[i]+j-gaps_human, context_human = context_human, context_orang = context_orang))
    }
  }
    }
}
} else {
    human <- tolower(strsplit(df$human_seq[i], "")[[1]])
  chimp <- tolower(strsplit(df$chimp_seq[i], "")[[1]])
  orang <- tolower(strsplit(df$orang_seq[i], "")[[1]])
  anc <- tolower(strsplit(df$anc_seq[i], "")[[1]])
print(human)
     if (df$Strand_1[i]=="+"){
for (j in 1:length(human)) {  
  if (j > 1 && j < length(human)) {
      print(c(j, human[j],chimp[j]))# Ensure j-1 and j+1 are valid indices
    if (human[j] != chimp[j]) {
      # Extract three-nucleotide contexts
      context_human <- paste(human[(j-1):(j+1)], collapse = "")
      context_chimp <- paste(chimp[(j-1):(j+1)], collapse = "")
        gaps_human <- sum(human[1:j]=="-")
      #print(c(j, df$human_start[i]+j,context_human,context_anc))
        # Add to data frame
      human_chimp_diffs <- rbind(human_chimp_diffs, data.frame(position = df$human_start[i]+j-gaps_human, context_human = context_human, context_chimp = context_chimp))
    }
    if (!is.na(human[j]) && !is.na(orang[j]) && human[j] != orang[j]) {
      # Extract three-nucleotide contexts
      context_human <- paste(human[(j-1):(j+1)], collapse = "")
      context_orang <- paste(orang[(j-1):(j+1)], collapse = "")
        gaps_human <- sum(human[1:j]=="-")
      # Add to data frame
      human_orang_diffs <- rbind(human_orang_diffs, data.frame(position = df$human_start[i]+j-gaps_human, context_human = context_human, context_orang = context_orang))
    }
  }
    }
         }
            if (df$Strand_1[i]=="-"){
    for (j in 1:length(human)) {  
  if (j > 1 && j < length(human)) {
      k = length(human)+1-j
      print(c(j, human[j],chimp[j]))# Ensure j-1 and j+1 are valid indices
    if (human[j] != chimp[j]) {
      # Extract three-nucleotide contexts
      context_human <- paste(human[(j-1):(j+1)], collapse = "")
      context_chimp <- paste(chimp[(j-1):(j+1)], collapse = "")
        gaps_human <- sum(human[1:j]=="-")
      #print(c(j, df$human_start[i]+j,context_human,context_anc))
      human_chimp_diffs <- rbind(human_chimp_diffs, data.frame(position = df$human_start[i]+k-gaps_human, context_human = context_human, context_chimp = context_chimp))
    }
    if (!is.na(human[j]) && !is.na(orang[j]) && human[j] != orang[j]) {
      # Extract three-nucleotide contexts
      context_human <- paste(human[(j-1):(j+1)], collapse = "")
      context_orang <- paste(orang[(j-1):(j+1)], collapse = "")
        gaps_human <- sum(human[1:j]=="-")
      # Add to data frame
      human_orang_diffs <- rbind(human_orang_diffs, data.frame(position = df$human_start[i]+j-gaps_human, context_human = context_human, context_orang = context_orang))
    }
  }
    }
}
    }
}
#write.csv(human_chimp_diffs, human_output_all)
#write.csv(human_orang_diffs, orang_output_all)


# Filter out rows where context contains "-" or "n"
human_chimp_diffs_new <- human_chimp_diffs[!grepl("-", human_chimp_diffs$context_human) & 
                                     !grepl("n", human_chimp_diffs$context_human) & 
                                     !grepl("-", human_chimp_diffs$context_chimp) & 
                                     !grepl("n", human_chimp_diffs$context_chimp), ]
human_orang_diffs_new <- human_orang_diffs[!grepl("-", human_orang_diffs$context_human) & 
                                     !grepl("n", human_orang_diffs$context_human) & 
                                     !grepl("-", human_orang_diffs$context_orang) & 
                                     !grepl("n", human_orang_diffs$context_orang), ]

human_chimp_diffs_new <- human_chimp_diffs_new %>% subset(substr(context_chimp, 1, 1)==substr(context_human, 1, 1) & substr(context_chimp, 3, 3)==substr(context_human, 3, 3))

human_orang_diffs_new <- human_orang_diffs_new %>% subset(substr(context_orang, 1, 1)==substr(context_human, 1, 1) & substr(context_orang, 3, 3)==substr(context_human, 3, 3))

human_chimp_diffs_new <- human_chimp_diffs_new %>%
  mutate(
    middle_base_change = paste0(substr(context_chimp, 2, 2), ">",substr(context_human, 2, 2))
  )
  
human_orang_diffs_new <- human_orang_diffs_new %>%
  mutate(
    middle_orangtohuman_base_change = paste0(substr(context_orang, 2, 2), ">",substr(context_human, 2, 2))
  )
write.csv(human_chimp_diffs_new, human_output)
write.csv(human_orang_diffs_new, orang_output)
write.csv(human_chimp_diffs, human_output_all)
write.csv(human_orang_diffs, orang_output_all)
