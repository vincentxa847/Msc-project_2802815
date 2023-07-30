## read the "SRRname.txt" in R, convert it to pair character

SRR_and_ERR_for_index = readr::read_lines("SRRname.txt")

SRR_and_ERR_for_index = strsplit(SRR_for_index, split = "|",fixed = TRUE)

# Convert the list into pair character
SRR_and_ERR_for_index = do.call(rbind, SRR_and_ERR_for_index)
colnames(SRR_and_ERR_for_index) = c("Information","Index")
saveRDS(SRR_and_ERR_for_index, file = "SRRandERRTable.rds")
