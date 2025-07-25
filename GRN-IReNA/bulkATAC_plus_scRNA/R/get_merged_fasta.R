# Date: 250715
get_merged_fasta <- function(footprints, fastadir, distance = 4) {
  # validInput(footprints,'footprints','df') # added by yd
  # validInput(fastadir,'fastadir','direxists') # added by yd
  # validInput(distance,'distance','numeric') # added by yd
  merged_footprints <- merge_footprints(footprints, ditance = distance)
  fasta <- getfasta(merged_footprints, fastadir = fastadir)
  return(fasta)
}


merge_footprints <- function(footprints, ditance = 4, revise = TRUE) {
  con1 <- footprints
  str1 <- con1[1,1]
  start1 <- con1[1,2]
  end1 <- con1[1,3]
  pva1 <- con1[1,4]
  no1 <- 1
  out1 <- c()
  for (i in 2:nrow(con1)) {
    str2 <- as.character(con1[i, ][1])
    start2 <- as.numeric(con1[i, ][2])
    end2 <- as.numeric(con1[i, ][3])
    pva2 <- as.numeric(con1[i, ][4])
    if (str2 == str1 & (start1 >= start2 & start1 <= end2 | end1 >= start2 &
                        end1 <= end2 | start1 <= start2 & end1 >= end2 | start2 >
                        end1 & (start2 - end1) <= ditance | start1 > end2 &
                        (start1 - end2) <= ditance)) {
      acc22 <- sort(c(start1, start2, end1, end2))
      start1 <- acc22[1]
      end1 <- acc22[4]
      pva1 <- pva1 + pva2
      no1 <- no1 + 1
    } else {
      pva2 <- pva1 / no1
      col1 <- paste(str1, start1, end1, pva2, sep = "\t")
      out1 <- c(out1, col1)
      str1 <- str2
      start1 <- start2
      end1 <- end2
      pva1 <- as.numeric(con1[i, ][4])
      no1 <- 1
    }
  }
  out2 <- as.data.frame(t(as.data.frame(strsplit(out1, "\t"))))
  colnames(out2) <- c("chr", "strat", "stop", "pvalue")
  out2[, 2] <- as.numeric(out2[, 2])
  out2[, 3] <- as.numeric(out2[, 3])
  if (revise == TRUE) {
    out2[, 2] <- out2[, 2] - 2
    out2[, 3] <- out2[, 3] + 2
    return(out2)
  } else {
    return(out2)
  }
}


getfasta <- function(merged_footprints, fastadir) {
  fasta <- Biostrings::readBStringSet(fastadir, format = "fasta", nrec = -1L,
                                      skip = 0L, seek.first.rec = FALSE,
                                      use.names = TRUE)
  fasta1 <- c()
  for (i in 1:nrow(merged_footprints)) {
    name <- paste0(">", merged_footprints[i, 1], ":", merged_footprints[i, 2],
                   "-", merged_footprints[i, 3])
    if (merged_footprints[i, 3] > length(fasta[[merged_footprints[i, 1]]])) {
      sequence <- toupper(as.character(fasta[[merged_footprints[i, 1]]][
        (merged_footprints[i, 2] + 1):length(fasta[[merged_footprints[i, 1]]])]))
      name <- paste0(">", merged_footprints[i, 1], ":", merged_footprints[i, 2],
                     "-", length(fasta[[merged_footprints[i, 1]]]))
    } else{
      sequence <- toupper(as.character(fasta[[merged_footprints[i, 1]]][
        (merged_footprints[i, 2] + 1):merged_footprints[i, 3]]))
    }
    fasta1 <- c(fasta1, name, sequence)
  }
  fasta1 <- as.data.frame(fasta1)
  return(fasta1)
}