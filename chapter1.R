# PatternCount(Text, Pattern)
#   count ← 0
#   for i ← 0 to |Text| − |Pattern|
#     if Text(i, |Pattern|) = Pattern
#       count ← count + 1
#   return count
PatternCount <- function(Text, Pattern) {
  count <- 0
  
  for (i in 1:(nchar(Text) - nchar(Pattern))) {
    if (substring(Text, i, i + nchar(Pattern) - 1) == Pattern) {
      count <- count + 1
    }
  }
  
  return(count)
}

# Text <- 'AGGTAGCTGCTTCTCGAGGTAGCTTTACGCTAGGTAGCCCATAGACCAGGTAGCTCAGGGAGGTAGCAGGTAGCCAAGGTAGCACGTCTCAGGTAGCGAGGTAGCGATAGGTAGCAGGTAGCGTAGGTAGCCAGGTAGCTAGGTAGCACACGTAGGTAGCTAGGTAGCATAAAATCCATCAGGTAGCCTGAGGTAGCCAGGTAGCAGGTAGCAGGTAGCTTCTGAGGTAGCTGAATGAGGTAGCAAGTGAGGTAGCTGGAGGTAGCAAGGTAGCAGGTAGCAGGTAGCTAAGGTAGCATAGGTAGCCCCAAAGGTAGCTAGCAGAGGTAGCACCGGGGGGTTTGTTAGGTAGCAAGCAGGTAGCTAGGTAGCAGGTAGCAGGTAGCTGCAGGTAGCGTAGGTAGCAGGTAGCATGTAAGGTAGCCTGAGGTAGCAGGTAGCTTGCGGAGGTAGCAGGTAGCAGGTAGCAGGTAGCAGGTAGCGAAGGTAGCATGAAATATGTCGTTGAAGGTAGCGGAAGGTAGCAATGCCAGGTAGCGCGAGGTAGCGAGGTAGCGGGAGGTAGCTAGGTAGCAGGTAGCAAATAGGTAGCAAGGTAGCTAGGTAGCTAGGTAGCCACTAGGTAGCAAGGTAGCCATTTTCAGGTAGCGGGAGGTAGCTGAACAAGGTAGCAGGTAGCAGGTAGCAGGTAGCGAGGTAGCAAGGTAGCGTAGGTAGCTCAGGTAGCGGTAGGTAGCGGATAAGGTAGCAGGTAGCTAGGTAGCTAGGTAGCGAGGTAGCCAAACAAGAGGTAGCTAGGTAGCTAGGAGGTAGCTAGGTAGCTTCTAAAGGTAGCAGGTAGCAGGTAGCAGGTAGCGAGGTAGCTATGAGGTAGCTAGGTAGCGCATCCACGTGATGAGGTAGCTAGGTAGCACTAGGTAGCTGTAAGGTAGCGATGAGGTAGCCATCGTGGAAGGTAGCAAAAGGTAGCAAGGTAGCCTCAGAGGTAGCGTAGGTAGCACGCAGGTAGCA'
# Pattern <- 'AGGTAGCAG'
# print(PatternCount(Text, Pattern))

# FrequentWords(Text, k)
#   FrequentPatterns ← an empty set
#   for i ← 0 to |Text| − k
#     Pattern ← the k-mer Text(i, k)
#     Count(i) ← PatternCount(Text, Pattern)
#   maxCount ← maximum value in array Count
#   for i ← 0 to |Text| − k
#     if Count(i) = maxCount
#       add Text(i, k) to FrequentPatterns
#   remove duplicates from FrequentPatterns
#   return FrequentPatterns
FrequentWords <- function(Text, k) {
  FrequentPatterns <- character()
  Count <- integer(nchar(Text) - k + 1)
  
  for (i in 1:length(Count)) {
    Pattern <- substring(Text, i, i + k - 1)
    Count[i] <- PatternCount(Text, Pattern)
  }
  
  maxCount <- max(Count)
  
  for (i in 1:length(Count)) {
    if (Count[i] == maxCount) {
      FrequentPatterns <- c(FrequentPatterns, substring(Text, i, i + k - 1))
    }
  }
  
  FrequentPatterns <- unique(FrequentPatterns)
  return(FrequentPatterns)
}

# Text <- 'GGATACTCCACGGCCATGGGGATACTCCCCGCCCTGAATCGTTGAGCTCCGCCCTGAATCGTTGAGCTACGATTGATTACGGCCATGGACGATTGATTACGGCCATGGACGATTGATTGGATACTCCCCGCCCTGAACCGCCCTGAAACGGCCATGGCCGCCCTGAAACGATTGATTGGATACTCCACGGCCATGGGGATACTCCGGATACTCCTCGTTGAGCTTCGTTGAGCTACGATTGATTACGGCCATGGTCGTTGAGCTTCGTTGAGCTACGGCCATGGTCGTTGAGCTGGATACTCCGGATACTCCGGATACTCCACGATTGATTACGATTGATTGGATACTCCGGATACTCCGGATACTCCCCGCCCTGAAACGATTGATTACGGCCATGGACGGCCATGGGGATACTCCTCGTTGAGCTGGATACTCCTCGTTGAGCTGGATACTCCCCGCCCTGAATCGTTGAGCTGGATACTCCACGATTGATTCCGCCCTGAACCGCCCTGAAACGGCCATGGACGGCCATGGACGATTGATTACGATTGATTGGATACTCCACGGCCATGGCCGCCCTGAAACGATTGATTGGATACTCCTCGTTGAGCTCCGCCCTGAAACGATTGATTTCGTTGAGCTCCGCCCTGAAGGATACTCCCCGCCCTGAATCGTTGAGCTCCGCCCTGAACCGCCCTGAAACGGCCATGGCCGCCCTGAAACGATTGATTACGGCCATGGACGATTGATTCCGCCCTGAATCGTTGAGCTCCGCCCTGAAACGGCCATGGTCGTTGAGCTTCGTTGAGCTACGGCCATGGCCGCCCTGAATCGTTGAGCTGGATACTCCACGATTGATTACGGCCATGGGGATACTCCCCGCCCTGAAGGATACTCCACGGCCATGGCCGCCCTGAATCGTTGAGCTCCGCCCTGAAACGGCCATGGCCGCCCTGAAACGGCCATGGCCGCCCTGAACCGCCCTGAA'
# k <- 13
# print(paste(FrequentWords(Text, k), collapse = ' '))

# AAGCAAAGGTGGG
PatternToNumber <- function(Pattern) {
  tokens <- rev(strsplit(Pattern, NULL)[[1]])
  Number <- sum(4 ^ (grep('C', tokens) - 1))
  Number <- Number + sum(2 * 4 ^ (grep('G', tokens) - 1))
  Number <- Number + sum(3 * 4 ^ (grep('T', tokens) - 1))
  return(Number)
}

# PatternToNumber('ATGCAA')

NumberToPattern <- function(Number, k) {
  Pattern <- 'Invalid parameters'
  if (Number >= 0 & Number < 4 ^ k) {
    powers <- k:1
    Pattern <- as.character((Number %% 4 ^ powers) %/% (4 ^ (powers - 1)))
    Pattern[Pattern == '0'] <- 'A'
    Pattern[Pattern == '1'] <- 'C'
    Pattern[Pattern == '2'] <- 'G'
    Pattern[Pattern == '3'] <- 'T'
    Pattern <- paste(Pattern, collapse = '')
  }
  return(Pattern)
}

# NumberToPattern(2362090, nchar('AAGCAAAGGTGGG'))

# ComputingFrequencies(Text , k)
#   for i <- 0 to 4k − 1
#       FrequencyArray(i) <- 0
#   for i <- 0 to |Text| − k
#       Pattern <- Text(i, k)
#       j <- PatternToNumber(Pattern)
#       FrequencyArray(j) <- FrequencyArray(j) + 1
#   return FrequencyArray
ComputingFrequencies <- function(Text , k) {
  FrequencyArray <- integer(4 ^ k)
  
  for (i in 1:(nchar(Text) - k + 1)) {
    Pattern <- substring(Text, i, i + k - 1)
    j <- PatternToNumber(Pattern)
    FrequencyArray[j + 1] <- FrequencyArray[j + 1] + 1
  }
  
  return(FrequencyArray)
}

# result <- ComputingFrequencies('CAAAGACAGCACGCTAGGGTGGATCACAGCAAGGAGTGTTGCAATTGAACCTCGGTTACATCGAAGACTTGACGCCGTGCTGATTCCCGATTGTGATTGTCTAGAAACCGGCAGAAAGTAGGGTTTATGATATAGCAAGAGCTCATTCTATAGCCTCTGGACATATCACCGACGGTGCCGTGATAAGCCCGTTGTCCTTCCACAGTAATTCAATGATTGCGAGGGAGTTATGTAAAGGTGCCATGTGTCTTACGAACTGATGCTGAGCTGAAGGGTAGATGCAATCCGCATGTCTCTAGGTGTCGGCAGGAACTCCAGCCAGTTATTATGACTAGTCCCCGAGGCTCATAGTGGTAGTCGACTAGTAATCTACAAGGTAAAACCAAATTCTATCTGCTCTATTACAGCCGCGTTTAGTCCGACTTTCTCCCGGACGCCGCTAGTGGGCATGCAGTCTCAGGAAGGTAATAGACCTGCAGCGGGACCGTACAATTTGCATCCATTGTTCTGGATGAGGGAGCGCGTCCTAGCGCAGGGATCCGCAATACAAGGATCACTCGGCGTGGTACAGGCCTCAGGGGTACCAGCTGTCCGAAGTTGAGCACCTGAGTGGCAGCCAATTGAACGAGACCAATAGCTACGGTCTGTTCTGACC', 6)
# write(result, ncolumns = length(result))

# FasterFrequentWords(Text , k)
#   FrequentPatterns <- an empty set
#   FrequencyArray <- ComputingFrequencies(Text, k)
#   maxCount <- maximal value in FrequencyArray
#   for i <-0 to 4k − 1
#       if FrequencyArray(i) = maxCount
#           Pattern <- NumberToPattern(i, k)
#           add Pattern to the set FrequentPatterns
#   return FrequentPatterns

FasterFrequentWords <- function(Text , k) {
  FrequencyArray <- ComputingFrequencies(Text, k)
  return(sapply(which(FrequencyArray == max(FrequencyArray)) - 1, NumberToPattern, k = k))
}

# FasterFrequentWords(toupper('atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaacctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgaccacggaaagatgatcaagagaggatgatttcttggccatatcgcaatgaatacttgtgacttgtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggattacgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttaggatagacggtttttcatcactgactagccaaagccttactctgcctgacatcgaccgtaaattgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaagatcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtttccttaaccctctattttttacggaagaatgatcaagctgctgctcttgatcatcgtttc'), 9)

# Reverse Complement Problem: Find the reverse complement of a DNA string.
# Input: A DNA string Pattern.
# Output: Pattern, the reverse complement of Pattern.
ReverseComplement <- function(Pattern) {
  Reverse <- rev(strsplit(Pattern, NULL)[[1]])
  Reverse <- ifelse(Reverse == 'A', 'T',
                    ifelse(Reverse == 'C', 'G',
                           ifelse(Reverse == 'G', 'C',
                                  'A')))
  Reverse <- paste(Reverse, collapse = '')
  return(Reverse)
}

# ReverseComplement('TCTTGATCA')

# Pattern Matching Problem: Find all occurrences of a pattern in a string.
# Input: Strings Pattern and Genome.
# Output: All starting positions in Genome where Pattern appears as a substring.
PatternMatching <- function(Pattern, Genome) {
  Matches <- integer()
  k <- nchar(Pattern)
  
  for (i in 1:(nchar(Genome) - k + 1)) {
    if (substring(Genome, i, i + k - 1) == Pattern) {
      Matches <- c(Matches, i - 1) 
    } 
  }
  
  return(Matches)
}

# fileName <- 'Thermotoga-petrophila.txt'
# Genome <- readChar(fileName, file.info(fileName)$size)
# result <- PatternMatching('CTTGATCAT', )
# write(result, ncolumns = length(result))

# Clump Finding Problem: Find patterns forming clumps in a string.
# Input: A string Genome, and integers k, L, and t.
# Output: All distinct k-mers forming (L, t)-clumps in Genome.
ClumpFinding <- function(Genome, k, L, t) {
  candidates <- sapply(which(ComputingFrequencies(Genome, k) >= t) - 1,
                       NumberToPattern, k = k)
  Patterns <- character()
  
  for (candidate in candidates) {
    indices <- gregexpr(candidate, Genome)[[1]]
    
    for (i in 1:(length(indices) - t + 1)) {
      if (indices[i + t - 1] <= indices[i] + L - k) {
        Patterns <- c(Patterns, candidate)
        break
      }
    }
  }
  
  return(Patterns)
}

# BetterClumpFinding(Genome, k, t, L)
#   FrequentPatterns ← an empty set
#   for i ←0 to 4k − 1
#     Clump(i) ← 0
#   Text ← Genome(0, L)
#   FrequencyArray ← ComputingFrequencies(Text, k)
#   for i ← 0 to 4k − 1
#     if FrequencyArray(i) ≥ t
#       Clump(i) ← 1
#   for i ← 1 to |Genome| − L
#     FirstPattern ← Genome(i − 1, k)
#     j ← PatternToNumber(FirstPattern)
#     FrequencyArray(j) ← FrequencyArray(j) − 1
#     LastPattern ← Genome(i + L − k, k)
#     j ← PatternToNumber(LastPattern)
#     FrequencyArray(j) ← FrequencyArray(j) + 1
#     if FrequencyArray(j) ≥ t
#       Clump(j) ← 1
#   for i ← 0 to 4k − 1
#     if Clump(i) = 1
#       Pattern ← NumberToPattern(i, k)
#       add Pattern to the set FrequentPatterns
#   return FrequentPatterns
BetterClumpFinding <- function(Genome, k, t, L) {
  FrequentPatterns <- character()
  Clump <- logical(4 ^ k)
  Text <- substring(Genome, 1, L)
  FrequencyArray <- ComputingFrequencies(Text, k)
  
  for (i in 1:(4 ^ k)) {
    if (FrequencyArray[i] >= t) {
      Clump[i] <- TRUE  
    }
  }

  for (i in 2:(nchar(Genome) - L + 1)) {
    print(i)
    FirstPattern <- substring(Genome, i - 1, i + k - 2)
    j <- PatternToNumber(FirstPattern)
    FrequencyArray[j] <- FrequencyArray[j] - 1
    LastPattern <- substring(Genome, i + L - k, i + L -1)
    j <- PatternToNumber(LastPattern)
    FrequencyArray[j] <- FrequencyArray[j] + 1
    
    if (FrequencyArray[j] >= t) {
      Clump[j] <- TRUE
    }
  }
  
  for (i in 1:(4 ^ k)) {
    if (Clump[i]) {
      Pattern <- NumberToPattern(i, k)
      FrequentPatterns <- c(FrequentPatterns, Pattern)
    }
  }
  
  return(FrequentPatterns)
}

fileName <- 'E-coli.txt'
Genome <- readChar(fileName, file.info(fileName)$size)
print(paste(BetterClumpFinding(Genome, 9, 500, 3), collapse = ' '))