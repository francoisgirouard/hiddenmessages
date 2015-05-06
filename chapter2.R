Skew <- function(Genome, i, PrecedingSkew = NULL) {
    if (i > 0) {
        value <- ifelse(is.null(PrecedingSkew), Skew(Genome, i - 1), PrecedingSkew)
        nucleotide <- substring(Genome, i, i)
        
        if (nucleotide == 'G') {
            value <- value + 1   
        } else if (nucleotide == 'C') {
            value <- value - 1
        }
    }
    
    return(value)
}

# Genome <- 'GAGCCACCGCGATA'
# Skews <- mapply(Skew, i = 0:nchar(Genome), MoreArgs = list(Genome = Genome))
# print(paste(as.character(Skews), collapse = ' '))
# print(length(Skews))

# Minimum Skew Problem: Find a position in a genome minimizing the skew.
# Input: A DNA string Genome.
# Output: All integer(s) i minimizing Skewi (Genome) among all values of i (from 0 to |Genome|).
MinimumSkew <- function(Genome) {
    skews <- numeric(nchar(Genome) + 1)
    skews[1] <- 0
    
    for (i in 1:length(skews)) {
        skews[i + 1] <- Skew(Genome, i, skews[i])
    }
    
    return(which(skews == min(skews)) - 1)
}

# fileName <- 'C:/Temp/dataset_7_6.txt'
# Genome <- readChar(fileName, file.info(fileName)$size)
# MinimumSkew(Genome)

# Hamming Distance Problem: Compute the g distance between two strings.
# Input: Two strings of equal length.
# Output: The g distance between these strings.
HammingDistance <- function(p, q) {
  first <- p
  second <- q
  
  if (length(first) == 1) {
    first <- strsplit(first, '')[[1]]
  }
  
  if (length(second) == 1) {
    second <- strsplit(second, '')[[1]]
  }
  
  return(sum(mapply(function(first, second) {
    return(first != second)
  }, first = first,
  second = second)))
}


# p <- 'ATCAAGCTAACGCTCGTTAGCATGTCATTTACGAGCCCGAGTCCTTAGAAATACACGGTCTATGCTGTGTTCGCACATTTTACATCCGTCCACGCGAGCCAACATGAGACGGATCGATTTATCCAGAGCAACTTGCTACCAATTGGCGTCGAGCTATGAAGTCGTGACGTACAATGTATGCAGAAGGGCACCGGGGGTATCCGAGCATACTATCTTGGGCACTGAATGAAGAACAGAACCTTGACCCGGAATCAACACACGCGGTAACGACGCAATCCTGTTCAAAAAAAGCCTGCCACAACAGTCTCAAACCACATACGGTGTAAGGACGTGTCGAATATGCAGACTGTAGATTGTAGGGCAGTACACTCTACAGTAACCGCTGGCACTATTACACACGTGCTCTACCGAGCTTGGATACACACCGTTAATCTCACTTAGGACTCGAAATGTCCCCTGGACTCTCGATGAGCTTAGGTAAGCACAGGTACCTAGGACATCTGCGCCTCCTAGGTGAGCCTCCCGAGAAAATGAGGGGTGATTTATGGGAAGTGAAGTACGGTGAAATAACTCGCAAGGGGGAAAACCGTTGAGACAGACTACCACCTTGACTGCTCCACGATGCGCGCACTTGACAGCCGCTCCTTCCGGACCCTGGACAGTTCTCACTCTGACACGTGATAAGTACGGCCATACGCAGTACGGAACTGCATAGGGCGGCATCGGTTTTCATGCATCGAGTATAGTTGTACCTATTGGCATACAATACTCCCTCGCTTAGCAATTGGCGTTAAATGTCGGGTCCACTGATGTGCACAAGGCTTACTCTCTACTTAGCTGAAAAATGTGCCTGTCCCTCTGCCATTTATACCCTTAAGAAACATCAACAGTTGGCCTGCGAGTCCTCGCCGAGCGCCGAGCCTACTTCTAGCTGCGAACGGTGTAGGCATTAAACAAAACGATGAACATCAATTGCAAGCCGAAATCCGCGGAGTTCGATTTTTTTGGATGCGGTAGATCCTACTCAACCCATATGACTGCAAAGTCCA'
# q <- 'ACCTATTGCAGCAACTCGTCCGCAGCCCCGGCTGTTTTCACACTAATGGATGCCACTCTGATGGACAGTTCTACGGCAGCAACCGATCAACTCCAACATTATTGGGCAGAGCAGAGACGATTCAGGAATGTACATTGAGCGGGACACGTTAGCTAGAATGAACGCAACTGGAGACATCCTCGAGTGCCGCTGCTAACTCGTAACTCCATGTGGTGTAAGATATTGGCCGCACGATGTCCCCTAACGGGAAAAAATGTTCCATCGGCATACCTCACCGATCCATCAAGCATAGCATGAGCAGCGAACGTCTGGAAGTTTGTCCTGTCACTGCTTCGTCTAAGAAGGGGTAGCGGGATTTGCTAATCTGAAGAACTTGTGCGCAAGGGTCCATGGTGACTATTAGGGAAAGCTAGGTCGTAAAATACACCTCGGGACAATTAATGGTGCGAACAAGGTTACCCCCGCTGTGGTGCTCTGAGGAGGAGCAGTCCGTAGGTATCACATCTCGCCTTCGCCGTCAAAAAGCGGTCATTCAGGAAGAGTCAGATTCGCCATGGGGCATCTACTACCTGAACACACCCTCGGACATGAGTCGCAAGGTGTGAAGGCCAGATCGTAGACCATTGCGGGTGTGACCACTTTGCGGAGGAGAGAATTCAGAAGAATGTGCCCAGATTCTGGCCGCTGAGCCGCCACTTAACAGCTCGAGAAATAAGCATGTGGGCTCTGCCCATGCTACTGCGCCTGCTTTTACACAAGCCATGCGGTAGTGCACTTCTCGAGATGAATTCGCGGAACCGTAGCGGTGTTCGGTCCCACTGCTAAATTACTTGACGTGTAGCCAATCGTTCAGGTTTTATCCGAGGGCTATAGTTAGAGAAGCCCGCGTACCGGTAGCTCCAAACCACGCGGCCAGCGGTCCTCATTTGTATAGTTGGATAAATTTACACGCCGTTCGAATTCGGTCGAGTTTGTGATAAGATGGCCAGTCAGGGAGGTAAAGCGGGGAAGTATAGCACTTTCCGACGAGATTGGTTCTGAAGAGGGTA'
# HammingDistance(p, q)

# Approximate Pattern Matching Problem: Find all approximate occurrences of a pattern in a string.
# Input: Strings Pattern and Text along with an integer d.
# Output: All starting positions where Pattern appears as a substring of Text with at most d mismatches.

ApproximatePatternMatching <- function(Pattern, Text, d) {
  HammingDistances <- numeric(nchar(Text) - nchar(Pattern) + 1)
  
  for (i in 1:length(HammingDistances)) {
    HammingDistances[i] <- HammingDistance(Pattern,
                                           substring(Text, i,
                                                     i + nchar(Pattern) - 1))
  }
  
  return(which(HammingDistances <= d) - 1)
}


# fileName <- 'dataset_9_4.txt'
# Text <- readChar(fileName, file.info(fileName)$size)
# result <- ApproximatePatternMatching('TGAAACGGAGA', Text, 6)
# write(result, ncolumns = length(result))
ApproximatePatternCount <- function(Text, Pattern, d = 0) {
  if (d == 0) {
    return(PatternCount(Text, Pattern))
  }
  
  count <- 0
  
  for (i in 1:(nchar(Text) - nchar(Pattern))) {
    if (HammingDistance(Pattern,
                        substring(Text, i, i + nchar(Pattern) - 1)) <= d) {
      count <- count + 1
    }
  }
  
  return(count)
}

ApproximatePatternCount('AACAAGCTGATAAACATTTAAAGAG', 'AAAAA', 2)
ApproximatePatternCount('TTTAGAGCCTTCAGAGG', 'GAGG', 2)

fileConn <- file('dataset_9_6.txt', 'rc')
Text <- trimws(readLines(fileConn, 1))
Pattern <- trimws(readLines(fileConn, 1))
d <- as.integer(trimws(readLines(fileConn, 1)))
close(fileConn)
ApproximatePatternCount(Text, Pattern, d)

# Frequent Words with Mismatches Problem: Find the most frequent k-mers with mismatches in a string.
# Input: A string Text as well as integers k and d. (You may assume k ≤ 12 and d ≤ 3.)
# Output: All most frequent k-mers with up to d mismatches in Text.
ComputingFrequenciesWithMismatches <- function(Text, k, d) {
  FrequencyArray <- integer(4 ^ k)
  alphabet <- c('A', 'C', 'G', 'T')
  
  for (i in 1:(nchar(Text) - k + 1)) {
    Pattern <- substring(Text, i, i + k - 1)
    splitPattern <- strsplit(Pattern, '')[[1]]
    mismatches <- splitPattern
    
    for (distance in 1:d) {
      mismatchesTemplate <- matrix(rep(splitPattern, times = 3 ^ distance), 3 ^ distance, k, TRUE)
      substitutionCandidates <- character(4 ^ distance)
      
      for (index in 1:distance) {
        substitutionCandidates <- cbind(substitutionCandidates,
                                        rep(alphabet,
                                            each = 4 ^ (distance - index),
                                            times = 4 ^ (index - 1)))
      }
      
      substitutionCandidates <- substitutionCandidates[, -1, drop = FALSE]
      
      for (indices in combn(1:k, distance, simplify = FALSE)) {
        currentMismatches <- mismatchesTemplate
        substitutes <- substitutionCandidates
        
        for (index in seq_along(indices)) {
          substitutes <- substitutes[substitutes[, index] != splitPattern[indices[index]],] 
        }
        
        currentMismatches[, indices] <- substitutes
        mismatches <- rbind(mismatches, currentMismatches)
      }
    }
    
    mismatches <- apply(mismatches, 1, function(row) {return(paste(row, collapse = ''))})
    mismatches <- as.numeric(lapply(mismatches, PatternToNumber))
    FrequencyArray[mismatches + 1] <- FrequencyArray[mismatches + 1] + 1
  }
  
  return(FrequencyArray)
}

FrequentWordsWithMismatches <- function(Text, k, d) {
  FrequencyArray <- ComputingFrequenciesWithMismatches(Text, k, d)
  return(sapply(which(FrequencyArray == max(FrequencyArray)) - 1, NumberToPattern, k = k))
}

# fileConn <- file('sample.txt', 'rc')
# Text <- readLines(fileConn, 1)
# otherParams <- as.integer(strsplit(readLines(fileConn, 1), ' ')[[1]])
# FrequentWordsWithMismatches(Text, otherParams[1], otherParams[2])
# close(fileConn)

# Frequent Words with Mismatches and Reverse Complements Problem: Find the most frequent k-mers (with mismatches and reverse complements) in a DNA string.
# Input: A DNA string Text as well as integers k and d.
# Output: All k-mers Pattern maximizing the sum Countd(Text, Pattern)+ Countd(Text, Pattern)
# over all possible k-mers.
FrequentWordsWithMismatchesAndReverseComplements <- function(Text, k, d) {
  FrequencyArray <- ComputingFrequenciesWithMismatches(Text, k, d)
  FrequencyArray <- FrequencyArray + ComputingFrequenciesWithMismatches(ReverseComplement(Text), k, d)
  return(sapply(which(FrequencyArray == max(FrequencyArray)) - 1, NumberToPattern, k = k))
}

# fileConn <- file('sample.txt', 'rc')
# Text <- readLines(fileConn, 1)
# otherParams <- as.integer(strsplit(readLines(fileConn, 1), ' ')[[1]])
# FrequentWordsWithMismatchesAndReverseComplements(Text, otherParams[1], otherParams[2])
# close(fileConn)

# fileConn <- file('Salmonella_enterica.txt', 'rc')
# readLines(fileConn, 1)
# Genome <- paste(readLines(fileConn), collapse = '')
# close(fileConn)
# MinimumSkew(Genome)
