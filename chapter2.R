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

# Hamming Distance Problem: Compute the Hamming distance between two strings.
# Input: Two strings of equal length.
# Output: The Hamming distance between these strings.
HammingDistance <- function(p, q) {
    return(sum(mapply(function(first, second) {return(first != second)},
                      first = strsplit(p, ''),
                      second = strsplit(q, ''))))
}

# p <- 'ATCAAGCTAACGCTCGTTAGCATGTCATTTACGAGCCCGAGTCCTTAGAAATACACGGTCTATGCTGTGTTCGCACATTTTACATCCGTCCACGCGAGCCAACATGAGACGGATCGATTTATCCAGAGCAACTTGCTACCAATTGGCGTCGAGCTATGAAGTCGTGACGTACAATGTATGCAGAAGGGCACCGGGGGTATCCGAGCATACTATCTTGGGCACTGAATGAAGAACAGAACCTTGACCCGGAATCAACACACGCGGTAACGACGCAATCCTGTTCAAAAAAAGCCTGCCACAACAGTCTCAAACCACATACGGTGTAAGGACGTGTCGAATATGCAGACTGTAGATTGTAGGGCAGTACACTCTACAGTAACCGCTGGCACTATTACACACGTGCTCTACCGAGCTTGGATACACACCGTTAATCTCACTTAGGACTCGAAATGTCCCCTGGACTCTCGATGAGCTTAGGTAAGCACAGGTACCTAGGACATCTGCGCCTCCTAGGTGAGCCTCCCGAGAAAATGAGGGGTGATTTATGGGAAGTGAAGTACGGTGAAATAACTCGCAAGGGGGAAAACCGTTGAGACAGACTACCACCTTGACTGCTCCACGATGCGCGCACTTGACAGCCGCTCCTTCCGGACCCTGGACAGTTCTCACTCTGACACGTGATAAGTACGGCCATACGCAGTACGGAACTGCATAGGGCGGCATCGGTTTTCATGCATCGAGTATAGTTGTACCTATTGGCATACAATACTCCCTCGCTTAGCAATTGGCGTTAAATGTCGGGTCCACTGATGTGCACAAGGCTTACTCTCTACTTAGCTGAAAAATGTGCCTGTCCCTCTGCCATTTATACCCTTAAGAAACATCAACAGTTGGCCTGCGAGTCCTCGCCGAGCGCCGAGCCTACTTCTAGCTGCGAACGGTGTAGGCATTAAACAAAACGATGAACATCAATTGCAAGCCGAAATCCGCGGAGTTCGATTTTTTTGGATGCGGTAGATCCTACTCAACCCATATGACTGCAAAGTCCA'
# q <- 'ACCTATTGCAGCAACTCGTCCGCAGCCCCGGCTGTTTTCACACTAATGGATGCCACTCTGATGGACAGTTCTACGGCAGCAACCGATCAACTCCAACATTATTGGGCAGAGCAGAGACGATTCAGGAATGTACATTGAGCGGGACACGTTAGCTAGAATGAACGCAACTGGAGACATCCTCGAGTGCCGCTGCTAACTCGTAACTCCATGTGGTGTAAGATATTGGCCGCACGATGTCCCCTAACGGGAAAAAATGTTCCATCGGCATACCTCACCGATCCATCAAGCATAGCATGAGCAGCGAACGTCTGGAAGTTTGTCCTGTCACTGCTTCGTCTAAGAAGGGGTAGCGGGATTTGCTAATCTGAAGAACTTGTGCGCAAGGGTCCATGGTGACTATTAGGGAAAGCTAGGTCGTAAAATACACCTCGGGACAATTAATGGTGCGAACAAGGTTACCCCCGCTGTGGTGCTCTGAGGAGGAGCAGTCCGTAGGTATCACATCTCGCCTTCGCCGTCAAAAAGCGGTCATTCAGGAAGAGTCAGATTCGCCATGGGGCATCTACTACCTGAACACACCCTCGGACATGAGTCGCAAGGTGTGAAGGCCAGATCGTAGACCATTGCGGGTGTGACCACTTTGCGGAGGAGAGAATTCAGAAGAATGTGCCCAGATTCTGGCCGCTGAGCCGCCACTTAACAGCTCGAGAAATAAGCATGTGGGCTCTGCCCATGCTACTGCGCCTGCTTTTACACAAGCCATGCGGTAGTGCACTTCTCGAGATGAATTCGCGGAACCGTAGCGGTGTTCGGTCCCACTGCTAAATTACTTGACGTGTAGCCAATCGTTCAGGTTTTATCCGAGGGCTATAGTTAGAGAAGCCCGCGTACCGGTAGCTCCAAACCACGCGGCCAGCGGTCCTCATTTGTATAGTTGGATAAATTTACACGCCGTTCGAATTCGGTCGAGTTTGTGATAAGATGGCCAGTCAGGGAGGTAAAGCGGGGAAGTATAGCACTTTCCGACGAGATTGGTTCTGAAGAGGGTA'
# HammingDistance(p, q)

# Approximate Pattern Matching Problem: Find all approximate occurrences of a pattern in a string.
# Input: Strings Pattern and Text along with an integer d.
# Output: All starting positions where Pattern appears as a substring of Text with at most d mismatches.

ApproximatePatternMatching <- function(Pattern, Text, d) {
  hammingDistances <- numeric(nchar(Text) - nchar(Pattern) + 1)
  
  for (i in 1:length(hammingDistances)) {
    hammingDistances[i] <- HammingDistance(Pattern,
                                           substring(Text, i,
                                                     i + nchar(Pattern) - 1))
  }
  
  return(which(hammingDistances <= d) - 1)
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
FrequentWordsWithMismatches <- function(Text, k, d) {
    
}

k <- 12
possibles <- character(4^k)
for (i in k:1) {
  possibles <- cbind(possibles, rep(c('A', 'C', 'G', 'T'),
                                    each = 4^(k-(k-i+1)),
                                    times = 4^(k-i)))
}
rep(c('A', 'C', 'G', 'T'), each = 4^(k-1), times = 4^(k-3))
rep(c('A', 'C', 'G', 'T'), each = 4^(k-2), times = 4^(3-2))
rep(c('A', 'C', 'G', 'T'), each = 4^(k-3), times = 4^(3-1))
