# Implanted Motif Problem: Find all (k, d)-motifs in a collection of strings.
# Input: A collection of strings Dna, and integers k and d.
# Output: All (k, d)-motifs in Dna.

# ImmediateNeighbors(Pattern)
#   Neighborhood <- set consisting of single string Pattern
#   for i = 1 to |Pattern|
#       symbol <- i-th nucleotide of Pattern
#       for each nucleotide x different from symbol
#           Neighbor <- Pattern with the i-th nucleotide substituted by x
#           add Neighbor to Neighborhood
#   return Neighborhood
ImmediateNeighbors <- function(Pattern) {
    Neighborhood <- character(0)
    
    for (i in 1:nchar(Pattern)) {
        symbol <- substring(Pattern, i, i)
        
        for (nucleotide in setdiff(c('A', 'C', 'G', 'T'), symbol)) {
            if (i > 1) {
                Neighbor <- substring(Pattern, 1, i - 1)
            } else {
                Neighbor <- ''
            }
            
            Neighbor <- paste(Neighbor, nucleotide, sep = '')
            
            if (i < nchar(Pattern)) {
                Neighbor <- paste(Neighbor, substring(Pattern, i + 1, nchar(Pattern)), sep = '')
            }
            
            Neighborhood <- c(Neighborhood, Neighbor)
        }
    }
    
    return(Neighborhood)
}

# Neighbors(Pattern, d)
#   if d = 0
#       return {Pattern}
#   if |Pattern| = 1 
#       return {A, C, G, T}
#   Neighborhood <- an empty set
#   SuffixNeighbors <- Neighbors(Suffix(Pattern), d)
#   for each string Text from SuffixNeighbors
#       if HammingDistance(Suffix(Pattern), Text) < d
#           for each nucleotide x
#               add x • Text to Neighborhood
#       else
#           add FirstSymbol(Pattern) • Text to Neighborhood
#   return Neighborhood
Suffix <- function(Pattern) {
    return(substring(Pattern, 2, nchar(Pattern)))
}

FirstSymbol <- function(Pattern) {
    return(substring(Pattern, 1, 1))
}

Neighbors <- function(Pattern, d) {
    if (d == 0) {
        Neighborhood <- Pattern
    } else if (nchar(Pattern) == 1 ) {
        Neighborhood <- c('A', 'C', 'G', 'T')
    } else {
        Neighborhood <- character(0)
        SuffixNeighbors <- Neighbors(Suffix(Pattern), d)
        
        for (Text in SuffixNeighbors) {
            if (HammingDistance(Suffix(Pattern), Text) < d) {
                for (nucleotide in c('A', 'C', 'G', 'T')) {
                    Neighborhood <- c(Neighborhood, paste(nucleotide, Text, sep = ''))
                }
            } else {
                Neighborhood <- c(Neighborhood, paste(FirstSymbol(Pattern), Text, sep = ''))
            }
        }
    }
    
    return(Neighborhood)
}

# result <- Neighbors('AGCCCGAAA', 3)
# write(result, ncolumns = 1)

# IterativeNeighbors(Pattern, d)
#   Neighborhood <- set consisting of single string Pattern
#   for j = 1 to d
#       for each string Pattern’ in Neighborhood
#           add ImmediateNeighbors(Pattern') to Neighborhood
#           remove duplicates from Neighborhood
#   return Neighborhood
IterativeNeighbors <- function(Pattern, d) {
    Neighborhood <- Pattern
    for (j in 1:d) {
        for (PatternPrime in Neighborhood) {
            Neighborhood <- c(Neighborhood, ImmediateNeighbors(PatternPrime))
            Neighborhood <- unique(Neighborhood)
        }
    }
    return(Neighborhood)
}

# MOTIFENUMERATION(Dna, k, d)
#   Patterns <- an empty set
#   for each k-mer Pattern in Dna
#       for each k-mer Pattern’ differing from Pattern by at most d mismatches
#           if Pattern' appears in each string from Dna with at most d mismatches
#               add Pattern' to Patterns
#   remove duplicates from Patterns
#   return Patterns
MotifEnumeration <- function(Dna, k, d) {
    Patterns <- character(0)
    
    for (string1 in Dna) {
        for (Pattern in sapply(1:(nchar(string1) - k + 1),
                               function(X, Text, k) {
                                   return(substring(Text, X, X + k - 1))
                               }, Tex
                               = string1, k = k)) {
            
            for (PatternPrime in Neighbors(Pattern, d)) {
                include <- TRUE
                
                for (string2 in Dna) {
                    include <- include & (ApproximatePatternCount(string2, PatternPrime, d) > 0)
                }
                
                if (include) {
                    Patterns <- c(Patterns, PatternPrime)   
                }
            }
        }
    }
    
    Patterns <- unique(Patterns)
    return(Patterns)
}

Dna <- c('AACTTGATCACATAGTCGCTGTACC',
         'CATCCCATCGGGGTAATCACCCGTG',
         'AGATATCAAACACCGCTCCATCCTC',
         'CAGAGGCGAGCAGGGTATTTGGATC',
         'CAACCCCCGACATGTCACCGGGTAA',
         'GCAACTCGGACGTCGCAAAGCCGGG')
# MotifEnumeration(Dna, 5, 2)


Motifs <- t(sapply(c('TCGGGGgTTTtt',
                     'cCGGtGAcTTaC',
                     'aCGGGGATTTtC',
                     'TtGGGGAcTTtt',
                     'aaGGGGAcTTCC',
                     'TtGGGGAcTTCC',
                     'TCGGGGATTcat',
                     'TCGGGGATTcCt',
                     'TaGGGGAacTaC',
                     'TCGGGtATaaCC'),
                   function(string) {
                       return(strsplit(string, '')[[1]])
                       }))

Motifs <- toupper(Motifs)

Score <- function(Motifs) {
    return(sum(apply(Motifs, 2, function(column) {
        return(sum(column != names(which.max(table(column)))))
        })))
}

Count <- function(Motifs) {
    return(t(sapply(c('A', 'C', 'G', 'T'), function(nucleotide, Motifs) {
        return(apply(Motifs, 2, function(column, nucleotide) {
            return(sum(column == nucleotide))
        }, nucleotide = nucleotide))
    }, Motifs = Motifs)))
}

Profile <- function(Motifs) {
    return(Count(Motifs) / nrow(Motifs))
}

Consensus <- function(Motifs) {
    return(paste(apply(Profile(Motifs), 2, function(column) {
        return(names(column)[which(column == max(column))])
    }), collapse = ''))
}

ColumnEntropy <- function(Motifs) {
    return(apply(Profile(Motifs), 2, function(column) {
        return(-sum(sapply(column, function(probability) {
            if (probability == 0) {
                return(0)
            } else {
                return(probability * log2(probability))
            }
        })))}))
}

Entropy <- function(Motifs) {
    return(sum(ColumnEntropy(Motifs)))    
}

# Score(Motifs)
# Count(Motifs)
# Profile(Motifs)
# Consensus(Motifs)
# ColumnEntropy(Motifs)
# Entropy(Motifs)

# Motif Finding Problem: Given a collection of strings, find a set of k-mers, one from each string, that minimizes the
#                        score of the resulting motif.
# Input: A collection of strings Dna and an integer k.
# Output: A collection Motifs of k-mers, one from each string in Dna, minimizing Score(Motifs) among all possible
#         choices of k-mers.
dMotifs <- function(Pattern, Motifs) {
    return(sum(apply(Motifs, 1, HammingDistance, q = Pattern)))
}

# d(Consensus(Motifs), Motifs)

dText <- function(Pattern, Text) {
  return(min(sapply(sapply(1:nchar(Pattern),
                           function(first,
                                    text,
                                    length) {
                             return(substring(text,
                                              first,
                                              first + length - 1))
                           },
                           text = text,
                           length = nchar(Pattern)),
                    HammingDistance,
                    q = Pattern)))
}

Motif <- function(Pattern, Text) {
  HammingDistances <- sapply(sapply(1:nchar(Pattern),
                                    function(first,
                                             text,
                                             length) {
                                      return(substring(text,
                                                       first,
                                                       first + length - 1))
                                    },
                                    text = text,
                                    length = nchar(Pattern)),
                             HammingDistance,
                             q = Pattern)
  return(names(HammingDistances[which.min(HammingDistances)]))  
}

dDna <- function(Pattern, Dna) {
  return(sum(sapply(Dna, function(Text, Pattern) {
    return(sum(dText(Pattern, Text)))
  })))
}

Motifs <- function(Pattern, Dna) {
  return(sapply(Dna, function(text, Pattern) {
    HammingDistances <- sapply(sapply(1:nchar(Pattern),
                                      function(first,
                                               text,
                                               length) {
                                        return(substring(text,
                                                         first,
                                                         first + length - 1))
                                      },
                                      text = text,
                                      length = nchar(Pattern)),
                               HammingDistance,
                               q = Pattern)
    return(strsplit(names(HammingDistances[which.min(HammingDistances)]),
                    '')[[1]])
  },
  Pattern = Pattern, USE.NAMES = FALSE))
}

# MEDIANSTRING(Dna, k)
#   distance <- ∞
#   for each k-mer Pattern from AA…AA to TT…TT
#     if distance > d(Pattern, Dna)
#       distance <- d(Pattern, Dna)
#       Median <- Pattern
#   return Median
MedianString <- function(Dna, k)
  distance <- Inf
  for each k-mer Pattern from AA…AA to TT…TT
    if distance > d(Pattern, Dna)
      distance <- d(Pattern, Dna)
      Median <- Pattern
  return Median