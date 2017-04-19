def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count += 1
    return count

def PatternMatching(Pattern, Genome):
    indices = [] # positions
    for i in range(len(Genome)-len(Pattern)+1):
        if Genome[i:i+len(Pattern)] == Pattern:
            indices.append(i)
    return indices

def CountDict(Text, k):
    Count = {}
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        Count[i] = PatternCount(Pattern, Text)
    return Count

def FrequentWords(Text, k):
    FrequentPatterns = []
    Count = CountDict(Text, k)
    m = max(Count.values())
    for i in Count:
        if Count[i] == m:
            FrequentPatterns.append(Text[i:i+k])
    FrequentPatternsNoDuplicates = RemoveDuplicates(FrequentPatterns)
    return FrequentPatternsNoDuplicates

def RemoveDuplicates(Words):
    clean_list = []
    for word in Words:
        if word not in clean_list:
            clean_list.append(word)
    return clean_list

def ReverseComplement(Word):
    reversed_strand = ReverseStrand(Word)
    complemented_reverse = ""
    dictionary = { 'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G' }
    for letter in reversed_strand:
        complemented_reverse += dictionary[letter]
    return complemented_reverse

def ReverseStrand(Strand):
    mirror_word = ""
    counter = len(Strand) - 1
    while counter >= 0:
        mirror_word += Strand[counter]
        counter -= 1
    return mirror_word

# Week 2 starts here

def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array

def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    array[0] = PatternCount(symbol, Genome[0:n//2])
    for i in range(1, n):
        array[i] = array[i-1]
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i] - 1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i] + 1
    return array

def Skew(Genome):
    skew = {}
    skew[0] = 0

    # First solution:
    #n = len(Genome)
    #for i in range(1, n+1):
    #    if Genome[i-1] == "G":
    #        skew[i] = skew[i-1] + 1
    #    elif Genome[i-1] == "C":
    #        skew[i] = skew[i-1] - 1
    #    else:
    #        skew[i] = skew[i-1]

    # Second solution:
    i = 1 # skew key
    for nucleotide in Genome:
        if nucleotide == "G":
            skew[i] = skew[i-1] + 1
        elif nucleotide == "C":
            skew[i] = skew[i-1] - 1
        else:
            skew[i] = skew[i-1]
        i += 1
    return skew

def MinSkew(Genome):
    skewArray = Skew(Genome)
    minValue = min(skewArray.values())
    minIndices = [] #positions
    for key in skewArray:
        if skewArray[key] == minValue:
            minIndices.append(key)
    return minIndices

def HammingDistance(p, q):
    # p and q are two k-mers
    # the total No. of mismatches is called Hamming distance
    HammingDistance = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            HammingDistance += 1
    return HammingDistance

def ApproximatePatternMatching(Pattern, Text, d):
    positions = []
    for i in range(len(Text)-len(Pattern)+1):
        substrEqualsPattern = Text[i:i+len(Pattern)] == Pattern
        substrAproximateToPattern = HammingDistance(Pattern, Text[i:i+len(Pattern)]) <= d
        if substrEqualsPattern or substrAproximateToPattern:
            positions.append(i)
    return positions

def ApproximatePatternCount(Pattern, Text, d):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        substrEqualsPattern = Text[i:i+len(Pattern)] == Pattern
        substrAproximateToPattern = HammingDistance(Pattern, Text[i:i+len(Pattern)]) <= d
        if substrEqualsPattern or substrAproximateToPattern:
            count += 1
    return count

#print SymbolArray("AAAAGGGG", "A")
#print FasterSymbolArray("AAAAGGGG", "A")
#print Skew("CATGGGCATCGGCCATACGCC")
#print Skew("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT")
#print MinSkew("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT")
#print HammingDistance("GGGCCGTTGGT", "GGACCGTTGAC") # 3
#print ApproximatePatternMatching(
#    "ATTCTGGA",
#    "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT",
#    3
#)

#print ApproximatePatternCount("GAGG", "TTTAGAGCCTTCAGAGG", 2)

#print HammingDistance(
#    "CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA",
#    "CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG"
#)

print MinSkew("GATACACTTCCCGAGTAGGTACTG")
