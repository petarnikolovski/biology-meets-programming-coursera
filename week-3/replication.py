# Week 1 starts here

def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count += 1
    return count

def PatternMatching(Pattern, Genome):
    indices = []
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

    i = 1
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
    minIndices = []
    for key in skewArray:
        if skewArray[key] == minValue:
            minIndices.append(key)
    return minIndices

def HammingDistance(p, q):
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

# Week 3 starts here
Motif1 = "TCGGGGGTTTTT"
Motif2 = "CCGGTGACTTAC"
Motif3 = "ACGGGGATTTTC"
Motif4 = "TTGGGGACTTTT"
Motif5 = "AAGGGGACTTCC"
Motif6 = "TTGGGGACTTCC"
Motif7 = "TCGGGGATTCAT"
Motif8 = "TCGGGGATTCCT"
Motif9 = "TAGGGGAACTAC"
Motif10 = "TCGGGTATAACC"

Motifs = [Motif1, Motif2, Motif3, Motif4, Motif5, Motif6, Motif7, Motif8, Motif9, Motif10]

def Score(Motifs):
    score = 0
    count = Count(Motifs)
    for j in range(len(count['A'])):
        temp = [count['A'][j], count['C'][j], count['G'][j], count['T'][j]]
        # Find me a maximum in a columns, and return me it's index
        # then delete that element in the list
        # this doesn't account if there are two same max values
        del temp[temp.index(max(temp))]
        score += sum(temp)
    return score

def Profile(Motifs):
    t = len(Motifs)
    count = Count(Motifs)
    profile = {}
    for nucleotide in count:
        profile[nucleotide] = [float(num) / t for num in count[nucleotide]]
    return profile

# Count matrix has four rows: A C G T, and k columns
def Count(Motifs):
    count = {
        'A' : [],
        'C' : [],
        'G' : [],
        'T' : []
    }
    for j in range(len(Motifs[0])):
        A = 0
        C = 0
        G = 0
        T = 0
        for i in range(len(Motifs)):
            if Motifs[i][j] == 'A':
                A += 1
            elif Motifs[i][j] == 'C':
                C += 1
            elif Motifs[i][j] == 'G':
                G += 1
            elif Motifs[i][j] == 'T':
                T += 1
            else:
                print "Not a valid nucleotide."
                return {}
        count['A'].append(A)
        count['C'].append(C)
        count['G'].append(G)
        count['T'].append(T)
    return count

def Consensus(Motifs):
    consensus = []
    count = Count(Motifs)
    for j in range(len(count['A'])):
        temp = [count['A'][j], count['C'][j], count['G'][j], count['T'][j]]
        idx = temp.index(max(temp))
        if idx == 0:
            consensus.append('A')
        elif idx == 1:
            consensus.append('C')
        elif idx == 2:
            consensus.append('G')
        elif idx == 3:
            consensus.append('T')
        else:
            print "Wrong index value. Index out of bounds."
            return []
    return consensus

print Profile(Motifs)
