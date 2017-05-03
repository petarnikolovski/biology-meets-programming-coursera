from random import randint
from random import uniform

def Profile(Motifs):
    t = len(Motifs)
    count = Count(Motifs)
    profile = {}
    for nucleotide in count:
        profile[nucleotide] = [float(num) / t for num in count[nucleotide]]
    return profile

def Consensus(Motifs):
    consensus = ""
    k = len(Motifs[0])
    count = CountWithPseudocounts(Motifs)
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

def Score(Motifs):
    consensus = Consensus(Motifs)
    count = 0
    k = len(Motifs[0])
    t = len(Motifs)
    for j in range(k):
        for i in range(t):
            if Motifs[i][j] != consensus[j]:
                count += 1
    return count

def Pr(Text, Profile):
    running_product = 1
    for index, nucleotide in enumerate(Text):
        running_product *= Profile[nucleotide][index]
    return running_product

def ProfileMostProbablePattern(Text, k, Profile):
    ProbabilityList = []
    for i in range(len(Text)-k+1):
        ProbabilityList.append(Pr(Text[i:i+k], Profile))
    MostProbableAt = ProbabilityList.index(max(ProbabilityList))
    return Text[MostProbableAt:MostProbableAt+k]

def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs

    return BestMotifs

def CountWithPseudocounts(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

def ProfileWithPseudocounts(Motifs):
    t = len(Motifs) + 4
    count = CountWithPseudocounts(Motifs)
    profile = {}
    for nucleotide in count:
        profile[nucleotide] = [float(num) / t for num in count[nucleotide]]
    return profile

def Motifs(Profile, Dna):
    k = len(Profile['A'])
    t = len(Dna)
    motifs = []
    for i in range(t):
        motifs.append(ProfileMostProbablePattern((Dna[i]), k, Profile))
    return motifs

def RandomMotifs(Dna, k, t):
    motifs = []
    for i in range(t):
        idx = randint(0,len(Dna[i])-k)
        motifs.append(Dna[i][idx:idx+k])
    return motifs

def RandomizedMotifSearch(Dna, k, t):
    M = [
        'GTC',
        'CCC',
        'ATA',
        'GCT'
    ]
    BestMotifs = M

    Profile = ProfileWithPseudocounts(M)
    M = Motifs(Profile, Dna)
    BestMotifs = M

    return BestMotifs

dna_mat = [
    'ATGAGGTC',
    'GCCCTAGA',
    'AAATAGAT',
    'TTGTGCTA'
]

print RandomizedMotifSearch(dna_mat, 3, 4)
