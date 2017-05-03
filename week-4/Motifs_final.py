from math import log
from random import randint
from random import uniform

def Count(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

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
    count = Count(Motifs)
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

def ConsensusWithPseudocounts(Motifs):
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

def ScoreWithPseudocounts(Motifs):
    consensus = ConsensusWithPseudocounts(Motifs)
    count = 0
    k = len(Motifs[0])
    t = len(Motifs)
    for j in range(k):
        for i in range(t):
            if Motifs[i][j] != consensus[j]:
                count += 1
    return count

def Entropy(Motif):
    profile = Profile(Motif)
    entropy = 0
    k = len(Motif[0])
    for j in range(k):
        for symbol in "ACGT":
            if profile[symbol][j] != 0:
                entropy += -profile[symbol][j] * log(profile[symbol][j], 2)
    return entropy

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

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs

    return BestMotifs

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
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs

def Normalize(Probabilities):
    running_sum = 0
    for kmer in Probabilities:
        running_sum += Probabilities[kmer]
    for kmer in Probabilities:
        Probabilities[kmer] /= running_sum
    return Probabilities

def WeightedDie(Probabilities):
    toss = uniform(0, 1)
    start = 0
    for key in Probabilities:
        if start < toss < start+Probabilities[key]:
            return key
        start += Probabilities[key]

def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}

    for i in range(n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)

    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

def GibbsSampler(Dna, k, t, N):
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = Motifs
    for j in range (1, N+1):
        i = randint(0, t-1)
        Profile = ProfileWithPseudocounts([motif for idx, motif in enumerate(Motifs) if idx != i])
        Motifs[i] = ProfileGeneratedString(Dna[i], Profile, k)
        if ScoreWithPseudocounts(Motifs) < ScoreWithPseudocounts(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

dna_mat = [
    'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
    'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
    'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
    'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
    'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'
]

print GibbsSampler(dna_mat, 8, 5, 100)

print "------------- QUIZ -----------------"
print " Question 6. "
print
print Normalize({'A' : 0.22, 'B' : 0.54, 'C' : 0.58, 'D' : 0.36, 'E' : 0.3})
print
print " Question 6. EOF"
print "------------- QUIZ -----------------"
