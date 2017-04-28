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

def Count(Motifs):
    count = {}
    # Create a list of zeroes
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    # Populate cout "matrix" with values
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

from math import log

def Entropy(Motif):
    profile = Profile(Motif)
    entropy = 0
    k = len(Motif[0])
    for j in range(k):
        for symbol in "ACGT":
            if profile[symbol][j] != 0:
                entropy += -profile[symbol][j] * log(profile[symbol][j], 2)
    return entropy

# Find a probability of a certain Motif [Text] being the consensus string in
# a Profile matrix of a certain Motif matrix (i.e. collection of motifs).
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

#def ProfileMostProbablePattern(Text, k, Profile):
#    most_prob = Text[0:k]
#    p_max = Pr(Text[0:k], Profile)
#    for i in range(1, len(Text) - k + 1):
#         if Pr(Text[i:i+k], Profile) > p_max:
#                p_max = Pr(Text[i:i+k], Profile)
#                most_prob = Text[i:i+k]
#    return most_prob

# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    # uzmi mi od svih DNA kojih je ukupno t, svaki prvi k-mer [0:k]
    # proglasi mi ove k-mere da su najbolji, dok ne nadjem bolje
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])
    # predji mi preko celog isecka DNA
    for i in range(n-k+1):
        Motifs = []
        # uzmi mi iz DNA i do i+k-ti isecak i stavi mi u motive
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            # sracunaj mi profil prvi_put/iznova
            P = Profile(Motifs[0:j])
            # sad mi uzmi sledeci Dna i za unteo k, na osnovu
            # dosadasnjeg profila, nadji mi najverovatniji motiv
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs

    return BestMotifs


Text = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
Prrrofile = {
    'A' : [0.2, 0.2, 0.3, 0.2, 0.3],
    'C' : [0.4, 0.3, 0.1, 0.5, 0.1],
    'G' : [0.3, 0.3, 0.5, 0.2, 0.4],
    'T' : [0.1, 0.2, 0.1, 0.1, 0.2]
}

#print ProfileMostProbablePattern(Text, 5, Prrrofile)

Dna1 = "GGCGTTCAGGCA"
Dna2 = "AAGAATCAGTCA"
Dna3 = "CAAGGAGTTCGC"
Dna4 = "CACGTCAATCAC"
Dna5 = "CAATAATATTCG"

Dna = [Dna1, Dna2, Dna3, Dna4, Dna5]

# Protein translator
def ProteinTranslator(RNA):
    GeneticCode = {
        'H' : ['CAU', 'CAC'],
        'Q' : ['CAA', 'CAG'],
        'P' : ['CCU', 'CCC', 'CCA', 'CCG'],
        'R' : ['CGU', 'CGC', 'CGA', 'CGG'],
        'L' : ['CUU', 'CUC', 'CUA', 'CUG'],
        'D' : [],
        'E' : [],
        'A' : [],
        'G' : [],
        'V' : [],

    }



#print GreedyMotifSearch(Dna, 3, 5)
"""
#  dormancy survival regulator
DosR = open("DosR.txt", "r")
DosR_DNAs = DosR.read().splitlines()
DosR.close()

best_motifs = GreedyMotifSearch(DosR_DNAs, 15, len(DosR_DNAs))
score_motif = Score(best_motifs)

print best_motifs
print
print score_motif
"""
print "--- QUIZ WEEK 3 ---"
print "Question 5"
prof = {
    'A': [0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
    'C': [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
    'G': [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
    'T': [0.3, 0.1, 0.0, 0.4, 0.5, 0.0],
}

print Pr("AAGTTC", prof)
print

print Pr("GAGCTA", prof)
print

print "--- QUIZ WEEK 3 ---"
