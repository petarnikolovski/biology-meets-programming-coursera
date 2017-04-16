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

def complement(Nucleotide):
    comp = '' # output variable
    return comp

oriC_vibrio_cholerae = "ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC"
k = 10
#k_mer_4 = FrequentWords("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4)
#k_mer_2 = FrequentWords("GATCCAGATCCCCATAC", 2)

#print k_mer_4
#print k_mer_2

print FrequentWords(oriC_vibrio_cholerae, k)
print ReverseComplement("AAAACCCGGT")
print PatternMatching("ATAT", "GATATATGCATATACTT")

vibrio_cholerae_genome = open("vibrio_cholerae_genome.txt", "r")
print PatternMatching("CTTGATCAT", vibrio_cholerae_genome.readline())

vibrio_cholerae_genome.close()

thermotoga_petrophila = "AACTCTATACCTCCTTTTTGTCGAATTTGTGTGATTTATAGAGAAAATCTTATTAACTGAAACTAAAATGGTAGGTTTGGTGGTAGGTTTTGTGTACATTTTGTAGTATCTGATTTTTAATTACATACCGTATATTGTATTAAATTGACGAACAATTGCATGGAATTGAATATATGCAAAACAAACCTACCACCAAACTCTGTATTGACCATTTTAGGACAACTTCAGGGTGGTAGGTTTCTGAAGCTCTCATCAATAGACTATTTTAGTCTTTACAAACAATATTACCGTTCAGATTCAAGATTCTACAACGCTGTTTTAATGGGCGTTGCAGAAAACTTACCACCTAAAATCCAGTATCCAAGCCGATTTCAGAGAAACCTACCACTTACCTACCACTTACCTACCACCCGGGTGGTAAGTTGCAGACATTATTAAAAACCTCATCAGAAGCTTGTTCAAAAATTTCAATACTCGAAACCTACCACCTGCGTCCCCTATTATTTACTACTACTAATAATAGCAGTATAATTGATCTGA"
print PatternMatching("CTTGATCAT", thermotoga_petrophila)
print PatternMatching("ATGATCAAG", thermotoga_petrophila)

print PatternCount("CTTGATCAT", thermotoga_petrophila)
print PatternCount("ATGATCAAG", thermotoga_petrophila)

print PatternCount("AAA", "GACCATCAAAACTGATAAACTACTTAAAAATCAGT")

print ReverseComplement("TTGTGTC")

print "---------------------------------"
print
print FrequentWords("CGCCTAAATAGCCTCGCGGAGCCTTATGTCATACTCGTCCT", 3)
print
print "---------------------------------"
