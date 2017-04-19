def CountNucleotideOccurence(Genome, Nucleotide):
    ExtGenome = ExtendedGenome(Genome)
    NucleotideCount = []
    for i in range(0, len(Genome)):
        counter = 0
        for nukeleotide in ExtGenome[i:i+len(Genome)//2]:
            if nukeleotide == Nucleotide:
                counter += 1
        NucleotideCount.append(counter)
        #print NucleotideCount[i], ExtGenome[i:i+len(Genome)//2]
    return NucleotideCount

def ExtendedGenome(Genome):
    return Genome + Genome[0:len(Genome)//2]

Genome = "AAAAGGGG"
Nucleotide = "A"

print CountNucleotideOccurence(Genome, Nucleotide)

Genome2 = "CTGCTTCGCCCGCCGGACCGGCCTCGTGATGGGGT"
Nucleotide2 = "C"

print CountNucleotideOccurence(Genome2, Nucleotide2)
