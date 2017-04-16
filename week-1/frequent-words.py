def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count += 1
    return count

def CountDict(Text, k):
    Count = {}
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        Count[i] = PatternCount(Pattern, Text)
    return Count

#def max(list):
#    m = list[0]
#    for item in list:
#        if item > m:
#            m = item
#    return m

def FrequentWords(Text, k):
    FrequentPatterns = []
    Count = CountDict(Text, k)
    m = max(Count.values())
    for i in Count:
        # If the item on the i-th place is equal to max, append it to FreqPtrns
        if Count[i] == m:
            FrequentPatterns.append(Text[i:i+k])
    return FrequentPatterns

def RemoveDuplicates(Words):
    clean_list = []
    for word in Words:
        if word not in clean_list:
            clean_list.append(word)
    return clean_list

k_mer_4 = FrequentWords("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4)
k_mer_2 = FrequentWords("GATCCAGATCCCCATAC", 2)

print k_mer_4
print k_mer_2

print RemoveDuplicates(k_mer_4)
print RemoveDuplicates(k_mer_2)

#import sys
#lines = sys.stdin.read().splitlines()
#print(CountDict(lines[1],int(lines[0])))
