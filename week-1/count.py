def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count += 1
    return count

# Count[i] stores PatternCount[Pattern, Text] for Pattern=Text[i:i+k]
def Count(Text, k):
    count_array = []
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        count_array.append(PatternCount(Pattern, Text))
    return count_array

print Count("CGATATATCCATAG", 3)
