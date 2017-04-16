def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count += 1
    return count

text = "CGATATATCCATAG"
pattern = "ATA"
count = PatternCount(pattern, text)
print "The %s appears in the given pattern %d times" % (pattern, count)
