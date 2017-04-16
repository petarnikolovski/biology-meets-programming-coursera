def PatternCount(Pattern, Text):
    counter = 0

    end = len(Pattern)
    start_index = 0
    end_index = len(Text)
    while end_index <= end:
        if Pattern[start_index:end_index] == Text:
            counter += 1

        print Pattern[start_index:end_index], start_index, end_index
        start_index +=1
        end_index += 1

    return counter

pattern = "CGATATATCCATAG"
text = "ATA"
count = PatternCount(pattern, text)
print "The %s appears in the given pattern %d times" % (text, count)
