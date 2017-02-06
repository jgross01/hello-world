def remove_duplicates(list):
    nodupes = []
    for i in list:
       if i not in nodupes:
           nodupes.append(i)
    return nodupes

def CountDict(Text, k):
    Count = {}
    # Output: CountDict(Text, k) Dict with counts for each k-mer of specific size
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        Count[i]= PatternCount(Text,Pattern)
    return Count

def PatternCount(Text, Pattern):
    #Outputs the number of times input pattern occers in text
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count

def FrequentWords(Text, k):
    #Outputs maximally occing k-mers in input text, without duplicates
    FrequentPatterns = [] # output variable
    Count=CountDict(Text, k)
    m = max(Count.values())
    for i in Count:
        if Count[i]==m:
            FrequentPatterns.append(Text[i:i+k])
    FrequentPatternsNoDupes = remove_duplicates(FrequentPatterns)
    return FrequentPatternsNoDupes








def SymbolToNumber(ACGT):
#Convert single base to int A=0, C=1, G=2, T=3.
    if ACGT == "A":
        return 0
    if ACGT == "C":
        return 1
    if ACGT == "G":
        return 2
    if ACGT == "T":
        return 3

def NumberToSymbol(r):
#Convert single int to base A=0, C=1, G=2, T=3.
    if r == 0:
        return 'A'
    if r == 1:
        return 'C'
    if r == 2:
        return 'G'
    if r == 3:
        return 'T'

def PatternToNumber(Pattern):
#Convert k-mer to lexicographical int, depends on SymbolToNumber
    if len(Pattern)==0:
        return 0
    symbol = Pattern[len(Pattern)-1]
    prefix = Pattern[:len(Pattern)-1]
    return 4*PatternToNumber(prefix)+SymbolToNumber(symbol)

def NumberToPattern(index,k):
#Convert int representing a k-mer to the DNA code for that k-mer.  Depends on NumberToSymbol
    if k==1:
        return NumberToSymbol(index)
    prefixIndex = index // 4
    r = index % 4
    symbol = NumberToSymbol(r)
    prefixPattern = NumberToPattern(prefixIndex, (k-1))
    return prefixPattern + symbol

def ComputeFrequencyArray(Text, k):
# Compute frequency array of k-mers for string of text.  Require PatternToNumber ==> SymbolToNumber.
    FrequencyArray = [0]*(4**k)
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        j = PatternToNumber(Pattern)
        FrequencyArray[j]= FrequencyArray[j] +1
    return FrequencyArray


def ReverseCompliment(Pattern):
    # Converts A<-->T and C<-->G and reverses order
    PatternRC = Pattern
    l = len(Pattern)
    i=0
    while i<l:
        if Pattern[i] == 'A':
            PatternRC = PatternRC[:l-i-1] + 'T' + PatternRC[l-i:]
        if Pattern[i] == 'C':
            PatternRC = PatternRC[:l-i-1] + 'G' + PatternRC[l-i:]
        if Pattern[i] == 'G':
            PatternRC = PatternRC[:l-i-1] + 'C' + PatternRC[l-i:]
        if Pattern[i] == 'T':
            PatternRC = PatternRC[:l-i-1] + 'A' + PatternRC[l-i:]
        i = i + 1
    return PatternRC

def PatternFinder(Text, Pattern):
    #Outputs the index positions of pattern in text
    found = []
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            found.append(i)
    if found == []:
        found = 'k-mer not found'
    return found

def ClumpFinder(Genome, k, t, L):
# Finds k-mers clumped together in a genome. Returns k-mers if full string of lenth k shows up at least t time in a sliding window L.
# Requires ComputeFrequencyArray ===> PatternToNumber ===> SymbolToNumber
    FrequentPatterns=[]
    Clump = [0]*(4**k)
    Window = Genome[:L]
    FrequencyArray = ComputeFrequencyArray(Window, k)
# Now that frequencies in intial window are computed, look for any initial clumps
    for i in range(4**k):
        if FrequencyArray[i] >= t:
            Clump[i] = 1
# Then remove first k-mer and add next k-mer and cheeck to see if addition creates a new clump
    for i in range(1, len(Genome)-L+1):
        FirstPattern = Genome[i-1:i+k-1]
        index = PatternToNumber(FirstPattern)
        FrequencyArray[index]= FrequencyArray[index] - 1
        LastPattern = Genome[i+L-k: i+L]
        index = PatternToNumber(LastPattern)
        FrequencyArray[index]= FrequencyArray[index] + 1
        if FrequencyArray[index] >= t:
            Clump[index] = 1
# Now convert the computed clumped k-mers to strings from indices and return them in a list
    for i in range(4**k):
        if Clump[i] == 1:
            FrequentPatterns.append(NumberToPattern(i, k))
    return FrequentPatterns




def MinSkew(Genome):
# Calculates a skew array for Genome in put of Total G - Total C at every position
# Then finds all potisions that minimize G-C skew in order to predict location of Ori
    L = len(Genome)
    SkewArray = [0]* (L+1)
    for i in range(L):
        if Genome[i] == "C":
            SkewArray[i+1]=SkewArray[i] - 1
        elif Genome[i] == "G":
            SkewArray[i+1]=SkewArray[i] + 1
        else:
            SkewArray[i+1] = SkewArray[i]
    MyMin = min(SkewArray)
    return  [ i   for i in range(L+1) if SkewArray[i] == MyMin]


def FormatList(List):
# Formating list of numerbers for pint
#use with 'print(" ".join(FormatList(XXXXXXX)))'
    string = [""] * len(List)
    for i in range(len(List)):
        string[i]= str(List[i])

    return string


def HammingDist(k1, k2):
#Computes Hamming Distance or # of mismatches between to k-mers of same size
    if len(k1) != len(k2):
        return "k-mers must be same size"
    else:
        ham = 0
        for i in range(len(k1)):
            if k1[i] != k2[i]:
                ham = ham +1
        return ham




def ApproxMatches(Pattern, Genome, d):
    #Output the positions of k-mers with at most d mismatches from input k-mer in Genome
    MatchPositions = []
    for i in range(len(Genome)-len(Pattern)+1):
        if HammingDist(Genome[i:i+len(Pattern)], Pattern) <= d:
           MatchPositions.append(i) #Fill in positions or fill list with positions
    return MatchPositions

def ApproxPatternCount(Text, Pattern, d):
    #Outputs the number of times input pattern occers in text
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDist(Text[i:i+len(Pattern)],Pattern) <= d:
            count = count+1
    return count

def ImmediateNeighbours(Pattern):
	#Outputs all strings with a Hamming distance of 1 from input string
	#Requires HammingDistance and NumberToSymbol
	neighbours = [Pattern]
	for i in range(len(Pattern)):
		for j in range(4):
			if NumberToSymbol(j) != Pattern[i]:
				neighbours.append(Pattern[:i]+NumberToSymbol(j)+Pattern[i+1:])
	return neighbours

def NeighboursIterative(Pattern, d):
	#Outputs the list of k-mers with up to d mismatches with input text pattern iteratively
	# Requires ImmediateNeighbours <===> HammingDistance <===> NumberToSymbol, remove_duplicates
	neighbours = [Pattern]
	for i in range(d):
		for j in range(len(neighbours)):
			nearneighbours= ImmediateNeighbours(neighbours[j])
			for x in nearneighbours:
				neighbours.append(x)
	neighbours = remove_duplicates(neighbours)
	return neighbours


def NeighboursRecursive(Pattern, d):
    #Outputs the list of k-mers with up to d mismatches with input text pattern recursively
	#Requires HammingDistance and NumberToSymbol
    if d == 0:
        return Pattern
    if len(Pattern) == 1:
        return ["A","C","G","T"]
    neighbours = []
    SuffixNeighbours = Neighbours(Pattern[1:],d)
    for i in SuffixNeighbours:
        if HammingDist(i, Pattern[1:]) < d:
            for j in range(4):
                neighbours.append(NumberToSymbol(j)+i)
        else:
            neighbours.append(Pattern[0] + i)
    return neighbours

def SymbolArray(Genome, symbol):
	#Outputs a dictionary "Symbol Array" with contains keys corresponding to position in extended genome and values that count the occurences of the base at that position in a sliding window
	# Requires PatternCount
	array = {}
	n = len(Genome)
	ExtendedGenome = Genome + Genome[0:n//2]
	for i in range(n):
		array[i]= PatternCount(ExtendedGenome[i:i+(n//2)],symbol)
	return array
