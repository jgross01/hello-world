#File for assembly algorithms
import random

def Composition(text, k):
    kmers = []
    for i in range(len(text)-k+1):
        kmers.append(text[i:i+k])
    return kmers

def SpelledFromPath(path):
    k = len(path[0])
    genome = path[0]
    for i in range(len(path)-1):
        genome = genome + path[i+1][k-1]
    return genome

def Prefix(kmer):
    return kmer[0:len(kmer)-1]

def Suffix(kmer):
    return kmer[1:len(kmer)]

def OverlapGraph(kmers):
    #takes a list of kmers and outputs an adjacency list
    # Requires Prefix and Suffix
    adjacency={}
    for i in range(len(kmers)):
        adjacency[kmers[i]] = []
        for j in range(len(kmers)):
            if Suffix(kmers[i]) == Prefix(kmers[j]):
                adjacency[kmers[i]] += [kmers[j]]
    return adjacency

def PathGraphk(Text,k):
    Path = [[0] * (len(Text) - (k -1)+ 2) for i in range((len(Text) - (k-1) + 2))]
    for i in range(len(Text)-k+2):
        Path[0][i+1]= Text[i:i+k-1]
        Path[i+1][0]= Text[i:i+k-1]
        if i > 0 and i < len(Text)-k+2:
            Path[i][i+1]+=1
    return Path

def DeBruijnk(Text,k):
    path = PathGraphk(Text, k)
    i=0
    # First glue like kmers together
    while i < (len(path)-1):
        j = 0
        while j < (len(path)-(i+1)):
            if i+j+2 == len(path):
                break
            elif path[0][i+1]== path [0][i+1+j+1]:
                for c in range(len(path)-1):
                    path[c+1][i+1] += path[c+1][i+1+j+1]
                for c in path:
                    del c[i+1+j+1]
                for r in range(len(path[0])-1):
                    path[i+1][r + 1] += path[i+1+j+1][r + 1]
                del path[i+1+j+1]
            j+=1
        i+=1
    #Then format output as graph
    for i in range(len(path)-1):
        match = []
        for j in range(len(path[0])-1):
            if path[i+1][j+1] > 0:
                count = path[i+1][j+1]
                while count >0:
                    match.append(path[0][j+1])
                    count-=1
        if match != []:
            print (path[i+1][0] + " -> " + str(','.join([str(i) for i in match])))

def remove_duplicates(list):
    nodupes = []
    for i in list:
       if i not in nodupes:
           nodupes.append(i)
    return nodupes

def DeBruijn(Patterns):
    kmers = []
    for i in Patterns:
        kmers.append(Prefix(i))
        kmers.append(Suffix(i))
    kmers = remove_duplicates(kmers)
    Path = [[0] * (len(kmers)) for i in range((len(kmers)))]
    for i in range(len(Patterns)):
        ith = kmers.index(Prefix(Patterns[i]))
        jth = kmers.index(Suffix(Patterns[i]))
        Path[ith][jth] += 1
    adjlist = []
    # Then format output as graph
    for i in range(len(Path)):
        match = []
        for j in range(len(Path)):
            if Path[i][j] > 0:
                count = Path[i][j]
                while count > 0:
                    match.append(kmers[j])
                    count -= 1

        if match != []:
            adjlist.append(kmers[i] + " -> " + str(','.join([str(i) for i in match])))
    return adjlist

def EulerianFormatt(adjlist):
    # First convert adj list to numerical structure from "["#### -> ###,###",...]\n" format
    indiv = adjlist
    adjdict = {}
    for i in indiv:
        adjdict[i.split(' -> ')[0]] = (i.split(' -> ')[1]).split(',')
        # print (adjdict)
    return adjdict

def EulerianCycle(adjlist, startnode = 'null'):
    #outputs one possible Eulerian cycle from the adjacency list of a directed graph (should be balanced and highly connected)
    # First convert adj list to numerical structure from "["#### -> ###,###",...]\n" format
    if type(adjlist) is not dict:
        adjdict = EulerianFormatt(adjlist)
    else:
        adjdict = adjlist
    #Run through a cycle once
    #Move from first node
    if startnode == 'null':
        startnode = random.choice(list(adjdict.keys()))

    currentnode = startnode
    path = [startnode]
    nextnode = random.choice(adjdict[currentnode])
    adjdict[currentnode].remove(nextnode)
    path += [nextnode]
    currentnode = nextnode
    #iterate until we return to first node
    while currentnode != startnode:
        nextnode = random.choice(adjdict[currentnode])
        adjdict[currentnode].remove(nextnode)
        path += [nextnode]
        currentnode = nextnode
    #print("First path: " + str(path))
    # while some edges are unexplored
    while bool([i for i in adjdict.values() if i != []]):
        #find a new startnode with edges
        newstart = path[0]
        j=0
        while adjdict[newstart]== []:
            j+=1
            newstart = path[j]
        #print ("Newstart: " + str(newstart))
        # Form newpath by moving through a cycle and then more random walking
        currentnode = newstart
        newpath = [] #will inset this path into old path so don't count newstart twice
        nextnode = random.choice(adjdict[currentnode])
        adjdict[currentnode].remove(nextnode)
        newpath += [nextnode]
        currentnode = nextnode
        # iterate until we return to newstart
        while currentnode != newstart:
            nextnode = random.choice(adjdict[currentnode])
            adjdict[currentnode].remove(nextnode)
            newpath += [nextnode]
            currentnode = nextnode
        path[(path.index(newstart)+1):(path.index(newstart)+1)] =  newpath
        #print ("Newpath: " + str(newpath))
        #print ("Path: " + str(path))
    #print (adjdict)
    #print ("Nodes touched: " + str(len(path)))
    return path

def EulerianPath(adjlist):
    #Converts an adjacency list of a directed graph with an Eulerian Path into one with an Eulerian Cycle and finds one such cycle
    adjdict = EulerianFormatt(adjlist)
    indegree = {}
    outdegree = {}
    for i in adjdict.keys():
        outdegree[i] = len(adjdict[i])
        for j in adjdict[i]:
            if j in indegree.keys():
                indegree [j] += 1
            else:
                indegree[j] = 1
    balance = {}
    for i in outdegree.keys():
        if i in indegree.keys():
            balance[i] = outdegree[i] - indegree[i]
            del indegree[i]
        else:
            balance[i] = outdegree[i]
    for i in indegree.keys():
        balance[i] = indegree[i]*-1
    for node, degree in balance.items():
        if degree == 1:
            start = node
        if degree == -1:
            end = node

    if end not in adjdict.keys():
        adjdict[end] = [start]
    else:
        adjdict[end] += [start]
    #adjdict now has a Eulerian cycle, so find one
    path = EulerianCycle(adjdict, start)
    #rotate list & remove the atificial edge
    path = path[0:len(path)-1]
    end_index = 0
    if path[len(path)-1] != end:
        for i in range(len(path)):
            if path[i] == end:
                end_index = i
        path = path[end_index+1:] + path[:end_index+1]
    return path

def StringReconstruction(Patterns, k):
    path = (EulerianPath(DeBruijn(Patterns)))
    genome = path[0]
    for i in range(len(path)-1):
        genome +=(path[i+1][k-2])
    return genome

def StringReconstructionCycle(Patterns, k):
    path = (EulerianCycle(DeBruijn(Patterns)))
    genome = path[0]
    for i in range(len(path)-1):
        genome +=(path[i+1][k-2])
    return genome

def SymbolToNumber(N):
#Convert single base to int A=0, C=1, G=2, T=3.
	if N == "A":
		return 0
	if N == "C":
		return 1
	if N == "G":
		return 2
	if N == "T":
		return 3
	else:
		return "Not a Base"

def NumberToBin(r):
#Convert single int to base A=0, C=1, G=2, T=3.
	if r == 0:
		return '0'
	if r == 1:
		return '1'

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

def NumberToBinar(index,k):
#Convert int representing a k-mer to the DNA code for that k-mer.  Depends on NumberToSymbol
	if k==1:
		return NumberToBin(index)
	prefixIndex = index // 2
	r = index % 2
	symbol = NumberToBin(r)
	prefixPattern = NumberToBinar(prefixIndex, (k-1))
	return prefixPattern + symbol
#with open('input.txt', 'r') as myfile: patterns = myfile.read().split('\n')

#adjlist =["0 -> 2","1 -> 3","2 -> 1","3 -> 0,4","6 -> 3,7","7 -> 8","8 -> 9","9 -> 6"]
#patterns = ["CTTA", "ACCA", "TACC", "GGCT", "GCTT", "TTAC"]
k= 9
patterns = []
for i in range(2**k):
    patterns.append(NumberToBinar(i, k))

output = StringReconstructionCycle(patterns, k)
print (output[:len(output)-(k-1)])



#print (str(output[0]), end = '')
#for i in range(len(output)-1):
   #print  (("->" + str(output[i+1])), end='')



    #with open('input.txt', 'r') as myfile: kmers = myfile.read().split('\n')

