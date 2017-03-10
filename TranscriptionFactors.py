#Print items in a list with spaces inbetween
#print(' '.join([str(i) for i in function(x))]))
import math
import random

def remove_duplicates(list):
    nodupes = []
    for i in range(len(list)):
       if list[i] not in nodupes:
           nodupes = nodupes + [list[i]]
    return nodupes

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

def ImmediateNeighbours(Pattern):
    #Outputs all strings with a Hamming distance of 1 from input string
	#Requires HammingDistance and NumberToSymbol
    neighbours = [Pattern]
    for i in range(len(Pattern)):
        for j in range(4):
            if NumberToSymbol(j) != Pattern[i]:
                neighbours = neighbours+ [(Pattern[:i]+NumberToSymbol(j)+Pattern[i+1:])]
    return neighbours

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

def NeighboursRecursive(Pattern, d):
    #Outputs the list of k-mers with up to d mismatches with input text pattern recursively
	#Requires HammingDistance and NumberToSymbol
    if d == 0:
        return Pattern
    if len(Pattern) == 1:
        return ["A","C","G","T"]
    neighbours = []
    SuffixNeighbours = NeighboursRecursive(Pattern[1:],d)
    for i in range(len(SuffixNeighbours)):
        if HammingDist(SuffixNeighbours[i], Pattern[1:]) < d:
            for j in range(4):
                neighbours = neighbours +[(NumberToSymbol(j)+SuffixNeighbours[i])]
        else:
            neighbours = neighbours + [Pattern[0] + SuffixNeighbours[i]]
    return neighbours

def MotifEnumerationBruteForce(DNA, k, d):
# Counts k-mers with d mismatches that occur in all entries of a list of strings, DNA
# Requires NeighboursRecursive ==> HammingDistance ====> NumberToSymbol
    Patterns = []
    for i in range(len(DNA[0])-k+1):
        neighbourhood = NeighboursRecursive(DNA[0][i:i+k], d)
        for j in range(len(neighbourhood)):
            JinS = False
            for s in range(1,len(DNA)):
                for p in range(len(DNA[s])-k+1):
                    if HammingDist(DNA[s][p:p+k], neighbourhood[j]) <= d:
                        JinS = True
                        break
                if JinS == True:
                    JinS = False
                elif JinS == False:
                    break
                if s == len(DNA) - 1:
                    Patterns = Patterns + [neighbourhood[j]]
    Patterns = remove_duplicates(Patterns)
    return Patterns

def Count(Motifs):
    # Takes a list of strings of same size and returns a count matrix as a dictionary of lists.
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol]= []
        for j in range(k):
            count[symbol] += [0]
    #Count matrix is initialized
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            count[Motifs[i][j]][j]+=1
    return count

def CountLaPlace(Motifs):
    # Takes a list of strings of same size and returns a count matrix as a dictionary of lists.
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol]= []
        for j in range(k):
            count[symbol] += [1]
    #Count matrix is initialized
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            count[Motifs[i][j]][j]+=1
    return count

def Profile(Motifs):
    # Output profile matrix for input list of strings (frequencies of each nucleotide at each position)
    # Requires Count
    count = Count(Motifs)
    for i in count:
        for j in range(len(count[i])):
            count[i][j] /= float(len(Motifs))
    return count

def ProfileLaPlace(Motifs):
    # Output profile matrix for input list of strings (frequencies of each nucleotide at each position)
    # Requires Count
    count = CountLaPlace(Motifs)
    for i in count:
        for j in range(len(count[i])):
            count[i][j] /= float(len(Motifs)+4.0)
    return count

def Consensus(Motifs):
    # Returns consensus sequence for input list of strings of uniform size
    #Requires Count
    count = Count(Motifs)
    k = len(Motifs[0])
    consensus =""
    for j in range(k):
        m=0
        frequentsymbol = ""
        for i in "ACGT":
            if count[i][j] > m:
                m = count[i][j]
                frequentsymbol = i
        consensus += frequentsymbol
    return consensus

def Score(Motifs):
#Compute score of consensus sequence for input list of strings of uniform lenth
#Requires Consensus ===> Count
    consensus = Consensus(Motifs)
    score = 0
    for j in range(len(Motifs[0])):
        for i in "ACGT":
            if i != consensus[j]:
                score += Count(Motifs)[i][j]
    return score

def Pr(Text, Profile):
# Takes a string and a profile matrix and outputs the probability that the string was generated from that profile
    probability = 1
    for i in range(len(Text)):
        probability *= Profile[Text[i]][i]
    return probability

def ProfileMostProbablePattern(Text, k , Profile):
# Takes a string of text, an int k and a profile matrix 4 x k
# Outputs the most probable k-mer in text to have been generated by profile matrix
# Requires Pr
    maxprob = 0
    probpattern = ""
    for i in range(len(Text)-k+1):
        prob = Pr(Text[i:i+k],Profile)
        if prob > maxprob:
            maxprob = prob
            probpattern = Text[i:i+k]
    if maxprob == 0:
        return Text[0:k]
    else:
        return probpattern

def GreedyMotifSearch(Dna,k,t):
# Finds a set of low scoring k-mer Motifs in a list of t strings of uniform lenth
#Requires Score ===> Consensus ===> Count, ProfileMostProbablePattern ===> Pr
    BestMotifs=[]
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1,t):
            profile = Profile(Motifs)
            Motifs.append(ProfileMostProbablePattern(Dna[j],k,profile))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

def d(Pattern, Text):
    k =len(Pattern)
    min = HammingDist(Pattern, Text[0:k])
    for i in range(1, len(Text)- k+1):
        dist = HammingDist(Pattern, Text[i:i+k])
        if dist < min:
            min = dist
    return min

def dList(Pattern, Dna):
    dist = 0
    for i in range(len(Dna)):
        dist += d(Pattern, Dna[i])
    return dist

def MedianString(Dna, k):
    distance = math.inf
    dict ={}
    Median = ""
    for i in range(4**k):
        bestd = dList(NumberToPattern(i,k),Dna)
        if distance > bestd:
            distance = bestd
            Median = NumberToPattern(i,k)
            dict[NumberToPattern(i,k)]=bestd
    return Median

def GreedyMotifSearchLaPlace(Dna,k,t):
# Finds a set of low scoring k-mer Motifs in a list of t strings of uniform lenth
#Requires Score ===> Consensus ===> Count, ProfileMostProbablePattern ===> Pr
    BestMotifs=[]
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1,t):
            profile = ProfileLaPlace(Motifs)
            Motifs.append(ProfileMostProbablePattern(Dna[j],k,profile))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

def Motifs(profile, Dna):
    #for a given profile matrix (4 x k)and list of strings, returns a list of k-mers, one from each string, which are Profile-most probable matches
    k = len(profile['A'])
    motifs = []
    for i in range(len(Dna)):
        motifs.append(ProfileMostProbablePattern(Dna[i],k,profile))
    return motifs

def RandomMotifs(Dna, k, t):
    motifs = []
    for i in range(t):
        randposition = random.randint(0,len(Dna[i])-k)
        motifs.append(Dna[i][randposition:(randposition+k)])
    return motifs

def RandomizedMotifSearch(Dna, k, t, N):
    # Randomly generates a set of motifs and then uses their profile to iteratively generate new better scoring motifs until a worse motif stops iteration
    motifs = RandomMotifs(Dna, k, t)
    Nbestmotifs = motifs
    for n in range(N):
        motifs = RandomMotifs(Dna, k, t)
        bestmotifs = motifs
        while True:
            profile = ProfileLaPlace(motifs)
            motifs = Motifs(profile, Dna)
            if Score(motifs) < Score(bestmotifs):
                bestmotifs = motifs
            else:
                break
        if Score(bestmotifs) < Score(Nbestmotifs):
            Nbestmotifs = bestmotifs
    return Nbestmotifs

def Normalize(Probabilities):
    #Takes a dict whos keys are kmers and whos values are probabilies generated with Pr
    #Normalizes probabilites so that the add to 1
    total = 0
    normalized = {}
    for keys in Probabilities:
        total += Probabilities[keys]
    for keys in Probabilities:
        normalized[keys] = Probabilities[keys]/total
    return normalized

def WeightedDie(Probabilities):
    tally = 0
    roll = random.uniform(0,1)
    for keys in Probabilities:
        tally += Probabilities[keys]
        if roll < tally:
            return keys

def ProfileGeneratedString(Text, profile, k):
    probabilities = {}
    for i in range(len(Text)-k+1):
        probabilities[Text[i:i+k]]= Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

def GibbsSampler(Dna, k, t, N):
    motifs = RandomMotifs(Dna, k ,t)
    bestmotifs = motifs
    for n in range(N):
        Roll1 = random.randint(0, t-1)
        profile = ProfileLaPlace(motifs[0:Roll1]+motifs[Roll1+1:])
        motifs[Roll1]= ProfileGeneratedString(Dna[Roll1], profile, k)
        if Score(motifs) < Score(bestmotifs):
            bestmotifs = motifs
    return bestmotifs



#DNA string list input
'''
DNA = [
    "ACGTAAAGAGCCATCTGAGCTATCTCAAGCGGTCTAGGACAA",
    "GAATACAATGGGTTCAGTGGCTACGCCATCTCATTACGTGAC",
    "CGGGCCTTACGCGCTATCTGGGAGCCGGTTGTGCTCCCCAGG",
    "GCGATCCGAAGGCCTGGGTGTTTGGCATGCACTTTCGCCTGT",
    "CGAGCAGGATCAGTCCAAGCGATCCCGCGTTCACACATCTGG",
    "GCAATCATTATGCAACCGAATACAGCTCATGGAAAAAAGGTT",
    "GTGCAAGCTATCGATTGGCGTGTAACTCTTGTGGATAAGTCA",
    "TCATTCTTTGCAGAAGAAGCCATCATCTTTTCATATGCTTGG",
    "AAGTTGTCCTAGTGACGTCGGTTGCGGAGTGCCATCATTTAT",
    "GAAGCAGCTATTCCTTTATATTGTCTGAACGCAATCCATTAG"
]
'''

#Profile matrix input
"""
Profile = {
    'A': [0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
    'C': [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
    'G': [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
    'T': [0.3, 0.1, 0.0, 0.4, 0.5, 0.0]
}
"""

#Benchmark "Subtle Implated Motif" input
"""
Benchmark tests time for a challenging "subtle motifs" solution (find patterns with mismatches of lenth k that appear in all strings in the list
Benchmark = ["TCCGAACAACGAGTAGGCGTACTCACCGGCATGGCCGGATACACCGACCATCGCGGACGAGAAAGGGAGGGCTGAAATACAGACAGCGTACTGTATTAAGCAGAAACGAGAGGAGACAGATCTCATCCCTGGTGTGGTGGAACTGGAGGACTCGCCTCGGTGTGAGTCGTAAGGTGACCGACGATGAAATGCAAGTTCCAACGGCCAACAGCGCGTCAACAACAATGCGCACGTGTCGTAAACTGACGTGAGGTCCCCTTATAGCCCATGAAGAACTTTGACTCGCCTCCGGTAGCCGCTAGTTTTATGCGTGAATGTGCGTTATGCCAACTCAAATGTCTCGCAAGTCAATGAATCAACCGCGATTCTTTATTAACTTCATATCAGGCTAACAAGGACAGACAGCAACAAAGTTCTGCAAAACTGTTCCGGTCTCATCCCTAACTCTCTAACTGATAACAGTCTAACTTGCACCAAGAGTCGCTCGATCGACCAAAGAAATTACCGCCGCCTTGCAGTTCCGATGCCTGGAGTCCCCCTCCGTGTGAAGGTGATAAACCATTTGTCCAACAATGTTAGACAATAAACCACGTAAAAGGC",
             "GGAACTAGCTTAAAAAATAGCAGGGTGTGCCTGATCCTTCCGGTGTTTAAGTAGAAGGCAGGACGGACAGAGTTCCATCCACAGAAGCATAGTTTGATCGTATTGGCGACAGGCTGATGCGAAGCTCCGCTCCAAACGAGAGAGATAAATGCATGCGGTTTGGCCTAAGGCGGGGGGGCAACCCGGCTTATCAATTAGCTAGCCTTGCTTTGGAACAAGGGCCAAGCGGGAGGTAAACTCTTCAGCCCGGGTGTCCCAGTAGCGCGATTTGGTGCTAGCCAGGTTTCGATCAAAAGGGGCTCTTGCAACGCTCTCTTCTAAAAATAAAATGCAATTAGTTGGCAGGGTTGATGAGTGTCGAATCCTTGCAAGCGAGATTCTTCCATGCAGTTGCAGCGGGGCGAGGCCAAAGAGCTCAGCTAGCTTGGGGACTCGCGCCCTGCTTATTCACCCTCGGTGCAACTAATCCTTACCGTGAATTTGGTAGATGTCCAAGCATTGTTCTTTATTAATGCACGTGTTAATAGGGGTAGACTATTCCCGTCCCGGCCTACGGTGTCAAAATCAAGTAGGCGCCGAGATTATTCTTGCATGCCTGTA",
             "TGCGAACTAGTTTTCGCAACTTAACGTAGCGCGTGGGGCGTCCCTAGTGGCTCTGTCAAAGCAATTTGGTTCGTTTAGCTGTTATAGTTTTGGATCACAGCGAATAGAGTTAGTCTTGCTAGTCCGTAAATCAACGGACCGCGTCCCATTAAACACAGCTTCGTCGAGTCTATGACGCTCATACTCTACCATGACCGCGCCGGGACGACCGCCAACTCATAAATGAACGCCTAATAGAACCGAAAAGGGTCGGCGGCACAAAACTCCGGAACGTGGTCTGGGTTAACAAAGGCGCGATGATATTGTTCGTAGATCCCTGTTGGACTCTCCAACAAGTTTCCCGGAGGACTCGAGGTTCCAGGCCGAGTAAATAAAAGTTTTCTCGGGGTGGTGCCGGAAGGCGGGAAGTGGTGGTTAGGACAGATAATGACGAAAACAATGGATCGTGGAAGAGATCGCCCAGAGGTTCGATAGGATGTTACGCTACTTGTGTTCGAGGGGGAGACGGTTTCTACCTAGGCGGGTACCACAAAGCTGTTCTCTATTCTGGAAATTATGTACTCTGTTACTTGAATAAAATAAACAGCGGGGTACGCGGAT",
             "ATCCTGACTACGGCGGTTTTCGTCTTGGGTAGGCACGGAGCTAGAGTATACACGGCAGCTCGTAGGGGGTCGATGCGTCTCGATTAGCTCGTTCCTATAGCTCAGCGATATCCCCGGGTTAAGAAGATTGCTCTCGTTACGCACTAGCCTCCGACTCGCGGGGCGTAACCAGTACAGTAAAAGACGCTAGAATCGACGCTTTCGCATAGTAGTCATTTAGAACCCGGGCTTAAACGATCGTACTTGATACACCCCGGGAGATGTGGATACCATTAAGTTAACCAGATCTATATGCGACCAGTCCTGCAGTAAAGATTGGCTGTCTTGGACTTGTATGCAAGCATAATCAGGGCAGAGGCAGTGGTCCGTTGCCTGAGGACGTCAAGAGTTCTCAGTCTAAAGTATTCCGGGGAAATAGTTAGTTGGCATAAGTCCGCCAAAGATCGCAGATGGTTAGTAGGTAACACTGGGGCCCTCCAGCTTAAGCCAAGCTAACTACGCTCAAGCAGGCTTTTTTTTTATGTTGAACAGAAAAAAAGGGGTTTTCACGCACACTTAGCCCTTTCTACGTAAGAGTCATTCTCAATACTGATGTCAGGA",
             "AATTATACATAGGTGTGACTCTATGCTCGGCTATGGAAATAAGGTTCGCGCCGACACCTATGAAGAATTGTCACCCATGTTTTTGTGTCTATCAGCTTTGAGTGAGATTTGGTTTTCACGGGAGAAAGAGGATGTTCTCTGCGTGCGGACTCCTGAGACTTTGCTGAATGATGATGTAGCGGATCCACGAGGAACTGAGGTCCCGCAGCTCCGAGACAGGTGCTGATGCTTTGGCAACGATTTGAGGGCACAATTCCCGAGTACCTAGGATGGTATTCTGTATTGATTGGTTTTTGAAATGTGCTTGATTCGAACCAAGCGAGCAATTGACAAACGCTGTGCCTAGGTATACCTAAAATAAAACTGCGACAGTTGATCAAACATAAAGTAGAGGGGGTCCAAGTATCCATGAGTGATGCTTAGCACACCCTGCTCCCTGGACTTTTGGATTACCCCCTTCTAGCTTGCTTCTAGCTCAAGCTAAGACCTACCCCAATAAGAGGTAGCTAAGAACGGGGTCTGGGCAGTCATCAACGCCCGTGATCGTAAATCGGTCGTCCCACCGCACTCGCCGCGAATTACGAATAGCCATAGATGAGC",
             "TCTAAAAATGGGGCGGCCAGTGAATAAAGCCTGCGCGTATTCGTAGCTGTTTACTCGGGAGACCGGCGCCCGAACAGCGCCCTGCCTAACGCCAGCTTACACCGATAGACGAACACGGTTGGGCTGATATACGTCGAACCTGCCTAACCTTAATACTTTCCCTAGTCAGAAGTTGGCCCGAACTTAAGCGTTCGAATGTAGGAGGACTATGAGAGCAAAGCGCGCGCCCGGTCATTTGCACAGAATTCACGTATGTAGTGTAGAGGCGAGACGGGTTTGTCGCGTACACTGCAGACCCAACAGTTTTACGGCAACACAATATCCGTCCAGCCGTAATACGAGCGCAAAGCACGTAGGGTCATCTGGCTAAAGAATTAGGCGCCACTCATTTTGACGGAGAGCGCTTTGCGATCAGATCAGTGGAGTCCAGATTTGATTGTAACTCACTTACCGCACGGCAACAACGCTCATTCCCGCTAATGTATGAGGTACAGGTTGCACTGGTCAGTTTAATGAAGGTCATAGAACACGGGTTTACGTGAATGCGTGTCGCCATCCTCGGCCGAAATGATGAGTTGCCAGGACCGATCTGGCGCCAGC",
             "AGGTAAGGCTCGCCTCTACATCTCCGTACAAACTATCAGACGTAAAGAAAGCTGGAGGATTGCCAGCGAAAAGTACATAACACAAAGAACAAAAAGAGAAGGGGGTACGGGCTATTCGATCTAGATGGAGGCTAGGCAATAGAAGTTCGATCATCCATGGTAACTAGATATATGCTGAGAGCAAACGATCCCTAGTACCGCCTGTGTTATATGCCACCAATCTTTCTTCAGTTAGAAACCTCATTGTCGGGCGACACCAGGTCGATTCAAGAGGCGAGAGCCCATATGCTCGACCTATGGCGTGAACGCTAAGCGGCTGGAGCAAGAGAGGTGTATCCAACGACGGTTTTGAATTTACAATTCAGCCCACTGATATAAGCTGTATGGACTGACTCTGGAGGGACGCGCTGATATCTAAGGGCTTCGCGTACTAGGGTCACTACGGAAGCCATCGGCACTGTGCATCTTACAAAACGGACGTCCTTGACGGCCCTATGACCTTAGCACAAACGAATTGATGACCGAATGTACAGTACTTTGTGCTGGCTGAGCACTCCCTACCACGATCCGGCCAGCCGATCTGCGTCGAGGCTGCCACGC",
             "AAGCTCAGCTAACTAGGCGTGAATAATAACGGAACACCTTAGGTAATGTTGGGGTCCTTACCACCATTTTACGTGGATCTCTAGACGGGCAGCACAAGCAGACGCTCAACGTAGTAATGCCAAGAAGAGATCTACTCCTGTGTTCACTTACATATATTCCCACTCAGAACCGCGTCTTCTGAACTGAGGAAGAAGTTACACTAACTGCACGAGATACCGGATCTGCACCTAGCCTGCTAGGCGTGGCACACGTAGCGCACCCTCACGGCTGCAATGGAATTTGCACAAAAACCAGCGCGTGGCGGATATTCCTCGTTTACAGAGTGGGTTGGAACATCCGGCGGTCCCGAGAGAACCGTCTTTCCGGTCGCCCATTTTATCAAAGATTGCAGTCTACTTGCCCGTATTCCTTGAGATGATTCGAAGGTCGAAATCGTAGCACATGGCTAACAATCCTGTTATTTATGCAGTAGCCGCGCCGCTTAGACGGCTTACCCCCGATATAGGGGAGCCCACCAGCTATGCCCTGGAAGGGACGATAAATAGCGTTGTGATTTATGATACCTTCACCAGCTTCGTACGTGCATAGAAAAGGAAGGG",
             "TTCGTATCTTTCTCGGCGCCCTGATTCCAGTGATGGATTGTGAGGTCACTTCAAGTGAGATGTGTATTCCCAGCCAATCTATCCGTGTTAACTGATCCTAAACAGAGTGTGCCCAGATTAATGGGAACCCCAGTGTCAAGCGGGCCCTTAACACGGCCTGGTTAGATTCGTTTTAAGTGGGTCCTCTAACTCCTAACATTTTGACTTAAGGGTTTAACCGCTGACAGGCAGTAGCAACGGCTGTAGGGGAACACGAGGTTTTTTAATAAGTCTTGCAGTTTCATGCGGTTCTCACCAGAACGTTATAATCGCGAGTGCCCCGCTCAGGAATAGGATCAATGACGATTCTTATATCTCCGGAATTATGGTTACAGCTTCGTCAACGGCCTCAGGGTCGGGTTTTAAGCGGGGCCCGTTATCCAGAAATTACCCTGCATGGCAGGTCTACGCTAAAGTCCGAGCAAGAAAAAAGAGAGGAGTTTGCCCACGTGCCGCACACCCGGAGCTAGTCAGCATTGGTCTTCGAGAGATGCTCGCTGGACTCGGTTCATCTACTCGATCTAATTTTATGGCCGCCAACCATCAAAACGTATGACCTAA",
             "TCTCATCCGTAGATTTAGTCCGGAGTGTTGAACAGCCCTCGGAGGTGCTACTAGCAATCACGAGATGCTAACGAGGAATATTTGGGATAGACGGTTCCTTCATGTTGTTCTGGGTACGCACTGCCGGCGAGTACCCCAGTGCCGAAACCGGTAAGAGTAAGTTCCTTAGGTTACGAGATTCCAGGCTTTTTGGGTAAGCGAGACCTACCCACTTGTTGCATCTACCCGTGTCTGTCAATCGCTGACTAGAACTGGTATCACGAGAAGAGAAACTTTCGATCTGTGCCCCATCAGTACCGAAGTTTGGTATAAATCGATGTGATATCCAAGACATGGAATAGCTTTCGCTCTTACGAGAGCATATGAAGGTTGCAACTAATTACCTATCTGATGTACGAAATTCAAGCTAAAGGGGGGTCAATCTCGTCCGAGTGCGACGGGGCAATAGCCCGGTACGATCTCCCATTTTCCCTTCCGGTACTCTACTGCTTTGGCGGGGTCGAGTTATCCGTGCGAACATTCAACCACCTCTGAGAACGGGGCCATAATGAACTGTGATCTTGATTCTACCTAAACACGCAGGACCAAAGCCTTCGCCGA"
             ]
"""

