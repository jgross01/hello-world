import collections

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

def SymbolToNumberRNA(N):
#Convert single base to int A=0, C=1, G=2, T=3.
    if N == "A":
        return 0
    if N == "C":
        return 1
    if N == "G":
        return 2
    if N == "U":
        return 3
    else:
        return "Not a Base"

def NumberToSymbolRNA(r):
#Convert single int to base A=0, C=1, G=2, T=3.
    if r == 0:
        return 'A'
    if r == 1:
        return 'C'
    if r == 2:
        return 'G'
    if r == 3:
        return 'U'

def PatternToNumberRNA(Pattern):
#Convert k-mer to lexicographical int, depends on SymbolToNumber
    if len(Pattern)==0:
        return 0
    symbol = Pattern[len(Pattern)-1]
    prefix = Pattern[:len(Pattern)-1]
    return 4*PatternToNumberRNA(prefix)+SymbolToNumberRNA(symbol)

def NumberToPatternRNA(index,k):
#Convert int representing a k-mer to the DNA code for that k-mer.  Depends on NumberToSymbol
    if k==1:
        return NumberToSymbolRNA(index)
    prefixIndex = index // 4
    r = index % 4
    symbol = NumberToSymbolRNA(r)
    prefixPattern = NumberToPatternRNA(prefixIndex, (k-1))
    return prefixPattern + symbol

def MakeGeneticCode():
    codons = [0]*64
    aminoacids = ["K", "N","K","N","T","T","T","T","R","S","R","S","I","I","M","I","Q","H","Q","H", "P","P","P","P","R","R","R","R","L","L","L","L","E","D","E","D","A","A","A","A","G","G","G","G","V","V","V","V","*","Y","*","Y","S","S","S","S","*","C","W","C","L","F","L","F"]
    geneticcode = []
    for i in range(64):
        codons[i] = NumberToPatternRNA(i,3)
    for i in range(64):
        geneticcode[codons[i]]= aminoacids[i]

def TranslatePeptide(Pattern):
    geneticcode = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAU': 'N', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
                   'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGU': 'S', 'AUA': 'I', 'AUC': 'I', 'AUG': 'M', 'AUU': 'I',
                   'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAU': 'H', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
                   'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R', 'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
                   'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAU': 'D', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
                   'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G', 'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
                   'UAA': '*', 'UAC': 'Y', 'UAG': '*', 'UAU': 'Y', 'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
                   'UGA': '*', 'UGC': 'C', 'UGG': 'W', 'UGU': 'C', 'UUA': 'L', 'UUC': 'F', 'UUG': 'L', 'UUU': 'F'}
    Peptide = ""
    for i in range(0,len(Pattern),3):
        if i+2 <= len(Pattern)-1:
            codon = Pattern[i:i+3]
            if geneticcode[codon] == "*":
                return Peptide
            else:
                Peptide += geneticcode[codon]
        else:
            return Peptide

def TranslateString(Pattern):
    geneticcode = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAU': 'N', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
                   'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGU': 'S', 'AUA': 'I', 'AUC': 'I', 'AUG': 'M', 'AUU': 'I',
                   'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAU': 'H', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
                   'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R', 'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
                   'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAU': 'D', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
                   'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G', 'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
                   'UAA': '*', 'UAC': 'Y', 'UAG': '*', 'UAU': 'Y', 'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
                   'UGA': '*', 'UGC': 'C', 'UGG': 'W', 'UGU': 'C', 'UUA': 'L', 'UUC': 'F', 'UUG': 'L', 'UUU': 'F'}
    Peptide = ""
    for i in range(0,len(Pattern),3):
        if i+2 <= len(Pattern)-1:
            codon = Pattern[i:i+3]
            Peptide += geneticcode[codon]
    return Peptide

def TranscribeString(Pattern):
    return Pattern.replace('T','U')

def ReverseCompliment(Pattern):
    # Converts A<-->T and C<-->G and reverses order
    PatternRC = ""
    i=len(Pattern)
    while i>=1:
        if Pattern[i-1] == 'A':
            PatternRC += 'T'
        if Pattern[i-1] == 'C':
            PatternRC += 'G'
        if Pattern[i-1] == 'G':
            PatternRC += 'C'
        if Pattern[i-1] == 'T':
            PatternRC += 'A'
        i -= 1
    return PatternRC

def ReverseComplimentRNA(Pattern):
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

def EncodingStrings(Text,Peptide):
    Frame1RC = ReverseCompliment(Text)
    #print ("Reverse Complimented. '"'You suck!'"' ")
    encoders = []
    Frames = [Text,Text[1:], Text[2:], Frame1RC, Frame1RC[1:], Frame1RC[2:]]
    #print("Frames set")
    Translated = []
    for i in range(6):
        #print ("Length of frame: " + str(len(Frames[i])))
        Translated.append(TranslateString(TranscribeString(Frames[i])))
        index = 0
        #print ("Translated frame: " + str(i))
        positions = []
        while index < len(Translated[i]):
            index = Translated[i].find(Peptide, index)
            if index == -1:
                break
            positions.append(index)
            index += 1
        for pos in positions:
            if i<=2:
                encoders.append(Frames[i][pos*3:pos*3+(len(Peptide)*3)])
            if i>=3:
                encoders.append(ReverseCompliment(Frames[i][pos*3:pos*3 + (len(Peptide) * 3)]))
        #print("Positions found frame " + str(i)+ ": "+str(len(positions)))
    return encoders

def CyclospectrumAminoAcid(Peptide):
    int_mass_table = {'G': '57', 'A': '71', 'S': '87', 'P': '97', 'V': '99', 'T': '101', 'C': '103', 'I': '113',
                      'L': '113', 'N': '114', 'D': '115', 'K': '128', 'Q': '128', 'E': '129', 'M': '131', 'H': '137',
                      'F': '147', 'R': '156', 'Y': '163', 'W': '186'}
    ExtendedPeptide = Peptide + Peptide[:len(Peptide)-1]
    MassSpec = [0] #Start with no mass and full peptide mass
    PeptideMass=0
    for aa in range(len(Peptide)):
        PeptideMass += int(int_mass_table[Peptide[aa]])
    MassSpec.append(PeptideMass)
    for N in range(1,len(Peptide)):
        for i in range(len(Peptide)):
            SubPeptide = ExtendedPeptide[i:i+N]
            #print (SubPeptide)
            SubPeptideMass = 0
            for aa in range(len(SubPeptide)):
                SubPeptideMass += int(int_mass_table[SubPeptide[aa]])
            #print (SubPeptideMass)
            MassSpec.append(SubPeptideMass)
    MassSpec.sort()
    return MassSpec

def LinearSpectrumAminoAcid(Peptide):
    int_mass_table = {'G': '57', 'A': '71', 'S': '87', 'P': '97', 'V': '99', 'T': '101', 'C': '103', 'I': '113',
                      'L': '113', 'N': '114', 'D': '115', 'K': '128', 'Q': '128', 'E': '129', 'M': '131', 'H': '137',
                      'F': '147', 'R': '156', 'Y': '163', 'W': '186'}
    MassSpec = [0] #Start with no mass and full peptide mass
    PeptideMass=0
    for aa in range(len(Peptide)):
        PeptideMass += int(int_mass_table[Peptide[aa]])
    MassSpec.append(PeptideMass)
    for N in range(1,len(Peptide)):
        for i in range(len(Peptide)- N+1):
            SubPeptide = Peptide[i:i+N]
            #print (SubPeptide)
            SubPeptideMass = 0
            for aa in range(len(SubPeptide)):
                SubPeptideMass += int(int_mass_table[SubPeptide[aa]])
            #print (SubPeptideMass)
            MassSpec.append(SubPeptideMass)
    MassSpec.sort(key=int)
    return MassSpec

def CyclospectrumMass(Peptide):
    mod_int_mass_table = ['57', '71', '87','97', '99','101', '103','113',
                      '114', '115', '128','129', '131', '137',
                      '147','156','163',  '186']
    ExtendedPeptide = Peptide + Peptide[:len(Peptide)-1]
    MassSpec = ['0'] #Start with no mass and full peptide mass
    PeptideMass=0
    for aa in Peptide:
        PeptideMass += int(aa)
    MassSpec.append(str(PeptideMass))
    for N in range(1,len(Peptide)):
        for i in range(len(Peptide)):
            SubPeptide = ExtendedPeptide[i:i+N]
            #print (SubPeptide)
            SubPeptideMass = 0
            for aa in range(len(SubPeptide)):
                SubPeptideMass += int(SubPeptide[aa])
            #print (SubPeptideMass)
            MassSpec.append(str(SubPeptideMass))
    MassSpec.sort(key=int)
    return MassSpec

def Expand(Peptides):
    mod_int_mass_table = ['57', '71', '87', '97', '99', '101', '103', '113', '114', '115', '128', '129', '131', '137',
                          '147', '156', '163', '186']
    NewPeptides = []
    for i in Peptides:
        for aa in mod_int_mass_table:
            stub = list(i)
            stub.append(aa)
            NewPeptides.append(stub)
    return NewPeptides

def Mass(Peptide):
    mass = 0
    for i in Peptide:
        mass+=int(i)
    return str(mass)

def LinearSpectrum(Peptide):
    mod_int_mass_table = ['57', '71', '87', '97', '99', '101', '103', '113',
                          '114', '115', '128', '129', '131', '137',
                          '147', '156', '163', '186']
    ExtendedPeptide = Peptide
    MassSpec = ['0']  # Start with no mass and full peptide mass
    PeptideMass = 0
    MassSpec.append(Mass(Peptide))
    for N in range(1, len(Peptide)):
        for i in range(len(Peptide)- N):
            SubPeptide = Peptide[i:i + N]
            # print (SubPeptide)
            SubPeptideMass = Mass(SubPeptide)
            # print (SubPeptideMass)
            MassSpec.append(SubPeptideMass)
    MassSpec.sort(key=int)
    return MassSpec

def ParentMass(Spectrum):
    return max(Spectrum, key=int)

def CyclopeptideSequencing(Spectrum):
    #Each peptide is reprensented as a list of masses of each aa
    mod_int_mass_table = ['57', '71', '87','97', '99','101', '103','113',
                      '114', '115', '128','129', '131', '137',
                      '147','156','163',  '186']
    Peptides = [[]]
    while Peptides != []:
    #for num in range(2):
        Peptides = Expand(Peptides)
        ##print ("Now testing " + str(Peptides))
        for i in list(Peptides):
            #print("Testing spectrum of " + str(i) + ": " + str(LinearSpectrum(i)))
            if Mass(i) == ParentMass(Spectrum):
                ##print ("I found a match, is it cyclic?")
                ##print(i)
                if CyclospectrumMass(i) == Spectrum:
                        print(("-".join(map(str, i))), end=' ')
                Peptides.remove(i)
            else:
                Test = list(Spectrum)
                ##print (Test)
                for m in LinearSpectrum(i):
                    try:
                        #print ("Is " + str(m) + " present?")
                        Test.remove(m)
                        #print ("Yes")
                    except ValueError:
                        #print ("No")
                        Peptides.remove(i)
                        #print(str(i) + " is out! Which leaves:" + str(Peptides))
                        break

def CyclopeptideScoring(Peptide,Spectrum):
    #Counts how many times the masses in the theoretical spectrum of peptide appear in Spectrum, outputs an int score (higher is better)
    Theoretical = CyclospectrumMass(Peptide)
    Score = 0
    for i in Spectrum:
        if i in Theoretical:
            Score += 1
            Theoretical.remove(i)
            #print (Score)
    return Score

def LinearPeptideScoringAminoAcid(Peptide,Spectrum):
    # Counts how many times the masses in the theoretical spectrum of peptide appear in Spectrum, outputs an int score (higher is better)
    Theoretical = LinearSpectrumAminoAcid(Peptide)
    Score = 0
    for i in Spectrum:
        if int(i) in Theoretical:
            Score += 1
            Theoretical.remove(int(i))
    #print (Score)
    return Score

def LinearPeptideScoring(Peptide,Spectrum):
    # Counts how many times the masses in the theoretical spectrum of peptide appear in Spectrum, outputs an int score (higher is better)
    Theoretical = LinearSpectrum(Peptide)
    #print(Theoretical)
    Score = 0
    for i in Spectrum:
        if i in Theoretical:
            Score += 1
            Theoretical.remove(i)
    #print (Score)
    return Score

def TrimAminoAcid(Leaderboard, Spectrum, N):
    #N is number of top scores including ties (top 10 scores might have 20 entries, but only 10 discrete scores)
    Leaders= sorted(Leaderboard, key = lambda leader: -LinearPeptideScoringAminoAcid(leader,Spectrum))
    LinearScores = []
    for i in Leaders:
        LinearScores.append(LinearPeptideScoringAminoAcid(i, Spectrum))
        #print (LinearScores)
    for i in range(N, len(Leaderboard)):
        #print (str(LinearScores[i]) + '<' + str(LinearScores[N-1]))
        if LinearScores[i] < LinearScores[N-1]:
            #print (i)
            Leaders = Leaders[0:i]
            return Leaders
    return Leaders

def Trim(Leaderboard, Spectrum, N):
    #N is number of top scores including ties (top 10 scores might have 20 entries, but only 10 discrete scores)
    #print ("Let's trim!")
    Leaders= sorted(Leaderboard, key = lambda leader: -LinearPeptideScoring(leader,Spectrum))
    LinearScores = []
    #print (Leaders)
    for i in Leaders:
        LinearScores.append(LinearPeptideScoring(i, Spectrum))
        #print (LinearScores)
    for i in range(N, len(Leaderboard)):
        #print (str(LinearScores[i]) + '<' + str(LinearScores[N-1]))
        if LinearScores[i] < LinearScores[N-1]:
            #print (i)
            Leaders = Leaders[0:i]
            return Leaders
    return Leaders

def LeaderboardCyclopeptideSequencing(Spectrum, N):
    Leaderboard = [[]]
    LeaderPeptide = []
    while Leaderboard != []:
        Leaderboard = Expand(Leaderboard)
        #print ("Now testing " + str(Leaderboard))
        for i in list(Leaderboard): #list(leaderboard to make new instance so we can modify original list)
            #print("Score of " + str(i) + ": " + str(LinearPeptideScoring(i, Spectrum)))
            #print (Mass(i))
            #print (ParentMass(Spectrum))
            if int(Mass(i)) == int(ParentMass(Spectrum)):
                #print (str(i) + " is a potential high score")
                #print (CyclopeptideScoring(i, Spectrum))
                if CyclopeptideScoring(i, Spectrum)> CyclopeptideScoring(LeaderPeptide, Spectrum):
                    #print (str(i) + " is a new leader")
                    LeaderPeptide = i
            elif int(Mass(i)) > int(ParentMass(Spectrum)):
                #print ("Too Heavy!")
                Leaderboard.remove(i)
        #print (Leaderboard)
        Leaderboard = Trim(Leaderboard, Spectrum, N)
        #print ('New leaderboard:' + str(Leaderboard))
    return LeaderPeptide

def SpectralConvolution(Spectrum):
    Convolution = []
    for i in Spectrum:
        for j in Spectrum:
            if int(i) > int(j):
                Convolution.append(int(i)-int(j))
    print ("Convolution Calculated")
    return Convolution

def NewAminoAcidAlphabet(Convolution, M):
    Candidates = []
    for i in Convolution:
        if int(i)>= 57 and int(i)<=200:
            Candidates.append(i)
    #print (Candidates)
    Count = collections.Counter(Candidates)
    #print (list(Count.values())[0])
    TopM = Count.most_common(M)
    #print (TopM)
    mincount = list(map(min,(zip(*TopM))))[1]
    #print (mincount)
    for i in range(len(list(Count.values()))):
        if list(Count.values())[i] == mincount:
            TopM.append((list(Count.keys())[i],mincount))
    TopM = list(set(TopM))
    print ("New Alphabet Written")
    return list(list((zip(*TopM)))[0])

def RestrictedAlphaExpand(Alphabet,Peptides):
    NewPeptides = []
    for i in Peptides:
        for aa in Alphabet:
            stub = list(i)
            stub.append(str(aa))
            NewPeptides.append(stub)
    print ("Expanded")
    return NewPeptides

def ConvolutionLeaderboardCyclopeptideSequencing(N, M, Spectrum):
    Alphabet = NewAminoAcidAlphabet(SpectralConvolution(Spectrum),M)
    Leaderboard = [[]]
    LeaderPeptide = []
    while Leaderboard != []:
        Leaderboard = RestrictedAlphaExpand(Alphabet, Leaderboard)
        print ("Now testing " + str(Leaderboard))
        for i in list(Leaderboard): #list(leaderboard to make new instance so we can modify original list)
            #print("Score of " + str(i) + ": " + str(LinearPeptideScoring(i, Spectrum)))
            #print (Mass(i))
            #print (ParentMass(Spectrum))
            if int(Mass(i)) == int(ParentMass(Spectrum)):
                #print (str(i) + " is a potential high score")
                #print (CyclopeptideScoring(i, Spectrum))
                if CyclopeptideScoring(i, Spectrum)> CyclopeptideScoring(LeaderPeptide, Spectrum):
                    #print (str(i) + " is a new leader")
                    LeaderPeptide = i
            elif int(Mass(i)) > int(ParentMass(Spectrum)):
                #print ("Too Heavy!")
                Leaderboard.remove(i)
        #print (Leaderboard)
        Leaderboard = Trim(Leaderboard, Spectrum, N)
        #print ('New leaderboard:' + str(Leaderboard))
    return LeaderPeptide


#with open('input.txt', 'r') as myfile: spectrum = myfile.read().split(' ')


#print ((" ".join(map(str, LinearSpectrumAminoAcid(('LCVAFSNSPCTTDNKFCDRVQAHVNFKLEMRDPAYYEWMNGTVF'))), )))
#output = EncodingStrings(pattern, "QIQVLEG")
#print (len(output))
Spectrum = [57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228, 265, 285, 299, 307, 323, 356, 364, 394, 422, 493]
output = (ConvolutionLeaderboardCyclopeptideSequencing(60,20,Spectrum))
print (output)



#0 71 101 113 131 184 202 214 232 285 303 315 345 416
