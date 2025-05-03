# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough (wschloug)
# Group Members: Junmo, Mikey, Wenyu


########################################################################
# FastAreader
# provided by Dr. B
########################################################################
import sys

class FastAreader :

    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''

        self.fname = fname

    def doOpen (self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta (self):

        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence


########################################################################
# CommandLine
# provided by Dr. B
########################################################################

class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse,
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination exception Usage.

    attributes:
    all arguments received from the commandline using .add_argument will be
    available within the .args attribute of an object instantiated from CommandLine.
    For example, if thisCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then thisCommandLine.args.requiredbool will
    name that option.

    '''

    def __init__(self, inOpts=None) :
        '''
        CommandLine constructor.

        Implement a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does',
                                              epilog = 'Program epilog - some other stuff you feel compelled to say',
                                              add_help = True, #default is True
                                              prefix_chars = '-',
                                              usage = 'python randomizedMotifSearch.py -i int -k int -p float < input.fa > output.out'
                                              )

        # options as specified by assignment
        self.parser.add_argument('-i', type=int, action='store', help='number of iterations in integer format')
        self.parser.add_argument('-k', type=int, action='store', help='motif length in integer format')
        self.parser.add_argument('-p', type=float, action='store', help='pseudo-count in float format')

        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)


class Usage(Exception):
    '''
    Signal a Usage error, evoking a usage statement and eventual exit when raised.
    '''
    def __init__(self, msg):
        self.msg = msg

########################################################################
# The Bioinformatics Class --> randomizedMotifSearch
# Here is where most of your solution will end up.
#
########################################################################

import random  # needed for randint
import math  # needed for calculating entropy function (log2)

class consensusCRISPR:

    '''

    class consensusCRISPR contains methods that are used to calculate a consensus motif sequence and entropy score for
    a given set of sequences, which in our case are housed in a fastA file.

    build the thing that needs doing

    things i want // overarching structure ----------------------------------------------------------

    initialization of user inputs   -- dna sequence, kmer size, pseudocount size, run size
    generate random motifs          -- can use random.() python library
    generate counts from motifs     -- include pseudocounts
    generate profile from counts    -- probability distribution -- DON'T NEED THIS
    generate most probable motifs from profile distribution
    generate entropy score          -- use equation discussed in class, textbook, and on assignment webpage
    generate a consensus            -- can call main algorithm from this function
    main algorithm found in textbook -- describing randomized motif search



    '''
    def __init__(self, DNA, kmer, pseudoCount, runs):
        '''

        DNA: list of seq taken from fastAreader --> LIST
        kmer: kmer (motif) length --> INT
        pseudoCount: psuedo-count --> FLOAT
        runs: number of iterations that we will run RMF --> INT
        bestMotifs: stores the most probable motifs
        highScore: global int

        '''

        self.dna = DNA
        self.k = kmer
        self.psc = pseudoCount
        self.i = runs

        self.bestMotifs = []
        self.highScore = 1000000000000000000

    def randomMotifs(self):
        '''

        take in -- list of sequences
        produce -- randomly generated motifs of length k (kmer), specified by user

        https://www.w3schools.com/python/ref_random_randint.asp
        "The randint() method returns an integer number selected element from the specified range."

        '''

        randomMotifs = [] # randomMotifs is a list that will be filled with randomly generated motifs of length k

        for sequence in self.dna: # for each sequence
            startP = random.randint(0, len(sequence) - self.k) # create a random starting point (startP)
            randomMotifs.append(sequence[startP: startP + self.k]) # append a cut from startP up to startP + length of kmer

        return randomMotifs # return list randomMotifs --> we now have a bunch of random motifs (of length k) to use

    def counts(self, motifs):
        '''

        take in -- list of motifs
        produce -- dictionary containing counts of nucs in positions across motif

        '''

        # dictionary of lists -- simulating matrix
        countDictionary = {'A': [], 'T': [], 'G': [], 'C': []} # countDictionary will host all counts of ATGC

        # we utilize our pseudocount variable here to increment the counts
        inc = self.psc

        for i in range(self.k): # for i in range 0 up to length of kmer
            columnX = [column[i] for column in motifs]

            countDictionary['A'].append(columnX.count('A') + inc) # add pseudocount
            countDictionary['T'].append(columnX.count('T') + inc) # add pseudocount
            countDictionary['G'].append(columnX.count('G') + inc) # add pseudocount
            countDictionary['C'].append(columnX.count('C') + inc) # add pseudocount

        #print("count dictionary")
        #print(countDictionary)
        return countDictionary

    def motif(self, counts):
        '''

        take in -- profile matrix for set of motifs
        produce -- matrix of most probable motifs

        https://www.w3schools.com/python/ref_keyword_none.asp
        "The None keyword is used to define a null value, or no value at all.
        None is not the same as 0, False, or an empty string.
        None is a data type of its own (NoneType) and only None can be None."

        https://stackoverflow.com/questions/521674/initializing-a-list-to-a-known-number-of-elements-in-python

        '''

        # new motif matrix
        # getting error "list assignment index out of range" unless I use [None] instead of []
        # i try [1] and it also seems to work but for now I stick with none
        motifMatrix = [None]*len(self.dna) # none inside []? trying to create a list with an element per element in dnaseq

        for i in range(len(self.dna)):  # search across length of dna
            sequence = self.dna[i]
            highProb = 0  # initialize variable for use of highest probability for each motif

            for x in range(len(sequence) - self.k + 1):  # search across range of sequence minus kmer size (+1 necessary to retain proper range)
                currMotif = sequence[x: x + self.k]  # store current motif that we are viewing
                probability = 1  # set probability

                for y in range(self.k): # search across length of kmer
                    nuc = currMotif[y]  # for each nucleotide base in the current motif of length kmer
                    probability *= counts[nuc][y] # generate a probability of the nuc at its location using profile matrix

                if probability > highProb: # starts true (we set highProb to 0)
                    motifMatrix[i] = currMotif # keep this motif since its probability was higher than our highest probability!
                    highProb = probability # set our highest probability to our new highest probability

        #print("motif matrix")
        #print(motifMatrix)
        return motifMatrix

    def entropy(self, motif):
        '''

        takes in -- motif matrix, which then has counts() called on it
        produce -- entropy score for motifs

        '''

        profile = self.counts(motif) # call profile on motif-input
        score = 0 # base score temp var

        totalCount = sum([seq[0] for seq in profile.values()]) # we need a total count of all values in the sequence (i.e total count of ALL ATGC)

        for nuc in profile.keys(): # search across each nucleotide
            for pos in profile[nuc]: # search across each position in matrix
                if pos != 0: # we only look for positions with values > 0
                    score += -1*((pos/totalCount) * math.log2(pos/totalCount)) # accumulate score using summation equation

        #score = score * -1 # finish summation equation

        #print("score")
        #print(score)
        return score

    def consensus(self):
        '''

        produce -- consensus motif (from best scoring motifs from randomizedMotifSearch())

        https://www.w3schools.com/python/ref_func_max.asp

        '''

        self.randomizedMotifSearch() # call rms() to calculate best motifs
        consensus = ''

        #print(self.bestMotifs)
        for i in range(self.k): # search across kmer
            # creates a list with each position being the ith element from each motif (does this for each .bestMotif)
            # example: we have 10 5-mers, each starting with A
            # 1st list created would house 10 A's
            # 2nd list created would house the second element from each motif, etc. (we would have 5 lists for this ex.)
            nuc = [base[i] for base in self.bestMotifs]
            #print(nuc)
            # we can then find the base that occurs the most using max() on the list, which will determine the base
            # that should be added to our consensus sequence
            # example continued: if there were 10 A's in the 1st list, then we add A to our consensus sequence
            # if there were 6 A's and 4 T's in the 2nd list, then we add another A to our consensus sequence, etc.
            consensus += max(nuc, key=nuc.count)


        #print("consensus")
        #print(consensus)
        return consensus # return consensus sequence


    def randomizedMotifSearch(self):
        '''

        produce -- best scoring motifs

        pseudocode from textbook is as follows ----------------

        randomly select k-mers motifs
        bestMotifs <-- motifs
        while true:
            profile
            motifs

            if SCORE(motifs) < SCORE(bestMotifs)
                bestMotifs <-- motifs
            else
                return bestMotifs


        tips from dr b. -----
        "Check the way you are comparing final profile scores across trajectories.
        Remember that we want the best score and the associated consensus that came from that best score.
        This will not likely be from the last iteration (trajectory) that we ran"

        '''

        for i in range(self.i): # search over all iterations, -i=INT defined by user in commandLine
            bestMotifs = self.randomMotifs() # set bestMotifs var to the random motifs we calculated with randomMotifs()

            while True:
                count = self.counts(bestMotifs) # create profile from random motifs
                currMotif = self.motif(count) # calculate new motifs based on profile

                currMotifScore = self.entropy(currMotif) # calculate entropy for current motifs
                bestMotifScore = self.entropy(bestMotifs) # calculate entropy for previous motifs

                if currMotifScore < bestMotifScore: # compare scoring, if current motifscore is less than best motif score
                    bestMotifs = currMotif # set previous equal to current

                    if currMotifScore < self.highScore: # if current score is less than global variable highScore
                        self.highScore = currMotifScore # set global variable highScore equal to current score
                        self.bestMotifs = bestMotifs # set best motif equal to bestMotif (which has been set to current motif)

                else: # otherwise (we finished!)
                    break # exit the loop

        return self.bestMotifs # return our best motifs


########################################################################
# Main
# Here is the main program
# This is a place to instantiate objects, read and write some data for
# those objects and maybe implement some of the command line options.
########################################################################

def main(myCommandLine=None):
    '''
    Implement the Usage exception handler that can be raised from anywhere in process.

    '''


    if myCommandLine is None:
        myCommandLine = CommandLine()  # read options from the command line
    else :
        myCommandLine = CommandLine(myCommandLine) # interpret the list passed from the caller of main as the commandline.


    seqL = []
    myReader = FastAreader()
    for head, seq in myReader.readFasta():
        seqL.append(seq)

    finalConsensus = consensusCRISPR(DNA=seqL,
                                     kmer=myCommandLine.args.k,
                                     pseudoCount=myCommandLine.args.p,
                                     runs=myCommandLine.args.i)

    # i built it so that consensus() runs in tandem with randomizedMotifSearch() --
    # calling .consensus() will call randomizedMotifSearch(), which uses all other methods to return
    # the best motifs, which consensus then sorts through to find the max occurrence of each nucleotide base,
    # and then builds the final sequence. Additionally, highscore has been calculated, but we want it in a string
    # format, which is why the second print line has the str() function.
    print("Final sequence consensus: " + finalConsensus.consensus())
    print("Final entropy score: " + str(finalConsensus.highScore))




if __name__ == "__main__":
    main()


########################################################################
# inspection markdown
# done by: Junmo
# Neat and well visualized code. Initial motif is random, I can get a profiling score and the consensus, too.
# I like the way you iterate by making one more function.
# But I got an entropy score of 16 from your code in which I think it's a little bit higher.
#
# WS -- when only using one of the "".fa files, the entropy score IS higher, but when running it with all of the
# "".fa files combined, my entropy score is around 8.8, which is good when compared to others.
#
########################################################################
# Inspector: Michael Collins
# Instead of 'take in' and 'produce' in docstrings use 'input' and 'output', and include the respective datatypes
#
# In Motif mehtod, you can used probablility *= instead of probablility = probablilty *
#
# could use count based profile instead of probablity based profile since you don't calcualte realative entorpy very often
#
# WS -- changes were made based on recommendations, code was shortened over 10% by removing profile() method,
# with same result still being achieved
########################################################################
# inspection markdown
# done by: Wenyu Liang
# I guess it is not neccesary to check whether current motifscore is less than best motif score (in line #347) because when
# we use randomized motif search, the RE will keep becoming greater
# line 304 using key function in max is cool!


'''
    This method has been removed as it is completely unnecessary (recommended by dr. b) -- we can just use counts() 
    and when needed, we can divide by total count to get a decimal probability -- this happens in entropy().
    
    def profile(self, motifs):

        take in -- list matrix of motifs
        produce -- dictionary containing probability distribution of nucs in each position

        https://www.programiz.com/python-programming/methods/built-in/sum
        "The sum() function adds start and items of the given iterable from left to right."


        # dictionary of lists -- simulating matrix
        profileDictionary = {'A': [], 'T': [], 'G': [], 'C': []} # profileDictionary will host all counts of ATGC

        counts = self.counts(motifs) # call function counts() with same motifs

        totalCount = sum([seq[0] for seq in counts.values()]) # we need a total count of all values in the sequence (i.e total count of ALL ATGC)

        for key in counts.keys(): # for key in counts
            for i in range(len(counts[key])): # for i in range 0 up to length of counts of specified key
                profileDictionary[key].append(counts[key][i] / totalCount) # append to our dictionary the counts of the key at position i, divided by total count of keys

        #print("profile dictionary")
        #print(profileDictionary)
        return profileDictionary
    '''