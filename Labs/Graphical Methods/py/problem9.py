# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough (wschloug)
# Group Members: Junmo, Mikey, Wenyu

# Find the string spelled by a genome path.
# Given: A sequence of k-mers Pattern1, ... , PatternN such that the last k - 1 symbols of Patterni are equal to the first k - 1 symbols of Patterni+1 for i from 1 to n-1.
# Return: A string Text of length k+n-1 where the i-th k-mer in Text is equal to Patterni for all i.

class genomePath:
    '''
    contains methods used to reconstruct a string from its genome path
    '''
    def __init__(self, kPattern):
        '''
        initializes necessary vars

        '''

        self.pattern = kPattern
        self.gps = '' # genomePathString

    def seq(self):
        '''
        build a sequence from given patterns (nodes)
        everything except the last symbol is the same as the previous
        for reference (sample dataset + output):
        ['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT']
        ACCGA A
        ACCGAA G
        ACCGAAG C
        ACCGAAGC T
        '''

        self.gps += self.pattern[0] # think like node zero

        for pattern in range(1, len(self.pattern)): # for each kmer pattern in range 1 (start at 1 dont include first pattern) up to length of all patterns
            self.gps += self.pattern[pattern][-1] # add to string gps (genomePathString) the last element of each pattern
            #print(self.gps)

    def printSeq(self):
        '''
        quick and easy way to print the genomePathString (gps)
        calls seq to build the sequence
        then prints the sequence
        '''

        self.seq()
        print(self.gps)


def main():
    '''
    reads from input file and instantiates patterns into class genomePath, then calls printSeq
    which in turn calls the necessary methods to build a sequence from the given patterns,
    along with printing out the formatted string which has been constructed from its genome path

    '''

    p = []
    inputFile = open('rosalind_ba3b.txt', 'r') # open the input file for reading
    for line in inputFile: # for each line in input file, append it as an element of list 'p' (patterns)
        p.append(line.strip())

    #print(p)
    gPath = genomePath(p)
    gPath.printSeq()


if __name__ == '__main__':
    main()