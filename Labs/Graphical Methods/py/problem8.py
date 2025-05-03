# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough (wschloug)
# Group Members: Junmo, Mikey, Wenyu

# Given a string Text, its k-mer composition Compositionk(Text) is the collection of all k-mer substrings of Text (including repeated k-mers).
# For example, Composition3(TATGGGGTGC) = {ATG, GGG, GGG, GGT, GTG, TAT, TGC, TGG}
# Note that we have listed k-mers in lexicographic order (i.e., how they would appear in a dictionary) rather than in the order of their appearance in TATGGGGTGC.
# We have done this because the correct ordering of the reads is unknown when they are generated.


class kmerComposition:
    '''
    contains methods used to find k-mer composition of string
    '''

    def __init__(self, k, text):
        '''
        initializes necessary vars

        "given an integer k and a string text"
        '''

        self.k = k # int k
        self.text = text # string text
        self.kComp = [] # kmer composition list

    def comp(self):
        '''
        generates strings of length kmer,
        and returns a filled kComp (a dictionary which houses all kmers)
        '''

        for i in range(len(self.text) - self.k + 1): # for i in range 0 up to length of text - kmersize + 1
            # kmer = self.text[i:i+self.k]
            # self.kComp.append(kmer)
            self.kComp.append(self.text[i: i+self.k])  # concatenate two previous lines, which append a cut of size kmer to self.kComp

        return self.kComp # return the filled kComp dictionary


    def printComp(self):
        '''
        quick method to print all kmers that are housed in kComp
        does this by first calling comp() to fill out dictionary
        then prints all kmers in filled dictionary
        '''
        self.comp()
        for kmer in sorted(self.kComp): # for kmer in kComp dictionary
            print(kmer) # print da kmer


def main():
    '''
    reads from an input file and calls necessary methods after initializing with class kmerComposition
    currently unable to print to output file - still learning best way to do so
    '''

    count = 0
    k = 0 # line 0 of given dataset is the kmer size
    seq = '' # line 1 is the string

    inputFile = open('rosalind_ba3a.txt', 'r') # open the input file for reading
    for line in inputFile:
        if count == 0: # line0
            k = int(line.strip()) # get the kmer number (which will dictate what size kmer we are finding in the seq)
        else: # line1
            seq = line.strip() # get the sequence
        count += 1

    kComp = kmerComposition(k, seq) # initialize our vars
    kComp.printComp() # call printComp(), which also calls comp() to fill in kmer list




if __name__ == '__main__':
    main()