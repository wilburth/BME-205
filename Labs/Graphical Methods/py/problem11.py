# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough (wschloug)
# Group Members: Junmo, Mikey, Wenyu

# Given a genome Text, PathGraphk(Text) is the path consisting of |Text| - k + 1 edges,
# where the i-th edge of this path is labeled by the i-th k-mer in Text
# and the i-th node of the path is labeled by the i-th (k - 1)-mer in Text.
# The de Bruijn graph DeBruijnK(Text) is formed by gluing identically labeled nodes in PathGraphk(Text).


class deBruijn:
    '''
    contains methods used to construct the deBruijn graph from a string
    '''

    def __init__(self, k, text):
        '''
        initializes necessary vars
        '''

        self.k = k # houses kmer size, we will use (k-1)-mer
        self.text = text # houses input seq
        self.dbDict = {} # houses deBruijn dictionary

    def nodes(self):
        '''
        find (k-1)-mer, and for each (k-1)-mer next (k-1)-mer will be shifted once to the right
        all of these connecting (k-1)-mers are added to the deBruinDictionary (dbDict)

        "the i-th node of the path is labeled by the i-th (k - 1)-mer in Text"

        https://stackoverflow.com/questions/3199171/append-multiple-values-for-one-key-in-a-dictionary
        '''

        for i in range(len(self.text) - self.k + 1): # path consisting of |Text| - k + 1 edges
            node = self.text[i: i + self.k - 1] # node from seq cut from i to i plus length of (k-1)mer
            nodeNext = self.text[i + 1: i + self.k] # node from seq cut from i+1 (shift right once) to i plus length of (k-1)-mer --> when i was doing i + self.k -1 my output was only 2/3 characters . . .? --> im dumb its cuz i have to shift right once XD

            if node not in self.dbDict.keys(): # if node isn't yet seen in deBruijnDictionary
                self.dbDict[node] = [] # make it a key with value empty list

            self.dbDict[node].append(nodeNext) # at that key, append next node to the list


    def print(self):
        '''
        calls self.nodes() which will build the deBruijn dictionary
        deBruijn dictionary -- {'key' : ['value']} -- {'tail' : ['head']}
        then prints (key -> key's value) -- there are instances where a key could have more than 1 value, in which case the values are also printed in alphabetical order

        "the edges may be given in any order, as long as the tail of each edge is provided in column 1"

        https://www.geeksforgeeks.org/python-concatenate-dictionary-value-lists/
        https://www.learnbyexample.org/python-string-join-method/ --> join() on dictionary section is super useful
        '''

        self.nodes()

        for key in sorted(self.dbDict.keys()): # for each key in deBruijnDictionary (dbDict)
            print(key + ' -> ' + ', '.join(str(i) for i in sorted(self.dbDict[key]))) # only keys with multiple values separated by a comma will have the comma printed in output



def main():
    '''
    reads from input file and instantiates k-mer size and sequence into class deBruijn, then calls print()
    which in turn calls nodes() which builds the deBruijn dictionary, and then prints it in the specified
    format (tail->head)

    '''

    count = 0
    k = 0 # line 0 of given dataset is the kmer size
    seq = '' # line 1 is the string

    inputFile = open('rosalind_ba3d.txt', 'r') # open the input file for reading
    for line in inputFile:
        if count == 0: # line0
            k = int(line.strip()) # get the kmer number (which will dictate what size kmer we are finding in the seq)
        else: # line1
            seq = line.strip() # get the sequence
        count += 1

    db = deBruijn(k, seq)
    db.print()

if __name__ == '__main__':
    main()

