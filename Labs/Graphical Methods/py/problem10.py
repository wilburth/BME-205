# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough (wschloug)
# Group Members: Junmo, Mikey, Wenyu

# In this chapter, we use the terms prefix and suffix to refer to the first k − 1 nucleotides and last k − 1 nucleotides of a k-mer, respectively.
# Given an arbitrary collection of k-mers Patterns, we form a graph having a node for each k-mer in Patterns and connect k-mers Pattern and Pattern' by a directed edge if Suffix(Pattern) is equal to Prefix(Pattern').
# The resulting graph is called the overlap graph on these k-mers, denoted Overlap(Patterns).

class overlapGraph:
    '''
    contains methods used to form an adjacency list (overlap graph) from a collection of k-mers
    '''

    def __init__(self, kmers):
        '''
        initializes necessary vars
        '''

        self.kmer = kmers
        self.adjacencyList = []

    def suffix(self, pattern):
        '''
        last k − 1 nucleotides of k-mer
        '''

        return pattern[1:len(pattern)] # cut from 2nd element up to last element, not including 1st element (in this case element 0)

    def prefix(self, pattern):
        '''
        first k − 1 nucleotides of k-mer
        '''

        return pattern[:len(pattern) - 1] # cut from start up to but not including last element

    def graph(self):
        '''
        we form a graph having a node for each k-mer in Patterns
        and connect k-mers Pattern and Pattern' by a directed edge if Suffix(Pattern) is equal to Prefix(Pattern')
        '''

        for k1 in self.kmer: # kmer
            for k2 in self.kmer: # kmer'
                if self.suffix(k1) == self.prefix(k2): # if suffix(pattern) == prefix(pattern')
                    self.adjacencyList.append((k1, k2)) # add it to list

    def print(self):
        '''
        calls self.graph to build the adjacency list
        adjacency list is then sorted so that it's in alphabetical order
        adjacency list is then printed as -- (tail -> head), sorted in alphabetical order
        '''

        self.graph()
        self.adjacencyList.sort(key=lambda x: x[0]) # sort adj list in order of descending element 0 (so it starts w/ 'aaaaaaaaaaaaaaaaa')
        for i in self.adjacencyList: # for element in adj list
            print(i[0] + ' -> ' + i[1]) # print contents of adjacency list in the form of an adjacency list (node -> node)



def main():
    '''
    reads from input file and instantiates patterns into class overlapGraph, then calls print() which in turn calls
    graph(), which builds the graph from the given patterns (nodes), and finally prints a sorted list in the
    specified format (tail->head)

    rosalind prompt --
    "Given: A collection Patterns of k-mers."
    "Return: The overlap graph Overlap(Patterns), in the form of an adjacency list."
    '''

    p = []
    inputFile = open('rosalind_ba3c.txt', 'r') # open the input file for reading
    for line in inputFile: # for each line in input file, append it as an element of list 'p' (patterns)
        p.append(line.strip())

    oGraph = overlapGraph(p)
    oGraph.print()

if __name__ == '__main__':
    main()

