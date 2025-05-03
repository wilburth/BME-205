# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough (wschloug)
# Group Members: Junmo, Mikey, Wenyu

# Given an arbitrary collection of k-mers Patterns (where some k-mers may appear multiple times),
# we define CompositionGraph(Patterns) as a graph with |Patterns| isolated edges.
# Every edge is labeled by a k-mer from Patterns,
# and the starting and ending nodes of an edge are labeled by the prefix and suffix of the k-mer labeling that edge.
# We then define the de Bruijn graph of Patterns, denoted DeBruijn(Patterns),
# by gluing identically labeled nodes in CompositionGraph(Patterns), which yields the following algorithm.


class deBruijn:
    '''
    contains methods used to construct a deBruijn graph from a collection of k-mers
    '''

    def __init__(self, patterns):
        '''
        initializes necessary vars
        '''

        self.p = patterns

        self.edges = [] # list housing edges
        self.dbDict = {} # dictionary housing (k-1)-mers


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


    def nodes(self):
        '''
        creates the edges for each k-mer and adds them to self.edges, which is a list housing all edges
        in a format --> [(prefix(pattern), suffix(pattern)), (prefix(pattern), suffix(pattern)), etc. . . ]
        '''

        for pattern in self.p:
            self.edges.append((self.prefix(pattern), self.suffix(pattern))) # add prefix of pattern + suffix of pattern to the edge list --> each element in self.edges is (all but last element of pattern, all but first element of pattern)

    def connect(self):
        '''
        glue all nodes with identical labels, yielding the graph DeBruijn(Patterns)
        '''

        self.nodes() # call self.nodes to fill out self.edges with (prefix(pattern), suffix(pattern))

        for i in range(len(self.edges)): # search across entire list of edges
            node = self.edges[i][0] # look at incoming edge at specified pos i, set it to be curr node

            if node not in self.dbDict.keys(): # if curr node not yet seen in our deBruijn dictionary
                self.dbDict[node] = [] # make it a key in deBruijn dictionary with value empty list

                for edge in self.edges: # for all edges
                    if edge[0] == node: # if incoming edge == node (which has been set to incoming edge at pos i)
                        self.dbDict[node].append(edge[1]) # append to curr key node in deBruijn dictionary the outgoing edge


    def print(self):
        '''
        calls self.connect() to create the deBruijn graph,
        then prints the deBruijn graph (dictionary) in specified format -- (tail->head)

        https://www.geeksforgeeks.org/python-concatenate-dictionary-value-lists/
        https://www.learnbyexample.org/python-string-join-method/ --> join() on dictionary section is super useful
        '''

        self.connect()
        for key in sorted(self.dbDict.keys()):
            print(key + ' -> ' + ','.join(str(i) for i in sorted(self.dbDict[key])))

def main():
    '''


    '''
    p = []
    inputFile = open('rosalind_ba3e.txt', 'r') # open the input file for reading
    for line in inputFile: # for each line in input file, append it as an element of list 'p' (patterns)
        p.append(line.strip())

    db = deBruijn(p)
    db.print()

if __name__ == '__main__':
    main()