# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough (wschloug)
# Group Members: Junmo, Mikey, Wenyu

# Reconstruct a string from its k-mer composition.
# Given: An integer k followed by a list of k-mers Patterns.
# Return: A string Text with k-mer composition equal to Patterns.
# (If multiple answers exist, you may return any one.)




class reconstruct:
    '''
    contains methods used to reconstruct a string from its k-mer composition
    '''

    def __init__(self, pattern):
        '''

        '''

        self.p = pattern
        self.edges = []
        self.nodes = {} # this will be our node dictionary --> housing {key: [input node], [output node]}
        self.walkpath = [] # keeps track of visited
        self.seq = ''

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

    def edgeFind(self):
        '''
        from problem10 -- we form a graph having a node for each k-mer in Patterns
        and connect k-mers Pattern and Pattern' by a directed edge if Suffix(Pattern) is equal to Prefix(Pattern')
        '''

        for k1 in self.p:
            for k2 in self.p:
                if self.suffix(k1) == self.prefix(k2): # if suffix(pattern) == prefix(pattern')
                    self.edges.append((k1, k2)) # add it to list

    def nodeFind(self):
        '''
        uses filled edgeList (self.edges) and fills self.nodes, which is a dictionary that houses a node as its key,
        along with a list as its value. The list contains the input (incoming) node, and the output (outgoing) node
        for the key node.
        '''

        self.edgeFind()

        for edge in self.edges: # for edge in edge list
            node1 = edge[0] # starting-from
            node2 = edge[1] # going-to

            if node1 not in self.nodes.keys(): # if starting node not yet seen in keys
                self.nodes[node1] = ([], []) # ([in], [out]) make it a key with value (2 lists for input + output)

            self.nodes[node1][1].append(node2) # add the output node to the output node list for the starting-from node

            if node2 not in self.nodes.keys(): # if going to node not yet seen in keys
                self.nodes[node2] = ([], []) # ([in], [out]) make it a key with value (2 lists for input + output)

            self.nodes[node2][0].append(node1) # add the input node to the input node list for the going-to node

    def startNode(self):
        '''
        finds a start node -- this start node MUST follow the rule that (# of inputNode) - (# of outputNode) == 1
        returns that start node

        https://www.youtube.com/watch?v=8MpoO2zA2l4 --> video i could found for additional information on eularian path
        '''

        startNode = int() # initialize startnode var

        for node in self.nodes.keys(): # for key (node) in keys
            if len(self.nodes[node][1]) - len(self.nodes[node][0]) == 1: # if (# of inputNode) - (# of outputNode) == 1
                startNode = node # choose it as our starting node

        return startNode # return the starting node that we count

    def path(self, visited):
        '''
        recursive implementation of eularianCycle. uses a start node, which follows the rule (# of output) - (# of input) == 1,
        and from this start node we begin the recursive method. the visited (startnode) will have its outgoing nodes
        set as the 'visited edge', and while the length of this edge (technically a List with multiple elements) is > 0
        do some stuff. this 'stuff' follows the textbook, where we set a newVisit (which will be the first edge of the outgoing nodes)
        and then remove that newVisit from the original list (this allows us to recurse properly; ex 6->3,7 -- 3 gets removed but we will
        still recurse with 7 once we complete the walkpath of 3).

        algorithm based on ch3 p146 of textbook

        visited1 6
        [3, 7]
        visited1 3
        [0, 4]
        visited1 0
        [2]
        visited1 2
        [1]
        visited1 1
        [3]
        visited1 3
        [4]
        visited1 4
        []
        4 will be the first element inserted into walkpath, followed by 3, 1, 2, 0, 3, 6. at this point, we have another
        outgoing edge from node 6 (7) that will continue the recursion until exit.

        https://www.geeksforgeeks.org/eulerian-path-and-circuit/ --> needed a bit help for the recursive method, might try keeping
        track of nodes without the recursive method (think dr. b suggested this way?)
        '''

        visitedEdge = self.nodes[visited][1] # retrieves outgoing node from startNode

        while len(visitedEdge) > 0:     # while the length of visited node is greater than 0 (in this case, visited edge is a list containing elements)
            newVisit = visitedEdge[0]   # retrieves the first element in list of visitedEdge (ex. 3 -> 0,4  newVisit set as 0, not 4) -- is this supposed to be random choice???
            self.nodes[visited][1].remove(newVisit)     # remove the newVisit element from outgoing nodes from original call
            self.path(newVisit)     # then, recursively call path with newVisit -- this is our new 'starting point', and wont be reached again since we deleted it as an outgoing node

        self.walkpath.insert(0, visited) # insert into head of list walkpath the visited node --  this node will be the node with no exits

    def print(self):
        '''
        calls self.path on return of function self.startNode,
        then parses through self.walkpath (which contains the order in which we walked through the string elements)
        and adds to variable self.seq the first whole entry of the walkpath, followed by the last element of each entry,
        then prints this sequence

        for reference (sample dataset + output):
        ['GGCT', 'GCTT', 'CTTA', 'TTAC', 'TACC', 'ACCA']
        GGCT
        GGCT T
        GGCTT A
        GGCTTA C
        GGCTTAC C
        GGCTTACC A

        GGCTTACCA <--> final sequence


        '''

        self.path(self.startNode())

        for i in range(len(self.walkpath)):
            if i == 0: # starting element of walkpath
                self.seq += self.walkpath[i] # append starting element of walkpath to seq
            else:
                self.seq += self.walkpath[i][-1] # add last char of element in walkpath at location i to seq

        print(self.seq)


# File "C:/Users/wwsch/PycharmProjects/BME 205/Labs/Graphical Methods/py/problem14.py", line 92, in path
# while len(visitedEdge) > 0:     # while the length of visited node is greater than 0 (in this case, visited edge is a list containing elements)
# RecursionError: maximum recursion depth exceeded while calling a Python object
# to get around this, the minimum recursion limit i could use was 10^4
# https://stackoverflow.com/questions/3323001/what-is-the-maximum-recursion-depth-in-python-and-how-to-increase-it
import sys
sys.setrecursionlimit(10 ** 4)
def main():
    '''
    reads contents of input file, then instantiates pattern list 'p' into class reconstruct. Calls nodeFind which calls
    edgeFind to build the edge graph and then build the node graph. Then calls print(), which calls path() on the return
    of function startNode(), ultimately filling out the walkpath which is then combed through and added to self.seq
    in the same way problem 9 did -- the first whole entry of walkpath is added, followed by the last element of
    the proceeding entries in walkpath
    '''

    count = 0
    p = []
    k = 0 # unnecessary

    inputFile = open('rosalind_ba3h.txt', 'r') # open the input file for reading
    for line in inputFile:
        if count == 0:
            k = int(line.strip())
        else:
            p.append(line.strip())
        count += 1

    reconstruction = reconstruct(p)
    reconstruction.nodeFind()
    reconstruction.print()


if __name__ == '__main__':
    main()