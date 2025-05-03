# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough (wschloug)
# Group Members: Junmo, Mikey, Wenyu

# Find an Eulerian path in a graph.
# Given: A directed graph that contains an Eulerian path, where the graph is given in the form of an adjacency list.
# Return: An Eulerian path in this graph.

import re

class eulerian:
    '''
    contains methods used to find an Eulerian path in a graph (of nodes)
    '''

    def __init__(self, edges):
        '''
        initializes necessary vars
        '''

        self.edges = edges # takes in list of edges
        self.nodes = {} # this will be our node dictionary --> housing {key (node): [input node], [output node]}
        self.walkpath = [] # keeps track of visited

    def node(self):
        '''
        uses filled edgeList (self.edges) and fills self.nodes, which is a dictionary that houses a node as its key,
        along with a list as its value. The list contains the input (incoming) node, and the output (outgoing) node
        for the key node.

        '''

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
            if len(self.nodes[node][1]) - len(self.nodes[node][0]) == 1: # if (# of output) - (# of input) == 1, ex 6->3,7 --> 2-1 = 1
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


        #print('visited: ' + str(visited))

        visitedEdge = self.nodes[visited][1] # retrieves outgoing node from startNode

        #print('visitedEdge: ' + str(visitedEdge))
        #print('length of visited edge: ' + str(len(visitedEdge)))




        while len(visitedEdge) > 0:     # while the length of visited node is greater than 0 (in this case, visited edge is a list containing elements)
            newVisit = visitedEdge[0]   # retrieves the first element in list of visitedEdge (ex. 6 -> 3,7  newVisit set as 3, not 7)
            #print('newVisit: ' + str(newVisit))
            self.nodes[visited][1].remove(newVisit)     # remove the newVisit element from outgoing nodes from original call (ex 6->3,7 has 3 removed, BUT 7 is still part of the list)
            self.path(newVisit)     # then, recursively call path with newVisit -- this is our new 'starting point', and wont be reached again since we deleted it as an outgoing node

        self.walkpath.insert(0, visited) # insert into head of list walkpath the visited node --  this node will be the node with no exits





    def print(self):
        '''
        calls self.path with the return of function self.startNode,
        then prints the walkpath separated by '->' per problem specifications.
        Learned some new python along the way!
        "Without using loops: * symbol is use to print the list elements in a single line with space.
        To print all elements in new lines or separated by space use sep=”\n” or sep=”, ” respectively."

        https://www.geeksforgeeks.org/print-lists-in-python-4-different-ways/
        https://treyhunner.com/2018/10/asterisks-in-python-what-they-are-and-how-to-use-them/
        '''

        self.path(self.startNode()) # call path with startNode

        #for x in range(len(self.walkpath)):
            #print(self.walkpath[x])

        #print(*self.walkpath)

        print(*self.walkpath, sep='->')



import sys
sys.setrecursionlimit(10 ** 4)
def main():
    '''

    reads contents of input file and uses regualar expression '.split()' to cut the '->' from the input
    and appends to edgeList the (incomingNode, outgoingNode).
    Then instantiates edgeList into class eulerian, calls node() which fills dictionary nodes(), finally calls print()
    which in turn calls path() with the return of function startNode() to fill the walkpath with our visiting order
    and prints it in the specified output -- a single string, with startNode->node->node->node->. . . .->endNode

    https://docs.python.org/3/library/re.html

    "re.split(pattern, string, maxsplit=0, flags=0)
    Split string by the occurrences of pattern. If capturing parentheses are used in pattern,
    then the text of all groups in the pattern are also returned as part of the resulting list.
    If maxsplit is nonzero, at most maxsplit splits occur,
    and the remainder of the string is returned as the final element of the list."

    https://stackoverflow.com/questions/12683201/python-re-split-to-split-by-spaces-commas-and-periods-but-not-in-cases-like
    '''

    edgeList = []
    inputFile = open('rosalind_ba3g.txt', 'r') # open the input file for reading
    for line in inputFile:
        lineList = re.split(', | -> |,', line.strip()) # (ex. 0 -> 2 gets split at '->', becomes ['0', '2'] with 0 being incoming node and 2 being outgoing node)
        for i in lineList[1:]: # need lineList[1: to the end] -- this way we get ALL outgoing edges in the case that an incoming node has multiple outgoing nodes (ex. 6-> 3, 7)
            edgeList.append((int(lineList[0]), int(i)))

    eularianPath = eulerian(edgeList)
    eularianPath.node()
    eularianPath.print()



if __name__ == '__main__':
    main()