import copy
import re


class directedAcyclicGraph:
    '''
    class directedAcyclicGraph (DAG) contains the methods required to find the longest path between
    two nodes in an edge-weighted DAG. It also contains a nested class 'Node', which is used to keep
    track of node data outside of our node dictionary. I found that using this made it a lot easier
    to work with the data, and apply the various methods to it.

     ROSALIND PROMPT: Find a longest path between two nodes in an edge-weighted DAG.
        Given: An integer representing the source node of a graph, followed by an integer representing
        the sink node of the graph, followed by an edge-weighted graph. The graph is represented by a
        modified adjacency list in which the notation "0->1:7" indicates that
        an edge connects node 0 to node 1 with weight 7.
        Return: The length of a longest path in the graph, followed by a longest path.
        (If multiple longest paths exist, you may return any one.)


    '''

    class Node:
        '''
        https://www.datacamp.com/community/tutorials/inner-classes-python

        using a nested class allows us to initialize data within a node, and then reference it
        wherever we want within the program. this makes it a lot easier to keep track of node data,
        which keeps track of the incoming edges along with saving all of them to a list. additionally,
        we have path and score, which is used to score our nodes and then create a path from the scoring.

        '''

        def __init__(self, nodeData):
            '''
            initialization of data of Node
            keeps track of inEdges in integer format, along with the previous nodes in list format
            keeps track of path and score of node

            :param nodeData: a node
            '''
            self.nodeData = nodeData  # contains node data

            self.inEdges = 0  # plural edges, denotes how many nodes connect to this specified node

            self.previousN = []  # a list to save all of the incoming nodes to this specified node

            self.path = None
            self.score = None



    def __init__(self, source, sink, data):
        '''
        initialization of necessary vars

        :param source: start node
        :param sink: end node
        :param data: edgeList, cleaned from main()
        '''

        self.graph = {}
        self.nodes = {}

        # for our node dictionary, we want to make sure we include the source (start) and sink (end) nodes
        # in our initialization -- this way we dont have to do it later
        self.source = self.nodes[source] = self.Node(source)
        self.sink = self.nodes[sink] = self.Node(sink)

        self.graphCreate(data) # calls graphCreate on data, passing it as an edgeList given by problem



    def graphCreate(self, edgeList):
        '''
        graphCreate creates the adjacency matrix

        https://www.tutorialspoint.com/python/dictionary_setdefault.htm

        :param edgeList: input edgeList which has been cleaned in main() and sent to initialization and passed to graphCreate()
        '''

        # for each item in our edgeList, make sure the key value and weight are all integers
        # example edgeList [('6', '26', '32'), ('10', '39', '30')] --> [(key, value, weight), (key, value, weight)]
        # [(inNode, outNode, edgeWeight), (inNode, outNode, edgeWeight)]
        # building of graph takes place in this for loop
        for item in edgeList:
            key = int(item[0])
            value = int(item[1])
            weight = int(item[2])

            # if the key is in nodes dictionary, temp var set to address of key
            # else (its not in nodes dict) temp var set to address of that key and initialize that key in Node class
            if key in self.nodes:
                curr = self.nodes[key]
            else:
                curr = self.nodes[key] = self.Node(key)

            # if the value is in nodes dictionary, temp var set to address of value
            # else (its not in nodes dict) temp var set to address of that value and initialize that value in Node class
            if value in self.nodes:
                next = self.nodes[value]
            else:
                next = self.nodes[value] = self.Node(value)

            # "This method returns the key value available in the dictionary and if given key is
            # not available then it will return provided default value." -- setdefault(data, []).append(value, weight)
            self.graph.setdefault(self.nodes[key].nodeData, []).append((value, weight))

            # before setdefault {}
            # after setdefault {6: [(26, 32)]}
            # value: 26
            # before setdefault {6: [(26, 32)], 26: []}
            # after setdefault {6: [(26, 32)], 26: [], 10: [(39, 30)]}
            # value: 39
            # before setdefault {6: [(26, 32)], 26: [], 10: [(39, 30)], 39: []}
            # after setdefault {6: [(26, 32)], 26: [(28, 24)], 10: [(39, 30)], 39: []}
            # value: 28

            # if the value hasn't been seen in graph dictionary
            # update the graph with that value as a key, and a list as its value
            if value not in self.graph:
                self.graph.update({value: []})

            # add to inEdges of value (outgoing node has += 1 inedge)
            self.nodes[value].inEdges += 1
            # append to previousNode list of value the key
            # [6]
            # [10]
            # [26]
            # [3] . . .
            self.nodes[value].previousN.append(key)

        self.source.inEdges = 0 # declare source inEdges to be 0 (its our start node, only has outEdges)

    def topologicalOrder(self):
        '''
        topologicalOrder is based off the pseudoCode provided in class / in the textbook, and appends to
        list iterable to topological ordering of the input graph.


        :return: returns an iterable list that is in topological order
        '''

        # create a copy of graph dictionary
        # create a copy of nodes dictionary
        deepGraph = copy.deepcopy(self.graph)
        deepNodes = copy.deepcopy(self.nodes)

        candidates = [] # list to house our candidate nodes
        iterable = [] # list to house nodes to iterate over

        # for each node in sorted node dictionary
        # if the number of inEdges is equal to 0
        # append to candidate list the data of each node
        for node in sorted(self.nodes):
            if self.nodes[node].inEdges == 0:
                candidates.append(self.nodes[node].nodeData)

        # while length of candidate list is greater than 0 (there are still elements in it)
        # pop each element
        # add the popped element to iterable list
        while len(candidates) > 0:
            x = candidates.pop(0)
            iterable.append(x)

            # while the length of deepGraph at location [a] (popped candidate)
            # pop element of deepgraph at pos x
            # count down the inEdges
            # if the inEdges reach 0, add it to our candidate list
            while len(deepGraph[x]) > 0:
                y, yWeight = deepGraph[x].pop(0)
                deepNodes[y].inEdges -= 1

                if deepNodes[y].inEdges == 0:
                    candidates.append(y)

        # the list iterable is filled with all nodes in topological order
        return iterable



    def longestPath(self, iterable):
        '''
        longestPath finds the longest path from the topological ordering with list iterable.
        based off the algorithm described in class and in textbook

        :param iterable: list of nodes that we are to iterate through, they have been sorted in order of inEdges count
        :return: the score of the longest path, along with the longest path itself
        '''

        # for each node in our graph dictionary
        # set the score in node dictionary of that node to be arbitrary negative value
        # set source score to be 0
        for node in self.graph:
            self.nodes[node].score = -1000000000000000
        self.source.score = 0

        # for item in iterable list
        # if the inEdges is 0
        # add to path that value
        # continue along iterable
        for value in iterable:
            if self.nodes[value].inEdges == 0:
                self.nodes[value].path = [value]

            # else (inEdges > 0)
            # var max previous set to previous
            else:
                maxP = self.nodes[value].previousN[0]

                # for node in all previous
                for node in self.nodes[value].previousN:
                    # if score + calculated edgeWeight of node is greater than
                    # score + calculated edgeWeight of previous node
                    # set previous node to current
                    if self.nodes[node].score + self.edgeWeight(node, value) > \
                       self.nodes[maxP].score + self.edgeWeight(maxP, value):
                        maxP = node

                # then, get the edgeWeight of current
                valueW = self.edgeWeight(maxP, value)
                # the score is equal to the edgeWeight + the original score
                self.nodes[value].score = valueW + self.nodes[maxP].score
                # the path has the value from list iterable added to it
                self.nodes[value].path = self.nodes[maxP].path + [value] # needs to be in list[] format to be added to path

        # return the score of the sink node, along with the path it took to get there
        return self.sink.score, self.sink.path

    def edgeWeight(self, x, y):
        '''
        quick method to calculate edgeWeight between two nodes

        :param x: incoming node
        :param y: outgoing node
        :return: edgeWeight between specified nodes
        '''

        # for all edges in self.graph at position x  (key)
        # if the outgoing edge equals the specified call (y)
        # return the edgeWeight of x->y
        # ex. x=1 --> retrieve edge[0] (outgoing node) of node1 --> 3 ;;
        # y also called with 3 ;;
        # return 36 (the edgeWeight of x1->y3)
        for edge in self.graph[x]:
            if edge[0] == y:
                return edge[1]


def main(inFile=None, options=None):
    '''
    main() opens the specified inFile, cleans it and sends it to be initialized in class directedAcyclicGraph.
    It then calls topologicalOrder to create the graph (adjacency matrix), and finally calls longestPath with
    the iterable adj matrix to calculate the score and path of the sink node, which should give us the path that we want
    to take.


    :param inFile: inFile provided by rosalind in .txt format and is specified in call to main (bottom of program)
    :param options: NONE
    :return: the score and walkpath of the longest path in the DAG

    https://stackoverflow.com/questions/4998629/split-string-with-multiple-delimiters-in-python
    '''

    edgeList = []
    source = 0
    sink = 0
    count = 0

    f = open(inFile)
    for line in f:
        if count == 0:
            source = int(line[0:])
        if count == 1:
            sink = int(line[0:])
        if count >= 2:
            lineList = (re.split(':|->', line.strip()))
            edgeList.append((lineList[0], lineList[1], lineList[2])) # [(inEdge, outEdge, weight)] --> [('0', '1', '7'), ('0', '2', '4'), ('2', '3', '2'), ('1', '4', '1'), ('3', '4', '3')]

        count +=1


    graph = directedAcyclicGraph(source, sink, edgeList)
    iterable = graph.topologicalOrder() # returns iterable
    score, path = graph.longestPath(iterable)

    print("{}".format(score))
    walkpath = ''
    for i in range(len(path) - 1):
        walkpath += str(path[i]) + '->'
    walkpath += str(path[-1])
    print(walkpath)


if __name__ == "__main__":
    main(inFile='rosalind_ba5d.txt', options=[None]) # no options for DAG