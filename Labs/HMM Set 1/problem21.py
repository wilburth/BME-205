import math
class HMM:
    '''
    class HMM contains the necessary methods required to find a path that maximizes probability over all possible paths
    '''

    def __init__(self, x, transition, emission):
        '''
        initialize necessary vars
        '''
        self.x = x
        self.transition = transition
        self.emission = emission

    def decoder(self, x, transition, emission):
        '''
        method decoder calculates the path that maximizes probability over all possible paths
        it does this by using two different matrices, one for score and one to calculate the backtracking
        our backtracking matrix is filled with the transitions that occur between nodes, so we can use it
        starting from the end to find the node it came from

        a really helpful image that i drew out from my backtracking list
        shows how we start at end and then go through finding the nodes that connect to each-other
        with the path that will also have the best probability
        (note: this image is of the backtracking list, so to get the correct order we start from the end)

                https://imgur.com/a/ceGTvBJ
        '''

        # sample transition{'A': {'A': 0.641, 'B': 0.359}, 'B': {'A': 0.729, 'B': 0.271}}
        # sample emission {'A': {'x': 0.117, 'y': 0.691, 'z': 0.192}, 'B': {'x': 0.097, 'y': 0.42, 'z': 0.483}}

        final = ''
        n = len(x)
        states = list(transition.keys())
        score = [{state: 0 for state in states} for i in range(n)]
        backtrack = [{state: -1 for state in states} for i in range(n)]

        # initialize score 0
        for state in states:
            score[0][state] = math.log(emission[state][x[0]] / len(states))

        # starting from 1st state
        for i in range(1, n):
            for state in states:
                layer = []
                for previous in states:
                    # calculate weight (probability)
                    weight = math.log(transition[previous][state] * emission[state][x[i]])
                    t = weight + score[i-1][previous]
                    layer.append((t, previous))

                # max(layer)[0] obtains the higher (more positive) weight from layer list
                score[i][state] = max(layer)[0]
                # max(layer)[1] obtains the node associated with the higher (more positive) weight from layer list
                backtrack[i][state] = max(layer)[1]

        # get the last element, start from it
        current = max(score[n-1], key=lambda z: score[n-1][z])
        final += current

        # start at end, going to 0, step down 1
        for i in range(n-1, 0, -1):
            final += backtrack[i][current]
            current = backtrack[i][current]

        return final[::-1]


def main(inFile=None):
    '''
    reads from input file, creates transition + emission matrices used to calculate(find) a path that
    maximizes probability over all possible paths
    '''

    with open(inFile) as inFile:
        lines = inFile.readlines()
        x = lines[0].strip()
        alphabet = lines[2].strip().split()
        states = lines[4].strip().split()

        transition = {}
        emission = {}

        # again, using Mikeys way of stripping lines, super useful!
        for i in range(len(states)):
            name = lines[7 + i].strip()
            row = list(map(float, lines[7+i].split()[1:]))
            t = []
            for j in range(len(states)):
                t.append((states[j], row[j]))
            transition[name[0]] = {char: prob for char, prob in t}

        # again, using Mikeys way of stripping lines, super useful!
        for i in range(len(states)):
            name = lines[12 + i].strip()
            row = list(map(float, lines[12 + i].split()[1:]))
            t = []
            for j in range(len(alphabet)):
                t.append((alphabet[j], row[j]))
            emission[name[0]] = {char: prob for char, prob in t}

    hmm = HMM(x, transition, emission)
    print(hmm.decoder(x, transition, emission))

if __name__ == "__main__":
    main(inFile='rosalind_ba10c.txt')
    # output for sample file -- BBBBBACCCACCCACACACCCCCCACCCCACCCCACACCCACCCCACCCACACACCCABBBBBBBBACCCCACCCCACCCCCCACCCCACACCCACCCAC
