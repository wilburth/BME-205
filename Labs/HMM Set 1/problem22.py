class HMM:
    '''
    class HMM contains the necessary methods required to find the probability that the HMM emits x, where x is
    a string
    '''

    def __init__(self, x, transition, emission):
        '''
        initialize necessary vars
        '''
        self.x = x
        self.transition = transition
        self.emission = emission

    def forward(self, x, transition, emission):
        '''
        method forward calculates the probability of an HMM emitting a string x
        does this through the use of matrix "score", traversing the input string
        and getting the values associated with the emission from the string node,
        and calculating the probability that we went from one node to the next.
        Summing the sink nodes scores will provide us with the probability that HMM emitted string x

        example string: xzyyzzyzyy
        we calculate emission of x for start
        then begin at z->end

        another really useful photo that shows what calculations are taking place
        (note: this happens per letter in string -- in other words, for each character in the string
        we will get the weight from transition*emission, the multiply it from source node, then sum both the values
        (from state a and state b) and add it as the score at the specified state)

                https://imgur.com/a/D469sQO

        '''
        n = len(x)
        states = list(transition.keys())
        score = [{state: 0 for state in states} for i in range(n)]
        # initialize score for state 0
        for state in states:
            score[0][state] = emission[state][x[0]] / len(states)

        # starting from 1st state
        for i in range(1, n):
            for state in states:
                layer = []
                for previous in states:
                    # calculate weight
                    weight = transition[previous][state] * emission[state][x[i]]
                    t = weight * score[i-1][previous]
                    layer.append(t)
                score[i][state] = sum(layer)

        # create sink node which is sum of the sink node, equaling the probability that the HMM emits path X
        sink = sum(score[n-1].values())

        return sink


def main(inFile=None):
    '''
    reads from input file, creates transition + emission matrices, then calculates probability of the HMM
    emitting string x
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
            name = lines[11 + i].strip()
            row = list(map(float, lines[11 + i].split()[1:]))
            t = []
            for j in range(len(alphabet)):
                t.append((alphabet[j], row[j]))
            emission[name[0]] = {char: prob for char, prob in t}

    hmm = HMM(x, transition, emission)
    print(hmm.forward(x, transition, emission))

if __name__ == "__main__":
    main(inFile='rosalind_ba10d.txt')
    # output for sample file -- 1.969174404840454e-49
