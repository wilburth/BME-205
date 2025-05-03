class HMM:
    '''
    class HMM contains the methods required to calculate the conditional probability that
    a given string x will be emitted from the HMM, given the hidden path
    '''

    def __init__(self, emissionMatrix):
        '''
        initialize necessary vars
        '''

        self.emissionMatrix = emissionMatrix

    def outcomeProbability(self, path, x, emission):
        '''
        calculates the probability of the HMM emitting a given string x
        '''

        probability = 1
        for i in range(len(path)):
            probability *= emission[path[i]][x[i]]

        return probability


def main(inFile=None):
    '''
    reads from input file, initializes emission matrix into class HMM, then calls method to calculate
    probability of the HMM emitting given string x

    '''

    with open(inFile) as inFile:
        lines = inFile.readlines()
        x = lines[0].strip()
        alphabet = lines[2].strip().split()
        path = lines[4].strip()
        states = lines[6].strip().split()
        emit = {}

        # used Mikeys way of reading in -- was having a lot of trouble with
        # creating the actual emission matrix (for state in states couldn't manipulate properly)
        # but using the for i in range of the length of states and then the 9+i to get lines 11 and on
        # was really cool and easy way to create the rows that I wanted
        for i in range(len(states)):
            name = lines[9+i].strip()
            row = list(map(float, lines[9+i].split()[1:]))

            t = []
            for j in range(len(alphabet)):
               t.append((alphabet[j], row[j]))
            emit[name[0]] = {char: prob for char, prob in t}

    hmm = HMM(emit)
    print(hmm.outcomeProbability(path, x, emit))


if __name__ == "__main__":
    main(inFile='rosalind_ba10b.txt')
    # output for sample file -- 3.668370394918308e-30
