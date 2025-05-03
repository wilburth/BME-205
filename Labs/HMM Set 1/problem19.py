class HMM:
    '''
    class HMM contains the methods required to calculate the probability of a hidden path in an HMM
    '''

    def __init__(self, transitionMatrix):
        '''
        initialize necessary vars
        '''

        self.transMatrix = transitionMatrix

    def hiddenPathProb(self, path, transition):
        '''
        calculates the probability of the hidden path in an HMM
        '''

        probability = 1/len(transition)
        for i in range(len(path) - 1):
            x = path[i]
            y = path[i + 1]
            probability *= transition[x][y]

        return probability


def main(inFile=None):
    '''
    reads from input file, initializes transition matrix into class HMM + calls method to calculate the probability
    of the hidden path from the HMM transition matrix
    '''

    transition = {}
    with open(inFile) as inFile:
        lines = inFile.readlines()
        path = lines[0].strip()
        states = lines[2].strip().split()
        columns = lines[4].strip().split()

        # used Mikeys way of reading in -- was having a lot of trouble with
        # creating the actual transition matrix (for state in states couldn't manipulate properly)
        # but using the for i in range of the length of states and then the 5+i to get lines 7
        # was really cool and easy way to create the rows that I wanted
        for i in range(len(states)):
            name = lines[5+i].strip()
            row = list(map(float, lines[5+i].split()[1:]))
            t = []

            for j in range(len(columns)):
                t.append([columns[j], row[j]])
            transition[name[0]] = {name: prob for name, prob in t}

    hmm = HMM(transition)
    print(hmm.hiddenPathProb(path, transition))


if __name__ == "__main__":
    main(inFile='sample.txt')
    # output for sample file == 2.0777937919245553e-23
