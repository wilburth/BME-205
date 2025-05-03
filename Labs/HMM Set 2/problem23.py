class HMM:
    '''
    class HMM contains the methods required to calculate a matrix of transition probabilities
    and a matrix of emission probabilities that maximize the path over all possible matrices
    '''

    def __init__(self, x, alphabet, path, states, transition, emission):
        '''
        initialize necessary vars
        '''
        self.x = x
        self.alphabet = alphabet
        self.path = path
        self.states = states
        self.transition = transition
        self.emission = emission

    def param(self):
        '''
        def param uses the passed + initialized variables to calculate both the transition
        and emission matrices that house the probabilities of maximizing a path over all possible
        matrices

        '''
        # initialize emission state 0
        self.emission[self.path[0]][self.x[0]] += 1

        # search range of path starting from pos1
        for i in range(1, len(self.x)):
            previous = self.path[i-1]
            current = self.path[i]

            self.transition[previous][current] += 1
            self.emission[current][self.x[i]] += 1

        # print(self.transition)
        # {'A': {'A': 0, 'B': 4, 'C': 0}, 'B': {'A': 4, 'B': 1, 'C': 0}, 'C': {'A': 0, 'B': 0, 'C': 0}}
        # print(self.emission)
        # {'A': {'x': 1, 'y': 1, 'z': 2}, 'B': {'x': 3, 'y': 1, 'z': 2}, 'C': {'x': 0, 'y': 0, 'z': 0}}

        for state in self.states:
            tRow = self.transition[state]
            # {'A': 0, 'B': 4, 'C': 0}
            eRow = self.emission[state]
            # {'x': 1, 'y': 1, 'z': 2}

            # added due to divide by zero error
            # if there aren't any values in the row we +1 to each value at each key (exit case)
            # in the given sample, there were 3 states (a,b,c) but the path only contained a,b so we normalize
            if not any(tRow.values()):
                for key in tRow:
                    tRow[key] += 1
            if not any(eRow.values()):
                for key in eRow:
                    eRow[key] += 1

            # sum row values for transition + emission table
            tSum = sum(tRow.values())
            eSum = sum(eRow.values())

            # calculate probabilities w/ rounding to nearest 3 decimals
            # https://www.w3schools.com/python/ref_func_round.asp
            for key in tRow:
                tRow[key] = round(tRow[key] / tSum, 3)
            for key in eRow:
                eRow[key] = round(eRow[key] / eSum, 3)

        return self.transition, self.emission

def main(inFile=None):
    '''
    def main() reads from input file, parses variables, initializes them into class HMM
    and then calls method param(), which returns the transition and emission matrices
    these two matrices house probabilities that maximize the given path over all possible matrices

    '''

    with open(inFile) as inFile:
        lines = inFile.readlines()
        x = lines[0].strip()
        alphabet = lines[2].strip().split()
        path = lines[4].strip()
        states = lines[6].strip().split()

    alphabet = list(alphabet)
    # ['x', 'y', 'z']
    path = list(path)
    # ['B', 'B', 'A', 'B', 'A', 'B', 'A', 'B', 'A', 'B']

    transition = {state: {state: 0 for state in states} for state in states}
    # {'A': {'A': 0, 'B': 0, 'C': 0}, 'B': {'A': 0, 'B': 0, 'C': 0}, 'C': {'A': 0, 'B': 0, 'C': 0}}
    emission = {state: {char: 0 for char in alphabet} for state in states}
    # {'A': {'x': 0, 'y': 0, 'z': 0}, 'B': {'x': 0, 'y': 0, 'z': 0}, 'C': {'x': 0, 'y': 0, 'z': 0}}

    hmm = HMM(x, alphabet, path, states, transition, emission)
    hmm.param()

    # print transition matrix
    print('', *transition.keys())
    for key in transition:
        print(key, *transition[key].values(), sep='\t')
    print('-' * 8)
    # print emission matrix
    print('', *alphabet)
    for key in emission:
        print(key, *emission[key].values(), sep='\t')


if __name__ == "__main__":
    main(inFile='rosalind_ba10h.txt')
    # output for sample file --
    #  A B C D
    # A	0.406	0.219	0.094	0.281
    # B	0.318	0.182	0.273	0.227
    # C	0.292	0.208	0.333	0.167
    # D	0.238	0.286	0.333	0.143
    # --------
    #  x y z
    # A	0.312	0.375	0.312
    # B	0.364	0.455	0.182
    # C	0.625	0.292	0.083
    # D	0.409	0.364	0.227
