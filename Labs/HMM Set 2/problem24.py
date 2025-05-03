class HMM:
    '''
    class HMM contains the methods required to calculate the probability that the HMM was in state k
    at step i (for each state k and each step i)
    '''

    def __init__(self, x, alphabet, states, transition, emission):
        '''
        initialize necessary vars
        '''

        self.x = x
        self.alphabet = alphabet
        self.states = states
        self.transition = transition
        self.emission = emission

    def backwards(self, x, transition, emission, step, st):
        '''
        backwards works almost like forwards except we have to change ranges to fit the HMM
        utilizes steps and state to record probability of state at each step
        '''

        n = len(x)
        states = transition.keys()
        if step >= n-1:
            return 1
        # step + 1 -- need 9 not 10
        score = [{state: 0 for state in states} for i in range(step + 1, n)]
        # [{'A': 0, 'B': 0}, {'A': 0, 'B': 0}, {'A': 0, 'B': 0}, {'A': 0, 'B': 0}, {'A': 0, 'B': 0}, {'A': 0, 'B': 0}, {'A': 0, 'B': 0}, {'A': 0, 'B': 0}, {'A': 0, 'B': 0}]
        # [{'A': 0, 'B': 0}, {'A': 0, 'B': 0}, {'A': 0, 'B': 0}, {'A': 0, 'B': 0}, {'A': 0, 'B': 0}, {'A': 0, 'B': 0}, {'A': 0, 'B': 0}, {'A': 0, 'B': 0}] etc....

        for state in states:
            score[0][state] = transition[st][state] * emission[state][x[step+1]]

        # step + 2 -- need 8
        for i in range(step + 2, n):
            index = i - step - 1  # 1 2 3 4 5 6 7 8
            # for each state
            for state in states:
                layer = []
                # need previous
                for previous in states:
                    weight = transition[previous][state] * emission[state][x[i]]
                    t = weight * score[index-1][previous]
                    layer.append(t)
                score[index][state] = sum(layer)
        # at each step return sum of sink node
        return sum(score[-1].values())

    def forwards(self, x, transition, emission, step, st):
        '''
        forwards works like it does in problem 22, except we implement the use of steps so that we can record
        the probability of state at each step, along with the given state
        '''

        n = len(x)
        states = transition.keys()
        if step == 0:
            return emission[st][x[0]] / len(states)

        score = [{state: 0 for state in states} for i in range(step)]
        # [{'A': 0, 'B': 0}]
        # [{'A': 0, 'B': 0}, {'A': 0, 'B': 0}] etc . . .
        for state in states:
            score[0][state] = emission[state][x[0]] / len(states)
        # starting from pos1
        for i in range(1, step):
            # for each state
            for state in states:
                layer = []
                # need previous
                for previous in states:
                    weight = transition[previous][state] * emission[state][x[i]]
                    t = weight * score[i-1][previous]
                    layer.append(t)
                score[i][state] = sum(layer)
        # if step is less than len string
        if step < n:
            t = 0
            for previous in states:
                # calculate probability
                weight = transition[previous][st] * emission[st][x[step]]
                t += weight * score[step-1][previous]
            return t
        # otherwise
        else:
            # return sink node probability
            return sum(score[-1].values())

    def softDecoder(self):
        '''
        our decoder will work like it does in problem 21, except that we calculate based on
        forward*backward/forward(sink)
        '''
        n = len(self.x)
        states = list(self.transition.keys())
        p = {state: [0] * n for state in states}
        # utilize position across string x for steps
        # utilize state across states
        # for each position in range of path
        for pos in range(n):
            # for each state call forwards + backwards + forwards(sink)
            for state in states:
                forward = self.forwards(self.x, self.transition, self.emission, step=pos, st=state)
                s = self.forwards(self.x, self.transition, self.emission, step=n, st=state)
                backward = self.backwards(self.x, self.transition, self.emission, step=pos, st=state)
                # pr[state][pos] = forward*backward / forward(sink)
                p[state][pos] = round(forward*backward / s, 4)

        return p

def main(inFile=None):
    '''
    def main reads from inFile and parses input, creating the emission and transition matrices, along with passing them
    into class HMM and calling softDecoder to find the probability that the HMM was in state k at step i (for each state
    k and each step i)
    '''

    with open(inFile) as inFile:
        lines = inFile.readlines()
        x = lines[0].strip()
        alphabet = lines[2].strip().split()
        states = lines[4].strip().split()

        transName = lines[6].strip().split()
        # ['A', 'B']
        emitName = lines[10].strip().split()
        # ['x', 'y', 'z']

        transition = {}
        emission = {}

        for i in range(len(states)):
            name = lines[7+i].strip()
            row = list(map(float, lines[7+i].split()[1:]))
            # (0, '0.911') (1, '0.089')
            # (0, '0.228') (1, '0.772')
            t = []
            for j in range(len(states)):
                t.append((states[j], row[j]))
            transition[name[0]] = {char: prob for char, prob in t}

        for i in range(len(states)):
            name = lines[11 + i].strip()
            row = list(map(float, lines[11+i].split()[1:]))
            t = []
            for j in range(len(alphabet)):
                t.append((alphabet[j], row[j]))
            emission[name[0]] = {char: prob for char, prob in t}


    #print(transition)
    # {'A': {'A': 0.911, 'B': 0.089}, 'B': {'A': 0.228, 'B': 0.772}}
    #print(emission)
    # {'A': {'x': 0.356, 'y': 0.191, 'z': 0.453}, 'B': {'x': 0.04, 'y': 0.467, 'z': 0.493}}

    hmm = HMM(x, alphabet, states, transition, emission)
    probs = hmm.softDecoder()
    # print(probs)

    print(*list(probs.keys()))
    for pos in range(len(x)):
        print(*[probs[state][pos] for state in probs.keys()])

if __name__ == "__main__":
    main(inFile='rosalind_ba10j.txt')
    # output for sample file --
    # A B
    # 0.5168 0.4832
    # 0.2655 0.7345
    # 0.4343 0.5657
    # 0.4279 0.5721
    # 0.2672 0.7328
    # 0.4342 0.5658
    # 0.4302 0.5698
    # 0.4303 0.5697
    # 0.4303 0.5697
    # 0.4307 0.5693
