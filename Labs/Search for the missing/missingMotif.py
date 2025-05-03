# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough (wschloug)
# Group Members: Junmo, Mikey, Wenyu


# Write a proper style BME205 notebook that reads a fasta file from STDIN and ranks motifs based on how statistically underrepresented the specific motif is.
# We need only consider sequence motifs that are up to about 8 in length ( this should be a program option).
# In order to use the equation that we developed in class, the minimum size should be at least 3 ( this should be an option also).
# We will find that there are quite a few of these, so we need to specify a statistical cutoff.
# For this assignment, we will use Z-scores.
# These will be negative since we are looking for underrepresented sequences (this cutoff needs to be a program option too. Finally, we need to sort the output by z-score and print it out, maybe to a file, so let's use STDOUT for this.
# Remember that sequences that we get from DNA sequencing are from only one strand of DNA. This means that when we see a AAGGTT, this implies that AACCTT ( the reverse complement) is present on the other strand.
# We should count these sequences as equivalent.
# Your report should be sorted by Zscore within motif-Size, such that the longest motifs are at the beginning of the report and, within motifs of a given size, print the most significant motifs first.
# Print the k-mer pairs in alpha order. ( AAA:TTT and not TTT:AAA)

########################################################################
# FastAreader
# provided by Dr. B
########################################################################
import sys


class FastAreader:

    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''

        self.fname = fname

    def doOpen(self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):

        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence


########################################################################
# CommandLine
# provided by Dr. B
########################################################################

class CommandLine():
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse,
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination exception Usage.

    attributes:
    all arguments received from the commandline using .add_argument will be
    available within the .args attribute of an object instantiated from CommandLine.
    For example, if thisCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then thisCommandLine.args.requiredbool will
    name that option.

    '''

    def __init__(self, inOpts=None):
        '''
        CommandLine constructor.

        Implement a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(
            description='Program prolog - a brief description of what this thing does',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='python missingMotif.py -minMotif int -maxMotif int -cutoff int < input.fa > output.out'
        )

        # options as specified by assignment
        self.parser.add_argument('-minMotif', type=int, choices=range(3, 9), action='store', default=3,
                                 help='minimum motif size in integer format')
        self.parser.add_argument('-maxMotif', type=int, choices=range(3, 9), action='store', default=8,
                                 help='maximum motif size in integer format')
        self.parser.add_argument('-cutoff', type=int, action='store', default=0,
                                 help='z score cutoff in integer format')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


class Usage(Exception):
    '''
    Signal a Usage error, evoking a usage statement and eventual exit when raised.
    '''

    def __init__(self, msg):
        self.msg = msg


########################################################################
# Genome
#
########################################################################
class Genome:
    '''

    class genome contains methods that are used to find the z-score, expected value, and count of all motifs

    setup --------------
    initialization
    z-score
    expected value
    motifcount / kmer list

    way to reverse easily -- bme 160 method

    '''



    def __init__(self, seq, N, min, max, cutoff):
        '''

        seq -- input sequences taken from fastaReader
        N -- genome size (length of all sequences)
        min --
        max --
        cutoff --

        '''

        self.seq = seq
        self.min = min
        self.max = max
        self.cutoff = cutoff

        self.kmerCounts = dict() # keep track of kmer + times found in seq (along with counting for reverse compliment too)

        self.N = N  # genome size

        for sequence in self.seq: # for each sequence
            for k in range(1, self.max+1): # for each motif
                for i in range(0, len(sequence) - k + 1): # find all motifs of that size
                    kmer = ''.join(sequence[i:i + k]) # kmer cut from i to k

                    if self.reverseSeq(kmer) in self.kmerCounts: # check if the reverse compliment is in dictionary
                        self.kmerCounts[kmer] = self.kmerCounts[self.reverseSeq(kmer)] # if it is, set it to same val as self.kmerCounts[kmer] so that both go up in count
                        self.kmerCounts[kmer][0] += 1 # add one to the count of specified kmer (and also, the count of its reverseSeq)

                    elif kmer in self.kmerCounts: # else if its already in dictionary
                        self.kmerCounts[kmer][0] += 1 # add 1 to the count of specified kmer

                    else: # otherwise
                        self.kmerCounts[kmer] = [1] # add it as a key to dictionary

    def reverseSeq(self, sequence):
        '''

        used to find reverse compliment -- method taken from bme 160, very useful
        returns the reverse compliment

        https://www.w3schools.com/python/ref_string_translate.asp
        https://www.programiz.com/python-programming/methods/string/maketrans

        '''

        #(sequence.translate(str.maketrans("ATCG", "TAGC"))[::-1])
        return sequence.translate(str.maketrans("ATCG", "TAGC"))[::-1]

    def eval(self, kmer):
        '''

        using markovian model -- left + right size + middle kmer are used to calculate escore
        calculates and returns the expected value of a specific motif appearing

        '''
        tail = len(kmer) - 1 # like middle
        left = kmer[:-1] # left kmer
        right = kmer[1:] # right kmer

        # consider counts of left side + right side of motif
        numerator = self.kmerCounts[left][0] * self.kmerCounts[right][0]

        # consider counts of middle part of motif
        denominator = self.kmerCounts[kmer[1:tail]][0]

        # divide the counts of left + right side of motif by the counts of middle part of motif
        return numerator / denominator

    def zscore(self, kmer):
        '''


        For this assignment, we will use Z-scores.
        These will be negative since we are looking for underrepresented sequences (this cutoff needs to be a program option too)

        "A score that indicates how many standard deviations a value is above or below the mean."

        '''
        mean = self.eval(kmer)
        N = self.N # "n = genome size (approximately)"
        p = mean / N # Pr(K) = E(K)/N -- N being length of genome, and we use eval() method for E(K)

        standardDev = (N * p * (1 - p)) ** 0.5  # sqrt(np(1-p)

        if standardDev != 0: # if our stddev is 0 we can automatically count it out
            zscore = (self.kmerCounts[kmer][0] - mean) / standardDev # calculate z score
            return zscore
        else:
            return 1



    def kmerList(self):
        '''


        '''
        kmax = 8  # due to program requirements -- "We need only consider sequence motifs that are up to about 8 in length"
        kmerList = list()
        listcheck = list()

        # check all keys
        for kmer in self.kmerCounts.keys():
            # kmer between min and max
            # reverse not yet added (using listcheck to see if kmer has been added yet or not)
            # zscore below (but not including) cutoff
            if self.min <= len(kmer) <= self.max and \
                    self.reverseSeq(kmer) not in listcheck and \
                    self.zscore(kmer) < self.cutoff:

                        k = kmax - len(kmer) # max motif size
                        # sequence: reverse     count      Expect Zscore
                        kmerList.append([kmer,
                                         self.reverseSeq(kmer),
                                         self.kmerCounts[kmer][0],
                                         self.eval(kmer),
                                         self.zscore(kmer),
                                         k]) # motif size
                        listcheck.append(kmer)

        # per program requirements -- "sorted by Zscore within motif-Size"
        # such that the longest motifs are at the beginning of the report and,
        # within motifs of a given size, print the most significant motifs first.
        # Print the k-mer pairs in alpha order.
        # ( AAA:TTT and not TTT:AAA)
        # kmerList sorted by (motif-size, then by z-score)
        return sorted(kmerList, key=lambda kmer:(kmer[5], kmer[4]))


def main(myCommandLine=None):
    '''
    Implement the Usage exception handler that can be raised from anywhere in process.
    '''
    if myCommandLine is None:
        myCommandLine = CommandLine()  # read options from the command line
    else:
        myCommandLine = CommandLine(myCommandLine)  # interpret the list passed from the caller of main as the commandline.

    seqList = []
    N = 0 # will keep track of genome size by counting length of each sequence passed by FastAreader,

    myReader = FastAreader()
    for head, seq in myReader.readFasta():
        seqList.append(seq)
        N += len(seq)

    myMotifs = Genome(seqList,
                      N, # then passed as a var in myMotifs = Genome() initialization
                      myCommandLine.args.minMotif,
                      myCommandLine.args.maxMotif,
                      myCommandLine.args.cutoff)

    kmerList = myMotifs.kmerList()
    print("N: {}".format(N))
    print("sequence:reverse \t\t count \t\t Expect \t\t Zscore")
    for kmer in kmerList:  # prints each kmer with respective data (in formatted order -- done through kmerList())
        print('{0:8}:{1:8} \t\t {2:0d} \t\t {3:0.2f} \t\t {4:0.2f}'.format(kmer[0], # kmer
                                                                           kmer[1][::-1], # reverse compliment kmer [::-1] to reverse the sequence to fit the output parameters
                                                                           kmer[2], # count
                                                                           kmer[3], # eValue
                                                                           kmer[4])) # z-score


if __name__ == "__main__":
    main()
