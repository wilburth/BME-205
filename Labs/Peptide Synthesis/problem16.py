# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough (wschloug)
# Group Members: Junmo, Mikey, Wenyu


class peptide:
    '''
    class peptide contains all methods required to solve problem16 on rosalind.
    The prompt is as follows:

    Find substrings of a genome encoding a given amino acid sequence.
    Given: A DNA string Text and an amino acid string Peptide.
    Return: All substrings of Text encoding Peptide (if any such substrings exist).
    Sample:
    ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA
    MA

    '''

    '''
    aa_to_dna is a re-organized dna codon table -- it is still a dictionary,
    but it groups amino acid codons into lists.
    '''
    aa_to_dna = {'F': ['TTT', 'TTC'],
                 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
                 'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
                 'Y': ['TAT', 'TAC'],
                 'C': ['TGT', 'TGC'],
                 'W': ['TGG'],
                 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
                 'H': ['CAT', 'CAC'],
                 'Q': ['CAA', 'CAG'],
                 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                 'I': ['ATT', 'ATC', 'ATA'],
                 'M': ['ATG'],
                 'T': ['ACT', 'ACC', 'ACA', 'ACG'],
                 'N': ['AAT', 'AAC'],
                 'K': ['AAA', 'AAG'],
                 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
                 'A': ['GCT', 'GCC', 'GCA', 'GCG'],
                 'D': ['GAT', 'GAC'],
                 'E': ['GAA', 'GAG'],
                 'G': ['GGT', 'GGC', 'GGA', 'GGG']
                 }

    def __init__(self, sequence, peptide):
        '''
        declare necessary vars
        '''

        self.seq = sequence
        self.pep = peptide

    def reverseSeq(self, sequence):
        '''

        used to find reverse compliment -- method taken from bme 160, very useful
        returns the reverse compliment

        https://www.w3schools.com/python/ref_string_translate.asp
        https://www.programiz.com/python-programming/methods/string/maketrans

        '''

        #(sequence.translate(str.maketrans("ATCG", "TAGC"))[::-1])
        return sequence.translate(str.maketrans("ATCG", "TAGC"))[::-1]

    def genDNA(self, peptide):
        '''
        generates dna sequence from peptide sequence

        https://stackoverflow.com/questions/231767/what-does-the-yield-keyword-do
        '''

        seq = []

        # never starts true (out peptide will be at least 2 letters) but if peptide seq was hypothetically 1
        if len(peptide) == 1:
            # for codon in amino acid to dna dictionary matching peptide
            for codon in self.aa_to_dna[peptide]:
                # append it to seq list
                seq.append(codon)

        # otherwise (peptide is >1 letters)
        else:
            prefix = peptide[:-1] # cut all the way up to last element
            suffix = peptide[-1] # last element

            # for each codon that matches the suffix peptide
            for codon in self.aa_to_dna[suffix]:
                # for each dna generations of the prefix (remaining part of peptide)
                for p in self.genDNA(prefix): # call genDNA on prefix (peptide sequence minus last peptide)
                    # creates prefix + all suffix generations, but also recursively calls on prefix (since we
                    # can [will] have a prefix of greater than 1 element, as seen in given file from rosalind)
                    # for sample MA -- M only ATG, but A has 4 codon, we need to add all 4
                    seq.append(p + codon)

        # return seq list, which will contain every single formation of peptide sequence (not including reverseSeq)
        return seq

    def substring(self, sequence, peptide):
        '''
        generate substrings of dna sequence

        "We say that a DNA string Pattern encodes an amino acid string Peptide if the RNA string transcribed
        from either Pattern or its reverse complement Pattern translates into Peptide."
        '''

        n = len(sequence)
        k = len(peptide)*3

        encoding = []
        substrings = []

        # for each sequence in generated DNA from peptide
        for seq in self.genDNA(peptide):
            # we also want the reverse sequence, add them both to a list
            rs = self.reverseSeq(seq)
            encoding.append(seq)
            encoding.append(rs)
            # sample MA encoding: ['ATGGCT', 'AGCCAT', 'ATGGCC', 'GGCCAT', 'ATGGCA', 'TGCCAT', 'ATGGCG', 'CGCCAT']

        # range of (len(seq) - len(pep)*3 + 1) -- sample == (52 - 6 + 1)
        # range(n) also works, dont notice time difference in program runtime
        # ensures that we are calculating across proper range (peptide = codon = 3 nuc)
        for i in range(n):
            # cut substring from i up to i+6 (sample)
            substring = sequence[i:i+k]

            # if substring already seen in encoding
            if substring in encoding:
                # append it to substrings list
                substrings.append(substring)

        #return substrings list
        return substrings

    def print(self):
        '''
        quick method to print all substrings of a genome encoding a given amino acid sequence
        '''

        substrings = self.substring(self.seq, self.pep)

        # for substring in substrings list
        for substring in substrings:
            # print substring
            print(substring)


def main(inFile = None):
    '''
    main() opens input file and reads it, then instantiates sequence into class peptide and calls print()
    which calls the remaining methods to create all substrings of a genome encoding a given AA sequence

    in main(inFile = '') you can input your own file for testing, as long as format follows 2 lines, with
    first line being given AA sequence, and second line being the peptide

    '''

    count = 0
    seq = ''
    pep = ''

    f = open(inFile)
    for line in f:
        if count == 0:
            seq = line[0:].strip()
        if count == 1:
            pep = line[0:].strip()
        count += 1

    Pep = peptide(seq, pep)
    Pep.print()

if __name__ == "__main__":
    main(inFile='rosalind_ba4b.txt')