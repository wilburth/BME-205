# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough (wschloug)
# Group Members: Junmo, Mikey, Wenyu

class theoretical:
    '''
    class theoretical contains all methods required to solve rosalind problem17
    The prompt is as follows:

    Generate the theoretical spectrum of a cyclic peptide.
    Given: An amino acid string Peptide.
    Return: Cyclospectrum(Peptide).
    Sample:
    LEQN

    '''

    # http://rosalind.info/media/problems/ba4c/integer_mass_table.png
    aa_masses = {
        'G': 57, 'A': 71, 'S': 87, 'P': 97,
        'V': 99, 'T': 101, 'C': 103, 'I': 113,
        'L': 113, 'N': 114, 'D': 115, 'K': 128,
        'Q': 128, 'E': 129, 'M': 131, 'H': 137,
        'F': 147, 'R': 156, 'Y': 163, 'W': 186
    }

    def __init__(self, peptide):
        '''
        initialize necessary vars
        '''

        self.pep = peptide

    def cyclic(self, peptide):
        '''
        generate cyclic peptide from input peptide
        '''

        subPs = [''] # starting an single element string list ensures we have 0 slot
        n = len(peptide)

        # range of length of peptide (sample == 0 up to but not including 4)
        for i in range(0, n):
            # for range of length of peptide -1 (sample == 1 up to but not including 4)
            for j in range(1, n):
                # if i + j less than or equal to len of peptide (4)
                if i + j <= n:
                    # cut peptide from i to i + j, append it to sub peps list
                    subPs.append(peptide[i:i+j])
                # otherwise (i + j > n)
                else:
                    # to ensure that we also get "wraparound" (cyclic) peptide(s)
                    # example: LEQN will have 3 residual cuts -- QNL, NL, and NLE -- all other cuts are accounted for
                    residual = i + j - n
                    subPs.append(peptide[i:] + peptide[:residual])

        subPs.append(peptide)
        return subPs


    def cyclospectrum(self, peptide):
        '''
        generates theoretical spectrum (cyclospectrum) of cyclic peptide

        The theoretical spectrum of a cyclic peptide Peptide, denoted Cyclospectrum(Peptide), is the collection of
        all of the masses of its subpeptides, in addition to the mass 0 and the mass of the entire peptide.
        We will assume that the theoretical spectrum can contain duplicate elements, as is the case for
        "NQEL" (shown in Figure 2), where "NQ" and "EL" have the same mass.

        '''

        # generate cyclic peptide
        # begin spectrum list (of ints)
        subPs = self.cyclic(peptide)
        spectrum = []

        # for string in sub peptides
        for string in subPs:
            mass = 0
            # for each cut (character) in string
            for cut in string:
                # add to mass the matching mass of the AA
                # ensures that if string is LE, we get mass of both L and E added
                mass += self.aa_masses[cut]
            # append it to spectrum list
            spectrum.append(mass)

        return sorted(spectrum)


    def print(self):
        '''
        quick method to print all elements of spectrum
        '''

        spectrum = self.cyclospectrum(self.pep)
        print(*spectrum)

def main(inFile = None):
    '''
    main opens the input file and reads form it, then adds to string 'seq' and passes that to class theoretical
    class print which will print the cyclospectrum of the input string
    '''

    seq = ''
    f = open(inFile)
    for line in f:
        seq = line.strip()

    peptide = theoretical(seq)
    peptide.print()

if __name__ == "__main__":
    main(inFile='rosalind_ba4c.txt')