# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough (wschloug)
# Group Members: Junmo, Mikey, Wenyu


class idealMatch:
    '''
    class idealMatch contains all methods required to calculate the cyclic peptide whose theoretical spectrum matches
    the experimental spectrum
    The prompt is as follows:
    Given an ideal experimental spectrum, find a cyclic peptide whose theoretical spectrum matches the experimental spectrum.
    Given: A collection of (possibly repeated) integers Spectrum corresponding to an ideal experimental spectrum.
    Return: Every amino acid string Peptide such that Cyclospectrum(Peptide) = Spectrum (if such a string exists).
    Sample:
    0 113 128 186 241 299 314 427
    '''

    # mass list for aa
    masses = [
        57, 71, 87, 97, 99,
        101, 103, 113, 114, 115,
        128, 129, 131, 137, 147,
        156, 163, 186
    ]

    def __init__(self, spectrum):
        '''
        initialize necessary vars
        '''

        self.spectrum = spectrum

    def expand(self, peptides):
        '''
        define EXPAND(Peptides) as a new collection containing all possible extensions of peptides in Peptides
        by a single amino acid mass.
        '''

        # for peptide in copy of peptides (do this so we can remove)
        for peptide in peptides.copy():
            # remove the peptide from original
            peptides.remove(peptide)
            # for mass in masses dictionary
            for m in self.masses:
                # if peptide (true)
                if peptide:
                    # create new peptide which is a junction of original mass + new mass from mass list
                    # (does this for each entry)
                    new_pep = peptide + '-' + str(m)
                # otherwise (not peptide)
                else:
                    # new peptide is just set to a mass from mass list
                    new_pep = str(m)
                # add to peptide list whatever the new peptide was
                peptides.add(new_pep)

        return peptides


    def calcMass(self, peptide):
        '''
        quick method to calculate mass of peptide
        '''
        # map will create a map obj with items from peptide, split at '-' if conjunction
        # then, we can sum the items in the map obj to find the mass of the conjunction peptide
        mass = sum(map(int, peptide.split('-')))
        return mass

    def cycloSpec(self, peptide):
        '''
        method to calculate (sorted) cyclospectrum of input peptide
        utilizes map(), sum()

        sample example:
        [128, 113, 186]
        128
        113
        186
        sum1-2: 241
        sum2-3: 299
        residual: [186, 128]
        residual sum (3-1): 314
        '''

        cs = [0]
        # create map obj of integers of peptides split at '-',
        # then create a list from map obj
        comp = list(map(int, peptide.split('-')))
        n = len(comp)
        # search in range of 1 up to but not including n (sample == 1 up to but not including 3)
        for k in range(1, n):
            # search in range of 0 up to but not including n (sample == 0 up to but not including 3)
            for i in range(n):
                # if range1 + range2 less than total length of composition
                if i + k <= n:
                    # cut the comp from location i to location i+k and append it to our list
                    cs.append(sum(comp[i:i+k]))
                # else (range1 + range2 > total length of composition)
                else:
                    # residual set to i+k-n
                    r = i + k - n
                    cs.append(sum(comp[i:] + comp[:r]))

        cs.append(self.calcMass(peptide))
        cs.sort()
        return cs

    def linearSpec(self, peptide):
        '''
        similar to cycloSpec, method to calculate (sorted) linear spectrum of peptide
        utilizes map(), sum()

        sample example:
        comp cut: [186, 128]
        comp sum : 314
        '''

        cs = [0]
        # create map obj of integers of peptides split at '-',
        # then create a list from map obj
        comp = list(map(int, peptide.split('-')))
        n = len(comp)

        # search across range of len(comp)
        for k in range(n):
            # search across range of (len(comp) + 1)
            for i in range(n + 1):
                # add a sum of comp cut from i to i+k
                cs.append(sum(comp[i:i+k]))

        cs.append(self.calcMass(peptide))
        cs.sort()
        return cs

    def notConsistent(self, spectrum, peptide):
        '''
        method nonConsistent will return false if the subpeptide is consistent,
        and true is the subpeptide isn't consistent

        "The key to our new algorithm is that every linear subpeptide of a cyclic peptide Peptide is consistent with
        CYCLOSPECTRUM(Peptide). Thus, to solve the Cyclopeptide Sequencing Problem for Spectrum, we can safely ban
        all peptides that are inconsistent with Spectrum from the growing set Peptides, which powers the
        bounding step that we described above." -- textbook pg 196

        '''

        # generate linearSpec of peptide
        peptide_spectrum = self.linearSpec(peptide)
        # for each value in linear spectrum
        for val in peptide_spectrum:
            # if the val isn't in the spectrum, return true (meaning it wasn't consistent)
            if val not in spectrum:
                return True
        # otherwise return false (meaning is was consistent)
        return False

    def cyclopeptideSequencing(self, spectrum):
        '''
        pseudocode taken from texbook/rosalind prompt, methods created are based off pseudocode
        https://www.programiz.com/python-programming/methods/built-in/set

        from pseudocde, will know i want expand(), caluclate-mass(), cyclospectrum(), not consistent(), also linear()
        because it looks like not-consistent() will use the linear spectrum

        '''

        peptides = {''}
        matches = set() # use set so we can use add()

        # while peptides not empty
        while peptides:
            # call expand on peptides
            peptides = self.expand(peptides)
            # for each peptide in peptides (use copy)
            for peptide in peptides.copy():
                # if mass == parent mass (max(spectrum))
                if self.calcMass(peptide) == max(spectrum):
                    # if cyclospectrum of peptide == spectrum
                    if self.cycloSpec(peptide) == spectrum:
                        # add it as a match
                        matches.add(peptide)
                    # since we used a copy, its ok for us to remove form original
                    peptides.remove(peptide)
                # else if peptide not consistent with spectrum
                elif self.notConsistent(spectrum, peptide):
                    # since we used a copy, its ok for us to remove from original
                    peptides.remove(peptide)
        # return all matches
        return matches

    def print(self):
        '''
        quick method to print all elements returned form cyclopeptideSequencing() after calling it with spectrum
        '''

        print(*self.cyclopeptideSequencing(self.spectrum))


def main(inFile = None):
    '''
    learning how to use map() -- seems like a super useful function that i can use throughout
    this program, along with in future assignments (hopefully :P)
    https://www.geeksforgeeks.org/python-map-function/
    https://realpython.com/python-map-function/
    '''

    spectrum = list()
    f = open(inFile)
    for line in f:
        spectrum = list(map(int, line.split()))

    # for sample, spectrum being taken into idealMatch looks like this:
    # [0, 113, 128, 186, 241, 299, 314, 427]

    pep = idealMatch(spectrum)
    pep.print()




if __name__ == "__main__":
    main(inFile='rosalind_ba4e.txt')