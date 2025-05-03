# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough(wschloug)
# Group Members: N/A


def main():

    userInput = input("Input nucleotide sequence: ")
    string = str(userInput)

    nucDict = {"A":0, "T":0, "G":0, "C":0} # starting key values are all 0

    for nuc in string: # for nucleotide in string input
        if nuc in nucDict: # if nucleotide matches A, T, G, or C
            nucDict[nuc] += 1  # nuc matching key gets value +1 (everytime its counted)

    print(nucDict["A"], nucDict["C"], nucDict["G"], nucDict["T"]) # print in ACGT order --> specified by rosalind Q

if __name__ == "__main__":
    main()