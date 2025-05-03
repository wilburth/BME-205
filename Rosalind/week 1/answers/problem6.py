# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough(wschloug)
# Group Members: N/A


def main():
    f = open('rosalind_ini5.txt', 'r')  # open the input file with 'r'
    lineList = f.readlines()            # read the lines in input file
    for line in range(len(lineList)):   # for each line in the range of length of lines
        if (line+1)%2 is 0:             # technically the even lines we read are actually read as odd lines in code (hence the +1)
            print(lineList[line])       # print the even lines!


if __name__ == "__main__":
    main()