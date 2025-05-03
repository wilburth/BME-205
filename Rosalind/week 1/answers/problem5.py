# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough(wschloug)
# Group Members: N/A


def main():
    userInput = input("Input two positive integers with a < b < 10000:")
    sumOdd = 0 # temp variable
    a, b = userInput.split() # since we know its just 2 numbers, we can use .split()
    for i in range(int(a), int(b) + 1): # for i in range of integer(a) to integer(b)+1     (to count the end)
        if i%2 == 1: # if the integer divided by 2 has a remainder of 1 --> meaning its odd
            sumOdd += i # start summing them up!

    print(sumOdd) # print filled temp variable "sumOdd" to system output

if __name__ == "__main__":
    main()