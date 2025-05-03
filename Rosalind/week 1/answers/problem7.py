# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough(wschloug)
# Group Members: N/A

def main():
    userInput = str(input("Input a string <= 10000 letters: "))

    wordDict = {}
    words = userInput.split()
    for word in words:
        if word in wordDict:
            wordDict[word] += 1
        else:
            wordDict.update({word:1})

    for word in wordDict:
        print(word, wordDict[word])

if __name__ == "__main__":
    main()