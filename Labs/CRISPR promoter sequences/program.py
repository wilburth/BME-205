#!/usr/bin/env python3

########################################################################
# File:program.py
#  executable: program.py
# Purpose:scaffold file for python executables
#   stderr: errors and status
#   stdout:
#          
# Author: David Bernick
# History:      dlb 08/20/2011 Created
#               dlb 08/26/2021 minor updates
########################################################################

########################################################################
# CommandLine
########################################################################
class CommandLine() :
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
    
    def __init__(self, inOpts=None) :
        '''
        CommandLine constructor.
        
        Implement a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output' 
                                             )
        self.parser.add_argument('-o', '--option', action = 'store', help='foo help')
        self.parser.add_argument('-i', '--integer', type=int, choices= range(5, 10), action = 'store', help='help for a boolean option')
        self.parser.add_argument('-c', '--character', choices ='abcdef', action = 'store', help='help for a charcter option')
        self.parser.add_argument('-b', '--bool', action = 'store', nargs='?', const=True, default=False, help='boolean switch')
        self.parser.add_argument('-r', '--requiredBool', action = 'store', nargs='?', required=True,const=True, default=False, help='required boolean switch')
        self.parser.add_argument('-l', '--list', action = 'append', nargs='?', help='list help') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  

        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)
  

class Usage(Exception):
    '''
    Signal a Usage error, evoking a usage statement and eventual exit when raised.
    '''
    def __init__(self, msg):
        self.msg = msg


########################################################################
# The Bioinformatics Class
# Here is where most of your solution will end up.
#
########################################################################
class YourSpecificStuff:
    '''
    Build the thing that needs doing.
    '''
    def __init__(self, someParameter):
        self.data = someParameter
        
    def someMethod (self):
        """
        Do something to the data
        """
        pass
########################################################################
# Main
# Here is the main program
# This is a place to instantiate objects, read and write some data for
# those objects and maybe implement some of the command line options.
########################################################################

def main(myCommandLine=None):
    '''
    Implement the Usage exception handler that can be raised from anywhere in process.

    '''
    if myCommandLine is None:
        myCommandLine = CommandLine()  # read options from the command line
    else :
        myCommandLine = CommandLine(myCommandLine) # interpret the list passed from the caller of main as the commandline.

    try:
        
        print (myCommandLine.args) # print the parsed argument string .. as there is nothing better to do

        if myCommandLine.args.requiredBool:
            print ('requiredBool is', str(myCommandLine.args.requiredBool) )  ## this is just an example
        else :
            pass
        raise Usage ('testing')  # this is an example of how to raise a Usage exception and you can include some text that will get printed. Delete this if you dont need it

    except Usage as err:
       print (err.msg)

if __name__ == "__main__":
    main()
#   main(['-r'])  # this would make this program.py behave as written

# if you want to make use of a test commandLine, you could do it like this ( notice that it is a list of strings ):
#   main([ '--bool',
#          '--character=b',
#          '-i=10'])

