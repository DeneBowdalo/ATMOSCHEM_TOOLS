import fileinput
import sys

def read(filename):
    filename = filename + str('.fa')
    for line in fileinput.input(filename):  
        find_result = line.find('>')
        if find_result != -1: 
            headers = open(headers_file,'a')
            print(line,file=headers)
        else: 
            sequences = open(sequences_file,'a')
            print(line,file=sequences)   


filea = input('Input the name of the FASTA file to open:  ')
headers_file = input('What filename would you like to give to the headers?')    
sequences_file = input('What filename would you like to give to the sequences?')
 
while 1==1:
    try:
        read(filea)
        sys.exit()
    except IOError: 
        print('Invalid FASTA filename, please enter valid name')
        filea = input('Input the name of the FASTA file to open:  ')
        

