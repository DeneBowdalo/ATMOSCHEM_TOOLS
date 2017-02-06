import fileinput

count = {}
for key in ['A', 'C', 'G', 'T']: 
    count[key] = 0


for line in fileinput.input():
    count['A'] = count['A'] + line.count('A') 
    count['C'] = count['C'] + line.count('C') 
    count['G'] = count['G'] + line.count('G') 
    count['T'] = count['T'] + line.count('T')
   
    
print("Number of A's:", count['A'])
print("Number of C's:", count['C'])
print("Number of G's:", count['G'])
print("Number of T's:", count['T'])


