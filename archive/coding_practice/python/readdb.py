import fileinput

database = {}

for line in fileinput.input():
    line_list = line.split()
    database[line_list[0]] = line_list[1:] 
    
print(database.items())
