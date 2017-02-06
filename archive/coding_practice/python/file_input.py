def process(line):
    print(line, line)


import fileinput
for line in fileinput.input():  # 'line' includes the newline character
    process(line)
