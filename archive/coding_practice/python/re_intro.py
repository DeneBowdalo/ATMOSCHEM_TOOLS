import re

# make a regular expression object from a string
# the pattern is:
# an 'A' followed by 3 characters each of which must be either M, K or W
re_obj = re.compile('A[AWGV]{4}')

for line in open('globlins630.fa'):  # compact way of reading line
    j = re_obj.search(line)
    if re_obj.search(line):         # returns None if no match
        print(j.span(),line)
