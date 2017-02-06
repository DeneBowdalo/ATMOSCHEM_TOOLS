seq = 'acgtatgctagctatagcttaggatgctagtcgatcgatcgatcgatcaggg'

count = {}


for i in range(len(seq)-2):
    triple = seq[i:i+3]
    count[triple] = count.get(triple,0) + 1

for triple, total in count.items():
    print('{0} occurrences of {1}'.format(total,triple))
