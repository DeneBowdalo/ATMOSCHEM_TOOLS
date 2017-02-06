listed = [('Dene', 28),('Paul', 19),('Carmel',5)]

list1 = [x[0] for x in listed]
list2 = [x[1] for x in listed]

print(list1)
print(list2)

list3, list4 = zip(*listed)
print(list3)
print(list4)
