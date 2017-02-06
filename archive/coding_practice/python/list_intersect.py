list_1 =[]
list_2 =[]
list_3 =[]
def create_list(list_name):
 
    print("Please enter 4 Numbers for list\n")
    list_name.append(int(input("Please enter number: ")))
    list_name.append(int(input("Please enter number: ")))
    list_name.append(int(input("Please enter number: ")))
    list_name.append(int(input("Please enter number: ")))
    print("")

create_list(list_1)
create_list(list_2)

list_3.append(set(list_1) & set(list_2))        
print("Numbers that match between the 2 lists are: ", list_3)
