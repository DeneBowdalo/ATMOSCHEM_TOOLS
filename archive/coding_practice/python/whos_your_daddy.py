# Program to produce name of father from son's name

father_son = {}

print("Welcome\n\nThis program will alow you to retrive the name of the father, from the name of a son entered.")

choice = -1

while choice != 0:

    choice = int(input("Please choose an option.\n\nTo view the name of a father of a entered son, type 1.\nTo add a father/son pair, type 2\nTo replace a father/son pair, type 3\nTo delete a father son/pair, type 4.\nTo view all the first names of the sons in the system, type 5\nTo view all the first names of the fathers in the system, type 6\nTo view the whole system father/son pairs, type 7\n\nType 0 to exit.\n"))

    if choice == 1:
        son = input("Please enter the first name of a son:  ")
        print("The first name of the father of son:", son, "is:", father_son[son])
    elif choice == 2:
        son = input("Please enter the first name of a son:  ")
        father = input("Please enter the first name of a Father:  " )
        father_son[son] = father
        print("OK, the son:", son, "has been assigned the father:", father)

    elif choice == 3:
        son = input("Please input the first name of a son in the system:  ")
        if son in father_son:
            del father_son[son]
            new_son = input("What is the first name of the son you would like to replace the old name with?  ")
            new_father = input("What is the first name of the father, of the son?  ") 
            father_son[new_son] = new_father
            print("OK, the son:", son, "has been replaced with the new son:", new_son, "and the new father:", new_father)
        else:
            print("That son's name is not in the system.")
        
    elif choice == 4:
        son = input("Please input the first name of a son in the system:  ")
        if son in father_son:
            del father_son[son]
            print("The father/son pair originating from the son:", son, "have been deleted from the system.")
        else:
            print("That son's name is not in the system.")
    
    elif choice == 5:
        print("The son first names in the system are:", father_son.keys())
        

    elif choice == 6:
        print("The father first names in the system are:", father_son.values())

  
    elif choice == 7:
        print("The whole system father/son pairs are:", father_son.items())

input("Press tne enter key to exit")

