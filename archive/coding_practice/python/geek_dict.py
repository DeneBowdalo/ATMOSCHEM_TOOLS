#Geek Translator
#Demonstates using Dictionaries

geek = {"404": "clueless. From the web error message 404, meaning page not found.",
        "Googling": "searching the internet for background information on a person.",
        "Keyboard Plaque" : "the collection of debris found in computer keyboards.",
        "Link Rot" : "the process by which web page links become obselete.",
        "Percussive Maintenance" : "the act of striking an electronic device to make it work.",
        "Uninstalled" : "being fired. Especially popular during the dot-bomb era."} 

choice = None
while choice != "0":

    print(
    """

    Geek Translator

    0 - Quit
    1 - Look up a Geek Term
    2 - Add a Geek Term
    3 - Redefine a Geek Term
    4 - Delete a Geek Term
    5 - Print Full Dictionary
    """

    )


    choice = input("Choice: ")
    print()

#exit
    if choice == "0":
        print("Goodbye.")

# get a definition
    elif choice == "1":
        term = input("What term do you want me to translate?: ")
        if term in geek:
            definition = geek[term]
            print("\n", term, "means", definition)
        else:
            print("\nSorry , I don't know", term)

# add a term-definition pair
    elif choice == "2":
        term = input("What term do you want me to add?: ")
        if term not in geek:
            definition = input("\nWhat's the definiton?: ")
            geek[term] = definition
            print("\n", term, "has been added.")
        else:
            print("\nThat term already exists! Try redefining it.")

#redefine existing term
    elif choice == "3":
        term = input("What term do you want me to redefine?: ")
        if term in geek:
            definition = input("What's the new definition?: ")
            geek[term] = defintion
            print("\n", term, "has been redefined.")
        else:
            print("\nThat term doesn't exist! Try adding it.")

# delete a term-defintion pair
    elif choice == "4":
        term = input("What term do you want me to delete?: ")
        if term in geek:
            del geek[term]
            print("\nOkay, I deleted", term)
        else:
            print("\I can't do that!", term, "doesn't exist in the dictionary.")

# print entire dictionary
    elif choice == "5":
        print("\nEntire Dictionary:")
        for term in geek.keys():
            print("\nGeek Term: ", term, "\tMenaing:", geek[term]) 
        

# Some unknown choice
    else:
        print("\nSorry, but", choice, "isn't a valid choice.")

input("\n\nPress the enter key to exit.")


        
