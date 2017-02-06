scores = []

choice  = None
while choice != "0":
    print(
    """
    High Scores 2.0
 
    0 - Quit
    1 - List Scores
    2 - Add a Score
    """
    )

    choice = input("Choice: ")
    print()
    #exit
    if choice == "0":
        print("Goodbye.")

#Display high scores
    elif choice == "1":
        print("High Scores\n")
        print("NAME\tSCORE")
        for entry in scores:
            score, name = entry
            print(name, "\t", score)

#add a score

    elif choice == "2":
        name = input("What is the player's name?: ")
        score = int(input("What score did the player get?: "))
        entry = (score, name)
        scores.append(entry)
        scores.sort(reverse=True)
        scores = scores[:5]

#some unknown choice

    else:
        print("Sorry, but", chocie, "isn't a valid choice.")


    input("\n\nPress the enter key to exit.")
