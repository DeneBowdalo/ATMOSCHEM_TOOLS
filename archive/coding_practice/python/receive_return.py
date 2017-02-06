#receive and return
#demonstates paramters and return values

def display(message):
    print(message)

def give_me_five():
    five = 5
    return five

def ask_yes_no(question):
    """Ask a yes or no question."""
    response = None
    while response not in ("y","n"):
        response = input(question).lower()
    return response

#main
display("Here's a messgae for you.\n")

number = give_me_five()
print("Here's what I got from give_me_five():", number)

answer = ask_yes_no("\nPlease enter 'y' or 'n' : ")
print("Thanks for answering:", answer)

input("\n\nPress the enter key to exit")
