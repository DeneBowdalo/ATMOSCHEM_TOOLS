name = "Guillermo"

print("Guess my name, go on.")
print("I'll give you 3 tries")
print("No more, no less.")
print("And then you're locked out.")

guess = input("Go Ahead and Guess: ")

count = 1
while count < 3 and guess.lower() != name: 
    print("Gutted, try again.")
    guess = input("Try again, what is my name? ")
    count += 1
    
if guess.lower() != name: 
    print("No more guesses, chump. The name was", name)

else: 
    print("I don't know how you did it, but congratulations.")



