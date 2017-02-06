# Plays the guessing game higher or lower 
 
# This should actually be something that is semi random like the
# last digits of the time or something else, but that will have to
# wait till a later chapter.  (Extra Credit, modify it to be random
# after the Modules chapter)
from random import randrange 

number = randrange(0,101)

guess = -1
count = 0
 
print("Guess the number!")
while guess != number:
    guess = int(input("Is it... "))
    count += 1
    if guess == number:
        print("Hooray! You guessed it right!")
        if count >= 10:
            print("That was tedious!")
    elif guess < number:
        print("It's bigger...")
    elif guess > number:
        print("It's not so big.")
