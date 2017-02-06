#create program to pick random word and let player guess 

import random

WORD = ('dog','cat','me','animal','liverpool','bone','science','fraudulant')

word_rand = random.choice(WORD)
right_letters=[]
print("Welcome.\nThe challenge is to try guess a random word.\nYou have 5 chances to ask if a letter is in a word.\nAfter that you must guess.\nGood luck!")

for i in range(5):
    guess = input("Please choose a letter:  ")
    guess = guess.lower()
    if guess not in word_rand:
        print(guess, "is not in the word.")
    elif guess in word_rand:
        print(guess, "is in the word.")
        right_letters += guess

print("Your 5 guesses are up.\n Your must now guess the word.\n")

if right_letters:
    print("To remind you the letter/s:", right_letters, ",are in the word you are guessing")

else:
    print("No letters you guesses are in the word you are guessing.")

final_guess = input("Please guess what you think the word is:  ")
final_guess = final_guess.lower()
if final_guess == word_rand:
    print("Congratulations, you won!\nThe word was indeed:", word_rand)
else:
    print("You are wrong, you have lost on this occasion.\nThe word was:", word_rand)

print("Thank you for playing")

input("\nPress the enter key to exit\n")
