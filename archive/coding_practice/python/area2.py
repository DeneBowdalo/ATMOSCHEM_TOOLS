#! /usr/bin/python
#-*-coding: utf-8 -*-
# calculates a given rectangle area
print
def hello():
    print('Hello!')

def print_options():
    print("Options: ")
    print("R - Find the area of a rectangle")
    print("S - Find the area of a square")
    print("C - Find the area of a circle")
    print("Q - Quit the program")

def area_rectangle(width, height):
    return width * height

def area_square(length):
    return length**2

def area_circle(radius):
    return 3.14 * radius**2

def print_welcome(name):
    print('Welcome,', name)
 
def positive_input(prompt):
    number = float(input(prompt))
    while number <= 0:
        print('Must be a positive number')
        number = float(input(prompt))
    return number
 
name = input('Your Name: ')
hello()
print_welcome(name)


option = ""

while option != "q":
    print_options()
    option = input("Which shape's area do you wish to calculate? ")

    if option == "r":
        print
        print('To find the area of a rectangle,')
        print('enter the width and height below.')
        print
        w = positive_input('Width: ')
        h = positive_input('Height: ')
        print('Width =', w, 'Height =', h, 'so Area =', area_rectangle(w, h))

    elif option == "s":
        print('To find the area of a square,')
        print('enter the side length below.')
        print
        l = positive_input('Length: ')
        print('Length =', l, 'so Area =', area_square(l))

    elif option == "c":
        print('To find the area of a circle,')
        print('enter the radius below.')
        rad = positive_input('Radius: ')
        print('Radius =', rad, 'so Area =', area_circle(rad))

    elif option == "q":
        print('Goodbye!')
    
    else: 
        print('Unrecognised option.')
        print_options()
