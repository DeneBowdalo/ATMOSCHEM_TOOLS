#! /usr/bin/python
#-*-coding: utf-8 -*-
# converts temperature to Fahrenheit or Celsius
 
def print_options():
    print("Options:")
    print(" 'p' print options")
    print(" 'c' convert from Celsius")
    print(" 'f' convert from Fahrenheit")
    print(" 'q' quit the program")
 
def celsius_to_fahrenheit(c_temp):
    return 9.0 / 5.0 * c_temp + 32
 
def fahrenheit_to_celsius(f_temp):
    return (f_temp - 32.0) * 5.0 / 9.0
 
choice = "p"
while choice != "q":
    if choice == "c":
        temp = float(input("Celsius temperature: "))
        print("Fahrenheit:", celsius_to_fahrenheit(temp))
        choice = input("option: ")
    elif choice == "f":
        temp = float(input("Fahrenheit temperature: "))
        print("Celsius:", fahrenheit_to_celsius(temp))
        choice = input("option: ")
    elif choice == "p": #Alternatively choice != "q": so that print when anything unexpected inputed
        print_options()
        choice = input("option: ")
