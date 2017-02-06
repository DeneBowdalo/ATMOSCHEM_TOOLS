number = int(input("Input a number between one & a million - ")) # FIRST, set the initial value of the variable a to 0(zero).
while 0<number<1000000:    # While the value of the variable a is less than 10 do the following:
    if number > 500000:    # Increase the value of the variable a by 1, as in: a = a + 1! 
        number += 1
        print(number)     
    
    if number < 500000:
        number -= 1
        print(number)

