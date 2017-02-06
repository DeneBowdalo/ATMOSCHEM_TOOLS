number1 = int(input("Please enter number: "))
number2 = int(input("Please enter number: "))
number3 = int(input("Please enter number: "))
number4 = int(input("Please enter number: "))

limit = 21

score = number1+number2+number3+number4

if score > limit:
    print ("Your score", score, "exceeded the limit", limit)

else:
    print ("Your score", score, "did not exceed the limit", limit)
