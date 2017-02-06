## This program runs a test of knowledge
 
# First get the test questions
# Later this will be modified to use file io.

    # notice how the data is stored as a list of lists
questions = [["What color is the daytime sky on a clear day? ", "blue"],
            ["What is the answer to life, the universe and everything? ", "42"],
            ["What is a three letter word for mouse trap? ", "cat"]]
 
# This will test a single question
# it takes a single question in
# it returns True if the user typed the correct answer, otherwise False
 
def check_question(question_and_answer):
    # extract the question and the answer from the list
    # This function take a list with two elements, a question and an answer.  
    question = question_and_answer[0]   
    answer = question_and_answer[1]
    # give the question to the user
    given_answer = input(question)
    # compare the user's answer to the testers answer
    if answer == given_answer:
        print("Correct")
        return True
    else:
        print("Incorrect, correct was:", answer)
        return False
 
# This will run through all the questions
def run_test(questions):
    if len(questions) == 0:
        print("No questions were given.")
        # the return exits the function
        return
    index = 0
    right = 0
    while index < len(questions):
        # Check the question
        #Note that this is extracting a question and answer list from the lists of lists.
        if check_question(questions[index]): 
            right = right + 1
        # go to the next question
        index = index + 1
    # notice the order of the computation, first multiply, then divide
    print("You got", right * 100 / len(questions),\
           "% right out of", len(questions))
    return 
# now let's get the questions from the get_questions function, and
# send the returned list of lists as an argument to the run_test function.

option = ""

while option != 9:
    print("Test Version 1.1")
    print("Options:")
    print("1 - Take the Test")
    print("2 - View the list of Questions")
    print("9 - Quit") 
   
    option = int(input("Please input a number: "))
 
    if option == 1: 
        run_test(questions)

    elif option == 2:    
        q = 0
        while len(questions) > q:
            a = 0    
            print(questions[q][a])
            q += 1
   
    elif option == 9:
        print("Adios")

    else:
        print("Input valid number please")
        
        
