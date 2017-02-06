username = input("Set your username: ")
password = input("Set your password: ")

while 1 == 1:
    username_attempt = ""
    password_attempt = ""
    key = ""
    while (username_attempt != username) or (password_attempt != password):
        username_attempt = input("Username: ")
        password_attempt = input("Password: ")
    print("Welcome,", username,".","Type lock to lock.")
    while key != "lock":
        key = input("")
                   
