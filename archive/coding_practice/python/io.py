football_file = open("football.txt","w")
football_file.write("Arsenal\nAstonVilla\nChelsea\nEverton\nLiverpool\nManchester United\nManchester City\nNewcastle United\nTottenham Hotspur")
football_file.close()

football_file = open("football.txt","r")
read_football_file = football_file.read()
print(read_football_file)
football_file.close()

rugby_file = open("rugby.txt", "w") 
rugby_teams = ["\n\nBradford Bulls\n", "Leeds Rhinos\n", "Wigan Warriors\n", "Warrington Wolves\n", "St Helens\n", "Hull FC\n"]
rugby_file.writelines(rugby_teams)
rugby_file.close()

rugby_file = open("rugby.txt", "r")
print(rugby_file.read())
rugby_file.close()
