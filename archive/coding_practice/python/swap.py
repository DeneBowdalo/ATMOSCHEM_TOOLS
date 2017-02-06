n_of_colours =int(input("Please input a number of colours to put into list: "))
colours =[]

for i in range(n_of_colours):
    colours.append(input("Please input a colour:"))

for i in range(len(colours)):
    if colours[i] == "blue":
        if i > colours.index("red"):
            index_red = colours.index("red") 
            colours[i], colours[index_red] = colours[index_red], colours[i]
        
        if i > colours.index("white"):
            index_white = colours.index("white")
            colours[i], colours[index_white] = colours[index_white], colours[i]
                

        
    elif colours[i] == "red":
        if i > colours.index("white"):
            index_white2 = colours.index("white")        
            colours[i], colours[index_white2] = colours[index_white2], colours[i]
    
print(colours)
