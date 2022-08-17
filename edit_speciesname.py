
file = open("test.txt")
output = open("genera.txt", "w")
lines = file.readlines()
for line in lines:
    #name = "".join(filter(str.isalpha, line))
    genus = line.split(" ")[0]
    species = line.split(" ")[1]
    output.write(genus + "\n")



