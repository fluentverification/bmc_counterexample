import csv
input_file = open("./circuit/circuit_single_species_plus.txt", "r")
csv_file = "./csv/circuit/circuit_single_species_plus.csv"

csv_data = []
# temp = [None] * 3
temp = []

for line in input_file:
    if not (("==================================================" in line)):
        if "probability=" in line or "K = " in line:   
            temp.append(line[line.find("= ")+2 : line.find("\n")])
            # temp[1] = line[line.find("= ")+2 : line.find("\n")]
        if "number of transitions:" in line:
            temp.append(line[line.find(": ")+2 : line.find("\n")])
            # temp[0] = line[line.find(":")+2 : line.find("\n")]
        if "time=" in line and "total" not in line:
            temp.append(line[line.find("= ")+2 : line.find("\n")])
            # temp[2] = line[line.find("= ")+2 : line.find("\n")]


    else:
        csv_data.append(temp)
        temp = []
        # temp = [None] * 3


# csv_data.append(temp)
with open(csv_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(csv_data)

print("CSV file has been created successfully.")