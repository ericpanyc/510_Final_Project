import os
import csv

file_name_dict = dict()

with open("filename_conversion.csv") as fin:
    data = csv.reader(fin)
    for row in data:
        old_name = row[0]
        new_name = row[1]
        #print(type(old_name), type(new_name))
        file_name_dict[old_name] = new_name
        #print(new_name)

for filename in os.listdir("."):
    if filename in file_name_dict:
        #print(filename)
        new_name = file_name_dict[filename]
        os.rename(filename, new_name)

