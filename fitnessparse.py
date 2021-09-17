import sys
import os
import csv

wordCheck = "-fitness\n"


with open("ringprofiles.csv", 'w', newline='') as csv_file: 
    writer = csv.writer(csv_file)
    for dirpath, dirnames, files in os.walk('.'):
        for file_name in files:
            if file_name.endswith(".lint"):
                direc = os.path.join(dirpath, file_name)
                input_file = open(direc,"r")
                exper = dirpath.split('/')[-3:-1]
                results = [exper]
                for line in input_file:
                    spl = line.split(' ')
                    if wordCheck in spl:
                        results.append(float(spl[0]))
                writer.writerow(results)

        
