'''
A Quick and Dirty solution to the tedious problem of having to open each Gaussian logfile
on its own.
This Program parses all gauss.log files and extracts the line that contains the HF energy.
The following line is also extracted since the energy sometimes continues over a
linebreak.
'''



import os

searchstring = "SCF Done:  E" # Change this if you want to search for something different
# Ask the User what type of scan is to be evaluated:
scanType = "invalid"
while scanType != "1" and scanType != "2":
    scanType = input("Handelt es sich um einen 1D oder 2D Scan? (1/2): ")

# Make a list of the subdirectories
subdirs = [x for x in os.listdir('.') if os.path.isdir(x)]
subdirs.sort() # Sorting is not necessary for the program, but it simplifies readability later on.

outputstring = ""
outputcsv = ""
if scanType == "1":
    outputstring += "Evaluating 1D Scan...\n"
    for directory in subdirs:
        os.chdir(directory)
        print("Parsing " + directory)
        outputstring += "\nIn " + directory + " found:\n"
        temp = open("gauss.log", "a")
        temp.close()
        with open("gauss.log", "r") as log:
            for line in log:
                if searchstring in line:
                    outputstring += line
                    for single in line.split():
                        try:
                            outputcsv += str(float(single))
                            outputcsv += "\n"
                            break
                        except ValueError:
                            dummy = 0
        os.chdir('..')

if scanType == "2":
    outputcsv = "Remember that the matrix is transposed!\n"
    outputstring += "Evaluating 2D Scan...\n"
    for directory in subdirs:
        print("Directory " + directory + ":")
        outputstring += "\nDirectory " + directory + ":"
        os.chdir(directory)
        subsubdir = [x for x in os.listdir('.') if os.path.isdir(x)]
        subsubdir.sort()
        for otherdirectory in subsubdir:
            os.chdir(otherdirectory)
            print("\tParsing " + otherdirectory)
            outputstring += "\n\tIn " + otherdirectory + " found:\n"
            temp = open("gauss.log", "a")
            temp.close()
            with open("gauss.log", "r") as log:
                for line in log:
                    if searchstring in line:
                        outputstring += "\t" + line
                        for single in line.split():
                            try:
                                outputcsv += str(float(single))
                                outputcsv += ";"
                                break
                            except ValueError:
                                dummy = 0
            os.chdir('..')
        outputcsv += "\n"
        os.chdir('..')

outfile = open("Energy_Summary.txt", "w")
outfile.write(outputstring)
outfile.close()

csvfile = open("Energy_Summary.csv", "w")
csvfile.write(outputcsv)
csvfile.close()