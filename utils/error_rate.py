od = 0
dist = 0
tsplib = ""
time = ""
filename = raw_input()
print("Enter output filename: ")
output_file = raw_input()

with open(filename, "r" ) as f:
	for lines in f:
		line = lines.split(' ')
		if "Optimal_Distance" in line:
			od = line[len(line)-1]
		elif "Min_distance_found" in line:
			dist = line[len(line)-1]
		elif "TSP_Instance" in line:
			tsplib = line[len(line)-1]
			if "\n" in tsplib:
				tsplib = tsplib[:-1]
		elif "Time_taken" in line:
			time = line[len(line)-1]
import os.path
firstline = os.path.isfile(output_file)	

with open(output_file, "a+") as f:
	if not firstline:
		f.write("TSPInstance\tError Rate(%)\t\tDistance\t\tTime of Execution\n\n")
	error = "NA"
	if od!=0:
		error = 100*(float(dist) - float(od))/float(od)
	f.write(tsplib+"\t\t"+ "%.2f"%round(error,2) +"\t\t"+ str(dist)+"\t\t"+str(time)+"\n")
