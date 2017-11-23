import os
import subprocess

LibPath = '../tsplib/'
ExecPath = './2opt '+LibPath

tsplib_file = raw_input()
output_file = raw_input()

with open(tsplib_file,"r") as f:
	for line in f:
		f1 = line.split(':')
		if(f1[0][-1] == ' '):
			f1[0] = f1[0][:-1]
	
		print(f1)
		f1tsp = f1[0]+".tsp"
		f1dist = float(f1[1])
	#-------Read Optimal Distance and TSPLIB instance-------------------------#
	#---------Now execute Lib if path exists----------------------------------#

		if os.path.isfile(LibPath+f1tsp):
			process = subprocess.Popen(ExecPath+f1tsp, shell=True, stdout=subprocess.PIPE)
			process.wait()
			foundDist = 0
			timeTaken = 0
			for lines in process.stdout:
				if 'Min distance' in lines:
					foundDist = float(lines.split(' ')[-1])
				if 'Time' in lines:
					timeTaken = float(lines.split(' ')[-1])

			firstline = os.path.isfile(output_file)

			with open(output_file, "a+") as of:
			        if not firstline:
	        	       		of.write("TSPInstance\tError Rate(%)\t\tTime of Execution\n")
			        error = 0.0

			        if f1dist!=0:
			                error = 100*(float(foundDist) - float(f1dist))/float(f1dist)

				of.write(f1[0]+"\t\t"+ "%.2f"%round(error,2) +"\t\t"+ str(timeTaken)+"\n")

			of.close()

f.close()
