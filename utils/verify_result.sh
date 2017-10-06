#!/bin/bash

echo "Enter TSPLIB_Filename:";
read tsplib_file;
echo "Enter output filename:";
read log_file;

while read host; 
do
	file=`echo $host | cut -d":" -f1| cut -d" " -f1`;
	echo $file ;
	if [ -e ~/Downloads/tsp/TSPLIB\ instances/$file.tsp ]; 
	then
		./a.out ~/Downloads/tsp/TSPLIB\ instances/$file.tsp > $file.tmp; 
		echo "TSP_Instance : $file" > $file.ver;
		opt_dist=$(echo $host | rev | cut -d" " -f1 | rev);
		echo "Optimal_Distance : $opt_dist" >> $file.ver;
		time2=$(cat $file.tmp | grep -i "time" | rev | cut -d" " -f1 | rev);
		mindist=$(cat $file.tmp | grep -i "min dist" | rev | cut -d" " -f1 | rev);
		rm $file.tmp;
		echo "Min_distance_found : $mindist" >> $file.ver;
		echo "Time_taken : $time2" >> $file.ver;
		echo "$file.ver" > $file.name;
		echo $log_file >> $file.name;
		python error_rate.py < $file.name;
		rm $file.name;
		rm $file.ver;
	fi;
done < $tsplib_file;
