#!/bin/bash
# extract data from terachem job output files
for i in {1..1}
do
 	 echo $i
	 cd traj_$i
         cd step
	 cp ../../extract.sh .
         bash extract.sh
         cd ..
	 cd ..
done
