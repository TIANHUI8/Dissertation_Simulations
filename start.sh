#!/bin/bash

#Create folders to contain each dataset
for i in {1..240}; do
if [ -d "$i" ]; then
rm -rf  $i
fi
mkdir  $i

#Copy necessary files
cp Sim_a.R $i;
cp run.pbs $i;
cd $i;

#Begin simulation
qsub run.pbs;
cd ..
done;
