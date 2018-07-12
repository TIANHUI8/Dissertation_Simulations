#!/bin/bash

#Create folders to contain each dataset
for i in {1..240}; do
if [ -d "$i" ]; then
rm -rf  $i
fi
done;
