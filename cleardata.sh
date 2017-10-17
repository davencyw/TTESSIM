#!/bin/bash

#___________________________

#terminal notification
echo "______\n"
echo "cleardata.sh"
echo "______\n"

for file in $(find ./data/out/ -name "*.csv"); do
	echo "deleting " $file
done

rm data/out/*