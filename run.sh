#!/bin/bash

#variables and parameters
NCELLS=256
HEIGHT=33
DIAMETER=6
INITTEMP=294
OUTFOLDER="data/out/"

#___________________________

#terminal notification
echo "______\n"
echo "run.sh"
echo "______\n\n"

sh make.sh

#run
./build/bin/ttessim						\
						-N $NCELLS		\
						-h $HEIGHT		\
						-d $DIAMETER	\
						-t $INITTEMP	\
						-o $OUTFOLDER