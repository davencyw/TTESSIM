#!/bin/bash

#variables and parameters
OUTFOLDER="data/out/"

NCELLS=256
HEIGHT=33
DIAMETER=6
NUMCYCLES=100
TSTEPSPERC=10

TDS0=0.1
TDS1=0.1
TDS2=0.1
TDS3=0.1

INITTEMP=294
KF=0.1
KS=0.1
RHOS=0.1
RHOF=0.1
CF=0.1
CS=0.1
EPS=0.1

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
						-o $OUTFOLDER	\
						--tds0 $TDS0	\
						--tds1 $TDS1	\
						--tds2 $TDS2	\
						--tds3 $TDS3	\
						-c $NUMCYCLES 	\
						-b $TSTEPSPERC	\
						--kf $KF			\
						--ks $KS			\
						--rhof $RHOF		\
						--rhos $RHOS		\
						--cf $CF			\
						--cs $CS			\
						--epsilon $EPS	\