#!/bin/bash

#variables and parameters
OUTFOLDER="../data/out/"


NCELLS=512
HEIGHT=33
DIAMETER=6
NUMCYCLES=10
TSTEPSPERC=10
OPS=1

TDS0=0.1
TDS1=0.1
TDS2=0.1
TDS3=0.1

DELTAT=0.001

INITTEMP=294
KF=0.1
KS=0.1
RHOS=0.1
RHOF=0.1
CF=0.1
CS=0.1
EPS=0.1
UF=10
HF=1
HS=1

#___________________________

# sh make.sh

#terminal notification
# echo "______\n"
# echo "run.sh"
# echo "______\n\n"

#run
clear
../build/bin/ttessim					\
						--ops $OPS		\
						-N $NCELLS		\
						-h $HEIGHT		\
						-d $DIAMETER	\
						-t $INITTEMP	\
						--folder $OUTFOLDER	\
						--tds0 $TDS0	\
						--tds1 $TDS1	\
						--tds2 $TDS2	\
						--tds3 $TDS3	\
						-c $NUMCYCLES 	\
						-b $TSTEPSPERC	\
						--kf $KF		\
						--ks $KS		\
						--rhof $RHOF	\
						--rhos $RHOS	\
						--cf $CF		\
						--cs $CS		\
						--epsilon $EPS	\
						--uf $UF		\
						--dt $DELTAT 	\
						--hf $HF 	    \
						--hs $HS 	    \