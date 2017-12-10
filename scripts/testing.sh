#!/bin/bash

#variables and parameters
OUTFOLDER="../data/testing/"

#TODO(dave): SPECIFY TESTING PARAMETERS!!

NCELLS=1000
HEIGHT=1
DIAMETER=6
NUMCYCLES=100
TSTEPSPERC=10
OPS=1

TDS0=0.1
TDS1=0.1
TDS2=0.1
TDS3=0.1

DELTAT=1e-5

INITTEMP=773
KF=1
KS=1
RHOS=1
RHOF=1
CF=1
CS=1
EPS=0.4
UF=100
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
../build/bin/ttessim_testing			\
						--ops $OPS		\
						-N $NCELLS		\
						-h $HEIGHT		\
						-d $DIAMETER	\
						-t $INITTEMP	\
						-c $NUMCYCLES 	\
						-b $TSTEPSPERC	\
						--folder $OUTFOLDER	\
						--tds0 $TDS0	\
						--tds1 $TDS1	\
						--tds2 $TDS2	\
						--tds3 $TDS3	\
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


#TODO(Dave): WRITE PLOTTING SCRIPTS FOR ORDER VERIFICATION-STUDIES FOR FLUID AND SOLID