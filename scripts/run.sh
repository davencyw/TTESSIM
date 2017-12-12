#!/bin/bash

#variables and parameters
OUTFOLDER="../data/out/"


NCELLS=512
HEIGHT=5.968
DIAMETER=4
NUMCYCLES=10
TSTEPSPERC=0
OPS=0

TDS0=1
TDS1=1
TDS2=1
TDS3=1

DELTAT=0.001

INITTEMP=288.15
TEMPC=873
TEMPD=293

KF=0.52
KS=2.0
RHOS=2600
RHOF=1835.6
CF=1511.8
CS=900.0
EPS=0.4
UF=0.000108
HF=480
HS=480

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
						--tc $TEMPC		\
						--td $TEMPD		\
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