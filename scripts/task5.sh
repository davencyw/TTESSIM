#!/bin/bash

#variables and parameters
OUTFOLDER="../data/out/"

#TODO(dave): hf, hs, uf!!

NCELLS=512
HEIGHT=1
DIAMETER=1
NUMCYCLES=1
TSTEPSPERC=0
OPS=0

TDS0=5000
TDS1=0
TDS2=0
TDS3=0

DELTAT=0.00001

INITTEMP=293
TEMPC=873
TEMPD=293

KF=2.0
KS=0.52
RHOS=2600
RHOF=1835.6
CF=1511.8
CS=900.0
EPS=0.4
UF=0.00001734
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