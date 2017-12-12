#!/bin/bash

#variables and parameters
OUTFOLDER="../data/out/"


NCELLS=512
HEIGHT=1.492
DIAMETER=8
NUMCYCLES=10
TSTEPSPERC=0
OPS=0

TDS0=21600
TDS1=21600
TDS2=21600
TDS3=21600

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
UF=2.709e-5
HF=960

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
						--hs $HF 	    \