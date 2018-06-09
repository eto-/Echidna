#!/bin/bash

if [ "$#" -eq "0" ]
  then
	echo "This script expects at least one argument: run number"
	echo "Optional second argument : calibration run number"
	echo "If no calibration run number is entered (or calibration run number '0'): last uploaded calibration run will be used"
	exit
fi

if [ "$1" -lt "5000" -o "$1" -gt "999999" ]
  then
	echo "Please enter a valid run number"
	exit
fi 

if [ -z "$2" ]   #check if second argument is set
   then
	#No calibration run number entered ; Last uploaded calibration run will be used
	calib_run=0
elif [ "$2" -eq "0" ]
  then
  	#Last uploaded calibration run used
	calib_run=0
elif [ "$2" -lt "5000" -o "$2" -gt "999999" ]
   then
  	echo "Please enter a valid calibration run number"
	echo "If you enter no calibration run number at all (or '0'), last uploaded calibration run will be used"
	exit
else	#if calib run is set, non-zero and valid
	calib_run=$2
fi


#find the necessary information for the normal run
WP=`find /bxstorage/rawdata/ -name \*$1\* | sed 's/\/bxstorage\/rawdata\///' | sed 's/\/Run.*//' | sed 's/\//,"/'`

#take care of a bug which produces double notation in $WP (noticed for an extra long run))
length=`echo ${#WP}`
if [ "$length" -gt "12" ]
   then
   	echo "Warning : Found path has too many characters."
	WP="$(echo $WP | cut -c 1-12)"	#Cut away the double statement
	echo "Reducing path to 12 characters : $WP "
fi

#find the necessary information for the calibration run
tempfile_path="/home/production/Echidna/dbMon/calibration_runs_temp.txt"
ls /bxstorage/rootfiles/cycle_15/*/*/ancillary/*laser_calibrations_c*.root > $tempfile_path

COMM="root -l -b -q /home/production/Run_Coordinator/BadChPlot.C+($1,$WP\",\"$tempfile_path\",$calib_run)"
$COMM

rm $tempfile_path
