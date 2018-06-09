#!/bin/bash

if [  "$#" -ne "1" ]
  then
	echo "This script expects only one argument: run number"
	exit
fi

if [ "$1" -lt "5000" -o "$1" -gt "999999" ]
  then
	echo "Please enter a valid run number"
	exit
fi 

WP=`find /bxstorage/rawdata/ -name \*$1\* | sed 's/\/bxstorage\/rawdata\///' | sed 's/\/Run.*//' | sed 's/\//,"/'`

COMM="root -l -b -q /home/production/Run_Coordinator/Run_Coordinator/BadChPlot.C+($1,$WP\")"

$COMM
