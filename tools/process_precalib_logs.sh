#! /bin/sh

mkdir -p empty
mkdir -p ok
mkdir -p no_gray_shift
mkdir -p zero_rms
mkdir -p bad
mkdir -p no_laben

for i in *; do 
   grep -q empty\ file $i
   if [ $? -ne 1 ]; then
     mv $i empty/
     continue
   fi

   grep -q laben\ detector\ is\ off $i
   if [ $? -ne 1 ]; then 
     mv $i no_laben/
     continue
   fi

   grep -q no\ gray\ shift $i
   if [ $? -ne 1 ]; then 
     mv $i no_gray_shift/
     continue
   fi

   grep -q \+\-0\) $i
   if [ $? -ne 1 ]; then 
     mv $i zero_rms
     continue
   fi

   grep -q precalibrations\ failed $i
   if [ $? -ne 1 ]; then 
     mv $i bad/
     continue
   fi

   grep -q too\ bad $i
   if [ $? -ne 1 ]; then 
     mv $i bad/
     continue
   fi

   grep -q 0\ errors $i
   if [ $? -ne 1 ]; then 
     mv $i ok/; 
     continue
   fi
done
