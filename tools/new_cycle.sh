#!/bin/sh

new_cycle=$1;
echo "cycling to $new_cycle"

uncycled=`grep -r ClassDef . |grep -v CYCLE_NUMBER |grep -v tools`
if [ -z "$a" ]; then 
  echo "the following files may not conform to cycle classdef policy, please check them"
  grep -r ClassDef . |grep -v CYCLE_NUMBER |grep -v tools
fi

TMPFILE=`mktemp /tmp/cycling.XXXXXXXXXX` 
if [ ! -e "$TMPFILE" ]; then
  echo "cannot create temp file"; 
  exit 1
fi

files=`grep -r CYCLE_NUMBER . |grep -v tools |cut -f1 -d:|sort -u`

for i in $files; do
  if [ `basename $i .hh` = `basename $i` ]; then
    echo "$i is not a .hh file, ignoring";
  else
    sed s/define\ CYCLE_NUMBER.*\$/define\ CYCLE_NUMBER\ $new_cycle/ $i > $TMPFILE
    cat $TMPFILE > $i
  fi
done

rm $TMPFILE
    
